#!/usr/bin/env python3
# Copyright (c) 2026 Renjie Yu and Daniel Oron
# SPDX-License-Identifier: MIT
"""Paired model comparison on an existing, frozen gamma-sensitivity panel.

Nominal objectives are context, not out-of-sample performance. The nominal
assignment is evaluated with the same exact recourse verifier as production LBBD.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import pathlib
import shutil
import statistics
import subprocess
import sys

from run_sensitivity import atomic_csv, atomic_json, read_csv, read_instance, split_values, finite_real

ROOT = pathlib.Path(__file__).resolve().parents[1]
MODELS = ("static-robust", "nominal-assignment")
SCHEMA = "agv-model-comparison-v1"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def solution(path, data, row=None):
    tokens = path.read_text().split()
    n, upper, lower = map(int, tokens[:3])
    assignment = list(map(int, tokens[3:3+n]))
    if n != data["n"] or len(assignment) != n or any(k < 0 or k >= data["m"] for k in assignment):
        raise ValueError(f"invalid assignment in {path}")
    if not 0 <= lower <= upper < 1_000_000_000:
        raise ValueError(f"invalid bounds in {path}")
    if row and (lower != int(row["final_lower"]) or upper != int(row["final_upper"])):
        raise ValueError(f"source result and solution disagree: {path}")
    return lower, upper, assignment


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--reference-run", type=pathlib.Path, default=ROOT / "runs/sensitivity_gamma_n40")
    p.add_argument("--output", type=pathlib.Path, required=True)
    p.add_argument("--build", type=pathlib.Path, default=ROOT / "build")
    p.add_argument("--models", nargs="+", choices=MODELS, default=list(MODELS))
    p.add_argument("--gamma-values", nargs="+", help="existing budget levels, space- or comma-separated; default: all")
    p.add_argument("--replications", nargs="+", help="existing replication numbers; default: all")
    p.add_argument("--total-limit", type=finite_real, default=1800, help="per static robust optimization")
    p.add_argument("--evaluation-limit", type=finite_real, default=1800, help="per fixed nominal-assignment robust evaluation")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--resume", action="store_true")
    p.add_argument("--prepare-only", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--summarize-only", action="store_true")
    return p


def panel(args):
    reference = args.reference_run.resolve()
    source_design = json.loads((reference / "design.json").read_text())
    if source_design["study"] != "gamma" or source_design["configuration"] != "lbbd-hig":
        raise ValueError("reference must be a gamma-sensitivity campaign using the full lbbd-hig method")
    rows = read_csv(reference / "sensitivity_results.csv")
    if not rows:
        raise ValueError("reference campaign has no completed results")
    indexed = {}
    for row in rows:
        key = (float(row["gamma"]), int(row["replication"]))
        if row["configuration"] != "lbbd-hig" or key in indexed:
            raise ValueError("ambiguous/nonproduction results in reference panel")
        indexed[key] = row
    gammas = split_values(args.gamma_values, finite_real, "--gamma-values")
    reps = split_values(args.replications, int, "--replications")
    gammas = sorted({key[0] for key in indexed}) if gammas is None else sorted(gammas)
    reps = sorted({key[1] for key in indexed}) if reps is None else sorted(reps)
    if any(g < 0 for g in gammas) or any(r < 1 for r in reps):
        raise ValueError("budgets must be nonnegative and replications positive")
    if args.total_limit <= 0 or args.evaluation_limit <= 0 or not 0 <= args.seed <= 2_000_000_000:
        raise ValueError("invalid time limit or solver seed")
    if len(set(args.models)) != len(args.models):
        raise ValueError("duplicate model")
    frozen = json.loads((reference / "solver/binaries.json").read_text())
    for name in ("agv_hig", "agv_lbbd"):
        if sha(reference / "solver/bin" / name) != frozen[name]:
            raise ValueError(f"reference binary fingerprint mismatch: {name}")
    selected = []
    for replication in reps:
        if (0.0, replication) not in indexed:
            raise ValueError(f"missing Gamma=0 nominal-model solve for replication {replication}")
        nominal_row = indexed[(0.0, replication)]
        nominal_path = reference / "prepared_instances" / nominal_row["source_file"]
        nominal_data = read_instance(nominal_path)
        nominal_solver_instance = reference / "solver/instances" / f"{nominal_row['instance']}.agv"
        if sha(nominal_path) != nominal_row["instance_sha256"] or sha(nominal_solver_instance) != sha(nominal_path):
            raise ValueError("nominal reference is not the instance used by its saved solver run")
        nominal_sol = reference / "solver/cases" / nominal_row["instance"] / "lbbd-hig/lbbd.sol"
        nominal_lower, nominal_upper, _ = solution(nominal_sol, nominal_data, nominal_row)
        if nominal_lower != nominal_upper or nominal_row["optimal"] != "1":
            raise ValueError("nominal assignment must come from a certified optimal Gamma=0 solve")
        for gamma in gammas:
            if (gamma, replication) not in indexed:
                raise ValueError(f"missing reference case Gamma={gamma}, replication={replication}")
            row = indexed[(gamma, replication)]
            instance = reference / "prepared_instances" / row["source_file"]
            data = read_instance(instance)
            if data["gamma"] != gamma or any(data[key] != nominal_data[key] for key in data if key != "gamma"):
                raise ValueError(f"paired instances differ in more than Gamma: {instance}")
            if sha(instance) != row["instance_sha256"]:
                raise ValueError(f"reference instance fingerprint mismatch: {instance}")
            solver_instance = reference / "solver/instances" / f"{row['instance']}.agv"
            if sha(solver_instance) != sha(instance):
                raise ValueError(f"reference differs from the instance actually solved: {instance}")
            robust_sol = reference / "solver/cases" / row["instance"] / "lbbd-hig/lbbd.sol"
            lower, upper, robust_assignment = solution(robust_sol, data, row)
            _, _, nominal_assignment = solution(nominal_sol, nominal_data, nominal_row)
            # An interrupted tree can have a weaker dual bound than a solved
            # nominal model. That is legitimate; the nominal optimum is also
            # a valid robust LB and can tighten the imported interval.
            source_lower = lower
            lower = max(lower, nominal_lower)
            if lower > upper:
                raise ValueError("nominal optimum contradicts the robust reference upper bound")
            selected.append({
                "case": f"g{gamma:g}_r{replication:02d}", "gamma": gamma,
                "replication": replication, "n": data["n"], "m": data["m"], "alpha": data["alpha"],
                "nominal_optimum": nominal_upper, "two_stage_lower": lower, "two_stage_upper": upper,
                "two_stage_status": row["status"], "two_stage_seconds": float(row["wall_seconds"]),
                "two_stage_source_lower": source_lower,
                "nominal_vs_two_stage_reassigned_jobs": sum(k != h for k, h in zip(nominal_assignment, robust_assignment)),
                "source_case": row["instance"], "instance_source": str(instance),
                "nominal_solution_source": str(nominal_sol), "robust_solution_source": str(robust_sol),
                "instance_sha256": sha(instance), "nominal_solution_sha256": sha(nominal_sol),
                "robust_solution_sha256": sha(robust_sol),
            })
    selected.sort(key=lambda c: (c["gamma"], c["replication"]))
    design = {
        "schema": SCHEMA, "reference_run": str(reference), "models": sorted(args.models),
        "gamma_values": gammas, "replications": reps, "total_limit": args.total_limit,
        "evaluation_limit": args.evaluation_limit, "seed": args.seed,
        "reference_design": source_design, "reference_binaries": frozen, "cases": selected,
        "nominal_tie_rule": "reuse the recorded optimal assignment, with no selection using uncertain outcomes",
        "static_warm_start": "best certified static block construction of the nominal and two-stage assignments",
        "bound_credit": "reused two-stage lower bound is also valid for both restricted comparison models",
    }
    return design


def prepare(output, design, resume):
    target = output / "design.json"
    if target.exists():
        if not resume:
            raise ValueError("output exists; use --resume with identical settings or a new output directory")
        if json.loads(target.read_text()) != design:
            raise ValueError("resume settings or reference inputs changed; use a new output directory")
    elif output.exists() and any(output.iterdir()):
        raise ValueError("output is nonempty and has no model-comparison design")
    output.mkdir(parents=True, exist_ok=True)
    for case in design["cases"]:
        folder = output / "cases" / case["case"]
        folder.mkdir(parents=True, exist_ok=True)
        for source, name, key in (
            (case["instance_source"], "instance.agv", "instance_sha256"),
            (case["nominal_solution_source"], "nominal.sol", "nominal_solution_sha256"),
            (case["robust_solution_source"], "two_stage.sol", "robust_solution_sha256"),
        ):
            path = folder / name
            if path.exists():
                if sha(path) != case[key]:
                    raise ValueError(f"prepared input changed: {path}")
            else:
                shutil.copy2(source, path)
    atomic_json(target, design)


def completed_result(prefix, model, case, require_done=True):
    path = prefix.with_suffix(".csv")
    marker = prefix.parent / f"{model}_done.json"
    if require_done:
        if not marker.exists():
            return None
        for name, fingerprint in json.loads(marker.read_text()).items():
            if sha(prefix.parent / name) != fingerprint:
                raise ValueError(f"completed certificate changed: {prefix.parent / name}")
    if not path.exists():
        return None
    rows = read_csv(path)
    if len(rows) != 1 or rows[0]["model"] != model:
        raise ValueError(f"invalid result file {path}")
    row = rows[0]
    lower, upper, assignment = solution(prefix.with_suffix(".sol"), {"n": case["n"], "m": case["m"]})
    if (lower, upper) != (int(row["lower"]), int(row["upper"])) or row["optimal"] != str(int(lower == upper)):
        raise ValueError(f"result/certificate disagreement: {path}")
    if lower < case["two_stage_lower"]:
        raise ValueError(f"comparison lost its valid initial lower bound: {path}")
    if model == "nominal-assignment":
        _, _, nominal = solution(prefix.parent / "nominal.sol", {"n": case["n"], "m": case["m"]})
        if assignment != nominal:
            raise ValueError("nominal assignment was changed during evaluation")
    return row


def summarize(output, design):
    results, paired = [], []
    for case in design["cases"]:
        folder = output / "cases" / case["case"]
        for name, key in (("instance.agv", "instance_sha256"), ("nominal.sol", "nominal_solution_sha256"),
                          ("two_stage.sol", "robust_solution_sha256")):
            if sha(folder / name) != case[key]:
                raise ValueError(f"frozen input changed: {folder / name}")
        base = {k: case[k] for k in ("case", "n", "m", "gamma", "alpha", "replication", "nominal_optimum",
                                    "nominal_vs_two_stage_reassigned_jobs")}
        base.update({"two_stage_lower": case["two_stage_lower"], "two_stage_upper": case["two_stage_upper"]})
        comparison = dict(base)
        for model in MODELS:
            prefix = output / "cases" / case["case"] / model
            row = completed_result(prefix, model, case)
            if row is None:
                continue
            results.append({**base, **row})
            name = "static" if model == "static-robust" else "nominal_assignment"
            comparison[f"{name}_lower"] = int(row["lower"])
            comparison[f"{name}_upper"] = int(row["upper"])
            comparison[f"{name}_optimal"] = int(row["optimal"])
            # Interval for restricted optimum - two-stage optimum. Nonnegative
            # nesting is mathematical, but a positive lower endpoint is evidence.
            comparison[f"{name}_penalty_lower"] = max(0, int(row["lower"])-case["two_stage_upper"])
            comparison[f"{name}_penalty_upper"] = int(row["upper"])-case["two_stage_lower"]
            if comparison[f"{name}_penalty_upper"] < 0:
                raise ValueError("model nesting violated by recorded bounds")
            comparison[f"{name}_penalty_exact"] = int(
                int(row["lower"]) == int(row["upper"]) and case["two_stage_lower"] == case["two_stage_upper"])
            if case["gamma"] == 0 and not int(row["lower"]) <= case["nominal_optimum"] <= int(row["upper"]):
                raise ValueError("Gamma=0 model equivalence failed")
        paired.append(comparison)
    # Monotonicity is an invariant of these nested uncertainty sets. Overlapping
    # interrupted intervals are inconclusive and do not count as violations.
    for replication in design["replications"]:
        group = sorted((r for r in paired if r["replication"] == replication), key=lambda r: r["gamma"])
        for name in ("two_stage", "static", "nominal_assignment"):
            previous_lower = 0
            for row in group:
                if f"{name}_lower" not in row:
                    continue
                if row[f"{name}_upper"] < previous_lower:
                    raise ValueError(f"certified nonmonotonic {name} results for replication {replication}")
                previous_lower = max(previous_lower, row[f"{name}_lower"])
    if results:
        fields = list(dict.fromkeys(k for row in results for k in row))
        atomic_csv(output / "model_results.csv", fields, results)
    fields = list(dict.fromkeys(k for row in paired for k in row))
    atomic_csv(output / "paired_comparison.csv", fields, paired)
    grouped = []
    for gamma in design["gamma_values"]:
        group = [r for r in paired if r["gamma"] == gamma]
        summary = {"gamma": gamma, "replications": len(group),
                   "mean_nominal_optimum": statistics.mean(r["nominal_optimum"] for r in group),
                   "mean_two_stage_lower": statistics.mean(r["two_stage_lower"] for r in group),
                   "mean_two_stage_upper": statistics.mean(r["two_stage_upper"] for r in group)}
        for name in ("static", "nominal_assignment"):
            available = [r for r in group if f"{name}_lower" in r]
            summary[f"{name}_completed"] = len(available)
            summary[f"{name}_optimal"] = sum(r[f"{name}_optimal"] for r in available)
            summary[f"{name}_penalty_exact"] = sum(r[f"{name}_penalty_exact"] for r in available)
            summary[f"{name}_positive_penalty_certified"] = sum(r[f"{name}_penalty_lower"] > 0 for r in available)
            # Do not silently average a partial or easy-only sample.
            if len(available) == len(group):
                for suffix in ("lower", "upper", "penalty_lower", "penalty_upper"):
                    summary[f"mean_{name}_{suffix}"] = statistics.mean(r[f"{name}_{suffix}"] for r in available)
        grouped.append(summary)
    atomic_csv(output / "model_summary.csv", list(dict.fromkeys(k for r in grouped for k in r)), grouped)
    return len(results)


def main():
    args = parser().parse_args()
    if args.dry_run and (args.resume or args.prepare_only or args.summarize_only):
        raise ValueError("--dry-run cannot be combined with other run modes")
    if args.summarize_only:
        design = json.loads((args.output / "design.json").read_text())
        count = summarize(args.output, design)
        print(f"summarized {count} completed comparison solves")
        return
    design = panel(args)
    count = len(design["cases"])*len(design["models"])
    print(f"{len(design['replications'])} paired replications, {len(design['gamma_values'])} budgets; "
          f"{count} additional solves; nominal and full LBBD-HIG results reused", flush=True)
    if args.dry_run:
        return
    output = args.output.resolve()
    prepare(output, design, args.resume)
    if args.prepare_only:
        summarize(output, design)
        return
    build = args.build.resolve()
    frozen = output / "bin/agv_model_compare"
    if frozen.exists():
        manifest = json.loads((output / "binaries.json").read_text())
        if sha(frozen) != manifest["agv_model_compare"]:
            raise ValueError("frozen comparison binary changed")
    else:
        # GNU make splits target lists on whitespace, even when subprocess
        # passes BUILD as one argument. Use the relative build name so the
        # spaces in this project's absolute path never enter make's targets.
        make_build = os.path.relpath(build, ROOT)
        if any(character.isspace() for character in make_build):
            raise ValueError("--build must have a path relative to the project without whitespace")
        subprocess.run(["make", f"BUILD={make_build}", f"{make_build}/agv_model_compare"], cwd=ROOT, check=True)
        frozen.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(build / "agv_model_compare", frozen)
        atomic_json(output / "binaries.json", {"agv_model_compare": sha(frozen)})
        # Save the exact sources used for this model-comparison executable.
        sources = {str(p.relative_to(ROOT)): sha(p) for p in (ROOT / "src").rglob("*") if p.is_file()}
        sources["scripts/run_model_comparison.py"] = sha(pathlib.Path(__file__))
        atomic_json(output / "source_hashes.json", sources)
    env = os.environ.copy();env["LC_ALL"] = "C"
    for case in design["cases"]:
        folder = output / "cases" / case["case"]
        for model in design["models"]:
            prefix = folder / model
            if completed_result(prefix, model, case) is not None:
                continue
            limit = args.total_limit if model == "static-robust" else args.evaluation_limit
            command = [str(frozen), model, str(folder / "instance.agv"), str(limit), str(prefix),
                       str(args.seed), str(folder / "nominal.sol"), str(folder / "two_stage.sol"),
                       str(case["two_stage_lower"])]
            atomic_json(folder / f"{model}_command.json", command)
            print(f"{case['case']} {model} limit={limit:g}s", flush=True)
            # Results become visible to resume only after a successful process.
            for extension in (".csv", ".sol"):
                prefix.with_suffix(extension).unlink(missing_ok=True)
            with (folder / f"{model}.log").open("w") as log:
                done = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT, check=False)
            if done.returncode:
                prefix.with_suffix(".csv").unlink(missing_ok=True)
                raise ValueError(f"{model} failed; inspect {folder / (model+'.log')}")
            row = completed_result(prefix, model, case, require_done=False)
            if row is None:
                raise ValueError("successful process wrote no result")
            certificate_paths = [prefix.with_suffix(".csv"), prefix.with_suffix(".sol")]
            if model == "static-robust":
                certificate_paths.append(folder / f"{model}_schedule.csv")
            atomic_json(folder / f"{model}_done.json", {p.name: sha(p) for p in certificate_paths})
            summarize(output, design)
            print(f"  [{row['lower']}, {row['upper']}] {row['status']}", flush=True)
    summarize(output, design)
    print(f"results saved to {output}")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, KeyError, subprocess.CalledProcessError) as error:
        print(f"error: {error}", file=sys.stderr)
        sys.exit(1)
