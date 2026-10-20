#!/usr/bin/env python3
# Copyright (c) 2026 Renjie Yu and Daniel Oron
# SPDX-License-Identifier: MIT
"""Rerun selected published instances without regenerating their random data."""

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
METHODS = ("lbbd-basic", "lbbd-hig", "milp-ccg-basic", "milp-ccg-hig", "cp-basic", "cp-hig")
PANELS = {
    "gamma": ["sensitivity_gamma_n40"],
    "crossed": ["sensitivity_crossed_n40_gamma0.5", "sensitivity_crossed_n40_gamma2", "sensitivity_crossed_n40_gamma4"],
    "technology": ["sensitivity_technology_n40"],
}


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def selected(row, key, choices):
    return choices is None or float(row[key]) in choices


def run(command):
    environment = os.environ.copy()
    environment["LC_ALL"] = "C"
    subprocess.run([str(value) for value in command], cwd=ROOT, env=environment, check=True)


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("experiment", choices=("benchmark", "sensitivity", "models"))
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--build", type=Path, default=ROOT / "build")
    parser.add_argument("--jobs", nargs="+", type=int)
    parser.add_argument("--agvs", nargs="+", type=int)
    parser.add_argument("--gamma", nargs="+", type=float)
    parser.add_argument("--beta", nargs="+", type=float)
    parser.add_argument("--replications", nargs="+", type=int)
    parser.add_argument("--methods", nargs="+", choices=METHODS, default=["lbbd-hig"])
    parser.add_argument("--study", choices=PANELS, default="gamma")
    parser.add_argument("--deviations", nargs="+", type=float)
    parser.add_argument("--capacity-scales", nargs="+", type=float)
    parser.add_argument("--recharge-factors", nargs="+", type=float)
    parser.add_argument("--models", nargs="+", choices=("static-robust", "nominal-assignment"), default=["static-robust", "nominal-assignment"])
    parser.add_argument("--total-limit", type=float, default=1800)
    parser.add_argument("--hig-limit", type=float, default=60)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    if args.experiment == "models":
        args.hig_limit = 0
    if not math.isfinite(args.total_limit) or args.total_limit <= 0 or not math.isfinite(args.hig_limit) or not 0 <= args.hig_limit <= args.total_limit or args.seed < 0:
        parser.error("limits must be finite with 0 <= HIG <= total, and seed must be nonnegative")
    for key in ("methods", "models"):
        if len(getattr(args, key)) != len(set(getattr(args, key))):
            parser.error(f"duplicate --{key}")
    for key in ("gamma", "beta", "deviations", "capacity_scales", "recharge_factors"):
        if getattr(args, key) and any(not math.isfinite(v) for v in getattr(args, key)):
            parser.error(f"nonfinite --{key.replace('_', '-')}")
    rows = []
    if args.experiment == "benchmark":
        for row in read_csv(ROOT / "data/benchmark/specifications.csv"):
            row["input"] = str(ROOT / "data/benchmark/instances" / (row["case"] + ".agv"))
            row["file"] = row["case"] + ".agv"
            rows.append(row)
    elif args.experiment == "sensitivity":
        for panel in PANELS[args.study]:
            for row in read_csv(ROOT / "data" / panel / "sensitivity_specifications.csv"):
                row["input"] = str(ROOT / "data" / panel / "prepared_instances" / row["source_file"])
                row["file"] = panel + "__" + row["source_file"]
                rows.append(row)
    else:
        rows = read_csv(ROOT / "results/recorded/model_comparison_n30_m3/paired_comparison.csv")
    filters = [("n", args.jobs), ("m", args.agvs), ("gamma", args.gamma), ("replication", args.replications)]
    if args.experiment == "benchmark":
        filters.append(("tardiness_factor", args.beta))
    if args.experiment == "sensitivity":
        filters += [("deviation_level", args.deviations), ("capacity_scale", args.capacity_scales), ("recharge_factor", args.recharge_factors)]
    rows = [row for row in rows if all(selected(row, key, choices) for key, choices in filters)]
    if not rows:
        parser.error("no archived instances match the selected parameters")
    configurations = args.models if args.experiment == "models" else args.methods
    print(f"{len(rows)} instances, {len(rows) * len(configurations)} solver runs", flush=True)
    if args.dry_run:
        return
    output = args.output.resolve()
    if output.exists() and any(output.iterdir()) and not args.resume:
        parser.error("output is nonempty; use --resume with unchanged settings or choose a new directory")
    output.mkdir(parents=True, exist_ok=True)
    design = {"experiment": args.experiment, "configurations": configurations, "total_limit": args.total_limit,
              "hig_limit": args.hig_limit, "seed": args.seed, "cases": []}
    for row in rows:
        if args.experiment == "models":
            folder = ROOT / "data/model_comparison/cases" / row["case"]
            inputs = {name: digest(folder / name) for name in ("instance.agv", "nominal.sol", "two_stage.sol")}
            design["cases"].append({"case": row["case"], "inputs": inputs})
        else:
            design["cases"].append({"case": row["file"], "sha256": digest(Path(row["input"]))})
    manifest = output / "selection.json"
    if manifest.exists() and json.loads(manifest.read_text()) != design:
        parser.error("selected instances or settings differ from the saved run")
    write_json(manifest, design)
    build = args.build.resolve()
    if args.experiment == "models":
        if args.prepare_only:
            return
        binary_record = output / "model_binary.json"
        binary_hash = {"agv_model_compare": digest(build / "agv_model_compare")}
        if binary_record.exists() and json.loads(binary_record.read_text()) != binary_hash:
            raise ValueError("model-comparison executable changed during the run")
        write_json(binary_record, binary_hash)
        for row in rows:
            source = ROOT / "data/model_comparison/cases" / row["case"]
            target = output / row["case"]
            target.mkdir(exist_ok=True)
            bound = max(int((source / name).read_text().split()[2]) for name in ("nominal.sol", "two_stage.sol"))
            for model in args.models:
                prefix = target / model
                marker = target / (model + ".done.json")
                if marker.exists():
                    saved = json.loads(marker.read_text())
                    if any(digest(target / name) != checksum for name, checksum in saved.items()):
                        raise ValueError("completed model output was changed")
                    continue
                run([build / "agv_model_compare", model, source / "instance.agv", args.total_limit,
                     prefix, args.seed, source / "nominal.sol", source / "two_stage.sol", bound])
                result = read_csv(prefix.with_suffix(".csv"))[0]
                if int(result["lower"]) > int(result["upper"]):
                    raise ValueError("invalid comparison bounds")
                write_json(marker, {p.name: digest(p) for p in (prefix.with_suffix(".csv"), prefix.with_suffix(".sol"))})
        return
    instances = output / "selected_instances"
    instances.mkdir(exist_ok=True)
    expected = {row["file"] for row in rows}
    if {p.name for p in instances.glob("*.agv")} - expected:
        raise ValueError("unexpected files in selected_instances")
    for row in rows:
        destination = instances / row["file"]
        if destination.exists() and digest(destination) != digest(Path(row["input"])):
            raise ValueError("a selected input was changed")
        shutil.copy2(row["input"], destination)
    command = [sys.executable, ROOT / "scripts/run_experiment.py", "--instances", instances,
               "--output", output / "solver", "--build", build, "--configurations", *args.methods,
               "--total-limit", args.total_limit, "--hig-limit", args.hig_limit, "--seed", args.seed]
    if args.resume:
        command.append("--resume")
    if args.prepare_only:
        command.append("--prepare-only")
    run(command)
    result_path = output / "solver/results.csv"
    if result_path.exists():
        metadata = {row["file"]: row for row in rows}
        enriched = []
        for row in read_csv(result_path):
            extras = {key: value for key, value in metadata[row["source_file"]].items() if key not in ("input", "file", "case", "source_file")}
            enriched.append({**row, **extras})
        fields = list(dict.fromkeys(key for row in enriched for key in row))
        with (output / "results_with_factors.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(enriched)


if __name__ == "__main__":
    main()
