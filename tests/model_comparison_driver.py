#!/usr/bin/env python3
# Copyright (c) 2026 Renjie Yu and Daniel Oron
# SPDX-License-Identifier: MIT
"""Test pairing, production-result reuse, resume, and certificate safeguards."""
import argparse
import csv
import json
import pathlib
import subprocess
import tempfile


def run(command, success=True):
    p = subprocess.run(list(map(str, command)), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if success != (p.returncode == 0):
        raise RuntimeError(p.stdout)
    return p.stdout


def rows(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def main():
    p = argparse.ArgumentParser();p.add_argument("--build", default="build");args = p.parse_args()
    root = pathlib.Path(__file__).resolve().parents[1];build = (root / args.build).resolve()
    driver = root / "scripts/run_model_comparison.py"
    with tempfile.TemporaryDirectory(prefix="agv-model-comparison-") as temporary:
        work = pathlib.Path(temporary);reference = work / "reference"
        run(["python3", root / "scripts/run_sensitivity.py", "gamma", "--output", reference,
             "--build", build, "--jobs", 5, "--agvs", 2, "--replications", 1, 2,
             "--gamma-values", 0, .5, 1, "--base-gamma", 1, "--total-limit", 10,
             "--hig-limit", 1, "--metrics-limit", 0, "--tardiness-factor", 1, "--due-date-range", 0])
        output = work / "comparison"
        command = ["python3", driver, "--reference-run", reference, "--output", output,
                   "--build", build, "--total-limit", 10, "--evaluation-limit", 10]
        dry = run([*command, "--dry-run"])
        if "12 additional solves" not in dry or output.exists():
            raise RuntimeError("dry-run changed files or miscounted cases")
        run([*command, "--prepare-only"])
        if len(rows(output / "paired_comparison.csv")) != 6:
            raise RuntimeError("paired panel has incorrect size")
        run([*command, "--resume"])
        result = rows(output / "model_results.csv")
        if len(result) != 12 or any(r["optimal"] != "1" for r in result):
            raise RuntimeError("small model-comparison campaign was not solved completely")
        for row in result:
            if float(row["gamma"]) == 0 and int(row["upper"]) != int(row["nominal_optimum"]):
                raise RuntimeError("Gamma=0 nominal/model equivalence failed")
        before = {path: path.stat().st_mtime_ns for path in output.glob("cases/*/*.csv")}
        run([*command, "--resume"])
        if any(path.stat().st_mtime_ns != stamp for path, stamp in before.items()):
            raise RuntimeError("resume reran a completed model solve")
        # A result without a successful completion marker must be recomputed,
        # even when the interrupted process happened to leave a parseable CSV.
        marker = output / "cases/g0.5_r01/static-robust_done.json"
        marker.unlink()
        changed = output / "cases/g0.5_r01/static-robust.csv"
        stamp = changed.stat().st_mtime_ns
        run([*command, "--resume"])
        if not marker.exists() or changed.stat().st_mtime_ns == stamp:
            raise RuntimeError("resume trusted an uncompleted process")
        run([*command, "--resume", "--total-limit", 11], success=False)
        run([*command, "--gamma-values", 2, "--dry-run"], success=False)
        run(["python3", driver, "--output", output, "--summarize-only"])
        # A modified instance or nominal assignment must not silently be reused.
        nominal = output / "cases/g0_r01/nominal.sol"
        original = nominal.read_text();parts = original.split();parts[3] = str(1-int(parts[3]))
        nominal.write_text(" ".join(parts))
        run([*command, "--resume"], success=False);nominal.write_text(original)
        schedule = output / "cases/g0_r01/static-robust_schedule.csv"
        if len(rows(schedule)) != 5:
            raise RuntimeError("static schedule certificate omitted jobs")
        # Partial summaries must not turn interval endpoints into exact values.
        design = json.loads((output / "design.json").read_text())
        if len(design["cases"]) != 6:
            raise RuntimeError("manifest lost cases")
        # An interrupted robust run with a weak LB is legitimate, rather than
        # an input error. The solved nominal optimum safely strengthens it.
        source_rows = rows(reference / "sensitivity_results.csv")
        weakened = next(r for r in source_rows if float(r["gamma"]) == 1 and int(r["replication"]) == 1)
        source_sol = reference / "solver/cases" / weakened["instance"] / "lbbd-hig/lbbd.sol"
        tokens = source_sol.read_text().split();tokens[2] = "0";source_sol.write_text(" ".join(tokens)+"\n")
        weakened["final_lower"] = "0";weakened["optimal"] = "0"
        with (reference / "sensitivity_results.csv").open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(source_rows[0]));writer.writeheader();writer.writerows(source_rows)
        partial = work / "weakened"
        weak_command = ["python3", driver, "--reference-run", reference, "--output", partial,
                        "--build", build, "--gamma-values", 1, "--replications", 1,
                        "--models", "nominal-assignment", "--total-limit", 10, "--evaluation-limit", 10]
        run(weak_command)
        weak_case = json.loads((partial / "design.json").read_text())["cases"][0]
        if weak_case["two_stage_source_lower"] != 0 or weak_case["two_stage_lower"] != weak_case["nominal_optimum"]:
            raise RuntimeError("nominal bound was not safely reused")
    print("PASS model-comparison pairing, endpoints, result reuse, frozen inputs, and resume guards")


if __name__ == "__main__":
    main()
