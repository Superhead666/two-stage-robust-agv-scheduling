#!/usr/bin/env python3
# Copyright (c) 2026 Renjie Yu and Daniel Oron
# SPDX-License-Identifier: MIT
"""End-to-end checks for the flexible Part 1 experiment driver."""

import argparse
import csv
import json
import os
import pathlib
import subprocess
import tempfile


def run(command, *, expect_success=True):
    environment = os.environ.copy()
    environment["LC_ALL"] = "C"
    completed = subprocess.run(
        [str(value) for value in command],
        env=environment,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    if expect_success and completed.returncode:
        raise RuntimeError(completed.stdout)
    if not expect_success and not completed.returncode:
        raise RuntimeError("command unexpectedly succeeded")
    return completed.stdout


def read_rows(path):
    with path.open(newline="") as source:
        return list(csv.DictReader(source))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--build", default="build")
    args = parser.parse_args()
    root = pathlib.Path(__file__).resolve().parents[1]
    driver = root / "scripts" / "run_experiment.py"
    build = (root / args.build).resolve()
    with tempfile.TemporaryDirectory(prefix="agv-ejor-driver-") as temporary:
        work = pathlib.Path(temporary)
        dry = work / "dry"
        output = run([
            "python3", driver, "--output", dry,
            "--target-jobs", "4,5", "--agvs", 2, "--size-mode", "paired",
            "--gamma", 0, 1, "--deviation", "0.2,0.4",
            "--processing-max", 5, "--nominal-energy-max", 5,
            "--capacity-scale", 2, "--due-date-factor", 0.4,
            "--recharge-factor", 1, "--repetitions", 2,
            "--configurations", "lbbd-basic", "ccg-basic", "--dry-run",
        ])
        if "16 instances" not in output or "32 solver runs" not in output or dry.exists():
            raise RuntimeError(f"unexpected dry-run behavior: {output}")

        sentinel_dry = work / "sentinel-dry"
        output = run([
            "python3", driver, "--output", sentinel_dry,
            "--size-pairs", "40:3", "40:5", "--gamma", 1, 3,
            "--tardiness-factor", "0.30,0.50", "--due-date-range", 0.50,
            "--replications", 1, 2, "--configurations",
            "milp-ccg-basic", "milp-ccg-hig", "cp-basic", "cp-hig",
            "--total-limit", 1800, "--hig-limit", 60, "--dry-run",
        ])
        if "16 instances" not in output or "64 solver runs" not in output:
            raise RuntimeError(f"unexpected sentinel design: {output}")

        experiment = work / "panel"
        common = [
            "python3", driver, "--output", experiment, "--build", build,
            "--target-jobs", 4, "--agvs", 2, "--size-mode", "paired",
            "--gamma", 0, 1, "--deviation", 0.3,
            "--processing-max", 5, "--nominal-energy-max", 5,
            "--capacity-scale", 2, "--tardiness-factor", 0.4,
            "--due-date-range", 0.5,
            "--recharge-time", 3, "--repetitions", 1,
            "--generator-seed", 318903,
            "--configurations", "lbbd-basic", "lbbd-hig", "ccg-hig",
            "--total-limit", 10, "--hig-limit", 1, "--seed", 1,
        ]
        run([*common, "--prepare-only", "--resume"])
        specifications = read_rows(experiment / "specifications.csv")
        if len(specifications) != 2:
            raise RuntimeError("the prepared panel has the wrong number of cases")
        if {row["instance_seed"] for row in specifications} != {"318903"}:
            raise RuntimeError("common random numbers were not preserved across gamma levels")
        if {row["actual_recharge_time"] for row in specifications} != {"3"}:
            raise RuntimeError("the requested recharge time was not preserved")
        if set(json.loads((experiment / "binaries.json").read_text())) != {
            "agv_ccg", "agv_hig", "agv_initialize", "agv_lbbd", "run_case"
        }:
            raise RuntimeError("the frozen executable set is incorrect")

        run([*common, "--resume"])
        results = read_rows(experiment / "results.csv")
        if len(results) != 6:
            raise RuntimeError("the aggregate does not contain all case-configuration pairs")
        for case in {row["instance"] for row in results}:
            values = {
                int(row["final_upper"])
                for row in results
                if row["instance"] == case
            }
            if len(values) != 1:
                raise RuntimeError(f"methods disagree for {case}: {values}")
        if not (experiment / "method_summary.csv").is_file():
            raise RuntimeError("the method summary was not written")
        warm_rows = read_rows(experiment / "warm_start_cases.csv")
        if len(warm_rows) != 4 or any(
            not row["exact_progress_class"] for row in warm_rows
        ):
            raise RuntimeError("HIG progress attribution is incomplete")
        run([*common, "--resume"])
        if len(read_rows(experiment / "results.csv")) != 6:
            raise RuntimeError("a completed resume changed the aggregate")

        imported_source = work / "imported-source"
        imported_source.mkdir()
        source_instance = imported_source / "one.agv"
        source_instance.write_bytes((experiment / "instances" / "case_000001.agv").read_bytes())
        imported_run = work / "imported-run"
        imported = [
            "python3", driver, "--instances", imported_source,
            "--output", imported_run, "--build", build,
            "--configurations", "ccg-basic",
            "--total-limit", 10, "--hig-limit", 1, "--seed", 1,
            "--prepare-only",
        ]
        run(imported)
        if set(json.loads((imported_run / "binaries.json").read_text())) != {
            "agv_ccg", "run_case"
        }:
            raise RuntimeError("the imported mode froze unnecessary executables")
        if read_rows(imported_run / "specifications.csv")[0]["source_file"] != "one.agv":
            raise RuntimeError("the imported source identity was not recorded")
        run([*imported[:-1], "--prepare-only", "--resume"])

        changed = list(common)
        gamma_index = changed.index("--gamma") + 1
        changed[gamma_index:gamma_index + 2] = [0, 2]
        output = run([*changed, "--resume"], expect_success=False)
        if "do not match the saved design" not in output:
            raise RuntimeError(f"changed design was rejected for the wrong reason: {output}")
        print("PASS flexible panel generation, freezing, execution, and exact resume validation")


if __name__ == "__main__":
    main()
