#!/usr/bin/env python3
# Copyright (c) 2026 Renjie Yu and Daniel Oron
# SPDX-License-Identifier: MIT
"""End-to-end checks for the sensitivity experiment driver."""

import argparse
import csv
import hashlib
import json
import os
import pathlib
import subprocess
import tempfile


def run(command, *, expect_success=True):
    environment = os.environ.copy()
    environment["LC_ALL"] = "C"
    completed = subprocess.run(
        [str(value) for value in command], env=environment,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, check=False,
    )
    if expect_success and completed.returncode:
        raise RuntimeError(completed.stdout)
    if not expect_success and not completed.returncode:
        raise RuntimeError("command unexpectedly succeeded")
    return completed.stdout


def rows(path):
    with path.open(newline="") as source:
        return list(csv.DictReader(source))


def job_hash(path):
    lines = path.read_text().splitlines()
    return hashlib.sha256(("\n".join(lines[2:])+"\n").encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--build", default="build")
    args = parser.parse_args()
    root = pathlib.Path(__file__).resolve().parents[1]
    driver = root / "scripts" / "run_sensitivity.py"
    build = (root / args.build).resolve()

    with tempfile.TemporaryDirectory(prefix="agv-sensitivity-test-") as temporary:
        work = pathlib.Path(temporary)
        output = run([
            "python3", driver, "gamma", "--output", work / "dry",
            "--gamma-values", 0, .5, 1, 1.5, 2, 2.5, 3, 3.5, 4,
            "--dry-run",
        ])
        if "9 settings" not in output or "90 lbbd-hig runs" not in output:
            raise RuntimeError(f"unexpected gamma dry run: {output}")
        output = run([
            "python3", driver, "deviation", "--output", work / "deviation-dry",
            "--deviation-values", .1, .2, .3, .4, .5, "--base-gamma", 2,
            "--dry-run",
        ])
        if "5 settings" not in output or "50 lbbd-hig runs" not in output:
            raise RuntimeError(f"unexpected deviation dry run: {output}")

        technology = work / "technology"
        common = [
            "python3", driver, "technology", "--output", technology,
            "--build", build, "--jobs", 6, "--agvs", 2,
            "--replications", 1, 2, "--generator-seed", 419000,
            "--capacity-scales", 3, 4, "--recharge-factors", 1, 2,
            "--base-capacity-scale", 4, "--base-recharge-factor", 2,
            "--base-gamma", 1, "--processing-max", 5,
            "--nominal-energy-max", 5, "--total-limit", 10,
            "--hig-limit", 1, "--metrics-limit", 5, "--seed", 1,
        ]
        run([*common, "--prepare-only", "--resume"])
        specifications = rows(technology / "sensitivity_specifications.csv")
        if len(specifications) != 8:
            raise RuntimeError("technology panel has the wrong number of cases")
        for replication in (1, 2):
            group = [row for row in specifications if int(row["replication"]) == replication]
            hashes = {
                job_hash(technology / "prepared_instances" / row["source_file"])
                for row in group
            }
            recorded = {row["job_data_sha256"] for row in group}
            if len(hashes) != 1 or len(recorded) != 1:
                raise RuntimeError("job data or due dates changed across technology settings")
            alphas = {
                float(row["recharge_factor"]): int(row["actual_recharge_time"])
                for row in group
            }
            if alphas[2.0] < alphas[1.0]:
                raise RuntimeError("recharge duration does not follow its requested factor")
        protocol = json.loads((technology / "protocol_hashes.json").read_text())
        if len(protocol["prepared_instances"]) != 8:
            raise RuntimeError("the protocol does not cover every technology instance")
        run([*common, "--prepare-only", "--resume"])

        deviation = work / "deviation"
        run([
            "python3", driver, "deviation", "--output", deviation,
            "--build", build, "--jobs", 6, "--agvs", 2,
            "--replications", 1, "--generator-seed", 419900,
            "--deviation", .3, "--deviation-values", .1, .3, .5,
            "--base-gamma", 2, "--base-capacity-scale", 4,
            "--base-recharge-factor", 2, "--processing-max", 5,
            "--nominal-energy-max", 5, "--total-limit", 10,
            "--hig-limit", 1, "--metrics-limit", 5,
            "--prepare-only", "--resume",
        ])
        deviation_rows = rows(deviation / "sensitivity_specifications.csv")
        if len(deviation_rows) != 3:
            raise RuntimeError("deviation panel has the wrong number of cases")
        if len({row["fixed_data_sha256"] for row in deviation_rows}) != 1:
            raise RuntimeError("fixed job data changed across deviation levels")
        if len({row["job_data_sha256"] for row in deviation_rows}) != 3:
            raise RuntimeError("maximum deviations did not change across levels")
        variant_paths = [
            deviation / "prepared_instances" / row["source_file"]
            for row in deviation_rows
        ]
        capacity_rows = {path.read_text().splitlines()[1] for path in variant_paths}
        if len(capacity_rows) != 1:
            raise RuntimeError("battery capacity changed in the deviation study")

        gamma = work / "gamma"
        gamma_command = [
            "python3", driver, "gamma", "--output", gamma, "--build", build,
            "--jobs", 6, "--agvs", 2, "--replications", 1,
            "--generator-seed", 420000, "--gamma-values", 0, 2,
            "--base-gamma", 2, "--base-capacity-scale", 4,
            "--base-recharge-factor", 2, "--processing-max", 5,
            "--nominal-energy-max", 5, "--total-limit", 10,
            "--hig-limit", 1, "--metrics-limit", 5, "--seed", 1,
        ]
        run(gamma_command)
        results = rows(gamma / "sensitivity_results.csv")
        if len(results) != 2 or any(row["metrics_status"] != "complete" for row in results):
            raise RuntimeError("gamma results or operational metrics are incomplete")
        if len(rows(gamma / "sensitivity_summary.csv")) != 2:
            raise RuntimeError("gamma summary is incomplete")
        if len(rows(gamma / "sensitivity_pairs.csv")) != 2:
            raise RuntimeError("paired gamma summary is incomplete")
        run([*gamma_command, "--resume"])

        changed = list(gamma_command)
        changed[changed.index("--base-capacity-scale")+1] = 5
        output = run([*changed, "--resume"], expect_success=False)
        if "do not match the saved sensitivity design" not in output:
            raise RuntimeError(f"changed design was rejected for the wrong reason: {output}")

        print("PASS sensitivity panel freezing, common data, execution, metrics, and resume")


if __name__ == "__main__":
    main()
