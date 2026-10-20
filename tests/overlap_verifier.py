#!/usr/bin/env python3
# Copyright (c) 2026 Renjie Yu and Daniel Oron
# SPDX-License-Identifier: MIT
"""Check independent panel generation and byte-level overlap verification."""

import argparse
import csv
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


def rows(path):
    with path.open(newline="") as source:
        return list(csv.DictReader(source))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--build", default="build")
    args = parser.parse_args()
    root = pathlib.Path(__file__).resolve().parents[1]
    driver = root / "scripts" / "run_experiment.py"
    verifier = root / "scripts" / "verify_instance_overlap.py"
    build = (root / args.build).resolve()
    with tempfile.TemporaryDirectory(prefix="agv-overlap-") as temporary:
        work = pathlib.Path(temporary)
        baseline = work / "baseline"
        scalability = work / "scalability"
        common = [
            "--build", build,
            "--gamma", 0,
            "--tardiness-factor", 0.30,
            "--due-date-range", 0.50,
            "--deviation", 0.30,
            "--processing-max", 5,
            "--nominal-energy-max", 5,
            "--capacity-scale", 2,
            "--recharge-factor", 1,
            "--repetitions", 1,
            "--generator-seed", 431000,
            "--pair-seed-stride", 1000,
            "--configurations", "lbbd-hig",
            "--total-limit", 10,
            "--hig-limit", 1,
            "--seed", 1,
            "--prepare-only", "--resume",
        ]
        run([
            "python3", driver, "--output", baseline,
            "--size-pairs", "4:2", *common,
        ])
        run([
            "python3", driver, "--output", scalability,
            "--size-pairs", "4:2", "5:2", *common,
        ])
        output = run([
            "python3", verifier,
            "--baseline", baseline,
            "--scalability", scalability,
            "--jobs", 4,
        ])
        manifest = scalability / "overlap_manifest.csv"
        if "Verified 1 byte-identical overlap instances" not in output:
            raise RuntimeError(f"unexpected overlap report: {output}")
        if len(rows(manifest)) != 1:
            raise RuntimeError("overlap manifest has the wrong number of rows")

        overlap = scalability / "instances" / "case_000001.agv"
        overlap.write_text(overlap.read_text() + "\n")
        output = run([
            "python3", verifier,
            "--baseline", baseline,
            "--scalability", scalability,
            "--jobs", 4,
        ], expect_success=False)
        if "instance hash mismatch" not in output:
            raise RuntimeError(f"modified instance was rejected incorrectly: {output}")
        print("PASS independent panel generation and overlap hash verification")


if __name__ == "__main__":
    main()
