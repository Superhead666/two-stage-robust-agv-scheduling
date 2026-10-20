#!/usr/bin/env python3
# Copyright (c) 2026 Renjie Yu and Daniel Oron
# SPDX-License-Identifier: MIT
"""Small end-to-end equivalence check for the four Part 1 configurations."""

import argparse
import csv
import os
import pathlib
import subprocess
import tempfile


def run(command):
    environment = os.environ.copy()
    environment["LC_ALL"] = "C"
    subprocess.run([str(x) for x in command], check=True, env=environment)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--build", default="build")
    parser.add_argument("--cp", action="store_true")
    args = parser.parse_args()
    root = pathlib.Path(__file__).resolve().parents[1]
    build = (root / args.build).resolve()
    with tempfile.TemporaryDirectory(prefix="agv-ejor-smoke-") as temporary:
        work = pathlib.Path(temporary)
        instance = work / "tiny.agv"
        run([build / "agv_hig", "--generate", instance, "--n", 4, "--m", 2,
             "--seed", 318901, "--gamma", 1, "--variation", 0.4,
             "--capacity", 2, "--due", 0.4])
        common = ["--instance", instance, "--total-limit", 10, "--hig-limit", 1,
                  "--seed", 1, "--hig-binary", build / "agv_hig",
                  "--lbbd-binary", build / "agv_lbbd",
                  "--initialize-binary", build / "agv_initialize"]
        values = {}
        methods = ["lbbd-basic", "lbbd-hig", "ccg-basic", "ccg-hig"]
        if args.cp:
            methods += ["cp-basic", "cp-hig"]
        for method in methods:
            output = work / method
            ccg = build / ("agv_ccg_cp" if method.startswith("cp-") else "agv_ccg")
            run([build / "run_case", *common, "--ccg-binary", ccg,
                 "--method", method, "--output", output])
            with (output / "result.csv").open(newline="") as source:
                row = next(csv.DictReader(source))
            lower, upper = int(row["final_lower"]), int(row["final_upper"])
            if lower != upper:
                raise RuntimeError(f"{method} did not solve the tiny instance: [{lower},{upper}]")
            values[method] = upper
        if len(set(values.values())) != 1:
            raise RuntimeError(f"methods disagree: {values}")
        print(f"PASS {len(methods)} end-to-end configurations at optimum {next(iter(values.values()))}")


if __name__ == "__main__":
    main()
