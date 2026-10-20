#!/usr/bin/env python3
# Copyright (c) 2026 Renjie Yu and Daniel Oron
# SPDX-License-Identifier: MIT
"""Check archived input hashes, paired data, result counts, and bounds."""

import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RECORDED = ROOT / "results/recorded"


def read(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check(test, message):
    if not test:
        raise ValueError(message)


def main():
    provenance = json.loads((RECORDED / "provenance.json").read_text())
    for campaign, record in provenance.items():
        for name, checksum in record["recorded_csv_sha256"].items():
            check(sha(RECORDED / campaign / name) == checksum, f"recorded CSV changed: {campaign}/{name}")
    benchmark = read(ROOT / "data/benchmark/specifications.csv")
    check(len(benchmark) == 320, "benchmark must contain 320 instances")
    for row in benchmark:
        path = ROOT / "data/benchmark/instances" / (row["case"] + ".agv")
        header = path.read_text().splitlines()[0].split()
        check((int(header[0]), int(header[1]), float(header[3])) == (int(row["n"]), int(row["m"]), float(row["gamma"])), "benchmark metadata mismatch")
    for folder in sorted((ROOT / "data").glob("sensitivity_*")):
        specs = read(folder / "sensitivity_specifications.csv")
        for row in specs:
            check(sha(folder / "prepared_instances" / row["source_file"]) == row["instance_sha256"], "sensitivity instance fingerprint mismatch")
        for replication in {row["replication"] for row in specs}:
            check(len({row["fixed_data_sha256"] for row in specs if row["replication"] == replication}) == 1, "sensitivity pairing failed")
    expected = {"lbbd_scalability_formal": (320, {"lbbd-hig": 312}),
                "baseline_formal_milp": (320, {"ccg-basic": 56, "ccg-hig": 54}),
                "baseline_formal_cp": (320, {"cp-basic": 53, "cp-hig": 71})}
    totals = {}
    for campaign, (count, optima) in expected.items():
        rows = read(RECORDED / campaign / "results.csv")
        check(len(rows) == count, "algorithmic result count mismatch")
        for method, optimum_count in optima.items():
            group = [row for row in rows if row["method"] == method]
            check(sum(int(row["optimal"]) for row in group) == optimum_count, "optimality count mismatch")
            check(all(int(row["final_lower"]) <= int(row["final_upper"]) for row in group), "invalid recorded bounds")
            totals[method] = {"runs": len(group), "optimal": optimum_count}
    paired = read(RECORDED / "model_comparison_n30_m3/paired_comparison.csv")
    check(len(paired) == 90, "model comparison must contain 90 paired cases")
    check(sum(int(row["static_optimal"]) for row in paired) == 86, "static-model optimality count mismatch")
    for row in paired:
        check(float(row["two_stage_upper"]) <= float(row["static_upper"]), "static restriction violated")
        check(float(row["two_stage_upper"]) <= float(row["nominal_assignment_upper"]), "nominal restriction violated")
    print(json.dumps({"benchmark_instances": len(benchmark), "algorithmic_results": totals, "model_comparison_cases": len(paired), "status": "passed"}, indent=2))


if __name__ == "__main__":
    main()
