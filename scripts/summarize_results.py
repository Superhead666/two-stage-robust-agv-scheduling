#!/usr/bin/env python3
# Copyright (c) 2026 Renjie Yu and Daniel Oron
# SPDX-License-Identifier: MIT
"""Regenerate table data from recorded results without optimization software."""

import csv
from collections import defaultdict
from pathlib import Path
from statistics import mean, median

ROOT = Path(__file__).resolve().parents[1]
RECORDED = ROOT / "results/recorded"
OUTPUT = ROOT / "results/tables"


def read(name, file="results.csv"):
    with (RECORDED / name / file).open(newline="") as stream:
        return list(csv.DictReader(stream))


def write(name, rows):
    OUTPUT.mkdir(exist_ok=True)
    with (OUTPUT / name).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def key(row):
    return tuple(float(row[name]) for name in ("n", "m", "gamma", "tardiness_factor", "replication"))


def main():
    lbbd = read("lbbd_scalability_formal")
    milp = read("baseline_formal_milp")
    cp = read("baseline_formal_cp")
    pool = defaultdict(float)
    for row in lbbd + milp + cp:
        pool[key(row)] = max(pool[key(row)], float(row["final_lower"]))

    def summarize(rows):
        return {"runs": len(rows), "optimal": sum(int(row["optimal"]) for row in rows),
                "mean_time": mean(float(row["wall_seconds"]) if int(row["optimal"]) else 1800 for row in rows),
                "mean_deviation_percent": mean(100 * (float(row["final_upper"]) - pool[key(row)]) / max(1, float(row["final_upper"])) for row in rows),
                "mean_gap_percent": mean(100 * (float(row["final_upper"]) - float(row["final_lower"])) / max(1, float(row["final_upper"])) for row in rows)}

    grouped = defaultdict(list)
    for row in lbbd + milp + cp:
        if int(row["n"]) <= 30 and row["method"].endswith("hig"):
            grouped[(row["method"], *key(row)[:4])].append(row)
    write("algorithmic_comparison.csv", [{"method": group[0], "n": group[1], "m": group[2], "gamma": group[3], "beta": group[4], **summarize(rows)} for group, rows in sorted(grouped.items())])
    grouped = defaultdict(list)
    for row in milp + cp:
        grouped[(row["method"], int(row["n"]), int(row["m"]))].append(row)
    write("baseline_initialization.csv", [{"method": group[0], "n": group[1], "m": group[2], **summarize(rows)} for group, rows in sorted(grouped.items())])
    grouped = defaultdict(list)
    for row in lbbd:
        if int(row["n"]) >= 40:
            grouped[key(row)[:4]].append(row)
    write("scalability.csv", [{"n": group[0], "m": group[1], "gamma": group[2], "beta": group[3], **summarize(rows),
                             "mean_exact_seconds": mean(float(row["solver_seconds"]) for row in rows),
                             "median_exact_seconds": median(float(row["solver_seconds"]) for row in rows),
                             "mean_verification_seconds": mean(float(row["verification_seconds"]) for row in rows),
                             "mean_master_nodes": mean(float(row["master_nodes"]) for row in rows)} for group, rows in sorted(grouped.items())])
    write("unproved_instances.csv", [row for row in lbbd if not int(row["optimal"])])
    model = read("model_comparison_n30_m3", "paired_comparison.csv")
    grouped = defaultdict(list)
    for row in model:
        grouped[float(row["gamma"])].append(row)
    model_summary = []
    for gamma, rows in sorted(grouped.items()):
        result = {"gamma": gamma, "replications": len(rows), "mean_two_stage": mean(float(row["two_stage_upper"]) for row in rows)}
        for counterpart in ("static", "nominal_assignment"):
            for endpoint in ("lower", "upper"):
                result[f"mean_{counterpart}_{endpoint}"] = mean(float(row[f"{counterpart}_{endpoint}"]) for row in rows)
                result[f"mean_absolute_reduction_{counterpart}_{endpoint}"] = mean(float(row[f"{counterpart}_{endpoint}"]) - float(row["two_stage_upper"]) for row in rows)
                result[f"mean_percent_reduction_{counterpart}_{endpoint}"] = mean(100 * (float(row[f"{counterpart}_{endpoint}"]) - float(row["two_stage_upper"])) / max(1, float(row[f"{counterpart}_{endpoint}"])) for row in rows)
        model_summary.append(result)
    write("model_comparison.csv", model_summary)
    sensitivity = []
    for campaign in ("sensitivity_gamma_n40", "sensitivity_crossed_n40_gamma0.5", "sensitivity_crossed_n40_gamma2", "sensitivity_crossed_n40_gamma4", "sensitivity_technology_n40"):
        grouped = defaultdict(list)
        for row in read(campaign, "sensitivity_results.csv"):
            grouped[tuple(float(row[name]) for name in ("gamma", "deviation_level", "capacity_scale", "recharge_factor"))].append(row)
        for group, rows in sorted(grouped.items()):
            sensitivity.append({"study": campaign, "gamma": group[0], "deviation": group[1], "capacity_scale": group[2], "recharge_factor": group[3], "replications": len(rows), "optimal": sum(int(row["optimal"]) for row in rows), "mean_tardiness": mean(float(row["final_upper"]) for row in rows)})
    write("sensitivity.csv", sensitivity)
    print("Regenerated six table-data files in results/tables.")


if __name__ == "__main__":
    main()
