#!/usr/bin/env python3
# Copyright (c) 2026 Renjie Yu and Daniel Oron
# SPDX-License-Identifier: MIT
"""Prepare, run, and summarize the managerial sensitivity experiments."""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
import os
import pathlib
import shutil
import statistics
import subprocess
import sys
from typing import Any, Iterable


ROOT = pathlib.Path(__file__).resolve().parents[1]
SCHEMA = "agv-sensitivity-v1"
DEFAULT_GAMMA = [value / 2 for value in range(9)]
DEFAULT_DEVIATION = [0.10, 0.20, 0.30, 0.40, 0.50]
DEFAULT_CAPACITY = [3.0, 3.5, 4.0, 4.5, 5.0]
DEFAULT_RECHARGE = [1.0, 1.5, 2.0, 2.5, 3.0]
METRIC_FIELDS = [
    "metrics_status",
    "nominal_maximum_tardiness",
    "nominal_maximum_completion_time",
    "nominal_mean_completion_time",
    "nominal_total_recharges",
    "nominal_maximum_recharges_per_agv",
    "metrics_seconds",
]


def fail(message: str) -> None:
    raise ValueError(message)


def finite_real(text: str) -> float:
    value = float(text)
    if not math.isfinite(value):
        fail(f"real values must be finite, found {text}")
    return value


def split_values(values: Iterable[str] | None, convert, option: str):
    if values is None:
        return None
    pieces: list[str] = []
    for value in values:
        fields = value.split(",")
        if any(field == "" for field in fields):
            fail(f"{option} contains an empty value")
        pieces.extend(fields)
    try:
        result = [convert(value) for value in pieces]
    except (TypeError, ValueError) as error:
        fail(f"invalid value for {option}: {error}")
    if not result:
        fail(f"{option} requires at least one value")
    if len(result) != len(set(result)):
        fail(f"duplicate value in {option}")
    return result


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description=(
            "Run the uncertainty-budget or battery/recharge sensitivity study. "
            "All variants of a replication retain the same job data and due dates."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    result.add_argument("study", choices=("gamma", "deviation", "technology"))
    result.add_argument("--output", required=True, type=pathlib.Path)
    result.add_argument("--build", type=pathlib.Path, default=ROOT / "build")
    result.add_argument("--jobs", type=int, default=40)
    result.add_argument("--agvs", type=int, default=3)
    result.add_argument("--replications", nargs="+",
                        help="replication indices; default is 1,...,10")
    result.add_argument("--generator-seed", type=int, default=2_004_000)
    result.add_argument("--deviation", type=finite_real, default=0.30)
    result.add_argument("--processing-max", type=int, default=20)
    result.add_argument("--nominal-energy-max", type=int, default=20)
    result.add_argument("--tardiness-factor", type=finite_real, default=0.50)
    result.add_argument("--due-date-range", type=finite_real, default=0.50)
    result.add_argument("--base-gamma", type=finite_real, default=2.0,
                        help="reference budget and fixed budget in the technology study")
    result.add_argument("--base-capacity-scale", type=finite_real, default=4.0)
    result.add_argument("--base-recharge-factor", type=finite_real, default=2.0)
    result.add_argument("--gamma-values", nargs="+",
                        help="budget values in the gamma study")
    result.add_argument("--deviation-values", nargs="+",
                        help="maximum deviation rates in the deviation study")
    result.add_argument("--capacity-scales", nargs="+",
                        help="capacity levels in the technology study")
    result.add_argument("--recharge-factors", nargs="+",
                        help="recharge-duration levels in the technology study")
    result.add_argument("--technology-mode", choices=("crossed", "one-at-a-time"),
                        default="crossed")
    result.add_argument("--total-limit", type=finite_real, default=1800.0)
    result.add_argument("--hig-limit", type=finite_real, default=60.0)
    result.add_argument("--metrics-limit", type=finite_real, default=30.0,
                        help="postprocessing limit per assignment; set to 0 to omit schedule metrics")
    result.add_argument("--seed", type=int, default=1, help="solver seed")
    result.add_argument("--dry-run", action="store_true")
    result.add_argument("--prepare-only", action="store_true")
    result.add_argument("--resume", action="store_true")
    return result


def normalize(args: argparse.Namespace) -> dict[str, Any]:
    replications = split_values(args.replications, int, "--replications")
    if replications is None:
        replications = list(range(1, 11))
    gamma = split_values(args.gamma_values, finite_real, "--gamma-values")
    deviation_values = split_values(
        args.deviation_values, finite_real, "--deviation-values"
    )
    capacity = split_values(args.capacity_scales, finite_real, "--capacity-scales")
    recharge = split_values(args.recharge_factors, finite_real, "--recharge-factors")
    gamma = DEFAULT_GAMMA if gamma is None else gamma
    deviation_values = DEFAULT_DEVIATION if deviation_values is None else deviation_values
    capacity = DEFAULT_CAPACITY if capacity is None else capacity
    recharge = DEFAULT_RECHARGE if recharge is None else recharge

    if not 1 <= args.jobs <= 62 or not 1 <= args.agvs <= 32:
        fail("--jobs must lie in [1,62] and --agvs in [1,32]")
    if any(value < 1 for value in replications):
        fail("replication indices must be positive")
    if args.generator_seed < 0 or args.generator_seed + max(replications)-1 > 2_000_000_000:
        fail("generator seeds must lie in [0,2000000000]")
    if args.deviation < 0 or args.processing_max < 1 or args.nominal_energy_max < 1:
        fail("the reference deviation must be nonnegative and data maxima must be positive")
    if args.study == "deviation" and args.deviation <= 0:
        fail("the reference deviation must be positive in the deviation study")
    if not 0 <= args.tardiness_factor <= 1 or args.due_date_range < 0:
        fail("invalid due-date parameters")
    if not 0 <= args.base_gamma <= args.jobs:
        fail("--base-gamma must lie in [0,n]")
    if args.base_capacity_scale <= 0 or args.base_recharge_factor < 0:
        fail("base capacity must be positive and base recharge factor nonnegative")
    if any(value < 0 or value > args.jobs for value in gamma):
        fail("every gamma value must lie in [0,n]")
    if any(value < 0 for value in deviation_values):
        fail("deviation levels must be nonnegative")
    if any(value <= 0 for value in capacity) or any(value < 0 for value in recharge):
        fail("capacity scales must be positive and recharge factors nonnegative")
    if args.total_limit <= 0 or args.hig_limit <= 0 or args.hig_limit >= args.total_limit:
        fail("HIG time must be positive and smaller than the total time limit")
    if args.metrics_limit < 0 or args.seed < 0 or args.seed > 2_000_000_000:
        fail("invalid metrics limit or solver seed")
    if args.dry_run and (args.prepare_only or args.resume):
        fail("--dry-run cannot be combined with --prepare-only or --resume")
    if args.study == "gamma" and args.base_gamma not in gamma:
        fail("--base-gamma must occur in --gamma-values for paired summaries")
    if args.study == "deviation" and args.deviation not in deviation_values:
        fail("--deviation must occur in --deviation-values for paired summaries")
    if args.study == "technology":
        if args.base_capacity_scale not in capacity:
            fail("--base-capacity-scale must occur in --capacity-scales")
        if args.base_recharge_factor not in recharge:
            fail("--base-recharge-factor must occur in --recharge-factors")

    if args.study == "gamma":
        settings = [
            {"setting_index": index+1, "gamma": value,
             "deviation_level": args.deviation,
             "capacity_scale": args.base_capacity_scale,
             "recharge_factor": args.base_recharge_factor}
            for index, value in enumerate(gamma)
        ]
    elif args.study == "deviation":
        settings = [
            {"setting_index": index+1, "gamma": args.base_gamma,
             "deviation_level": value,
             "capacity_scale": args.base_capacity_scale,
             "recharge_factor": args.base_recharge_factor}
            for index, value in enumerate(deviation_values)
        ]
    elif args.technology_mode == "crossed":
        settings = [
            {"setting_index": index+1, "gamma": args.base_gamma,
             "deviation_level": args.deviation,
             "capacity_scale": cap, "recharge_factor": charge}
            for index, (cap, charge) in enumerate(itertools.product(capacity, recharge))
        ]
    else:
        pairs = [(cap, args.base_recharge_factor) for cap in capacity]
        pairs.extend(
            (args.base_capacity_scale, charge)
            for charge in recharge if charge != args.base_recharge_factor
        )
        settings = [
            {"setting_index": index+1, "gamma": args.base_gamma,
             "deviation_level": args.deviation,
             "capacity_scale": cap, "recharge_factor": charge}
            for index, (cap, charge) in enumerate(pairs)
        ]

    return {
        "schema": SCHEMA,
        "study": args.study,
        "jobs": args.jobs,
        "agvs": args.agvs,
        "replications": replications,
        "generator_seed": args.generator_seed,
        "deviation": args.deviation,
        "processing_max": args.processing_max,
        "nominal_energy_max": args.nominal_energy_max,
        "tardiness_factor": args.tardiness_factor,
        "due_date_range": args.due_date_range,
        "base_gamma": args.base_gamma,
        "base_capacity_scale": args.base_capacity_scale,
        "base_recharge_factor": args.base_recharge_factor,
        "gamma_values": gamma,
        "deviation_values": deviation_values,
        "capacity_scales": capacity,
        "recharge_factors": recharge,
        "technology_mode": args.technology_mode,
        "settings": settings,
        "configuration": "lbbd-hig",
        "total_limit": args.total_limit,
        "hig_limit": args.hig_limit,
        "metrics_limit": args.metrics_limit,
        "solver_seed": args.seed,
    }


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_text(path: pathlib.Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(value)
    os.replace(temporary, path)


def atomic_json(path: pathlib.Path, value: Any) -> None:
    atomic_text(path, json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def atomic_csv(path: pathlib.Path, fields: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", newline="") as destination:
        writer = csv.DictWriter(destination, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def read_csv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open(newline="") as source:
        return list(csv.DictReader(source))


def run(command: list[Any], *, capture: bool = False) -> str:
    environment = os.environ.copy()
    environment["LC_ALL"] = "C"
    completed = subprocess.run(
        [str(value) for value in command], env=environment,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.STDOUT if capture else None,
        text=True, check=False,
    )
    if completed.returncode:
        details = f"\n{completed.stdout}" if capture and completed.stdout else ""
        fail(f"command exited with code {completed.returncode}: {' '.join(map(str, command))}{details}")
    return completed.stdout or ""


def read_instance(path: pathlib.Path) -> dict[str, Any]:
    lines = path.read_text().splitlines()
    header = lines[0].split()
    if len(header) != 4:
        fail(f"invalid instance header in {path}")
    n, m, alpha = map(int, header[:3])
    gamma = finite_real(header[3])
    capacities = [finite_real(value) for value in lines[1].split()]
    if len(capacities) != m or len(lines) != n+2:
        fail(f"invalid instance dimensions in {path}")
    jobs = []
    for line in lines[2:]:
        fields = line.split()
        if len(fields) != 1+3*m:
            fail(f"invalid job row in {path}")
        due = int(fields[0])
        triples = []
        for k in range(m):
            triples.append((int(fields[1+3*k]), finite_real(fields[2+3*k]),
                            finite_real(fields[3+3*k])))
        jobs.append((due, triples))
    return {"n": n, "m": m, "alpha": alpha, "gamma": gamma,
            "capacities": capacities, "jobs": jobs}


def number(value: float) -> str:
    return format(value, ".17g")


def write_instance(path: pathlib.Path, data: dict[str, Any]) -> None:
    lines = [
        f"{data['n']} {data['m']} {data['alpha']} {number(data['gamma'])}",
        " ".join(number(value) for value in data["capacities"]),
    ]
    for due, triples in data["jobs"]:
        row = [str(due)]
        for processing, nominal, deviation in triples:
            row.extend((str(processing), number(nominal), number(deviation)))
        lines.append(" ".join(row))
    atomic_text(path, "\n".join(lines) + "\n")


def job_hash(data: dict[str, Any]) -> str:
    digest = hashlib.sha256()
    for due, triples in data["jobs"]:
        digest.update(str(due).encode())
        for processing, nominal, deviation in triples:
            digest.update(f"|{processing}|{number(nominal)}|{number(deviation)}".encode())
        digest.update(b"\n")
    return digest.hexdigest()


def fixed_data_hash(data: dict[str, Any]) -> str:
    """Hash data held fixed in every sensitivity study, excluding deviations."""
    digest = hashlib.sha256()
    for due, triples in data["jobs"]:
        digest.update(str(due).encode())
        for processing, nominal, _ in triples:
            digest.update(f"|{processing}|{number(nominal)}".encode())
        digest.update(b"\n")
    return digest.hexdigest()


def rescale_deviations(data: dict[str, Any], reference: float,
                       requested: float) -> list[tuple[int, list[tuple[int, float, float]]]]:
    ratio = requested/reference
    jobs = []
    for due, triples in data["jobs"]:
        jobs.append((due, [
            (processing, nominal, deviation*ratio)
            for processing, nominal, deviation in triples
        ]))
    return jobs


def technology_values(data: dict[str, Any], capacity_scale: float,
                      recharge_factor: float) -> tuple[int, list[float], int]:
    processing = [
        triple[0] for _, triples in data["jobs"] for triple in triples
    ]
    alpha = int(math.floor(recharge_factor*statistics.mean(processing)+0.5))
    capacities = []
    clipped = 0
    for k in range(data["m"]):
        upper = [triples[k][1]+triples[k][2] for _, triples in data["jobs"]]
        scaled = capacity_scale*statistics.mean(upper)
        maximum = max(upper)
        if maximum > scaled:
            clipped += 1
        capacities.append(max(maximum, scaled))
    return alpha, capacities, clipped


def freeze(path: pathlib.Path, destination: pathlib.Path) -> pathlib.Path:
    if not path.is_file():
        fail(f"required executable is missing: {path}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(path, destination)
    destination.chmod(destination.stat().st_mode | 0o111)
    return destination


def prepare(output: pathlib.Path, build: pathlib.Path, design: dict[str, Any]) -> None:
    binary = output / "bin" / "agv_hig"
    metrics = output / "bin" / "agv_schedule_metrics"
    reference = output / "reference_instances"
    prepared = output / "prepared_instances"
    reference.mkdir(parents=True)
    prepared.mkdir(parents=True)
    freeze(build / "agv_hig", binary)
    freeze(build / "agv_schedule_metrics", metrics)

    specifications: list[dict[str, Any]] = []
    reference_hashes: dict[int, str] = {}
    fixed_hashes: dict[int, str] = {}
    for replication in design["replications"]:
        seed = design["generator_seed"]+replication-1
        path = reference / f"rep_{replication:02d}.agv"
        run([
            binary, "--generate", path, "--n", design["jobs"], "--m", design["agvs"],
            "--seed", seed, "--gamma", design["base_gamma"],
            "--variation", design["deviation"], "--capacity", design["base_capacity_scale"],
            "--tardiness-factor", design["tardiness_factor"],
            "--due-date-range", design["due_date_range"],
            "--processing-max", design["processing_max"],
            "--nominal-energy-max", design["nominal_energy_max"],
            "--recharge-factor", design["base_recharge_factor"],
        ])
        base = read_instance(path)
        reference_hashes[replication] = job_hash(base)
        fixed_hashes[replication] = fixed_data_hash(base)
        for setting in design["settings"]:
            data = {
                "n": base["n"], "m": base["m"], "jobs": base["jobs"],
                "gamma": setting["gamma"], "alpha": base["alpha"],
                "capacities": list(base["capacities"]),
            }
            clipped = 0
            if design["study"] == "deviation":
                data["jobs"] = rescale_deviations(
                    base, design["deviation"], setting["deviation_level"]
                )
                for _, triples in data["jobs"]:
                    for k, (_, nominal, maximum_deviation) in enumerate(triples):
                        if nominal+maximum_deviation > data["capacities"][k]+1e-9:
                            fail(
                                "a requested deviation makes one job exceed the fixed battery "
                                "capacity; choose a narrower deviation range or a larger base capacity"
                            )
            elif design["study"] == "technology":
                data["alpha"], data["capacities"], clipped = technology_values(
                    base, setting["capacity_scale"], setting["recharge_factor"]
                )
            source_file = f"s{setting['setting_index']:02d}_r{replication:02d}.agv"
            variant = prepared / source_file
            write_instance(variant, data)
            specifications.append({
                "source_file": source_file,
                "study": design["study"],
                "setting_index": setting["setting_index"],
                "replication": replication,
                "instance_seed": seed,
                "n": base["n"], "m": base["m"],
                "gamma": setting["gamma"],
                "deviation": setting["deviation_level"],
                "deviation_level": setting["deviation_level"],
                "capacity_scale": setting["capacity_scale"],
                "recharge_factor": setting["recharge_factor"],
                "actual_recharge_time": data["alpha"],
                "minimum_capacity": min(data["capacities"]),
                "mean_capacity": statistics.mean(data["capacities"]),
                "maximum_capacity": max(data["capacities"]),
                "capacity_clipped_machines": clipped,
                "fixed_data_sha256": fixed_hashes[replication],
                "reference_job_data_sha256": reference_hashes[replication],
                "job_data_sha256": job_hash(data),
                "instance_sha256": sha256(variant),
            })

    fields = list(specifications[0])
    specification_path = output / "sensitivity_specifications.csv"
    atomic_csv(specification_path, fields, specifications)
    protocol = {
        "schema": SCHEMA,
        "generator_binary": sha256(binary),
        "metrics_binary": sha256(metrics),
        "sensitivity_specifications": sha256(specification_path),
        "reference_instances": {path.name: sha256(path) for path in sorted(reference.glob("*.agv"))},
        "prepared_instances": {path.name: sha256(path) for path in sorted(prepared.glob("*.agv"))},
    }
    atomic_json(output / "protocol_hashes.json", protocol)


def validate_prepared(output: pathlib.Path) -> None:
    protocol = json.loads((output / "protocol_hashes.json").read_text())
    checks = {
        "generator_binary": output / "bin" / "agv_hig",
        "metrics_binary": output / "bin" / "agv_schedule_metrics",
        "sensitivity_specifications": output / "sensitivity_specifications.csv",
    }
    for name, path in checks.items():
        if not path.is_file() or sha256(path) != protocol[name]:
            fail(f"frozen file changed: {path}")
    for group, directory in (("reference_instances", "reference_instances"),
                             ("prepared_instances", "prepared_instances")):
        actual = {
            path.name: sha256(path)
            for path in sorted((output / directory).glob("*.agv"))
        }
        if actual != protocol[group]:
            fail(f"frozen {group.replace('_', ' ')} changed")


def run_solver(output: pathlib.Path, build: pathlib.Path, design: dict[str, Any],
               prepare_only: bool) -> None:
    command: list[Any] = [
        sys.executable, ROOT / "scripts" / "run_experiment.py",
        "--instances", output / "prepared_instances",
        "--output", output / "solver", "--build", build,
        "--configurations", "lbbd-hig",
        "--total-limit", design["total_limit"],
        "--hig-limit", design["hig_limit"], "--seed", design["solver_seed"],
        "--resume",
    ]
    if prepare_only:
        command.append("--prepare-only")
    run(command)


def add_metrics(output: pathlib.Path, design: dict[str, Any]) -> list[dict[str, Any]]:
    solver = output / "solver"
    results = read_csv(solver / "results.csv")
    source_specs = {
        row["source_file"]: row for row in read_csv(output / "sensitivity_specifications.csv")
    }
    imported_specs = {
        row["case"]: row for row in read_csv(solver / "specifications.csv")
    }
    binary = output / "bin" / "agv_schedule_metrics"
    merged: list[dict[str, Any]] = []
    for result in results:
        case = result["instance"]
        source_file = imported_specs[case]["source_file"]
        specification = source_specs[source_file]
        row: dict[str, Any] = dict(result)
        row.update(specification)
        metrics_values = {field: "" for field in METRIC_FIELDS}
        if design["metrics_limit"] > 0:
            case_directory = solver / "cases" / case / "lbbd-hig"
            metrics_path = case_directory / "nominal_schedule_metrics.csv"
            if not metrics_path.is_file():
                run([
                    binary, solver / "instances" / f"{case}.agv",
                    case_directory / "lbbd.sol", metrics_path, design["metrics_limit"],
                ])
            metrics_row = read_csv(metrics_path)[0]
            metrics_values = {
                "metrics_status": metrics_row["status"],
                **{field: metrics_row.get(field, "") for field in METRIC_FIELDS[1:]},
            }
        row.update(metrics_values)
        merged.append(row)

    fields = list(results[0])
    for field in list(source_specs[next(iter(source_specs))]) + METRIC_FIELDS:
        if field not in fields:
            fields.append(field)
    atomic_csv(output / "sensitivity_results.csv", fields, merged)
    return merged


def values(rows: list[dict[str, Any]], field: str) -> list[float]:
    return [float(row[field]) for row in rows if row.get(field, "") not in ("", None)]


def mean(rows: list[dict[str, Any]], field: str) -> float | str:
    data = values(rows, field)
    return statistics.mean(data) if data else ""


def deviation(rows: list[dict[str, Any]], field: str) -> float | str:
    data = values(rows, field)
    return statistics.stdev(data) if len(data) > 1 else ""


def summarize(output: pathlib.Path, design: dict[str, Any],
              rows: list[dict[str, Any]]) -> None:
    if design["study"] == "gamma":
        group_fields = ["gamma"]
    elif design["study"] == "deviation":
        group_fields = ["deviation_level"]
    else:
        group_fields = ["capacity_scale", "recharge_factor"]
    grouped: dict[tuple[str, ...], list[dict[str, Any]]] = {}
    for row in rows:
        grouped.setdefault(tuple(row[field] for field in group_fields), []).append(row)

    summary = []
    for key, group in grouped.items():
        item: dict[str, Any] = dict(zip(group_fields, key))
        item.update({
            "runs": len(group),
            "optimal": sum(int(row["optimal"]) for row in group),
            "mean_final_lower": mean(group, "final_lower"),
            "mean_final_upper": mean(group, "final_upper"),
            "sd_final_upper": deviation(group, "final_upper"),
            "mean_absolute_gap": mean(group, "absolute_gap"),
            "mean_relative_gap": mean(group, "relative_gap"),
            "mean_wall_seconds": mean(group, "wall_seconds"),
            "median_wall_seconds": statistics.median(values(group, "wall_seconds")),
            "hig_found_final_incumbent": sum(int(row["hig_found_final_incumbent"]) for row in group),
            "schedule_metrics_complete": sum(row.get("metrics_status") == "complete" for row in group),
            "mean_nominal_maximum_tardiness": mean(group, "nominal_maximum_tardiness"),
            "mean_nominal_maximum_completion_time": mean(group, "nominal_maximum_completion_time"),
            "sd_nominal_maximum_completion_time": deviation(group, "nominal_maximum_completion_time"),
            "mean_nominal_mean_completion_time": mean(group, "nominal_mean_completion_time"),
            "mean_nominal_total_recharges": mean(group, "nominal_total_recharges"),
            "mean_actual_recharge_time": mean(group, "actual_recharge_time"),
            "mean_capacity": mean(group, "mean_capacity"),
        })
        summary.append(item)
    summary.sort(key=lambda row: tuple(float(row[field]) for field in group_fields))
    summary_fields = group_fields + [field for field in summary[0] if field not in group_fields]
    atomic_csv(output / "sensitivity_summary.csv", summary_fields, summary)

    if design["study"] == "gamma":
        reference_key = (number(design["base_gamma"]),)
    elif design["study"] == "deviation":
        reference_key = (number(design["deviation"]),)
    else:
        reference_key = (number(design["base_capacity_scale"]),
                         number(design["base_recharge_factor"]))
    # CSV conversion may use a shorter equivalent spelling than number().
    reference = {
        int(row["replication"]): row for row in rows
        if all(math.isclose(float(row[field]), float(target), rel_tol=0, abs_tol=1e-12)
               for field, target in zip(group_fields, reference_key))
    }
    if len(reference) != len(design["replications"]):
        fail("reference setting is incomplete in the sensitivity results")
    paired = []
    for key, group in grouped.items():
        pairs = [(reference[int(row["replication"])], row) for row in group]
        objective = [float(row["final_upper"])-float(base["final_upper"]) for base, row in pairs]
        makespan = [
            float(row["nominal_maximum_completion_time"])-float(base["nominal_maximum_completion_time"])
            for base, row in pairs
            if row.get("metrics_status") == "complete" and base.get("metrics_status") == "complete"
        ]
        recharge = [
            float(row["nominal_total_recharges"])-float(base["nominal_total_recharges"])
            for base, row in pairs
            if row.get("metrics_status") == "complete" and base.get("metrics_status") == "complete"
        ]
        item = dict(zip(group_fields, key))
        item.update({
            "pairs": len(pairs),
            "both_optimal": sum(int(base["optimal"]) and int(row["optimal"]) for base, row in pairs),
            "mean_final_upper_change_from_reference": statistics.mean(objective),
            "upper_better": sum(value < 0 for value in objective),
            "upper_tied": sum(value == 0 for value in objective),
            "upper_worse": sum(value > 0 for value in objective),
            "completion_pairs": len(makespan),
            "mean_nominal_maximum_completion_change_from_reference": (
                statistics.mean(makespan) if makespan else ""
            ),
            "mean_nominal_recharge_change_from_reference": (
                statistics.mean(recharge) if recharge else ""
            ),
        })
        paired.append(item)
    paired.sort(key=lambda row: tuple(float(row[field]) for field in group_fields))
    paired_fields = group_fields + [field for field in paired[0] if field not in group_fields]
    atomic_csv(output / "sensitivity_pairs.csv", paired_fields, paired)


def main() -> int:
    args = parser().parse_args()
    design = normalize(args)
    cases = len(design["settings"])*len(design["replications"])
    hours = cases*design["total_limit"]/3600
    print(
        f"Validated {len(design['settings'])} settings and {len(design['replications'])} "
        f"replications: {cases} lbbd-hig runs, at most {hours:g} solver hours."
    )
    if args.dry_run:
        return 0

    output = args.output.expanduser().resolve()
    build = args.build.expanduser().resolve()
    expected = json.dumps(design, indent=2, sort_keys=True, allow_nan=False) + "\n"
    design_path = output / "design.json"
    if output.exists():
        if not args.resume:
            fail(f"output directory already exists: {output}; use --resume or another path")
        if not design_path.is_file() or design_path.read_text() != expected:
            fail("resume parameters do not match the saved sensitivity design")
        validate_prepared(output)
    else:
        output.mkdir(parents=True)
        atomic_text(design_path, expected)
        prepare(output, build, design)

    run_solver(output, build, design, args.prepare_only)
    if args.prepare_only:
        print(f"Prepared {cases} frozen sensitivity instances: {output}")
        return 0
    rows = add_metrics(output, design)
    summarize(output, design, rows)
    print(
        f"Completed sensitivity study. Results: {output / 'sensitivity_results.csv'}; "
        f"summary: {output / 'sensitivity_summary.csv'}; paired changes: "
        f"{output / 'sensitivity_pairs.csv'}"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        raise SystemExit(1)
