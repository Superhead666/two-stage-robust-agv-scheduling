#!/usr/bin/env python3
# Copyright (c) 2026 Renjie Yu and Daniel Oron
# SPDX-License-Identifier: MIT
"""Prepare and run the reproducible Part 1 algorithm comparison."""

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
from dataclasses import asdict, dataclass
from typing import Any, Iterable


ROOT = pathlib.Path(__file__).resolve().parents[1]
SCHEMA = "ejor-part1-v5"


@dataclass(frozen=True)
class Configuration:
    method: str
    lbbd_binary: str | None = None
    ccg_binary: str | None = None


CONFIGURATIONS = {
    "lbbd-basic": Configuration("lbbd-basic", lbbd_binary="agv_lbbd"),
    "lbbd-hig": Configuration("lbbd-hig", lbbd_binary="agv_lbbd"),
    "ccg-basic": Configuration("ccg-basic", ccg_binary="agv_ccg"),
    "ccg-hig": Configuration("ccg-hig", ccg_binary="agv_ccg"),
    "milp-ccg-basic": Configuration("ccg-basic", ccg_binary="agv_ccg"),
    "milp-ccg-hig": Configuration("ccg-hig", ccg_binary="agv_ccg"),
    "cp-basic": Configuration("cp-basic", ccg_binary="agv_ccg_cp"),
    "cp-hig": Configuration("cp-hig", ccg_binary="agv_ccg_cp"),

}

DEFAULT_CONFIGURATIONS = ["lbbd-basic", "lbbd-hig", "ccg-basic", "ccg-hig"]


SPECIFICATION_FIELDS = [
    "case",
    "source_file",
    "n",
    "m",
    "gamma",
    "deviation",
    "processing_max",
    "nominal_energy_max",
    "capacity_scale",
    "due_date_factor",
    "tardiness_factor",
    "due_date_range",
    "recharge_factor",
    "requested_recharge_time",
    "actual_recharge_time",
    "replication",
    "instance_seed",
]

AGGREGATE_METADATA_FIELDS = [
    "configuration",
    "source_file",
    "deviation",
    "processing_max",
    "nominal_energy_max",
    "capacity_scale",
    "due_date_factor",
    "tardiness_factor",
    "due_date_range",
    "recharge_factor",
    "requested_recharge_time",
    "replication",
    "instance_seed",
]


def fail(message: str) -> None:
    raise ValueError(message)


def split_values(values: Iterable[str] | None, convert, option: str):
    if values is None:
        return None
    pieces = []
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
    return result


def finite_real(text: str) -> float:
    value = float(text)
    if not math.isfinite(value):
        fail(f"real values must be finite, found {text}")
    return value


def unique(values: list[Any], option: str) -> None:
    if len(values) != len(set(values)):
        fail(f"duplicate value in {option}")


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def json_text(value: Any) -> str:
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"


def atomic_text(path: pathlib.Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(content)
    os.replace(temporary, path)


def atomic_json(path: pathlib.Path, value: Any) -> None:
    atomic_text(path, json_text(value))


def atomic_csv(path: pathlib.Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", newline="") as destination:
        writer = csv.DictWriter(destination, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def load_json(path: pathlib.Path) -> Any:
    with path.open() as source:
        return json.load(source)


def run(command: list[pathlib.Path | str | int | float], log: pathlib.Path) -> None:
    log.parent.mkdir(parents=True, exist_ok=True)
    environment = os.environ.copy()
    environment["LC_ALL"] = "C"
    with log.open("w") as destination:
        completed = subprocess.run(
            [str(value) for value in command],
            stdout=destination,
            stderr=subprocess.STDOUT,
            env=environment,
            check=False,
        )
    if completed.returncode:
        fail(f"command exited with code {completed.returncode}; see {log}")


def read_header(path: pathlib.Path) -> tuple[int, int, int, float]:
    try:
        fields = path.read_text().splitlines()[0].split()
        if len(fields) != 4:
            fail("the first line must contain n, m, recharge time, and gamma")
        n, m, alpha, gamma = int(fields[0]), int(fields[1]), int(fields[2]), float(fields[3])
    except (IndexError, OSError, ValueError) as error:
        fail(f"invalid instance header in {path}: {error}")
    if n < 1 or n > 62 or m < 1 or m > 32 or alpha < 0 or not math.isfinite(gamma):
        fail(f"invalid instance header in {path}")
    if gamma < 0 or gamma > n:
        fail(f"uncertainty budget in {path} must lie in [0,n]")
    return n, m, alpha, gamma


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description=(
            "Run the Part 1 comparison on imported instances or on a generated "
            "factorial panel. List options accept spaces, commas, or both."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    result.add_argument("--instances", type=pathlib.Path,
                        help="directory of existing .agv files; omit to generate a panel")
    result.add_argument("--output", required=True, type=pathlib.Path)
    result.add_argument("--build", type=pathlib.Path, default=ROOT / "build",
                        help="directory containing newly built executables")
    result.add_argument("--configurations", nargs="+", default=DEFAULT_CONFIGURATIONS)
    result.add_argument("--total-limit", "--time-limit", dest="total_limit", type=finite_real,
                        default=1800.0)
    result.add_argument("--hig-limit", type=finite_real, default=60.0)
    result.add_argument("--seed", "--solver-seed", dest="solver_seed", type=int, default=1)

    result.add_argument("--target-jobs", nargs="+")
    result.add_argument("--agvs", nargs="+")
    result.add_argument("--size-mode", choices=("crossed", "paired"))
    result.add_argument(
        "--size-pairs", nargs="+", metavar="N:M",
        help="explicit job:AGV pairs, for example 40:3 40:5",
    )
    result.add_argument("--gamma", nargs="+")
    result.add_argument("--deviation", nargs="+")
    result.add_argument("--processing-max", nargs="+")
    result.add_argument("--nominal-energy-max", nargs="+")
    result.add_argument("--capacity-scale", nargs="+")
    dates = result.add_mutually_exclusive_group()
    dates.add_argument("--due-date-factor", nargs="+",
                       help="legacy multiplicative due-date factor")
    dates.add_argument("--tardiness-factor", nargs="+",
                       help="tardiness factors T for the reference-horizon design")
    result.add_argument("--due-date-range", type=finite_real, default=0.5,
                        help="relative due-date range R used with --tardiness-factor")
    recharge = result.add_mutually_exclusive_group()
    recharge.add_argument("--recharge-factor", nargs="+")
    recharge.add_argument("--recharge-time", nargs="+")
    replications = result.add_mutually_exclusive_group()
    replications.add_argument("--repetitions", type=int,
                              help="use replications 1 through this value")
    replications.add_argument("--replications", nargs="+",
                              help="explicit replication indices")
    result.add_argument("--generator-seed", type=int)
    result.add_argument("--pair-seed-stride", type=int, default=1000,
                        help="seed spacing between consecutive size--fleet pairs")
    result.add_argument("--pair-seed-offset", type=int, default=0,
                        help="number of preceding size--fleet pairs in a parent panel")

    result.add_argument("--dry-run", action="store_true")
    result.add_argument("--prepare-only", action="store_true")
    result.add_argument("--resume", action="store_true")
    return result


def normalize_configurations(raw: list[str]) -> list[str]:
    names = split_values(raw, str, "--configurations")
    assert names is not None
    unique(names, "--configurations")
    unknown = [name for name in names if name not in CONFIGURATIONS]
    if unknown:
        fail(
            "unknown configuration(s): " + ", ".join(unknown)
            + "; available values are " + ", ".join(CONFIGURATIONS)
        )
    return names


def validate_common(args: argparse.Namespace) -> list[str]:
    configurations = normalize_configurations(args.configurations)
    if args.total_limit <= 0:
        fail("--total-limit must be positive")
    if args.hig_limit <= 0 or args.hig_limit >= args.total_limit:
        fail("--hig-limit must be positive and smaller than --total-limit")
    if args.solver_seed < 0 or args.solver_seed > 2_000_000_000:
        fail("--seed must lie in [0,2000000000]")
    if args.dry_run and args.prepare_only:
        fail("--dry-run and --prepare-only cannot be combined")
    if args.dry_run and args.resume:
        fail("--dry-run and --resume cannot be combined")
    return configurations


def panel_arguments_present(args: argparse.Namespace) -> list[str]:
    mapping = {
        "target_jobs": "--target-jobs",
        "agvs": "--agvs",
        "size_mode": "--size-mode",
        "size_pairs": "--size-pairs",
        "gamma": "--gamma",
        "deviation": "--deviation",
        "processing_max": "--processing-max",
        "nominal_energy_max": "--nominal-energy-max",
        "capacity_scale": "--capacity-scale",
        "due_date_factor": "--due-date-factor",
        "tardiness_factor": "--tardiness-factor",
        "recharge_factor": "--recharge-factor",
        "recharge_time": "--recharge-time",
        "repetitions": "--repetitions",
        "replications": "--replications",
        "generator_seed": "--generator-seed",
    }
    return [option for name, option in mapping.items() if getattr(args, name) is not None]


def size_pairs(jobs: list[int], agvs: list[int], mode: str) -> list[tuple[int, int]]:
    if mode == "crossed":
        return list(itertools.product(jobs, agvs))
    if len(jobs) == len(agvs):
        return list(zip(jobs, agvs))
    if len(jobs) == 1:
        return [(jobs[0], m) for m in agvs]
    if len(agvs) == 1:
        return [(n, agvs[0]) for n in jobs]
    fail("paired size mode requires equal list lengths or one singleton list")


def explicit_size_pairs(values: list[str]) -> list[tuple[int, int]]:
    pairs: list[tuple[int, int]] = []
    for token in values:
        for piece in token.split(","):
            fields = piece.split(":")
            if len(fields) != 2:
                fail(f"invalid --size-pairs value {piece!r}; expected N:M")
            try:
                pairs.append((int(fields[0]), int(fields[1])))
            except ValueError:
                fail(f"invalid --size-pairs value {piece!r}; expected integers")
    unique(pairs, "--size-pairs")
    return pairs


def generated_design(args: argparse.Namespace) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    if args.size_pairs is not None:
        if args.target_jobs is not None or args.agvs is not None or args.size_mode is not None:
            fail("--size-pairs cannot be combined with --target-jobs, --agvs, or --size-mode")
        pairs = explicit_size_pairs(args.size_pairs)
        jobs = sorted({n for n, _ in pairs})
        agvs = sorted({m for _, m in pairs})
        mode = "explicit"
    else:
        jobs = split_values(args.target_jobs, int, "--target-jobs") or [30]
        agvs = split_values(args.agvs, int, "--agvs") or [3, 5]
        mode = args.size_mode or "crossed"
        pairs = size_pairs(jobs, agvs, mode)
    gamma = split_values(args.gamma, finite_real, "--gamma") or [2.0, 4.0]
    deviation = split_values(args.deviation, finite_real, "--deviation") or [0.3]
    processing = split_values(args.processing_max, int, "--processing-max") or [20]
    energy = split_values(args.nominal_energy_max, int, "--nominal-energy-max") or [20]
    capacity = split_values(args.capacity_scale, finite_real, "--capacity-scale") or [4.0]
    if args.tardiness_factor is not None:
        tardiness = split_values(
            args.tardiness_factor, finite_real, "--tardiness-factor"
        )
        due = None
    else:
        due = split_values(args.due_date_factor, finite_real, "--due-date-factor") or [0.3, 0.4]
        tardiness = None
    if args.recharge_time is not None:
        recharge_time = split_values(args.recharge_time, int, "--recharge-time")
        recharge_factor = None
    else:
        recharge_factor = split_values(args.recharge_factor, finite_real,
                                       "--recharge-factor") or [2.0]
        recharge_time = None
    if args.replications is not None:
        replication_values = split_values(args.replications, int, "--replications")
        assert replication_values is not None
    else:
        repetitions = 5 if args.repetitions is None else args.repetitions
        replication_values = list(range(1, repetitions + 1))
    generator_seed = 960000 if args.generator_seed is None else args.generator_seed

    for values, option in (
        (jobs, "--target-jobs"), (agvs, "--agvs"), (gamma, "--gamma"),
        (deviation, "--deviation"), (processing, "--processing-max"),
        (energy, "--nominal-energy-max"), (capacity, "--capacity-scale"),
    ):
        unique(values, option)
    if due is not None:
        unique(due, "--due-date-factor")
    if tardiness is not None:
        unique(tardiness, "--tardiness-factor")
    unique(replication_values, "--replications")
    if recharge_factor is not None:
        unique(recharge_factor, "--recharge-factor")
    if recharge_time is not None:
        unique(recharge_time, "--recharge-time")

    if any(n < 1 or n > 62 for n in jobs):
        fail("--target-jobs values must lie in [1,62]")
    if any(m < 1 or m > 32 for m in agvs):
        fail("--agvs values must lie in [1,32]")
    if any(value < 0 for value in gamma):
        fail("--gamma values must be nonnegative")
    if any(value < 0 for value in deviation):
        fail("--deviation values must be nonnegative")
    if any(value < 1 for value in processing):
        fail("--processing-max values must be positive")
    if any(value < 1 for value in energy):
        fail("--nominal-energy-max values must be positive")
    if any(value <= 0 for value in capacity):
        fail("--capacity-scale values must be positive")
    if due is not None and any(value < 0 for value in due):
        fail("--due-date-factor values must be nonnegative")
    if tardiness is not None and any(value < 0 or value > 1 for value in tardiness):
        fail("--tardiness-factor values must lie in [0,1]")
    if args.due_date_range < 0:
        fail("--due-date-range must be nonnegative")
    if recharge_factor is not None and any(value < 0 for value in recharge_factor):
        fail("--recharge-factor values must be nonnegative")
    if recharge_time is not None and any(value < 0 for value in recharge_time):
        fail("--recharge-time values must be nonnegative")
    if any(value < 1 for value in replication_values):
        fail("replication indices must be positive")
    if args.pair_seed_stride < 1 or max(replication_values) > args.pair_seed_stride:
        fail("--pair-seed-stride must be at least the largest replication index")
    if args.pair_seed_offset < 0:
        fail("--pair-seed-offset must be nonnegative")
    if generator_seed < 0 or generator_seed > 2_000_000_000:
        fail("--generator-seed must lie in [0,2000000000]")

    factor_values = recharge_time if recharge_time is not None else recharge_factor
    assert factor_values is not None
    date_values = tardiness if tardiness is not None else due
    assert date_values is not None
    specifications = []
    index = 0
    for size_index, (n, m) in enumerate(pairs):
        for values in itertools.product(
            gamma, deviation, processing, energy, capacity, date_values, factor_values,
            replication_values,
        ):
            g, dev, pmax, emax, cap, date_value, recharge_value, replication = values
            if g > n:
                fail(f"gamma {g} exceeds n={n}")
            instance_seed = (
                generator_seed + (args.pair_seed_offset+size_index) * args.pair_seed_stride + replication - 1
            )
            if instance_seed > 2_000_000_000:
                fail("a generated instance seed exceeds 2000000000")
            index += 1
            specifications.append({
                "case": f"case_{index:06d}",
                "source_file": "",
                "n": n,
                "m": m,
                "gamma": g,
                "deviation": dev,
                "processing_max": pmax,
                "nominal_energy_max": emax,
                "capacity_scale": cap,
                "due_date_factor": date_value if due is not None else "",
                "tardiness_factor": date_value if tardiness is not None else "",
                "due_date_range": args.due_date_range if tardiness is not None else "",
                "recharge_factor": recharge_value if recharge_factor is not None else "",
                "requested_recharge_time": recharge_value if recharge_time is not None else "",
                "actual_recharge_time": "",
                "replication": replication,
                "instance_seed": instance_seed,
            })

    panel = {
        "target_jobs": jobs,
        "agvs": agvs,
        "size_mode": mode,
        "size_pairs": [list(pair) for pair in pairs],
        "gamma": gamma,
        "deviation": deviation,
        "processing_max": processing,
        "nominal_energy_max": energy,
        "capacity_scale": capacity,
        "due_date_factor": due,
        "tardiness_factor": tardiness,
        "due_date_range": args.due_date_range if tardiness is not None else None,
        "recharge_factor": recharge_factor,
        "recharge_time": recharge_time,
        "replications": replication_values,
        "generator_seed": generator_seed,
        "pair_seed_stride": args.pair_seed_stride,
        "pair_seed_offset": args.pair_seed_offset,
        "common_random_numbers": "size_pair_and_replication",
    }
    return panel, specifications


def imported_design(args: argparse.Namespace) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    assert args.instances is not None
    if panel_arguments_present(args):
        fail("panel options cannot be combined with --instances")
    if args.pair_seed_offset:
        fail("--pair-seed-offset cannot be combined with --instances")
    source_root = args.instances.expanduser().resolve()
    if not source_root.is_dir():
        fail(f"instance directory does not exist: {source_root}")
    sources = sorted(path for path in source_root.rglob("*.agv") if path.is_file())
    if not sources:
        fail(f"no .agv files found in {source_root}")
    specifications = []
    inventory = []
    for index, source in enumerate(sources, 1):
        n, m, alpha, gamma = read_header(source)
        case = f"case_{index:06d}"
        relative = source.relative_to(source_root).as_posix()
        fingerprint = sha256(source)
        inventory.append({"case": case, "source_file": relative, "sha256": fingerprint})
        specifications.append({
            "case": case,
            "source_file": relative,
            "n": n,
            "m": m,
            "gamma": gamma,
            "deviation": "",
            "processing_max": "",
            "nominal_energy_max": "",
            "capacity_scale": "",
            "due_date_factor": "",
            "tardiness_factor": "",
            "due_date_range": "",
            "recharge_factor": "",
            "requested_recharge_time": "",
            "actual_recharge_time": alpha,
            "replication": "",
            "instance_seed": "",
        })
    return {"source_root": str(source_root), "files": inventory}, specifications


def make_design(args: argparse.Namespace, configurations: list[str]):
    if args.instances is None:
        source, specifications = generated_design(args)
        mode = "generated"
    else:
        source, specifications = imported_design(args)
        mode = "imported"
    design = {
        "schema": SCHEMA,
        "mode": mode,
        "source": source,
        "configurations": configurations,
        "total_limit": args.total_limit,
        "hig_limit": args.hig_limit,
        "solver_seed": args.solver_seed,
    }
    return design, specifications


def required_binary_names(configurations: list[str], generated: bool) -> list[str]:
    names = {"run_case"}
    if generated:
        names.add("agv_hig")
    for name in configurations:
        configuration = CONFIGURATIONS[name]
        if configuration.method.endswith("-hig"):
            names.add("agv_hig")
        if configuration.method == "lbbd-basic":
            names.add("agv_initialize")
        if configuration.lbbd_binary:
            names.add(configuration.lbbd_binary)
        if configuration.ccg_binary:
            names.add(configuration.ccg_binary)
    return sorted(names)


def freeze_binaries(build: pathlib.Path, output: pathlib.Path,
                    configurations: list[str], generated: bool) -> dict[str, str]:
    build = build.expanduser().resolve()
    destination = output / "bin"
    destination.mkdir(parents=True, exist_ok=True)
    fingerprints = {}
    for name in required_binary_names(configurations, generated):
        source = build / name
        if not source.is_file() or not os.access(source, os.X_OK):
            fail(f"missing executable {source}; build the requested configurations first")
        frozen = destination / name
        shutil.copy2(source, frozen)
        frozen.chmod(frozen.stat().st_mode | 0o111)
        fingerprints[name] = sha256(frozen)
    atomic_json(output / "binaries.json", fingerprints)
    return fingerprints


def validate_frozen_binaries(output: pathlib.Path, expected_names: list[str]) -> dict[str, str]:
    manifest_path = output / "binaries.json"
    if not manifest_path.is_file():
        fail(f"cannot resume because {manifest_path} is missing")
    saved = load_json(manifest_path)
    if sorted(saved) != sorted(expected_names):
        fail("the frozen executable set does not match the requested configurations")
    for name, expected in saved.items():
        path = output / "bin" / name
        if not path.is_file() or sha256(path) != expected:
            fail(f"frozen executable failed its SHA-256 check: {path}")
    return saved


def prepare_generated(output: pathlib.Path, specifications: list[dict[str, Any]],
                      regenerate: bool) -> None:
    generator = output / "bin" / "agv_hig"
    for specification in specifications:
        destination = output / "instances" / f"{specification['case']}.agv"
        if regenerate or not destination.is_file():
            temporary = destination.with_name(destination.name + ".tmp")
            temporary.unlink(missing_ok=True)
            command: list[pathlib.Path | str | int | float] = [
                generator, "--generate", temporary,
                "--n", specification["n"], "--m", specification["m"],
                "--seed", specification["instance_seed"],
                "--gamma", specification["gamma"],
                "--variation", specification["deviation"],
                "--capacity", specification["capacity_scale"],
                "--processing-max", specification["processing_max"],
                "--nominal-energy-max", specification["nominal_energy_max"],
            ]
            if specification["tardiness_factor"] != "":
                command += [
                    "--tardiness-factor", specification["tardiness_factor"],
                    "--due-date-range", specification["due_date_range"],
                ]
            else:
                command += ["--due", specification["due_date_factor"]]
            if specification["requested_recharge_time"] != "":
                command += ["--recharge-time", specification["requested_recharge_time"]]
            else:
                command += ["--recharge-factor", specification["recharge_factor"]]
            run(command, output / "instances" / f"{specification['case']}.log")
            os.replace(temporary, destination)
        n, m, alpha, gamma = read_header(destination)
        expected = (specification["n"], specification["m"], specification["gamma"])
        if (n, m, gamma) != expected:
            fail(f"generated header does not match the design for {specification['case']}")
        specification["actual_recharge_time"] = alpha


def prepare_imported(output: pathlib.Path, design: dict[str, Any],
                     specifications: list[dict[str, Any]], regenerate: bool) -> None:
    root = pathlib.Path(design["source"]["source_root"])
    inventory = {entry["case"]: entry for entry in design["source"]["files"]}
    for specification in specifications:
        entry = inventory[specification["case"]]
        source = root / entry["source_file"]
        if not source.is_file() or sha256(source) != entry["sha256"]:
            fail(f"source instance changed or is missing: {source}")
        destination = output / "instances" / f"{specification['case']}.agv"
        if regenerate or not destination.is_file():
            temporary = destination.with_name(destination.name + ".tmp")
            shutil.copy2(source, temporary)
            os.replace(temporary, destination)
        if sha256(destination) != entry["sha256"]:
            fail(f"frozen instance does not match its source: {destination}")


def instance_fingerprints(output: pathlib.Path,
                          specifications: list[dict[str, Any]]) -> dict[str, str]:
    return {
        f"{specification['case']}.agv": sha256(
            output / "instances" / f"{specification['case']}.agv"
        )
        for specification in specifications
    }


def verify_instance_fingerprints(output: pathlib.Path, expected: dict[str, str]) -> None:
    for name, fingerprint in expected.items():
        path = output / "instances" / name
        if not path.is_file() or sha256(path) != fingerprint:
            fail(f"frozen instance failed its SHA-256 check: {path}")


def prepare_instances(output: pathlib.Path, design: dict[str, Any],
                      specifications: list[dict[str, Any]], resume: bool) -> dict[str, str]:
    output.joinpath("instances").mkdir(parents=True, exist_ok=True)
    saved_path = output / "instances.json"
    if resume and saved_path.is_file():
        saved = load_json(saved_path)
        verify_instance_fingerprints(output, saved)
        if design["mode"] == "generated":
            for specification in specifications:
                _, _, alpha, _ = read_header(
                    output / "instances" / f"{specification['case']}.agv"
                )
                specification["actual_recharge_time"] = alpha
        else:
            prepare_imported(output, design, specifications, regenerate=False)
        return saved

    if design["mode"] == "generated":
        prepare_generated(output, specifications, regenerate=resume)
    else:
        prepare_imported(output, design, specifications, regenerate=resume)
    saved = instance_fingerprints(output, specifications)
    atomic_json(saved_path, saved)
    return saved


def write_protocol_files(output: pathlib.Path, design: dict[str, Any],
                         specifications: list[dict[str, Any]],
                         binaries: dict[str, str], instances: dict[str, str]) -> None:
    atomic_csv(output / "specifications.csv", SPECIFICATION_FIELDS, specifications)
    manifest = {
        "schema": SCHEMA,
        "design_sha256": sha256(output / "design.json"),
        "specifications_sha256": sha256(output / "specifications.csv"),
        "binaries": binaries,
        "instances": instances,
    }
    atomic_json(output / "manifest.json", manifest)


def validate_protocol_files(output: pathlib.Path, binaries: dict[str, str],
                            instances: dict[str, str]) -> None:
    manifest_path = output / "manifest.json"
    if not manifest_path.is_file():
        return
    manifest = load_json(manifest_path)
    if manifest.get("schema") != SCHEMA:
        fail("the saved manifest has an incompatible format")
    if manifest.get("design_sha256") != sha256(output / "design.json"):
        fail("the saved design failed its manifest check")
    if manifest.get("binaries") != binaries or manifest.get("instances") != instances:
        fail("the saved manifest does not match the frozen inputs")
    specifications = output / "specifications.csv"
    if (not specifications.is_file()
            or manifest.get("specifications_sha256") != sha256(specifications)):
        fail("the saved specifications failed their manifest check")


def command_for(output: pathlib.Path, case: str, configuration_name: str,
                args: argparse.Namespace) -> list[pathlib.Path | str | int | float]:
    configuration = CONFIGURATIONS[configuration_name]
    command: list[pathlib.Path | str | int | float] = [
        output / "bin" / "run_case",
        "--instance", output / "instances" / f"{case}.agv",
        "--output", output / "cases" / case / configuration_name,
        "--method", configuration.method,
        "--total-limit", args.total_limit,
        "--hig-limit", args.hig_limit,
        "--seed", args.solver_seed,
    ]
    if configuration.method.endswith("-hig"):
        command += ["--hig-binary", output / "bin" / "agv_hig"]
    if configuration.method == "lbbd-basic":
        command += ["--initialize-binary", output / "bin" / "agv_initialize"]
    if configuration.lbbd_binary:
        command += ["--lbbd-binary", output / "bin" / configuration.lbbd_binary]
    if configuration.ccg_binary:
        command += ["--ccg-binary", output / "bin" / configuration.ccg_binary]
    return command


def read_result(path: pathlib.Path) -> tuple[list[str], dict[str, str]]:
    try:
        with path.open(newline="") as source:
            reader = csv.DictReader(source)
            row = next(reader)
            if next(reader, None) is not None:
                fail(f"result file contains more than one row: {path}")
            if reader.fieldnames is None:
                fail(f"result file has no header: {path}")
            return reader.fieldnames, row
    except (OSError, StopIteration) as error:
        fail(f"invalid result file {path}: {error}")


def numeric_values(rows: list[dict[str, Any]], field: str) -> list[float]:
    return [float(row[field]) for row in rows if row.get(field, "") != ""]


def shifted_geometric_mean(values: list[float], shift: float = 1.0) -> float | str:
    if not values:
        return ""
    return math.exp(sum(math.log(value + shift) for value in values) / len(values)) - shift


def write_result_summaries(output: pathlib.Path, rows: list[dict[str, Any]]) -> None:
    summary_fields = [
        "configuration", "runs", "optimal", "hig_runs",
        "hig_found_final_incumbent", "hig_already_certified",
        "exact_proved_hig_incumbent", "incumbent_improved",
        "lower_bound_improved_only", "no_final_bound_improvement",
        "median_initial_gap", "median_final_gap", "median_relative_gap",
    ]
    summaries = []
    for configuration in sorted({str(row["configuration"]) for row in rows}):
        group = [row for row in rows if row["configuration"] == configuration]
        hig = [row for row in group if row["exact_progress_class"] != "not_applicable"]
        initial_gaps = numeric_values(hig, "initial_gap")
        summaries.append({
            "configuration": configuration,
            "runs": len(group),
            "optimal": sum(int(row["optimal"]) for row in group),
            "hig_runs": len(hig),
            "hig_found_final_incumbent": sum(
                int(row["hig_found_final_incumbent"]) for row in hig
            ),
            "hig_already_certified": sum(
                int(row["hig_already_certified"]) for row in hig
            ),
            "exact_proved_hig_incumbent": sum(
                int(row["exact_proved_hig_incumbent"]) for row in hig
            ),
            "incumbent_improved": sum(
                row["exact_progress_class"] == "incumbent_improved" for row in hig
            ),
            "lower_bound_improved_only": sum(
                row["exact_progress_class"] == "lower_bound_improved_only" for row in hig
            ),
            "no_final_bound_improvement": sum(
                int(row["exact_no_final_bound_improvement"]) for row in hig
            ),
            "median_initial_gap": statistics.median(initial_gaps) if initial_gaps else "",
            "median_final_gap": statistics.median(numeric_values(group, "absolute_gap")),
            "median_relative_gap": statistics.median(numeric_values(group, "relative_gap")),
        })
    atomic_csv(output / "method_summary.csv", summary_fields, summaries)

    warm_rows = [row for row in rows if row["exact_progress_class"] != "not_applicable"]
    warm_fields = [
        "instance", "configuration", "n", "m", "gamma", "tardiness_factor",
        "replication", "initial_lower", "initial_upper", "final_lower", "final_upper",
        "initial_gap", "absolute_gap", "upper_improvement", "lower_improvement",
        "gap_reduction", "hig_found_final_incumbent", "hig_already_certified",
        "exact_proved_hig_incumbent", "exact_no_final_bound_improvement",
        "exact_progress_class",
    ]
    atomic_csv(output / "warm_start_cases.csv", warm_fields, warm_rows)


def aggregate(output: pathlib.Path, specifications: list[dict[str, Any]],
              configurations: list[str]) -> int:
    metadata = {row["case"]: row for row in specifications}
    rows: list[dict[str, Any]] = []
    base_fields: list[str] | None = None
    for specification in specifications:
        case = specification["case"]
        for configuration_name in configurations:
            result = output / "cases" / case / configuration_name / "result.csv"
            if not result.is_file():
                continue
            fields, row = read_result(result)
            expected_method = CONFIGURATIONS[configuration_name].method
            if row.get("method") != expected_method:
                fail(f"method mismatch in {result}")
            if base_fields is None:
                base_fields = fields
            elif fields != base_fields:
                fail(f"result header mismatch in {result}")
            combined = dict(row)
            combined["configuration"] = configuration_name
            for field in AGGREGATE_METADATA_FIELDS[1:]:
                combined[field] = metadata[case][field]
            rows.append(combined)
    if base_fields is not None:
        atomic_csv(output / "results.csv", base_fields + AGGREGATE_METADATA_FIELDS, rows)
        write_result_summaries(output, rows)
    return len(rows)


def execute(output: pathlib.Path, specifications: list[dict[str, Any]],
            configurations: list[str], args: argparse.Namespace) -> None:
    output.joinpath("cases").mkdir(parents=True, exist_ok=True)
    completed = aggregate(output, specifications, configurations)
    total = len(specifications) * len(configurations)
    for index, specification in enumerate(specifications):
        order = configurations[index % len(configurations):] + configurations[:index % len(configurations)]
        for configuration_name in order:
            case = specification["case"]
            destination = output / "cases" / case / configuration_name
            result = destination / "result.csv"
            if result.is_file():
                read_result(result)
                continue
            if destination.exists():
                shutil.rmtree(destination)
            destination.mkdir(parents=True, exist_ok=True)
            print(
                f"Running {case} with {configuration_name} "
                f"({completed + 1}/{total})",
                flush=True,
            )
            try:
                run(
                    command_for(output, case, configuration_name, args),
                    destination / "runner.log",
                )
            finally:
                completed = aggregate(output, specifications, configurations)
    print(
        f"Completed {completed} solver runs. Results: {output / 'results.csv'}; "
        f"summary: {output / 'method_summary.csv'}; "
        f"HIG cases: {output / 'warm_start_cases.csv'}"
    )


def main() -> int:
    args = parser().parse_args()
    configurations = validate_common(args)
    design, specifications = make_design(args, configurations)
    cases = len(specifications)
    runs = cases * len(configurations)
    nominal_hours = runs * args.total_limit / 3600.0
    print(
        f"Validated {cases} instances, {len(configurations)} configurations, "
        f"and {runs} solver runs; maximum nominal solver time is "
        f"{nominal_hours:g} hours."
    )
    if args.dry_run:
        return 0

    output = args.output.expanduser().resolve()
    design_path = output / "design.json"
    expected_design = json_text(design)
    continuing = args.resume and design_path.is_file()
    if continuing:
        if design_path.read_text() != expected_design:
            fail("resume parameters or imported inputs do not match the saved design")
        binary_names = required_binary_names(configurations, design["mode"] == "generated")
        binaries = validate_frozen_binaries(output, binary_names)
    else:
        if output.exists():
            if args.resume:
                fail(f"cannot resume because {design_path} is missing")
            fail(f"output directory already exists: {output}; use --resume or another path")
        output.mkdir(parents=True)
        atomic_text(design_path, expected_design)
        binaries = freeze_binaries(
            args.build, output, configurations, design["mode"] == "generated"
        )

    instances = prepare_instances(output, design, specifications, continuing)
    validate_protocol_files(output, binaries, instances)
    write_protocol_files(output, design, specifications, binaries, instances)
    if args.prepare_only:
        print(f"Prepared {cases} frozen instances without starting solver runs: {output}")
        return 0
    execute(output, specifications, configurations, args)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        raise SystemExit(1)
