#!/usr/bin/env python3
# Copyright (c) 2026 Renjie Yu and Daniel Oron
# SPDX-License-Identifier: MIT
"""Verify that two prepared experiment panels use identical overlap instances."""

from __future__ import annotations

import argparse
import csv
import hashlib
import pathlib
import sys
from typing import Any


MATCH_FIELDS = (
    "n", "m", "gamma", "deviation", "processing_max",
    "nominal_energy_max", "capacity_scale", "due_date_factor",
    "tardiness_factor", "due_date_range", "recharge_factor",
    "requested_recharge_time", "actual_recharge_time", "replication",
    "instance_seed",
)


def fail(message: str) -> None:
    raise ValueError(message)


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_rows(path: pathlib.Path) -> list[dict[str, str]]:
    try:
        with path.open(newline="") as source:
            return list(csv.DictReader(source))
    except OSError as error:
        fail(f"cannot read {path}: {error}")


def key(row: dict[str, str]) -> tuple[str, ...]:
    try:
        return tuple(row[field] for field in MATCH_FIELDS)
    except KeyError as error:
        fail(f"specification file lacks field {error.args[0]}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compare prepared experiment instances by factors and SHA-256 hash."
    )
    parser.add_argument("--baseline", required=True, type=pathlib.Path)
    parser.add_argument("--scalability", required=True, type=pathlib.Path)
    parser.add_argument("--jobs", nargs="+", type=int, default=(20, 30))
    parser.add_argument("--output", type=pathlib.Path)
    parser.add_argument(
        "--allow-left-subset", action="store_true",
        help="allow the right panel to contain additional selected specifications",
    )
    args = parser.parse_args()
    baseline = args.baseline.expanduser().resolve()
    scalability = args.scalability.expanduser().resolve()
    selected = {str(value) for value in args.jobs}
    if any(value < 1 for value in args.jobs):
        fail("--jobs values must be positive")

    left = [
        row for row in read_rows(baseline / "specifications.csv")
        if row["n"] in selected
    ]
    right = [
        row for row in read_rows(scalability / "specifications.csv")
        if row["n"] in selected
    ]
    if not left or not right:
        fail("no selected overlap instances were found")
    right_by_key: dict[tuple[str, ...], dict[str, str]] = {}
    for row in right:
        identifier = key(row)
        if identifier in right_by_key:
            fail(f"duplicate scalability specification {identifier}")
        right_by_key[identifier] = row

    records: list[dict[str, Any]] = []
    left_keys = set()
    for row in left:
        identifier = key(row)
        if identifier in left_keys:
            fail(f"duplicate baseline specification {identifier}")
        left_keys.add(identifier)
        match = right_by_key.get(identifier)
        if match is None:
            fail(f"baseline specification has no scalability match {identifier}")
        left_path = baseline / "instances" / f"{row['case']}.agv"
        right_path = scalability / "instances" / f"{match['case']}.agv"
        left_hash, right_hash = sha256(left_path), sha256(right_path)
        if left_hash != right_hash:
            fail(f"instance hash mismatch for {identifier}")
        record: dict[str, Any] = {field: row[field] for field in MATCH_FIELDS}
        record.update({
            "baseline_case": row["case"],
            "scalability_case": match["case"],
            "sha256": left_hash,
        })
        records.append(record)
    extra = set(right_by_key) - left_keys
    if extra and not args.allow_left_subset:
        fail(f"scalability panel contains {len(extra)} unmatched selected specifications")

    output = (
        args.output.expanduser().resolve()
        if args.output is not None else scalability / "overlap_manifest.csv"
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".tmp")
    fields = list(MATCH_FIELDS) + ["baseline_case", "scalability_case", "sha256"]
    with temporary.open("w", newline="") as destination:
        writer = csv.DictWriter(destination, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)
    temporary.replace(output)
    print(f"Verified {len(records)} byte-identical overlap instances. Manifest: {output}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        raise SystemExit(1)
