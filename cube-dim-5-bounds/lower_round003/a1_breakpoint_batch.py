#!/usr/bin/env python3
"""Run one explicit batch of eight untested a1 endpoint crossings."""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import glob
import hashlib
import json
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from a1_breakpoint_search import endpoint_forms, evaluate, lift
sys.path.insert(0, str(HERE.parents[0] / "lower"))
from certify_rectangular import field


def reconstruct(F, sparse):
    coeffs = [Q(0)] * 25
    for key, value in sparse.items():
        i, j = map(int, key.split(","))
        coeffs[5*i + j] = Q(value)
    return F(a=tuple(coeffs))


def prior_coefficients():
    seen = set()
    for name in glob.glob(str(HERE / "a1_breakpoint*_results.json")):
        data = json.loads(Path(name).read_text())
        for row in data.get("candidates", []):
            sparse = row.get("a1_coefficients")
            if sparse:
                seen.add(tuple(sorted(sparse.items())))
    return seen


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch", type=int, required=True)
    args = parser.parse_args()
    if args.batch < 2:
        raise ValueError("batch 1 is the separately locked initial D result")
    started = time.monotonic()
    F, u, v, bounds, _, _ = field()
    center = 1 + u ** 3 / 2
    left, right = center - lift(F, Q(1, 8)), center + lift(F, Q(1, 8))
    forms = endpoint_forms(F, u, v)
    values = {left.a: left, right.a: right}
    for b1, m1 in forms:
        for b2, m2 in forms:
            if m1 != m2:
                point = (b2 - b1) / (m1 - m2)
                if left <= point <= right:
                    values[point.a] = point
    def sparse(point):
        return tuple(sorted((f"{i//5},{i%5}", str(c)) for i, c in enumerate(point.a) if c))
    used = prior_coefficients()
    def distance(point):
        d = point - center if center <= point else center - point
        return bounds(d.a)[0], point
    remaining = [point for point in values.values() if point != center and sparse(point) not in used]
    chosen = sorted(remaining, key=distance)[:8]
    lock_path = HERE / "lock_c_a2_equals_a3.json"
    lock = json.loads(lock_path.read_text())
    upper = Q(lock["exact_interval"][1])
    rows = []
    for point in chosen:
        row = evaluate(F, u, v, bounds, point, time.monotonic() + 120)
        row["strictly_beats_C_exact_upper"] = Q(row["ratio_interval"][0]) > upper
        rows.append(row)
    out = {
        "status": "completed_exact_a1_breakpoint_continuation_batch",
        "batch": args.batch,
        "total_exact_breakpoints_in_window": len(values),
        "previously_evaluated_count": len(used),
        "remaining_before_batch": len(remaining),
        "selected_count": len(chosen),
        "selection_rule": "exact nearest unrun values to C; ties exact ascending",
        "candidates": rows,
        "promoted_count": sum(row["strictly_beats_C_exact_upper"] for row in rows),
        "elapsed_seconds": time.monotonic() - started,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "C_lock_sha256": hashlib.sha256(lock_path.read_bytes()).hexdigest(),
        "scope": "Finite full-family endpoint batch; only the listed points are excluded in later batches."
    }
    raw = json.dumps(out, indent=2) + "\n"
    output = HERE / f"a1_breakpoint_batch_{args.batch}_results.json"
    if output.exists():
        raise FileExistsError(output)
    output.write_text(raw)
    print(raw, end="")


if __name__ == "__main__":
    main()
