#!/usr/bin/env python3
"""Exact completion check for the finite registered a1 breakpoint window."""
from __future__ import annotations

from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from a1_breakpoint_search import endpoint_forms, lift
sys.path.insert(0, str(HERE.parents[0] / "lower"))
from certify_rectangular import field


def key(row):
    return tuple(sorted(row["a1_coefficients"].items()))


def main():
    F, u, v, _, _, _ = field()
    center = 1 + u ** 3 / 2
    left, right = center - lift(F, Q(1, 8)), center + lift(F, Q(1, 8))
    values = {left.a: left, right.a: right}
    forms = endpoint_forms(F, u, v)
    for b1, m1 in forms:
        for b2, m2 in forms:
            if m1 != m2:
                x = (b2 - b1) / (m1 - m2)
                if left <= x <= right:
                    values[x.a] = x
    # Do not glob our own exhaustion result.  Only the initial selected batch
    # and explicitly numbered continuation batches are producer inputs.
    files = [HERE / "a1_breakpoint_search_results.json"]
    files += sorted(HERE.glob("a1_breakpoint_batch_[0-9]*_results.json"),
                    key=lambda p: int(p.name.split("_")[3]))
    rows = []
    for name in files:
        rows.extend(json.loads(name.read_text())["candidates"])
    distinct = {key(row): row for row in rows}
    assert len(values) == 83
    assert len(rows) == len(distinct) == 82
    lock_path = HERE / "lock_c_a2_equals_a3.json"
    lock = json.loads(lock_path.read_text())
    lower = Q(lock["exact_interval"][0])
    assert all(lower > Q(row["ratio_interval"][1]) for row in distinct.values())
    # Each raw crossing except center must have been evaluated exactly once.
    expected = set()
    for x in values.values():
        if x != center:
            expected.add(tuple(sorted((f"{i//5},{i%5}", str(c)) for i, c in enumerate(x.a) if c)))
    assert set(distinct) == expected
    largest_other = max(distinct.values(), key=lambda row: Q(row["ratio_interval"][1]))
    output_path = Path(__file__).with_name("a1_breakpoint_exhaustion_results.json")
    old_result_hash = (hashlib.sha256(output_path.read_bytes()).hexdigest()
                       if output_path.exists() else None)
    out = {
        "status": "PASS_complete_finite_a1_breakpoint_window",
        "raw_affine_endpoint_form_count": len(forms),
        "exact_breakpoint_count": len(values),
        "C_center_count": 1,
        "other_breakpoints_evaluated_once": len(distinct),
        "C_exact_lower_exceeds_every_other_exact_upper": True,
        "largest_other_display_only": largest_other["display_only"],
        "C_display_only": lock["display_only"],
        "window": "a1 in [1+u^3/2-1/8,1+u^3/2+1/8], a2=a3=a4=a5=(1+u)/2",
        "conclusion_scope": "By the separately stated full-family endpoint-order piecewise-affine lemma, C is the strict maximizer of this finite witness-union ratio in this one-coordinate window.  This is not a global lattice optimum or a new C5 statement by itself.",
        "result_files": [name.name for name in files],
        "previous_result_sha256_before_replay_fix": old_result_hash,
        "C_lock_sha256": hashlib.sha256(lock_path.read_bytes()).hexdigest(),
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    }
    raw = json.dumps(out, indent=2) + "\n"
    output_path.write_text(raw)
    print(raw, end="")


if __name__ == "__main__":
    main()
