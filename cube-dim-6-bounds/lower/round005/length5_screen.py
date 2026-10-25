#!/usr/bin/env python3
"""Finite numeric screen for up to 11,429 new length-five R6 witnesses."""
from functools import lru_cache
from itertools import product
from math import prod
from pathlib import Path
import hashlib
import json
import time

TOL = 2e-12
NEW_ELL_CAP = 11429
WALL_SECONDS = 20 * 60
HERE = Path(__file__).resolve()
MAIN_LEDGER = Path("/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/results/bounds_ledger.json")


def spacing():
    u, v = 2 ** (1 / 6), 3 ** (1 / 6)
    return (1 + u**4 / 2, (u**3 + u**4 * v**3) / 4,
            *((1 + u) / 2 for _ in range(4)))


def witness_box(a, ell):
    radius = prod(ell) ** (1 / 6) / 2
    box = []
    for aj, k in zip(a, ell):
        if k % 2:
            lo, hi = 0.0, min(aj / 2, radius - (k - 1) * aj / 2)
            start = -(k // 2)
        else:
            lo, hi = max(0.0, k * aj / 2 - radius), aj / 2
            start = 1 - k // 2
        if hi <= lo + TOL:
            return None
        for x in (lo, hi):
            for atom in (start * aj, (start + k - 1) * aj):
                if abs(x - atom) > radius + 4e-12:
                    raise AssertionError((ell, aj, k, x, atom, radius))
        box += [lo, hi]
    return tuple(box)


def contained(inner, outer):
    return all(outer[2*j] <= inner[2*j] + TOL and inner[2*j+1] <= outer[2*j+1] + TOL
               for j in range(6))


def maximal(records):
    unique = {box: ell for ell, box in records}
    keep = []
    for box, ell in unique.items():
        if any(other != box and contained(box, other) for other in unique):
            continue
        keep.append((ell, box))
    return keep


def box_union(B):
    @lru_cache(None)
    def rec(X):
        if not X:
            return 0.0
        if len(X[0]) == 2:
            ints = sorted(X)
            lo, hi = ints[0]
            result = 0.0
            for c, d in ints[1:]:
                if c > hi + TOL:
                    result += hi - lo
                    lo, hi = c, d
                elif d > hi:
                    hi = d
            return result + hi - lo
        result = 0.0
        cuts = sorted({x for b in X for x in b[:2]})
        for lo, hi in zip(cuts, cuts[1:]):
            if hi <= lo + TOL:
                continue
            child = tuple(sorted({b[2:] for b in X if b[0] <= lo + TOL and hi <= b[1] + TOL}))
            result += (hi - lo) * rec(child)
        return result
    return rec(tuple(sorted(B)))


def main():
    started = time.monotonic()
    ledger_bytes = MAIN_LEDGER.read_bytes()
    a = spacing()
    old = [(ell, witness_box(a, ell)) for ell in product(range(1, 5), repeat=6)]
    old = [(ell, box) for ell, box in old if box is not None]
    old_keep = maximal(old)
    new_all = [ell for ell in product(range(1, 6), repeat=6) if 5 in ell]
    new_all.sort(key=lambda e: (-prod(e), e))
    selected, omitted = new_all[:NEW_ELL_CAP], new_all[NEW_ELL_CAP:]
    out = {"status": "running", "scope": "numeric discovery; no C6 claim",
           "spacing": "a2=(u^3+u^4*v^3)/4", "ledger_sha256": hashlib.sha256(ledger_bytes).hexdigest(),
           "frozen_complete_upper": json.loads(ledger_bytes)["lower"]["exact_interval"][1],
           "full_new_index_count": len(new_all), "new_index_cap": NEW_ELL_CAP,
           "ordering": "decreasing product(ell), lexicographic tie break",
           "untested_new_index_count": len(omitted), "untested_new_index_sha256": hashlib.sha256(repr(omitted).encode()).hexdigest(),
           "old_reference_raw_nonempty": len(old), "old_reference_maximal": len(old_keep)}
    new = []
    for n, ell in enumerate(selected, 1):
        if time.monotonic() - started > WALL_SECONDS:
            out.update(status="wall_time_stop", screened_new_indices=n-1)
            break
        box = witness_box(a, ell)
        if box is not None:
            new.append((ell, box))
    else:
        out["status"] = "screen_completed"
        out["screened_new_indices"] = len(selected)
    direct_contained = [(ell, box) for ell, box in new if any(contained(box, oldbox) for _, oldbox in old_keep)]
    not_direct = [(ell, box) for ell, box in new if not any(contained(box, oldbox) for _, oldbox in old_keep)]
    out.update(new_raw_nonempty=len(new), new_directly_contained_by_old_maximal=len(direct_contained),
               new_not_directly_contained=len(not_direct),
               new_not_directly_contained_ells=[list(ell) for ell, _ in not_direct])
    # Only the non-contained new boxes can affect the old union; deduplicate and
    # remove any boxes contained by another member of the joint collection.
    joint_keep = maximal(old_keep + not_direct)
    out["joint_maximal_boxes"] = len(joint_keep)
    old_ratio = 64 * box_union(tuple(box for _, box in old_keep))
    joint_ratio = 64 * box_union(tuple(box for _, box in joint_keep))
    out.update(old_reference_ratio_float=old_ratio, joint_ratio_float=joint_ratio,
               numeric_gain_over_old_reference=joint_ratio-old_ratio,
               requires_exact_promotion=joint_ratio > 2.0669238260108809 + 1e-7,
               exact_direct_containment_ready=(len(not_direct) == 0),
               elapsed_seconds=time.monotonic()-started,
               code_sha256=hashlib.sha256(HERE.read_bytes()).hexdigest())
    HERE.with_name("length5_screen_results.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
