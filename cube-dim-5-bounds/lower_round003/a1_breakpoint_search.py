#!/usr/bin/env python3
"""Finite exact a1-breakpoint search using the full 1024-block endpoint family."""
from __future__ import annotations

from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import prod
from pathlib import Path
import hashlib
import json
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "lower"))
from certify_rectangular import field

NODE_LIMIT = 200_000


def lift(F, q: Q):
    return F(q.numerator) / q.denominator


def root_value(n: int, u, v):
    e2 = e3 = 0
    while n % 2 == 0:
        n //= 2; e2 += 1
    while n % 3 == 0:
        n //= 3; e3 += 1
    assert n == 1
    return u ** e2 * v ** e3


def endpoint_forms(F, u, v):
    """Return all raw affine endpoint forms b+m*a1 from the 1024 family."""
    forms = {(F(0).a, Q(0)), (F(0).a, Q(1, 2))}
    for lengths in product(range(1, 5), repeat=5):
        ell = lengths[0]
        radius = root_value(prod(lengths), u, v) / 2
        if ell % 2:
            forms.add((radius.a, Q(0)))
            forms.add((radius.a, -Q(ell - 1, 2)))
        else:
            forms.add(((-radius).a, Q(ell, 2)))
    return tuple((F(a=basis), slope) for basis, slope in forms)


def evaluate(F, u, v, bounds, a1, deadline):
    s = (1 + u) / 2
    spacing = (a1, s, s, s, s)
    boxes = []
    for lengths in product(range(1, 5), repeat=5):
        radius, box = root_value(prod(lengths), u, v) / 2, []
        for aa, ell in zip(spacing, lengths):
            if ell % 2:
                lo, hi = F(0), min(aa / 2, radius - (ell - 1) * aa / 2)
            else:
                lo, hi = max(F(0), ell * aa / 2 - radius), aa / 2
            if hi <= lo:
                break
            box.extend((lo, hi))
        else:
            boxes.append((lengths, tuple(box)))
    unique = {box: lengths for lengths, box in boxes}
    kept = [(lengths, box) for box, lengths in unique.items()
            if not any(other != box and all(other[2*j] <= box[2*j] and box[2*j+1] <= other[2*j+1]
                                             for j in range(5)) for other in unique)]
    nodes = 0
    @lru_cache(None)
    def sweep(items):
        nonlocal nodes
        nodes += 1
        if nodes > NODE_LIMIT:
            raise RuntimeError("recursive node checkpoint reached")
        if time.monotonic() > deadline:
            raise TimeoutError("120-second soft checkpoint reached")
        if not items:
            return F(0)
        if len(items[0]) == 2:
            left = right = None; out = F(0)
            for lo, hi in sorted(items):
                if left is None:
                    left, right = lo, hi
                elif right < lo:
                    out += right - left; left, right = lo, hi
                elif right < hi:
                    right = hi
            return out + right - left
        cuts = sorted({z for item in items for z in item[:2]})
        out = F(0)
        for lo, hi in zip(cuts, cuts[1:]):
            sub = tuple(sorted({item[2:] for item in items if item[0] <= lo and hi <= item[1]}))
            out += (hi - lo) * sweep(sub)
        return out
    retained = tuple(sorted(box for _, box in kept))
    recursive = sweep(retained)
    cells = []
    for j in range(5):
        cuts = sorted({box[2*j+k] for box in retained for k in (0, 1)})
        axis = []
        for lo, hi in zip(cuts, cuts[1:]):
            mask = sum(1 << idx for idx, box in enumerate(retained)
                       if box[2*j] <= lo and hi <= box[2*j+1])
            if mask:
                axis.append((hi - lo, mask))
        cells.append(axis)
    def grid(j=0, mask=None):
        if j == 5:
            return F(1)
        out = F(0)
        for length, nxt in cells[j]:
            use = nxt if mask is None else mask & nxt
            if use:
                out += length * grid(j + 1, use)
        return out
    assert grid() == recursive
    ratio = 32 * recursive
    lo, hi = bounds(ratio.a)
    return {"a1_coefficients": {f"{i//5},{i%5}": str(c) for i, c in enumerate(a1.a) if c},
            "ratio_interval": [str(lo), str(hi)], "display_only": [float(lo), float(hi)],
            "retained_witness_count": len(kept), "coordinate_segment_counts": [len(x) for x in cells],
            "union_engines_agree": True, "recursive_nodes": nodes,
            "ratio_coefficients": {f"{i//5},{i%5}": str(c) for i, c in enumerate(ratio.a) if c}}


def main():
    started = time.monotonic()
    F, u, v, bounds, roots, den = field()
    center = 1 + u ** 3 / 2
    left, right = center - lift(F, Q(1, 8)), center + lift(F, Q(1, 8))
    forms = endpoint_forms(F, u, v)
    candidates = {left.a: left, right.a: right}
    for b1, m1 in forms:
        for b2, m2 in forms:
            if m1 == m2:
                continue
            point = (b2 - b1) / (m1 - m2)
            if left <= point <= right:
                candidates[point.a] = point
    # a1=center is C itself, and is already locked: it is recorded but not rerun.
    all_points = sorted(candidates.values())
    unrun = [point for point in all_points if point != center]
    def exact_distance_key(point):
        distance = point - center if center <= point else center - point
        # Distance is nonnegative, so its isolated lower endpoint orders it
        # exactly; different endpoint crossings cannot be ranked from floats.
        return bounds(distance.a)[0], point
    ordered = sorted(unrun, key=exact_distance_key)
    selected = ordered if len(ordered) <= 8 else ordered[:8]
    reference = json.loads(Path(__file__).with_name("lock_c_a2_equals_a3.json").read_text())
    reference_upper = Q(reference["exact_interval"][1])
    rows = []
    for point in selected:
        row = evaluate(F, u, v, bounds, point, time.monotonic() + 120)
        row["strictly_beats_C_exact_upper"] = Q(row["ratio_interval"][0]) > reference_upper
        rows.append(row)
    result = {"status": "completed_exact_full_family_a1_breakpoint_batch",
              "window": {"center_coefficients": {f"{i//5},{i%5}": str(c) for i, c in enumerate(center.a) if c},
                         "radius": "1/8"},
              "raw_endpoint_affine_form_count": len(forms),
              "exact_breakpoint_count_including_window_endpoints": len(all_points),
              "selection_rule": "all if <=8 unrun; else nearest eight to center by exact absolute difference, then exact value",
              "selected_count": len(selected),
              "selected_a1_coefficients": [{f"{i//5},{i%5}": str(c) for i, c in enumerate(point.a) if c} for point in selected],
              "candidates": rows,
              "promoted_count": sum(row["strictly_beats_C_exact_upper"] for row in rows),
              "elapsed_seconds": time.monotonic() - started,
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "reference_lock_sha256": hashlib.sha256(Path(__file__).with_name("lock_c_a2_equals_a3.json").read_bytes()).hexdigest(),
              "field": "Q(u,v), u^5=2, v^5=3",
              "scope": "The raw endpoint enumeration protects this finite search from fixed-42-box assumptions; it is not a global lattice optimum proof."}
    raw = json.dumps(result, indent=2) + "\n"
    Path(__file__).with_name("a1_breakpoint_search_results.json").write_text(raw)
    print(raw, end="")


if __name__ == "__main__":
    main()
