#!/usr/bin/env python3
"""Exact finite common-tail perturbation batch around the certified C lattice."""
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

LIMIT = 200_000


def lift(F, q):
    return F(q.numerator) / q.denominator


def root_value(n, u, v):
    e2 = e3 = 0
    while n % 2 == 0:
        n //= 2; e2 += 1
    while n % 3 == 0:
        n //= 3; e3 += 1
    assert n == 1
    return u ** e2 * v ** e3


def evaluate(F, u, v, bounds, delta, deadline):
    a1 = 1 + u ** 3 / 2
    tail = (1 + u) / 2 + lift(F, delta)
    spacing = (a1, tail, tail, tail, tail)
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
        if nodes > LIMIT:
            raise RuntimeError("recursive node checkpoint reached")
        if time.monotonic() > deadline:
            raise TimeoutError("120-second soft checkpoint reached")
        if not items:
            return F(0)
        if len(items[0]) == 2:
            left = right = None; answer = F(0)
            for lo, hi in sorted(items):
                if left is None:
                    left, right = lo, hi
                elif right < lo:
                    answer += right - left; left, right = lo, hi
                elif right < hi:
                    right = hi
            return answer + right - left
        cuts = sorted({z for item in items for z in item[:2]})
        answer = F(0)
        for lo, hi in zip(cuts, cuts[1:]):
            sub = tuple(sorted({item[2:] for item in items if item[0] <= lo and hi <= item[1]}))
            answer += (hi - lo) * sweep(sub)
        return answer
    retained = tuple(sorted(box for _, box in kept))
    recursive = sweep(retained)
    axes = []
    for j in range(5):
        cuts = sorted({box[2*j+k] for box in retained for k in (0, 1)})
        axis = []
        for lo, hi in zip(cuts, cuts[1:]):
            mask = sum(1 << index for index, box in enumerate(retained)
                       if box[2*j] <= lo and hi <= box[2*j+1])
            if mask:
                axis.append((hi - lo, mask))
        axes.append(axis)
    def grid(j=0, mask=None):
        if j == 5:
            return F(1)
        answer = F(0)
        for length, nxt in axes[j]:
            chosen = nxt if mask is None else mask & nxt
            if chosen:
                answer += length * grid(j + 1, chosen)
        return answer
    assert grid() == recursive
    ratio = 32 * recursive
    lo, hi = bounds(ratio.a)
    return {"delta": str(delta), "ratio_interval": [str(lo), str(hi)],
            "display_only": [float(lo), float(hi)], "retained_witness_count": len(kept),
            "coordinate_segment_counts": [len(axis) for axis in axes],
            "union_engines_agree": True, "recursive_nodes": nodes,
            "ratio_coefficients": {f"{i//5},{i%5}": str(c) for i, c in enumerate(ratio.a) if c}}


def main():
    started = time.monotonic()
    F, u, v, bounds, _, _ = field()
    lock_path = ROOT / "lower_round003" / "lock_c_a2_equals_a3.json"
    lock = json.loads(lock_path.read_text())
    upper = Q(lock["exact_interval"][1])
    deltas = tuple(Q(n, 256) for n in (-4, -3, -2, -1, 1, 2, 3, 4))
    rows = []
    for delta in deltas:
        row = evaluate(F, u, v, bounds, delta, time.monotonic() + 120)
        row["strictly_beats_C_exact_upper"] = Q(row["ratio_interval"][0]) > upper
        rows.append(row)
    result = {"status": "completed_exact_common_tail_batch", "candidate_count": len(rows),
              "fixed_a1": "1+u^3/2", "tail_base": "(1+u)/2", "candidates": rows,
              "promoted_count": sum(row["strictly_beats_C_exact_upper"] for row in rows),
              "elapsed_seconds": time.monotonic() - started,
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "C_lock_sha256": hashlib.sha256(lock_path.read_bytes()).hexdigest(),
              "scope": "Exact finite ratios at eight registered points only; no common-tail polynomial optimization claim."}
    raw = json.dumps(result, indent=2) + "\n"
    Path(__file__).with_name("common_tail_batch_results.json").write_text(raw)
    print(raw, end="")


if __name__ == "__main__":
    main()
