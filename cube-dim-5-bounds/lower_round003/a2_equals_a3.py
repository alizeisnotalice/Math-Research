#!/usr/bin/env python3
"""Exact finite witness test at the algebraic spacing breakpoint a2=a3."""
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


def main():
    started = time.monotonic()
    deadline = started + 120
    F, u, v, bounds, roots, denominator = field()
    # a2 equals a3, so its delta over the round-001 a2 is exactly stated.
    delta2 = (2 * u - u ** 4 * v) / 4
    a = (1 + u ** 3 / 2, (1 + u) / 2, (1 + u) / 2,
         (1 + u) / 2, (1 + u) / 2)

    def root(n: int):
        e2 = e3 = 0
        while n % 2 == 0:
            n //= 2; e2 += 1
        while n % 3 == 0:
            n //= 3; e3 += 1
        assert n == 1
        return u ** e2 * v ** e3

    boxes = []
    for lengths in product(range(1, 5), repeat=5):
        radius, box = root(prod(lengths)) / 2, []
        for spacing, ell in zip(a, lengths):
            if ell % 2:
                lo, hi = F(0), min(spacing / 2, radius - (ell - 1) * spacing / 2)
            else:
                lo, hi = max(F(0), ell * spacing / 2 - radius), spacing / 2
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
        cuts = sorted({x for item in items for x in item[:2]})
        answer = F(0)
        for lo, hi in zip(cuts, cuts[1:]):
            sub = tuple(sorted({item[2:] for item in items if item[0] <= lo and hi <= item[1]}))
            answer += (hi - lo) * sweep(sub)
        return answer

    retained = tuple(sorted(box for _, box in kept))
    union = sweep(retained)
    segments = []
    for j in range(5):
        cuts = sorted({box[2*j+k] for box in retained for k in (0, 1)})
        rows = []
        for lo, hi in zip(cuts, cuts[1:]):
            mask = sum((1 << index) for index, box in enumerate(retained)
                       if box[2*j] <= lo and hi <= box[2*j+1])
            if mask:
                rows.append((hi - lo, mask))
        segments.append(rows)

    def grid(j=0, mask=None):
        if j == 5:
            return F(1)
        answer = F(0)
        for length, next_mask in segments[j]:
            overlap = next_mask if mask is None else mask & next_mask
            if overlap:
                answer += length * grid(j + 1, overlap)
        return answer

    assert grid() == union
    ratio = 32 * union
    lo, hi = bounds(ratio.a)
    previous = json.loads(Path(__file__).with_name("lock_b_1_32.json").read_text())
    strict = Q(lo) > Q(previous["exact_interval"][1])
    output = {
        "status": "completed_exact_algebraic_breakpoint_test",
        "spacing": ["1+u^3/2", "(1+u)/2", "(1+u)/2", "(1+u)/2", "(1+u)/2"],
        "delta2_expression": "(2*u-u^4*v)/4",
        "ratio_interval": [str(lo), str(hi)],
        "display_only": [float(lo), float(hi)],
        "strictly_beats_003_B_exact_upper": strict,
        "retained_witness_count": len(kept),
        "coordinate_segment_counts": [len(s) for s in segments],
        "union_engines_agree": True,
        "recursive_nodes": nodes,
        "coefficient_basis": "u^i v^j, u^5=2, v^5=3",
        "coefficients": {f"{i//5},{i%5}": str(c) for i, c in enumerate(ratio.a) if c},
        "root_isolation_denominator": str(denominator),
        "root_lower_numerators": roots,
        "elapsed_seconds": time.monotonic() - started,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "locked_reference_sha256": hashlib.sha256(Path(__file__).with_name("lock_b_1_32.json").read_bytes()).hexdigest(),
        "scope": "Exact finite ratio only; a strict improvement would still require independent union and transfer review."
    }
    raw = json.dumps(output, indent=2) + "\n"
    Path(__file__).with_name("a2_equals_a3_results.json").write_text(raw)
    print(raw, end="")


if __name__ == "__main__":
    main()
