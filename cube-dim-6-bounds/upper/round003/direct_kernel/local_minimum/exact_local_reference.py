#!/usr/bin/env python3
"""Exact fmpq/Arb reference for one integer-gated local-search point."""
from __future__ import annotations

import hashlib
import itertools
import json
import sys
from pathlib import Path
from time import perf_counter

from flint import arb, ctx, fmpq

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent


def parse_fan():
    tokens = iter((PARENT / "fan_integer_input.txt").read_text().split())
    count = int(next(tokens)); cones = []
    for _ in range(count):
        determinant = int(next(tokens)); rays = []
        for _ in range(6):
            q = int(next(tokens)); rays.append((q, tuple(int(next(tokens)) for _ in range(6))))
        cones.append((determinant, rays))
    return cones


def dd5(nodes):
    nodes = sorted(nodes)
    assert all(nodes[i] != nodes[i - 1] for i in range(1, 6))
    table = [[fmpq(0) for _ in range(6)] for _ in range(6)]
    for i, node in enumerate(nodes): table[i][0] = fmpq(1 if node > 0 else -1)
    for order in range(1, 6):
        for i in range(order, 6): table[i][order] = (table[i][order - 1] - table[i - 1][order - 1]) / (nodes[i] - nodes[i - order])
    return table[5][5]


def main() -> None:
    index = int(sys.argv[1]) if len(sys.argv) == 2 else 0
    values = iter((HERE / "local_points.txt").read_text().split())
    denominator, count = int(next(values)), int(next(values))
    points = [tuple(int(next(values)) for _ in range(6)) for _ in range(count)]
    assert 0 <= index < count
    cones = parse_fan()
    started = perf_counter(); total = fmpq(0); terms = 0
    for perm in itertools.permutations(range(6)):
        for mask in range(64):
            y = tuple((-1 if (mask >> i) & 1 else 1) * points[index][perm[i]] for i in range(6))
            for determinant, rays in cones:
                nodes = []; alpha = fmpq(determinant)
                for q, ray in rays:
                    numerator = sum(y[i] * ray[i] for i in range(6))
                    assert numerator != 0
                    nodes.append(fmpq(numerator, denominator * q)); alpha /= q
                total += alpha * dd5(nodes); terms += 1
    q_over_pi5 = total / 128
    ctx.prec = 384
    interval = (arb(int(q_over_pi5.numerator)) / int(q_over_pi5.denominator)) / arb.pi()**5
    payload = {
        "status": "PASS_EXACT_FMPQ_LOCAL_REFERENCE", "index": index,
        "point_numerators": list(points[index]), "denominator": denominator, "terms": terms,
        "Q_over_pi5_numerator": str(q_over_pi5.numerator), "Q_over_pi5_denominator": str(q_over_pi5.denominator),
        "G_arb_384": interval.str(80, more=True), "elapsed_seconds": perf_counter() - started,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "fan_input_sha256": hashlib.sha256((PARENT / "fan_integer_input.txt").read_bytes()).hexdigest(),
        "points_input_sha256": hashlib.sha256((HERE / "local_points.txt").read_bytes()).hexdigest(),
    }
    (HERE / "exact_local_reference_results.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({k: payload[k] for k in payload if not k.startswith("Q_")}, indent=2))


if __name__ == "__main__": main()
