#!/usr/bin/env python3
"""One Arb calibration of the complete confluent H_epsilon cone sum.

This is a point certificate only.  Nodes are sorted and equated with integer
rational arithmetic before Arb evaluates atan and confluent jets.
"""
from __future__ import annotations

from functools import cmp_to_key
import hashlib
import itertools
import json
import sys
from pathlib import Path
from time import perf_counter

from flint import acb, arb, ctx

HERE = Path(__file__).resolve().parent


def parse_fan():
    tokens = iter((HERE / "../../round003/direct_kernel/fan_integer_input.txt").read_text().split())
    count = int(next(tokens)); cones = []
    for _ in range(count):
        determinant = int(next(tokens)); rays = []
        for _ in range(6):
            q = int(next(tokens)); rays.append((q, tuple(int(next(tokens)) for _ in range(6))))
        cones.append((determinant, rays))
    return cones


def cmp_node(a, b): return (a[0] * b[1] > b[0] * a[1]) - (a[0] * b[1] < b[0] * a[1])
def equal(a, b): return a[0] * b[1] == b[0] * a[1]
def aval(a): return arb(a[0]) / a[1]


def jet(x, epsilon, order):
    if order == 0: return (x / epsilon).atan()
    ii = acb(0, 1)
    first = ii**order / acb(epsilon, x)**order
    second = (-ii)**order / acb(epsilon, -x)**order
    result = ((1 if order & 1 else -1) * (first - second) / (2 * order * ii))
    return result.real


def hermite(nodes, epsilon):
    z = sorted(nodes, key=cmp_to_key(cmp_node))
    table = [[arb(0) for _ in range(6)] for _ in range(6)]
    values = [aval(a) for a in z]
    for i in range(6): table[i][0] = jet(values[i], epsilon, 0)
    for order in range(1, 6):
        for i in range(order, 6):
            table[i][order] = jet(values[i], epsilon, order) if equal(z[i], z[i-order]) else (table[i][order-1] - table[i-1][order-1]) / (values[i] - values[i-order])
    return table[5][5]


def main():
    point_index = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    epsilon_index = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    tokens = iter((HERE / "points.txt").read_text().split()); den, count = int(next(tokens)), int(next(tokens))
    points = [tuple(int(next(tokens)) for _ in range(6)) for _ in range(count)]
    epst = iter((HERE / "epsilons.txt").read_text().split()); epscount = int(next(epst)); epss = [(int(next(epst)), int(next(epst))) for _ in range(epscount)]
    assert 0 <= point_index < count and 0 <= epsilon_index < epscount
    cones = parse_fan(); ctx.prec = 256
    epsilon = arb(epss[epsilon_index][0]) / epss[epsilon_index][1]
    started = perf_counter(); total = arb(0); equal_events = 0; terms = 0
    for perm in itertools.permutations(range(6)):
        for mask in range(64):
            y = tuple((-1 if (mask >> i) & 1 else 1) * points[point_index][perm[i]] for i in range(6))
            for determinant, rays in cones:
                nodes = []; alpha = arb(determinant)
                for q, ray in rays:
                    nodes.append((sum(y[i] * ray[i] for i in range(6)), den * q)); alpha /= q
                ordered = sorted(nodes, key=cmp_to_key(cmp_node))
                equal_events += sum(equal(ordered[i-1], ordered[i]) for i in range(1, 6))
                total += alpha * hermite(nodes, epsilon); terms += 1
    result = total / (2 * arb.pi())**6
    payload = {
        "status": "PASS_ARB256_POINT_H_EPSILON", "point_index": point_index,
        "point_numerators": list(points[point_index]), "denominator": den,
        "epsilon": f"{epss[epsilon_index][0]}/{epss[epsilon_index][1]}", "terms": terms,
        "equal_adjacent_events": equal_events, "H_epsilon_arb_256": result.str(80, more=True),
        "elapsed_seconds": perf_counter()-started, "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "points_sha256": hashlib.sha256((HERE / "points.txt").read_bytes()).hexdigest(),
        "epsilons_sha256": hashlib.sha256((HERE / "epsilons.txt").read_bytes()).hexdigest(),
        "fan_sha256": hashlib.sha256((HERE / "../../round003/direct_kernel/fan_integer_input.txt").read_bytes()).hexdigest(),
    }
    (HERE / "arb_reference_results.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__": main()
