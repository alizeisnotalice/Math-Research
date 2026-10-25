#!/usr/bin/env python3
"""High-precision calibration of the six-node truncated Green identity.

The finite examples test a stable confluent Newton divided-difference routine
for F_e(b)=atan(b/e).  They are calibration evidence only: no cone cover or
continuous Green certificate is asserted here.
"""
from __future__ import annotations

from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
from time import perf_counter
import json
import mpmath as mp


def partitions(total: int, largest: int | None = None):
    if total == 0:
        yield ()
        return
    if largest is None or largest > total:
        largest = total
    for first in range(largest, 0, -1):
        for rest in partitions(total - first, first):
            yield (first,) + rest


def derivative_f(x, epsilon, order):
    return mp.diff(lambda z: mp.atan(z / epsilon), x, order)


def hermite_dd(nodes, epsilon):
    """Confluent DD in sorted-node Newton-table form."""
    z = sorted(mp.mpf(v) for v in nodes)
    n = len(z)
    table = [[mp.mpf('0') for _ in range(n)] for _ in range(n)]
    for i, x in enumerate(z):
        table[i][0] = mp.atan(x / epsilon)
    for order in range(1, n):
        for i in range(order, n):
            if z[i] == z[i - order]:
                table[i][order] = derivative_f(z[i], epsilon, order) / mp.factorial(order)
            else:
                table[i][order] = (table[i][order - 1] - table[i - 1][order - 1]) / (z[i] - z[i - order])
    return table[-1][-1]


def direct_integral(nodes, epsilon):
    def integrand(t):
        product = mp.mpc(1)
        for b in nodes:
            product /= t - 1j * b
        return mp.re(product)
    return mp.quad(integrand, [epsilon, 1, mp.inf])


def nodes_for(partition):
    centers = [-5, -3, -1, 0, 2, 4]
    return [centers[i] for i, multiplicity in enumerate(partition) for _ in range(multiplicity)]


def main() -> None:
    started = perf_counter()
    mp.mp.dps = 90
    epsilon = mp.mpf(3) / 2
    rows = []
    for multiplicities in partitions(6):
        nodes = nodes_for(multiplicities)
        dd = hermite_dd(nodes, epsilon)
        integral = direct_integral(nodes, epsilon)
        rows.append({
            "multiplicities": list(multiplicities),
            "nodes": list(map(str, nodes)),
            "dd": mp.nstr(dd, 45),
            "integral_re_product": mp.nstr(integral, 45),
            "absolute_residual": mp.nstr(abs(dd - integral), 8),
        })
    zero_anchor = hermite_dd([0] * 6, epsilon)
    zero_exact = 1 / (5 * epsilon ** 5)
    nonzero_anchor = hermite_dd([2] * 6, epsilon)
    nonzero_exact = mp.re(1 / (5 * (epsilon - 2j) ** 5))
    payload = {
        "status": "numerical_calibration_only",
        "epsilon": "3/2",
        "formula_tested": "integral_e^infty Re prod_j(t-i b_j)^-1 dt = DD5[atan(b/e); b1,...,b6]",
        "partition_types_tested": len(rows),
        "types": rows,
        "analytic_anchors": {
            "all_zero": {"computed": mp.nstr(zero_anchor, 45), "exact": mp.nstr(zero_exact, 45), "residual": mp.nstr(abs(zero_anchor-zero_exact), 8)},
            "all_two": {"computed": mp.nstr(nonzero_anchor, 45), "exact": mp.nstr(nonzero_exact, 45), "residual": mp.nstr(abs(nonzero_anchor-nonzero_exact), 8)}
        },
        "elapsed_seconds": perf_counter() - started,
        "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    Path(__file__).with_name("hermite6_calibration.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
