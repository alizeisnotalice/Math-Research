#!/usr/bin/env python3
"""Exact Laurent cancellation for selected non-S2 body-frame Green points."""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

from flint import arb, ctx

HERE = Path(__file__).resolve().parent
FAN = HERE.parent / "fan" / "exact_cells_results.json"
GROUP = 2**6 * math.factorial(6)
ZERO = (Fraction(), Fraction(), Fraction(), Fraction())  # pi, e^-1, e^-3, e^-5


def parse_fraction(text):
    return Fraction(text)


def add(left, right):
    return tuple(x + y for x, y in zip(left, right, strict=True))


def subtract(left, right):
    return tuple(x - y for x, y in zip(left, right, strict=True))


def divide(vector, scalar):
    return tuple(x / scalar for x in vector)


def cones():
    payload = json.loads(FAN.read_text())
    result = []
    for cell in payload["cells"]:
        rays = [tuple(parse_fraction(x) for x in ray) for ray in cell["rays_a0_equals_one"]]
        rates = [parse_fraction(x) for x in cell["rates"]]
        alpha = parse_fraction(cell["determinant"]) / math.prod(rates)
        result.append((rays, rates, alpha))
    return result, parse_fraction(payload["full_exp_minus_symbol_integral"])


CONES, FULL_PEAK_INTEGRAL = cones()


def orbit(point):
    points = set()
    for permutation in itertools.permutations(range(6)):
        reordered = tuple(point[index] for index in permutation)
        for signs in itertools.product((-1, 1), repeat=6):
            points.add(tuple(Fraction(sign) * x for sign, x in zip(signs, reordered, strict=True)))
    return sorted(points)


def aggregate(point):
    points = orbit(point)
    assert GROUP % len(points) == 0
    stabilizer = GROUP // len(points)
    coefficients = defaultdict(Fraction)
    for transformed in points:
        for rays, rates, alpha in CONES:
            nodes = tuple(sorted(sum((x * r for x, r in zip(transformed, ray, strict=True)), Fraction()) / q for ray, q in zip(rays, rates, strict=True)))
            coefficients[nodes] += stabilizer * alpha
    assert sum(coefficients.values(), Fraction()) == FULL_PEAK_INTEGRAL
    return coefficients, len(points), stabilizer


def jet(node, derivative_order):
    """F^(m)(node)/m! Laurent coefficients in pi,e^-1,e^-3,e^-5."""
    if node:
        return (Fraction(1 if node > 0 else -1, 2), Fraction(), Fraction(), Fraction()) if derivative_order == 0 else ZERO
    if derivative_order == 1:
        return (Fraction(), Fraction(1), Fraction(), Fraction())
    if derivative_order == 3:
        return (Fraction(), Fraction(), Fraction(-1, 3), Fraction())
    if derivative_order == 5:
        return (Fraction(), Fraction(), Fraction(), Fraction(1, 5))
    return ZERO


def dd5_laurent(nodes):
    table = [[ZERO for _ in range(6)] for _ in range(6)]
    for i, node in enumerate(nodes):
        table[i][0] = jet(node, 0)
    for order in range(1, 6):
        for i in range(order, 6):
            if nodes[i] == nodes[i - order]:
                table[i][order] = jet(nodes[i], order)
            else:
                table[i][order] = divide(subtract(table[i][order - 1], table[i - 1][order - 1]), nodes[i] - nodes[i - order])
    return table[-1][-1]


def point(k, a):
    return (Fraction(1),) + (a,) * k + (Fraction(),) * (5 - k)


def string(value):
    return f"{value.numerator}/{value.denominator}" if value.denominator != 1 else str(value.numerator)


def source_hash():
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def evaluate(k, a):
    x = point(k, a)
    coefficients, orbit_size, stabilizer = aggregate(x)
    total = ZERO
    zero_node_tuples = 0
    for nodes, coefficient in coefficients.items():
        if Fraction() in nodes:
            zero_node_tuples += 1
        total = add(total, tuple(coefficient * value for value in dd5_laurent(nodes)))
    pi_coefficient, e1, e3, e5 = total
    if any((e1, e3, e5)):
        return {"status": "FAIL_PRINCIPAL_PART", "point": [string(v) for v in x], "orbit_size": orbit_size, "stabilizer": stabilizer, "merged_exact_node_tuples": len(coefficients), "zero_node_tuples": zero_node_tuples, "coefficients": {"pi": string(pi_coefficient), "epsilon^-1": string(e1), "epsilon^-3": string(e3), "epsilon^-5": string(e5)}}
    q_over_pi5 = pi_coefficient / 64
    ctx.prec = 384
    value = (arb(q_over_pi5.numerator) / q_over_pi5.denominator) / arb.pi()**5
    return {"status": "PASS_EXACT_GREEN_POINT", "point": [string(v) for v in x], "orbit_size": orbit_size, "stabilizer": stabilizer, "raw_group_terms": len(CONES) * GROUP, "merged_exact_node_tuples": len(coefficients), "zero_node_tuples": zero_node_tuples, "coefficients": {"pi": string(pi_coefficient), "epsilon^-1": "0", "epsilon^-3": "0", "epsilon^-5": "0"}, "G_exact": f"({string(q_over_pi5)})/pi^5", "G_arb_384": value.str(80, more=True)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--two", action="store_true", help="also evaluate k=4,a=1/2 after the primary point passes")
    args = parser.parse_args()
    rows = [evaluate(4, Fraction(3, 4))]
    if rows[0]["status"] == "PASS_EXACT_GREEN_POINT" and args.two:
        rows.append(evaluate(4, Fraction(1, 2)))
    payload = {"scope": "finite exact Laurent cancellation at selected non-S2 points; no cube minimum, face integral, or C6 claim", "fan_input_sha256": hashlib.sha256(FAN.read_bytes()).hexdigest(), "full_peak_integral": string(FULL_PEAK_INTEGRAL), "rows": rows, "source_sha256": source_hash()}
    (HERE / "laurent_points_results.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
