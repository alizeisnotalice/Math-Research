#!/usr/bin/env python3
"""Exact-node Arb evaluation of selected six-dimensional truncated potentials.

The compressed fan from ../fan is expanded only through the signed-permutation
orbit of each anchor.  It never materializes the full 967,680 cones.  Nodes
and all aggregation keys are Fraction objects; floating point is never used
to decide equality, multiplicity, or cancellation.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

from flint import arb, ctx

HERE = Path(__file__).resolve().parent
FAN = HERE.parent / "fan" / "exact_cells_results.json"
GROUP_SIZE = 2**6 * math.factorial(6)


def frac(text: str) -> Fraction:
    return Fraction(text)


def arbq(value: Fraction) -> arb:
    return arb(value.numerator) / value.denominator


def read_base_cones():
    source = json.loads(FAN.read_text())
    cones = []
    for cell in source["cells"]:
        rays = [tuple(frac(x) for x in ray) for ray in cell["rays_a0_equals_one"]]
        rates = [frac(x) for x in cell["rates"]]
        coefficient = frac(cell["determinant"]) / math.prod(rates)
        cones.append((rays, rates, coefficient))
    return cones, frac(source["full_exp_minus_symbol_integral"])


BASE_CONES, FULL_EXP_INTEGRAL = read_base_cones()


def signed_permutation_orbit(point: tuple[Fraction, ...]):
    values = set()
    for permutation in itertools.permutations(range(6)):
        moved = tuple(point[index] for index in permutation)
        for signs in itertools.product((-1, 1), repeat=6):
            values.add(tuple(Fraction(sign) * coordinate for sign, coordinate in zip(signs, moved, strict=True)))
    return sorted(values)


def aggregate_nodes(point: tuple[Fraction, ...]):
    orbit = signed_permutation_orbit(point)
    assert GROUP_SIZE % len(orbit) == 0
    stabilizer = GROUP_SIZE // len(orbit)
    coefficients: dict[tuple[Fraction, ...], Fraction] = defaultdict(Fraction)
    for transformed_point in orbit:
        for rays, rates, alpha in BASE_CONES:
            nodes = tuple(sorted(sum((coordinate * ray_coordinate for coordinate, ray_coordinate in zip(transformed_point, ray, strict=True)), Fraction()) / rate for ray, rate in zip(rays, rates, strict=True)))
            coefficients[nodes] += stabilizer * alpha
    return coefficients, len(orbit), stabilizer


@lru_cache(maxsize=None)
def atan_derivative(node: Fraction, epsilon: Fraction, order: int):
    x, e = arbq(node), arbq(epsilon)
    if order == 0:
        return (x / e).atan()
    d = e * e + x * x
    if order == 1:
        return e / d
    if order == 2:
        return -2 * e * x / d**2
    if order == 3:
        return 2 * e * (3 * x * x - e * e) / d**3
    if order == 4:
        return 24 * e * x * (e * e - x * x) / d**4
    if order == 5:
        return 24 * e * (e**4 - 10 * e * e * x * x + 5 * x**4) / d**5
    raise ValueError(order)


def hermite_dd(nodes: tuple[Fraction, ...], epsilon: Fraction):
    table = [[arb(0) for _ in range(6)] for _ in range(6)]
    for i, node in enumerate(nodes):
        table[i][0] = atan_derivative(node, epsilon, 0)
    for order in range(1, 6):
        for i in range(order, 6):
            if nodes[i] == nodes[i - order]:
                table[i][order] = atan_derivative(nodes[i], epsilon, order) / math.factorial(order)
            else:
                table[i][order] = (table[i][order - 1] - table[i - 1][order - 1]) / arbq(nodes[i] - nodes[i - order])
    return table[-1][-1]


def h_epsilon(point, epsilon):
    coefficients, orbit_size, stabilizer = aggregate_nodes(point)
    total = arb(0)
    for nodes, coefficient in coefficients.items():
        total += arbq(coefficient) * hermite_dd(nodes, epsilon)
    answer = total / (2 * arb.pi())**6
    return answer, {"orbit_size": orbit_size, "stabilizer": stabilizer, "raw_group_terms": len(BASE_CONES) * GROUP_SIZE, "merged_exact_node_tuples": len(coefficients)}


def anchor(k: int, a: Fraction):
    return (Fraction(1),) + (a,) * k + (Fraction(),) * (5 - k)


def render(value: arb):
    return value.str(70, more=True)


def calibrate_origin(epsilon: Fraction):
    actual, metadata = h_epsilon((Fraction(),) * 6, epsilon)
    expected = arbq(FULL_EXP_INTEGRAL) / (2 * arb.pi())**6 / (5 * arbq(epsilon)**5)
    difference = actual - expected
    if 0 not in difference:
        raise AssertionError("origin calibration interval excludes zero")
    return {"epsilon": str(epsilon), "actual": render(actual), "expected": render(expected), "difference": render(difference), **metadata}


def source_hash():
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", choices=("pilot", "all"), default="pilot")
    parser.add_argument("--precision", type=int, default=384)
    args = parser.parse_args()
    ctx.prec = args.precision
    epsilons = (Fraction(1, 1000), Fraction(1, 2000))
    cases = [(1, Fraction(1, 4))] if args.phase == "pilot" else [(k, a) for k in range(1, 5) for a in (Fraction(1, 4), Fraction(1, 2), Fraction(3, 4))]
    payload = {
        "status": "arb_truncated_anchor_evaluation",
        "scope": "finite anchor values only; no g lower bound, face integral, or C6 upper bound",
        "formula": "H_epsilon=(2*pi)^-6 sum alpha DD5 atan(b/epsilon)",
        "precision_bits": args.precision,
        "phase": args.phase,
        "fan_input_sha256": hashlib.sha256(FAN.read_bytes()).hexdigest(),
        "fan_base_cones": len(BASE_CONES),
        "full_signed_permutation_group_size": GROUP_SIZE,
        "origin_calibration": [calibrate_origin(epsilon) for epsilon in epsilons],
        "anchors": [],
        "source_sha256": source_hash(),
    }
    for k, a in cases:
        point = anchor(k, a)
        row = {"k": k, "a": str(a), "point": [str(value) for value in point], "values": []}
        for epsilon in epsilons:
            value, metadata = h_epsilon(point, epsilon)
            row["values"].append({"epsilon": str(epsilon), "interval": render(value), **metadata})
        payload["anchors"].append(row)
    destination = HERE / ("truncated_anchor_pilot_results.json" if args.phase == "pilot" else "truncated_anchor_results.json")
    destination.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
