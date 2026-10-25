#!/usr/bin/env python3
"""Independent exact audit of the seven-cone sorted body-frame fan."""
from fractions import Fraction
from itertools import permutations, product
from math import factorial, prod
from pathlib import Path
import hashlib
import json
import sympy as sp

CONES = (
    ((1,1,0,0,0),(1,0,0,0,0),(4,1,1,1,1),(3,1,1,1,0),(2,1,1,0,0)),
    ((1,1,0,0,0),(4,1,1,1,1),(2,1,1,1,1),(3,1,1,1,0),(2,1,1,0,0)),
    ((1,1,1,1,0),(1,1,0,0,0),(2,1,1,1,1),(3,1,1,1,0),(2,1,1,0,0)),
    ((1,1,1,1,0),(2,2,2,1,1),(1,1,0,0,0),(2,1,1,1,1),(2,1,1,0,0)),
    ((1,1,1,1,0),(2,2,2,1,1),(1,1,0,0,0),(3,3,2,2,2),(2,1,1,1,1)),
    ((1,1,1,1,1),(1,1,1,1,0),(2,2,2,1,1),(3,3,2,2,2),(2,1,1,1,1)),
    ((1,1,1,1,0),(1,1,1,0,0),(2,2,2,1,1),(1,1,0,0,0),(2,1,1,0,0)),
)
CHAMBER = ((1,0,0,0,0),(1,1,0,0,0),(1,1,1,0,0),(1,1,1,1,0),(1,1,1,1,1))


def det_columns(columns):
    return sum(
        (1 if sum(permutation[i] > permutation[j] for i in range(5) for j in range(i + 1, 5)) % 2 == 0 else -1)
        * prod(columns[j][permutation[j]] for j in range(5))
        for permutation in permutations(range(5))
    )


def functional(signs, ray):
    return ray[0] + sum(sign * coordinate for sign, coordinate in zip(signs, ray[1:]))


def normalized_section_volume(rays):
    return Fraction(abs(det_columns(rays)), prod(sum(ray) for ray in rays))


def rate(ray):
    return sum(abs(functional(signs, ray)) for signs in product((-1, 1), repeat=4))


def divided_difference(values, nodes):
    values = list(values)
    nodes = list(nodes)
    for order in range(1, len(nodes)):
        values = [(values[j + 1] - values[j]) / (nodes[j + order] - nodes[j]) for j in range(len(values) - 1)]
    return values[0]


def dd4_sign_checks():
    # For f_e(b)=log(e^2+b^2)/2 and five distinct real nodes,
    # d_e f_e[b0,...,b4] = Re product(e-i*bj)^-1.
    for epsilon in (Fraction(1), Fraction(7), Fraction(11, 3)):
        nodes = (Fraction(0), Fraction(1), Fraction(2), Fraction(3), Fraction(5))
        derivative_values = [epsilon / (epsilon * epsilon + node * node) for node in nodes]
        left = divided_difference(derivative_values, nodes)
        right = sp.re(prod(1 / (sp.Rational(epsilon.numerator, epsilon.denominator) - sp.I * sp.Rational(node.numerator, node.denominator)) for node in nodes))
        assert sp.Rational(left.numerator, left.denominator) == sp.simplify(right)
    # Hermite case: five coincident nodes.  The fourth b derivative divided by 4!
    # of d_e f equals the same real product.
    b, e = sp.symbols("b e", real=True, positive=True)
    hermite = sp.diff(e / (e**2 + b**2), b, 4) / factorial(4)
    target = sp.re((e - sp.I * b) ** -5)
    assert sp.simplify(hermite - target) == 0


def main():
    body = tuple(product((-1, 1), repeat=4))
    nontrivial = tuple(signs for signs in body if signs != (1, 1, 1, 1))
    signatures, determinants, rates = [], [], []
    for cone in CONES:
        assert all(all(ray[i] >= ray[i + 1] >= 0 for i in range(4)) for ray in cone)
        determinant = abs(det_columns(cone))
        assert determinant > 0
        determinants.append(determinant)
        rates.append(tuple(rate(ray) for ray in cone))
        signature = []
        for signs in nontrivial:
            values = [functional(signs, ray) for ray in cone]
            assert min(values) >= 0 or max(values) <= 0
            aggregate = functional(signs, tuple(map(sum, zip(*cone))))
            assert aggregate != 0
            signature.append(1 if aggregate > 0 else -1)
        signatures.append(tuple(signature))
    assert len(set(signatures)) == 7
    cone_volume_sum = sum((normalized_section_volume(cone) for cone in CONES), Fraction(0))
    chamber_volume = normalized_section_volume(CHAMBER)
    assert cone_volume_sum == chamber_volume == Fraction(1, 120)

    # At x=0, signed-permutation images contribute equally.  There are 2^5*5!
    # images, and (2*pi)^-5 times this factor reduces to 120/pi^5.
    base_sum = sum((Fraction(det, prod(ray_rates)) for det, ray_rates in zip(determinants, rates)), Fraction(0))
    peak_numerator = 120 * base_sum
    assert peak_numerator == Fraction(35509, 764411904)
    dd4_sign_checks()
    report = {
        "status": "PASS independent exact seven-cone and DD4 audit",
        "base_determinants": determinants,
        "base_rates": [list(row) for row in rates],
        "distinct_open_sign_cells": len(set(signatures)),
        "normalized_cone_volume_sum": str(cone_volume_sum),
        "normalized_chamber_volume": str(chamber_volume),
        "full_cone_count": 7 * 2**5 * factorial(5),
        "peak": "p_1(0)=35509/(764411904*pi^5)",
        "peak_coefficient": str(peak_numerator),
        "dd4_formula": "H_epsilon=-(2*pi)^-5 sum_C |det R_C|/prod(q_Cj) f_epsilon[b_C1,...,b_C5]",
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    Path(__file__).with_name("004_body_fan_independent_results.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
