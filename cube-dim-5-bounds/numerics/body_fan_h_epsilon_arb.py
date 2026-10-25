#!/usr/bin/env python3
"""Arb evaluation of the 5D truncated Green potential from Hermite DDs."""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from collections import Counter
from fractions import Fraction
from pathlib import Path

from flint import acb, arb, ctx, fmpq

from body_fan_arb_anchor import exact_terms


def ball(value: Fraction | int) -> arb:
    return arb(value) if isinstance(value, int) else arb(fmpq(value.numerator, value.denominator))


def f_value(node: Fraction, epsilon: Fraction) -> arb:
    return ball(node * node + epsilon * epsilon).log() / 2


def f_derivative_over_factorial(node: Fraction, epsilon: Fraction, order: int) -> arb:
    if order == 0:
        return f_value(node, epsilon)
    z_minus, z_plus = acb(ball(node), -ball(epsilon)), acb(ball(node), ball(epsilon))
    return (ball(Fraction((-1) ** (order - 1), 2 * order))
            * (z_minus ** (-order) + z_plus ** (-order))).real


def hermite_dd(nodes: tuple[Fraction, ...], epsilon: Fraction) -> arb:
    ordered, table = tuple(sorted(nodes)), []
    for right, node in enumerate(ordered):
        row = [f_value(node, epsilon)]
        for order in range(1, right + 1):
            left = ordered[right - order]
            row.append(f_derivative_over_factorial(node, epsilon, order) if node == left
                       else (row[order - 1] - table[right - 1][order - 1]) / ball(node - left))
        table.append(row)
    return table[-1][-1]


def truncated_green(point: tuple[Fraction, ...], epsilon: Fraction, terms) -> tuple[arb, Counter]:
    total, multiplicities = arb(0), Counter()
    for rays, determinant, rates in terms:
        nodes = tuple(sum(Fraction(coordinate) * value for coordinate, value in zip(ray, point)) / rate
                      for ray, rate in zip(rays, rates))
        multiplicities[tuple(sorted(Counter(nodes).values(), reverse=True))] += 1
        alpha = Fraction(determinant)
        for rate in rates:
            alpha /= rate
        total -= ball(alpha) * hermite_dd(nodes, epsilon)
    return total / (2 * acb.pi()) ** 5, multiplicities


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    ctx.prec = 256
    point, epsilons = ((Fraction(1), Fraction(1, 7), Fraction(2, 7), Fraction(3, 7), Fraction(4, 7)),
                       (Fraction(1, 32), Fraction(1, 64), Fraction(1, 128)))
    started, terms, values, pattern_counts = time.monotonic(), tuple(exact_terms()), {}, None
    for epsilon in epsilons:
        value, pattern_counts = truncated_green(point, epsilon, terms)
        values[str(epsilon)] = str(value)
    result = {"status": "arb_truncated_green_anchor_evaluation", "precision_bits": 256,
              "anchor": [str(value) for value in point], "epsilons": [str(value) for value in epsilons],
              "terms": len(terms), "formula": "H_epsilon=-(2*pi)^-5 sum_C alpha_C DD4[f_epsilon](b_C1,...,b_C5)",
              "values": values, "hermite_multiplicity_patterns": {str(key): value for key, value in pattern_counts.items()},
              "elapsed_seconds": time.monotonic() - started,
              "scope": "Finite exact-cone / Arb evaluation at one generic point only; it is not a global Green or C_5 certificate.",
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "anchor_source_sha256": hashlib.sha256(Path(__file__).with_name("body_fan_arb_anchor.py").read_bytes()).hexdigest()}
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
