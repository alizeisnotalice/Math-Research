#!/usr/bin/env python3
"""Generic-point epsilon=0 Green evaluations for the 5D body Fourier fan."""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from collections import Counter
from fractions import Fraction
from itertools import combinations
from pathlib import Path

import sympy as sp
from flint import acb, arb, ctx

from body_fan_arb_anchor import exact_terms
from body_fan_geometry import directions
from body_fan_h_epsilon_arb import ball, hermite_dd, truncated_green


DENOMINATOR = 1009
TAILS = ((936, 673, 68, 917), (201, 693, 737, 586), (901, 455, 74, 837), (926, 419, 373, 965), (148, 785, 707, 948), (139, 53, 348, 452), (399, 641, 776, 78), (931, 367, 597, 638))


def f0_dd(nodes: tuple[Fraction, ...]) -> arb:
    ordered = tuple(sorted(nodes))
    table: list[list[arb]] = []
    for right, node in enumerate(ordered):
        row = [ball(abs(node)).log()]
        for order in range(1, right + 1):
            left = ordered[right - order]
            if node == left:
                coefficient = Fraction((-1) ** (order - 1), order) / (node ** order)
                row.append(ball(coefficient))
            else:
                row.append((row[order - 1] - table[right - 1][order - 1]) / ball(node - left))
        table.append(row)
    return table[-1][-1]


def point_from_tail(tail: tuple[int, ...]) -> tuple[Fraction, ...]:
    return (Fraction(1),) + tuple(Fraction(value, DENOMINATOR) for value in tail)


def precheck(point: tuple[Fraction, ...], terms) -> dict:
    zero_terms = 0
    for rays, _, _ in terms:
        if any(sum(Fraction(coordinate) * value for coordinate, value in zip(ray, point)) == 0 for ray in rays):
            zero_terms += 1
    s2 = False
    for first, second in combinations(directions(), 2):
        if sp.Matrix.hstack(sp.Matrix(first), sp.Matrix(second), sp.Matrix(point)).rank() <= 2:
            s2 = True
            break
    return {"zero_projection_cones": zero_terms, "belongs_to_body_S2": s2}


def green(point: tuple[Fraction, ...], terms) -> tuple[arb, Counter]:
    total = arb(0)
    patterns: Counter = Counter()
    for rays, determinant, rates in terms:
        nodes = tuple(sum(Fraction(coordinate) * value for coordinate, value in zip(ray, point)) / rate for ray, rate in zip(rays, rates))
        patterns[tuple(sorted(Counter(nodes).values(), reverse=True))] += 1
        alpha = Fraction(determinant, 1)
        for rate in rates:
            alpha /= rate
        total -= ball(alpha) * f0_dd(nodes)
    return total / (2 * acb.pi()) ** 5, patterns


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--limit", type=int, default=8)
    args = parser.parse_args()
    ctx.prec = 256
    started = time.monotonic()
    terms = tuple(exact_terms())
    if args.output.exists():
        output = json.loads(args.output.read_text(encoding="utf-8"))
    else:
        output = {"status": "in_progress", "precision_bits": 256, "denominator": DENOMINATOR, "epsilon_for_loss": "1/128", "terms": len(terms), "points": [], "formula": "G=-(2*pi)^-5 sum_C alpha_C DD4[log|.|](b_C1,...,b_C5)", "scope": "Generic point values only; no global lower bound or face integral certificate."}
    completed = {tuple(item["tail_numerators"]) for item in output["points"]}
    for tail in TAILS:
        if tail in completed:
            continue
        if len(output["points"]) >= args.limit or time.monotonic() - started > 120:
            break
        point = point_from_tail(tail)
        checks = precheck(point, terms)
        assert checks["zero_projection_cones"] == 0 and not checks["belongs_to_body_S2"]
        green_value, patterns = green(point, terms)
        truncated_value, _ = truncated_green(point, Fraction(1, 128), terms)
        output["points"].append({"tail_numerators": list(tail), "point": [str(value) for value in point], "precheck": checks, "G": str(green_value), "H_1_128": str(truncated_value), "G_minus_H_1_128": str(green_value - truncated_value), "multiplicity_patterns": {str(key): value for key, value in patterns.items()}, "elapsed_seconds_at_save": time.monotonic() - started})
        output["source_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    output["status"] = "complete" if len(output["points"]) == len(TAILS) else "checkpoint"
    output["elapsed_seconds"] = time.monotonic() - started
    output["source_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
