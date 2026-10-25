#!/usr/bin/env python3
"""256-bit Arb evaluation of the 26,880-cone heat-kernel formula at one anchor."""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import time
from fractions import Fraction
from pathlib import Path

from flint import acb, arb, ctx, fmpq

from body_fan_geometry import base_terms


def exact_terms():
    for permutation in itertools.permutations(range(5)):
        for signs in itertools.product((-1, 1), repeat=5):
            for cone, determinant, rates in base_terms():
                yield (tuple(tuple(signs[row] * ray[permutation[row]] for row in range(5))
                             for ray in cone), determinant, rates)


def q(value: Fraction) -> arb:
    return arb(fmpq(value.numerator, value.denominator))


def heat(point: tuple[Fraction, ...], time_value: Fraction) -> arb:
    total = acb(0)
    for rays, determinant, rates in exact_terms():
        term = acb(determinant)
        for ray, rate in zip(rays, rates):
            dot = sum(coordinate * value for coordinate, value in zip(ray, point))
            term /= acb(q(time_value * rate), -q(dot))
        total += term
    return (total / (2 * acb.pi()) ** 5).real


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    ctx.prec = 256
    point = (Fraction(1), Fraction(1, 7), Fraction(2, 7), Fraction(3, 7), Fraction(4, 7))
    started = time.monotonic()
    values = {label: str(heat(point, value)) for label, value in
              (("1/2", Fraction(1, 2)), ("1", Fraction(1)), ("2", Fraction(2)))}
    direct_half = heat(point, Fraction(1, 2))
    scaled_half = q(Fraction(1, 2) ** -5) * heat(tuple(value * 2 for value in point), Fraction(1))
    direct_double = heat(point, Fraction(2))
    scaled_double = q(Fraction(2) ** -5) * heat(tuple(value / 2 for value in point), Fraction(1))
    result = {"status": "arb_anchor_evaluation", "precision_bits": 256,
              "anchor": [str(value) for value in point], "terms": 26880, "values": values,
              "homogeneity_differences": {"t=1/2": str(direct_half - scaled_half),
                                           "t=2": str(direct_double - scaled_double)},
              "elapsed_seconds": time.monotonic() - started,
              "scope": "A finite Fourier representation evaluation only; no global bound or integral certificate.",
              "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "geometry_source_sha256": hashlib.sha256(Path(__file__).with_name("body_fan_geometry.py").read_bytes()).hexdigest()}
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
