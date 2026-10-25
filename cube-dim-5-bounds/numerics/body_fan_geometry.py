#!/usr/bin/env python3
"""Exact common geometry for the five-dimensional body Fourier fan.

This is intentionally independent of the former numerical checker: it holds
only immutable integer fan data and an exact sorted-chamber coverage check.
Consumers must import this module, rather than importing a script which also
performs floating-point heat-kernel work.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from fractions import Fraction
from functools import lru_cache
from math import prod
from pathlib import Path

import sympy as sp


BASE_CONES = (
    ((1, 1, 0, 0, 0), (1, 0, 0, 0, 0), (4, 1, 1, 1, 1), (3, 1, 1, 1, 0), (2, 1, 1, 0, 0)),
    ((1, 1, 0, 0, 0), (4, 1, 1, 1, 1), (2, 1, 1, 1, 1), (3, 1, 1, 1, 0), (2, 1, 1, 0, 0)),
    ((1, 1, 1, 1, 0), (1, 1, 0, 0, 0), (2, 1, 1, 1, 1), (3, 1, 1, 1, 0), (2, 1, 1, 0, 0)),
    ((1, 1, 1, 1, 0), (2, 2, 2, 1, 1), (1, 1, 0, 0, 0), (2, 1, 1, 1, 1), (2, 1, 1, 0, 0)),
    ((1, 1, 1, 1, 0), (2, 2, 2, 1, 1), (1, 1, 0, 0, 0), (3, 3, 2, 2, 2), (2, 1, 1, 1, 1)),
    ((1, 1, 1, 1, 1), (1, 1, 1, 1, 0), (2, 2, 2, 1, 1), (3, 3, 2, 2, 2), (2, 1, 1, 1, 1)),
    ((1, 1, 1, 1, 0), (1, 1, 1, 0, 0), (2, 2, 2, 1, 1), (1, 1, 0, 0, 0), (2, 1, 1, 0, 0)),
)
CHAMBER = ((1, 0, 0, 0, 0), (1, 1, 0, 0, 0), (1, 1, 1, 0, 0),
           (1, 1, 1, 1, 0), (1, 1, 1, 1, 1))


def directions() -> tuple[tuple[int, ...], ...]:
    return tuple((1,) + signs for signs in itertools.product((-1, 1), repeat=4))


def determinant(rays: tuple[tuple[int, ...], ...]) -> int:
    return abs(int(sp.Matrix.hstack(*(sp.Matrix(ray) for ray in rays)).det()))


def direction_rate(ray: tuple[int, ...]) -> int:
    return sum(abs(sum(direction[k] * ray[k] for k in range(5))) for direction in directions())


@lru_cache(maxsize=1)
def base_terms():
    return tuple((cone, determinant(cone), tuple(direction_rate(ray) for ray in cone))
                 for cone in BASE_CONES)


def _form(signs, ray):
    return ray[0] + sum(sign * coordinate for sign, coordinate in zip(signs, ray[1:]))


def _normalized_volume(rays):
    return Fraction(determinant(rays), prod(sum(ray) for ray in rays))


def sorted_chamber_certificate():
    # These 15 hyperplanes induce the body-direction cells in the open chamber.
    forms = tuple(signs for signs in itertools.product((-1, 1), repeat=4)
                  if signs != (1, 1, 1, 1))
    signatures = []
    for cone in BASE_CONES:
        assert all(all(ray[i] >= ray[i + 1] >= 0 for i in range(4)) for ray in cone)
        assert determinant(cone) > 0
        signature = []
        for signs in forms:
            values = [_form(signs, ray) for ray in cone]
            assert min(values) >= 0 or max(values) <= 0
            interior_test = _form(signs, tuple(map(sum, zip(*cone))))
            assert interior_test != 0
            signature.append(1 if interior_test > 0 else -1)
        signatures.append(tuple(signature))
    total = sum((_normalized_volume(cone) for cone in BASE_CONES), Fraction(0))
    chamber = _normalized_volume(CHAMBER)
    assert len(set(signatures)) == len(BASE_CONES)
    assert total == chamber == Fraction(1, 120)
    return {
        "status": "exact_a.e._coverage_of_sorted_chamber",
        "sorted_chamber_cones": len(BASE_CONES),
        "distinct_open_sign_cells": len(set(signatures)),
        "normalized_volume_sum": str(total),
        "sorted_chamber_normalized_volume": str(chamber),
        "base_determinants": [item[1] for item in base_terms()],
        "base_rates": [list(item[2]) for item in base_terms()],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = sorted_chamber_certificate()
    output["source_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
