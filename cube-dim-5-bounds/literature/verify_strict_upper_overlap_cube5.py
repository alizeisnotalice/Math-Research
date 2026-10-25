#!/usr/bin/env python3
"""Exact audit of the 32-fold l-infinity Besicovitch configuration in R^5."""

from fractions import Fraction
from itertools import product
import json
from pathlib import Path


DIMENSION = 5
RADIUS = Fraction(1)
OFFSET = Fraction(3, 4)
SIGNS = tuple(product((-1, 1), repeat=DIMENSION))
CENTERS = tuple(tuple(OFFSET * sign for sign in signs) for signs in SIGNS)


def linf_distance(left, right):
    return max(abs(x - y) for x, y in zip(left, right))


pair_distances = [
    linf_distance(CENTERS[i], CENTERS[j])
    for i in range(len(CENTERS))
    for j in range(i + 1, len(CENTERS))
]
assert len(CENTERS) == 2**DIMENSION
assert all(linf_distance(center, (Fraction(0),) * DIMENSION) < RADIUS for center in CENTERS)
assert min(pair_distances) == 2 * OFFSET > RADIUS

# Intersection over both signs in every coordinate: [-1/4,1/4]^5.
intersection_half_side = RADIUS - OFFSET
intersection_volume = (2 * intersection_half_side) ** DIMENSION
assert intersection_half_side == Fraction(1, 4)
assert intersection_volume == Fraction(1, 32)

certificate = {
    "dimension": DIMENSION,
    "number_of_cubes": len(CENTERS),
    "radius": str(RADIUS),
    "center_coordinate_magnitudes": str(OFFSET),
    "common_point": ["0"] * DIMENSION,
    "center_distance_to_common_point": str(OFFSET),
    "minimum_pairwise_linf_center_distance": str(min(pair_distances)),
    "besicovitch_condition": "each pairwise center distance is > radius, so no center lies in another equal-radius closed cube",
    "common_open_overlap": "(-1/4,1/4)^5",
    "common_overlap_volume": str(intersection_volume),
    "conclusion": "the multiplicity 32 occurs on a set of positive Lebesgue measure",
}

output = Path(__file__).with_name("strict_upper_overlap_cube5_certificate.json")
output.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
print(json.dumps(certificate, indent=2, sort_keys=True))
