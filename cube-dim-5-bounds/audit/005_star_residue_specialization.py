#!/usr/bin/env python3
"""Exact generic-specialization check of n=5 Fourier-fan star residues.

This is a check, not the proof of cancellation.  It evaluates the rational
star coefficient at x=(1,B,...,B^4), B=10000, after summing all 26880 signed
permutation images of the seven base cones.  At a simple t-pole b=x.r/q,
the local cone contribution is

  (|det R|/prod q_j) / prod_{k != j}(b_j-b_k).

The factor (2*pi)^-5 is intentionally removed.  B is large enough to
separate the finite normalized-ray set in this fan; the program asserts both
within-cone simplicity and global normalized-ray separation at this x.
"""
from __future__ import annotations

from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path


BASE_CONES = (
    ((1, 1, 0, 0, 0), (1, 0, 0, 0, 0), (4, 1, 1, 1, 1), (3, 1, 1, 1, 0), (2, 1, 1, 0, 0)),
    ((1, 1, 0, 0, 0), (4, 1, 1, 1, 1), (2, 1, 1, 1, 1), (3, 1, 1, 1, 0), (2, 1, 1, 0, 0)),
    ((1, 1, 1, 1, 0), (1, 1, 0, 0, 0), (2, 1, 1, 1, 1), (3, 1, 1, 1, 0), (2, 1, 1, 0, 0)),
    ((1, 1, 1, 1, 0), (2, 2, 2, 1, 1), (1, 1, 0, 0, 0), (2, 1, 1, 1, 1), (2, 1, 1, 0, 0)),
    ((1, 1, 1, 1, 0), (2, 2, 2, 1, 1), (1, 1, 0, 0, 0), (3, 3, 2, 2, 2), (2, 1, 1, 1, 1)),
    ((1, 1, 1, 1, 1), (1, 1, 1, 1, 0), (2, 2, 2, 1, 1), (3, 3, 2, 2, 2), (2, 1, 1, 1, 1)),
    ((1, 1, 1, 1, 0), (1, 1, 1, 0, 0), (2, 2, 2, 1, 1), (1, 1, 0, 0, 0), (2, 1, 1, 0, 0)),
)
BASE_DETERMINANTS = (1, 2, 2, 2, 2, 1, 1)

# Exact classification from 005_vertex_rays_exact.py, invariant under signed
# coordinate permutations.  Keys are primitive rays with first nonzero entry
# positive, as every base ray already is.
TRUE_VERTEX_TYPES = {
    (1, 1, 0, 0, 0),
    (1, 1, 1, 1, 0),
    (2, 1, 1, 1, 1),
}


def directions() -> tuple[tuple[int, ...], ...]:
    return tuple((1,) + signs for signs in itertools.product((-1, 1), repeat=4))


def dot(left: tuple[int, ...], right: tuple[int, ...]) -> int:
    return sum(a * b for a, b in zip(left, right))


def q(ray: tuple[int, ...]) -> int:
    return sum(abs(dot(direction, ray)) for direction in directions())


def transformed(ray: tuple[int, ...], permutation: tuple[int, ...], signs: tuple[int, ...]) -> tuple[int, ...]:
    # Same convention as the producing body_fan_check.py: output i takes
    # sign_i times input permutation_i.
    return tuple(signs[i] * ray[permutation[i]] for i in range(5))


def main() -> None:
    B = 10000
    x = (1, B, B**2, B**3, B**4)
    residue_by_p: dict[tuple[Fraction, ...], Fraction] = {}
    kind_by_p: dict[tuple[Fraction, ...], bool] = {}
    b_by_p: dict[tuple[Fraction, ...], Fraction] = {}
    cone_count = 0
    contribution_count = 0
    for permutation in itertools.permutations(range(5)):
        for signs in itertools.product((-1, 1), repeat=5):
            for cone, determinant in zip(BASE_CONES, BASE_DETERMINANTS):
                rays = tuple(transformed(ray, permutation, signs) for ray in cone)
                rates = tuple(q(ray) for ray in rays)
                b = tuple(Fraction(dot(x, ray), rate) for ray, rate in zip(rays, rates))
                assert len(set(b)) == 5, (permutation, signs, cone, b)
                alpha = Fraction(determinant, 1)
                for rate in rates:
                    alpha /= rate
                cone_count += 1
                for j, ray in enumerate(rays):
                    p = tuple(Fraction(value, rates[j]) for value in ray)
                    coefficient = alpha
                    for k in range(5):
                        if k != j:
                            coefficient /= b[j] - b[k]
                    residue_by_p[p] = residue_by_p.get(p, Fraction()) + coefficient
                    expected_kind = cone[j] in TRUE_VERTEX_TYPES
                    if p in kind_by_p:
                        assert kind_by_p[p] == expected_kind
                    kind_by_p[p] = expected_kind
                    if p in b_by_p:
                        assert b_by_p[p] == b[j]
                    b_by_p[p] = b[j]
                    contribution_count += 1
    assert cone_count == 7 * 120 * 32
    assert len(set(b_by_p.values())) == len(b_by_p), "specialization failed to separate normalized rays"
    artificial = [p for p, true_vertex in kind_by_p.items() if not true_vertex]
    genuine = [p for p, true_vertex in kind_by_p.items() if true_vertex]
    zero_artificial = [p for p in artificial if residue_by_p[p] == 0]
    nonzero_genuine = [p for p in genuine if residue_by_p[p] != 0]
    assert len(zero_artificial) == len(artificial)
    result = {
        "status": "PASS_exact_single_generic_specialization",
        "x": list(x),
        "base_cones": 7,
        "full_cones": cone_count,
        "simple_pole_contributions": contribution_count,
        "distinct_normalized_rays": len(kind_by_p),
        "true_vertex_normalized_rays": len(genuine),
        "artificial_normalized_rays": len(artificial),
        "artificial_star_coefficients_exactly_zero": len(zero_artificial),
        "true_vertex_stars_nonzero_at_this_specialization": len(nonzero_genuine),
        "normalization": "reported residue excludes the common (2*pi)^-5 Fourier factor",
        "limitation": "This exact specialization is a consistency check only. Universal cancellation follows from the vertex-only Fourier-Laplace identity in 005_nonvertex_star_cancellation.md, not from one x value.",
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    target = Path(__file__).with_name("005_star_residue_specialization_results.json")
    target.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
