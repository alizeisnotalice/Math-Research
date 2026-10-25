#!/usr/bin/env python3
"""Exact sorted-chamber cells for the six-dimensional 32-body frame.

The chamber section is a0=1>=a1>=...>=a5>=0.  Only the twelve body
functionals that strictly change sign on this section are cut.  All candidate
vertices are obtained by exact intersections of five of the 18 boundary
hyperplanes (six chamber facets plus twelve crossing hyperplanes).  This lets
the script classify every one of the 2^12 sign choices without a floating LP.

This file establishes cells and rational interior witnesses.  It deliberately
does not yet assert a simplex triangulation or a Fourier representation.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DIM = 5


def solve(rows: list[list[Fraction]], rhs: list[Fraction]):
    """Exact Gauss elimination for a square system, returning None if singular."""
    a = [row[:] + [b] for row, b in zip(rows, rhs, strict=True)]
    for col in range(DIM):
        pivot = next((r for r in range(col, DIM) if a[r][col]), None)
        if pivot is None:
            return None
        a[col], a[pivot] = a[pivot], a[col]
        p = a[col][col]
        a[col] = [value / p for value in a[col]]
        for row in range(DIM):
            if row == col:
                continue
            factor = a[row][col]
            if factor:
                a[row] = [u - factor * v for u, v in zip(a[row], a[col], strict=True)]
    return tuple(a[row][-1] for row in range(DIM))


def affine_value(affine, point):
    constant, coefficients = affine
    return constant + sum((coefficient * value for coefficient, value in zip(coefficients, point, strict=True)), Fraction())


def body_forms():
    # The 32 forms are a0+eps1*a1+...+eps5*a5 with a0=1.
    all_forms = [(Fraction(1), tuple(Fraction(sign) for sign in signs)) for signs in itertools.product((-1, 1), repeat=5)]
    vertices = [tuple(Fraction(1 if i < k else 0) for i in range(DIM)) for k in range(DIM + 1)]
    return [form for form in all_forms if min(affine_value(form, vertex) for vertex in vertices) < 0 < max(affine_value(form, vertex) for vertex in vertices)]


CHAMBER = (
    ("a1=1", (Fraction(1), (Fraction(-1), 0, 0, 0, 0))),
    ("a1=a2", (Fraction(0), (Fraction(1), Fraction(-1), 0, 0, 0))),
    ("a2=a3", (Fraction(0), (0, Fraction(1), Fraction(-1), 0, 0))),
    ("a3=a4", (Fraction(0), (0, 0, Fraction(1), Fraction(-1), 0))),
    ("a4=a5", (Fraction(0), (0, 0, 0, Fraction(1), Fraction(-1)))),
    ("a5=0", (Fraction(0), (0, 0, 0, 0, Fraction(1)))),
)


def chamber_ok(point):
    return all(affine_value(affine, point) >= 0 for _, affine in CHAMBER)


def encode_fraction(value):
    return f"{value.numerator}/{value.denominator}" if value.denominator != 1 else str(value.numerator)


def point_json(point):
    return [encode_fraction(value) for value in point]


def determinant(matrix):
    """Exact determinant of a square Fraction matrix."""
    a = [list(row) for row in matrix]
    sign, result = 1, Fraction(1)
    for column in range(len(a)):
        pivot = next((row for row in range(column, len(a)) if a[row][column]), None)
        if pivot is None:
            return Fraction()
        if pivot != column:
            a[column], a[pivot] = a[pivot], a[column]
            sign *= -1
        p = a[column][column]
        result *= p
        for row in range(column + 1, len(a)):
            factor = a[row][column] / p
            for j in range(column, len(a)):
                a[row][j] -= factor * a[column][j]
    return sign * result


def all_body_forms():
    return [(Fraction(1), tuple(Fraction(sign) for sign in signs)) for signs in itertools.product((-1, 1), repeat=5)]


def rate(ray):
    return sum(abs(ray[0] + sum(sign * coordinate for sign, coordinate in zip(signs, ray[1:], strict=True))) for signs in itertools.product((-1, 1), repeat=5))


def source_hash():
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def main():
    crossing = body_forms()
    assert len(crossing) == 12
    boundaries = list(CHAMBER) + [(f"body_{i}", affine) for i, affine in enumerate(crossing)]
    vertices = set()
    singular_systems = 0
    for chosen in itertools.combinations(range(len(boundaries)), DIM):
        equations = [boundaries[index][1] for index in chosen]
        point = solve([list(affine[1]) for affine in equations], [-affine[0] for affine in equations])
        if point is None:
            singular_systems += 1
        elif chamber_ok(point):
            vertices.add(point)
    vertices = sorted(vertices)

    # A vertex is compatible with a sign word when each nonzero body form has
    # that sign.  Iterate its zero-bit completions, avoiding 4096*|V| scans.
    compatible: list[list[int]] = [[] for _ in range(1 << len(crossing))]
    for index, vertex in enumerate(vertices):
        positive = 0
        free = 0
        for bit, form in enumerate(crossing):
            value = affine_value(form, vertex)
            if value > 0:
                positive |= 1 << bit
            elif value == 0:
                free |= 1 << bit
        subset = free
        while True:
            compatible[positive | subset].append(index)
            if subset == 0:
                break
            subset = (subset - 1) & free

    cells = []
    for word, indices in enumerate(compatible):
        if not indices:
            continue
        denominator = Fraction(len(indices))
        center = tuple(sum((vertices[index][j] for index in indices), Fraction()) / denominator for j in range(DIM))
        if not chamber_ok(center):
            continue
        # The average of every vertex of a bounded closure is in its relative
        # interior.  Strictness for every chamber/body inequality exactly tests
        # whether this closure has a full-dimensional open sign cell.
        if not all(affine_value(affine, center) > 0 for _, affine in CHAMBER):
            continue
        if not all((1 if word & (1 << bit) else -1) * affine_value(form, center) > 0 for bit, form in enumerate(crossing)):
            continue
        assert len(indices) == DIM + 1, "A non-simplex cell requires recursive exact triangulation."
        rays = [(Fraction(1),) + vertices[index] for index in indices]
        det = abs(determinant([[ray[row] for ray in rays] for row in range(DIM + 1)]))
        assert det > 0
        all_signs_constant = all(
            min(affine_value(form, vertex) for vertex in (vertices[index] for index in indices)) >= 0
            or max(affine_value(form, vertex) for vertex in (vertices[index] for index in indices)) <= 0
            for form in all_body_forms()
        )
        assert all_signs_constant
        rates = [rate(ray) for ray in rays]
        cells.append({"word": word, "signs": [1 if word & (1 << bit) else -1 for bit in range(len(crossing))], "vertex_indices": indices, "interior": center, "rays": rays, "determinant": det, "rates": rates})

    section_volume = sum((cell["determinant"] / math.factorial(DIM) for cell in cells), Fraction())
    chamber_volume = Fraction(1, math.factorial(DIM))
    assert section_volume == chamber_volume
    # On a simplicial cone, int exp(-psi) equals det(R)/prod q_j.  The signed
    # permutation group has 2^6*6! chamber images with disjoint interiors.
    base_exp_integral = sum((cell["determinant"] / math.prod(cell["rates"]) for cell in cells), Fraction())
    group_size = 2**6 * math.factorial(6)
    full_exp_integral = group_size * base_exp_integral
    spectral_ball_volume = full_exp_integral / math.factorial(6)

    output = {
        "scope": "exact rational chamber-cell enumeration; not yet a triangulation or spectral certificate",
        "section": "a0=1>=a1>=a2>=a3>=a4>=a5>=0",
        "chamber_facets": [name for name, _ in CHAMBER],
        "crossing_forms": [{"constant": encode_fraction(form[0]), "coefficients": [encode_fraction(value) for value in form[1]]} for form in crossing],
        "boundary_hyperplanes": len(boundaries),
        "five_hyperplane_systems": math_combination(len(boundaries), DIM),
        "singular_systems": singular_systems,
        "unique_chamber_arrangement_vertices": len(vertices),
        "vertices": [point_json(point) for point in vertices],
        "all_sign_words_examined": 1 << len(crossing),
        "full_dimensional_cells_exact": len(cells),
        "cells": [{"word": cell["word"], "signs": cell["signs"], "vertex_indices": cell["vertex_indices"], "strict_interior": point_json(cell["interior"]), "rays_a0_equals_one": [point_json(ray) for ray in cell["rays"]], "determinant": encode_fraction(cell["determinant"]), "rates": [encode_fraction(value) for value in cell["rates"]], "all_32_body_form_signs_constant": True} for cell in cells],
        "coverage_reason": "Every non-boundary point has one of 2^12 sign words.  A nonempty full-dimensional bounded closure has a vertex formed by five independent listed boundaries.  Every such vertex was enumerated exactly, and strict rational barycenters classify the open cells.",
        "simplex_triangulation": "Each of the 21 cells has exactly 6 affinely independent rational vertices, so the listed cell is itself its unique one-simplex triangulation.",
        "section_volume_sum": encode_fraction(section_volume),
        "sorted_chamber_section_volume": encode_fraction(chamber_volume),
        "base_exp_minus_symbol_integral": encode_fraction(base_exp_integral),
        "signed_permutation_group_size": group_size,
        "full_exp_minus_symbol_integral": encode_fraction(full_exp_integral),
        "spectral_unit_ball_volume": encode_fraction(spectral_ball_volume),
        "heat_kernel_peak": f"({encode_fraction(full_exp_integral)})/(2*pi)^6",
        "next_step": "use these compressed exact cones to build and independently audit Fourier/Green evaluations; global g and face-integral certification remain open",
        "source_sha256": source_hash(),
    }
    (ROOT / "exact_cells_results.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({key: output[key] for key in ("boundary_hyperplanes", "five_hyperplane_systems", "singular_systems", "unique_chamber_arrangement_vertices", "all_sign_words_examined", "full_dimensional_cells_exact", "source_sha256")}, indent=2))


def math_combination(n, k):
    from math import comb
    return comb(n, k)


if __name__ == "__main__":
    main()
