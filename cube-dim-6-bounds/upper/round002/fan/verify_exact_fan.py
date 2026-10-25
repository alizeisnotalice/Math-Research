#!/usr/bin/env python3
"""Independent exact audit of the compressed six-dimensional body fan.

It reads the saved certificate data but does not import the producer.  It
rebuilds all 18-boundary intersections, reclassifies every sign word, and
checks rays, weak signs, simplex volumes, and the peak integral.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RESULT = ROOT / "exact_cells_results.json"
N = 5


def fraction(value):
    return Fraction(value)


def solve_square(rows, rhs):
    data = [list(row) + [right] for row, right in zip(rows, rhs, strict=True)]
    for column in range(N):
        pivot = next((r for r in range(column, N) if data[r][column] != 0), None)
        if pivot is None:
            return None
        data[column], data[pivot] = data[pivot], data[column]
        lead = data[column][column]
        for j in range(column, N + 1):
            data[column][j] /= lead
        for r in range(N):
            if r != column:
                multiplier = data[r][column]
                for j in range(column, N + 1):
                    data[r][j] -= multiplier * data[column][j]
    return tuple(data[i][N] for i in range(N))


def value(affine, point):
    return affine[0] + sum((a * x for a, x in zip(affine[1], point, strict=True)), Fraction())


def det(matrix):
    rows = [list(row) for row in matrix]
    sign = 1
    answer = Fraction(1)
    for column in range(len(rows)):
        pivot = next((r for r in range(column, len(rows)) if rows[r][column]), None)
        if pivot is None:
            return Fraction()
        if pivot != column:
            rows[column], rows[pivot] = rows[pivot], rows[column]
            sign = -sign
        lead = rows[column][column]
        answer *= lead
        for r in range(column + 1, len(rows)):
            multiplier = rows[r][column] / lead
            for j in range(column, len(rows)):
                rows[r][j] -= multiplier * rows[column][j]
    return sign * answer


def main():
    certificate = json.loads(RESULT.read_text())
    forms = [(fraction(item["constant"]), tuple(fraction(x) for x in item["coefficients"])) for item in certificate["crossing_forms"]]
    chamber = (
        (Fraction(1), (Fraction(-1), 0, 0, 0, 0)),
        (Fraction(0), (Fraction(1), Fraction(-1), 0, 0, 0)),
        (Fraction(0), (0, Fraction(1), Fraction(-1), 0, 0)),
        (Fraction(0), (0, 0, Fraction(1), Fraction(-1), 0)),
        (Fraction(0), (0, 0, 0, Fraction(1), Fraction(-1))),
        (Fraction(0), (0, 0, 0, 0, Fraction(1))),
    )
    boundaries = chamber + tuple(forms)
    rebuilt = set()
    for subset in itertools.combinations(boundaries, N):
        candidate = solve_square([list(item[1]) for item in subset], [-item[0] for item in subset])
        if candidate is not None and all(value(item, candidate) >= 0 for item in chamber):
            rebuilt.add(candidate)
    stored = {tuple(fraction(x) for x in point) for point in certificate["vertices"]}
    assert rebuilt == stored
    vertices = sorted(rebuilt)
    matching = [[] for _ in range(1 << len(forms))]
    for index, vertex in enumerate(vertices):
        pos = free = 0
        for bit, form in enumerate(forms):
            signed = value(form, vertex)
            if signed > 0:
                pos |= 1 << bit
            elif signed == 0:
                free |= 1 << bit
        sub = free
        while True:
            matching[pos | sub].append(index)
            if sub == 0:
                break
            sub = (sub - 1) & free
    full_words = []
    for word, indices in enumerate(matching):
        if not indices:
            continue
        barycenter = tuple(sum((vertices[i][j] for i in indices), Fraction()) / len(indices) for j in range(N))
        if all(value(item, barycenter) > 0 for item in chamber) and all((1 if word & (1 << bit) else -1) * value(form, barycenter) > 0 for bit, form in enumerate(forms)):
            full_words.append(word)
    expected_words = [cell["word"] for cell in certificate["cells"]]
    assert full_words == expected_words
    total_section_volume = Fraction()
    base_integral = Fraction()
    all_forms = [(Fraction(1), tuple(Fraction(sign) for sign in signs)) for signs in itertools.product((-1, 1), repeat=5)]
    for cell in certificate["cells"]:
        rays = [tuple(fraction(x) for x in ray) for ray in cell["rays_a0_equals_one"]]
        matrix = [[ray[i] for ray in rays] for i in range(6)]
        determinant = abs(det(matrix))
        assert determinant == fraction(cell["determinant"])
        assert determinant > 0
        assert all(min(value(form, ray[1:]) for ray in rays) >= 0 or max(value(form, ray[1:]) for ray in rays) <= 0 for form in all_forms)
        rates = []
        for ray in rays:
            q = sum(abs(ray[0] + sum(sign * coordinate for sign, coordinate in zip(signs, ray[1:], strict=True))) for signs in itertools.product((-1, 1), repeat=5))
            rates.append(q)
        assert rates == [fraction(item) for item in cell["rates"]]
        total_section_volume += determinant / math.factorial(5)
        base_integral += determinant / math.prod(rates)
    assert total_section_volume == Fraction(1, 120)
    assert base_integral == fraction(certificate["base_exp_minus_symbol_integral"])
    full = base_integral * 2**6 * math.factorial(6)
    assert full == fraction(certificate["full_exp_minus_symbol_integral"])
    output = {
        "status": "PASS",
        "reconstructed_chamber_vertices": len(rebuilt),
        "all_sign_words_reclassified": 1 << len(forms),
        "full_cells": len(full_words),
        "simplex_section_volume_sum": str(total_section_volume),
        "base_exp_minus_symbol_integral": str(base_integral),
        "full_exp_minus_symbol_integral": str(full),
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "input_sha256": hashlib.sha256(RESULT.read_bytes()).hexdigest(),
        "scope": "finite exact fan audit only; no Green minimum, face integral, or weak-type upper bound",
    }
    (ROOT / "verify_exact_fan_results.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
