#!/usr/bin/env python3
"""Independent exact audit of the n=6 sorted-section body fan.

The 12 crossing forms, every chamber vertex, every compatible vertex set, and
every simplex are reconstructed from the 32 directions.  The producer's JSON
is opened only after this construction, to compare its stored cells against
the independently bound word-to-vertex map.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
from math import comb, factorial, prod
from pathlib import Path
import json


N = 5
SOURCE = Path(__file__).resolve().parents[2] / "upper/round002/fan"
CERTIFICATE = SOURCE / "exact_cells_results.json"
ROOT = Path(__file__).resolve().parent


def value(form, point):
    return form[0] + sum((a * x for a, x in zip(form[1], point, strict=True)), Fraction())


def solve(rows, rhs):
    a = [list(row) + [b] for row, b in zip(rows, rhs, strict=True)]
    for col in range(N):
        pivot = next((r for r in range(col, N) if a[r][col]), None)
        if pivot is None:
            return None
        a[col], a[pivot] = a[pivot], a[col]
        lead = a[col][col]
        a[col] = [z / lead for z in a[col]]
        for r in range(N):
            if r != col and a[r][col]:
                q = a[r][col]
                a[r] = [u - q * v for u, v in zip(a[r], a[col], strict=True)]
    return tuple(a[i][-1] for i in range(N))


def determinant(rows):
    a = [list(row) for row in rows]
    out, sign = Fraction(1), 1
    for col in range(len(a)):
        pivot = next((r for r in range(col, len(a)) if a[r][col]), None)
        if pivot is None:
            return Fraction()
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            sign *= -1
        lead = a[col][col]
        out *= lead
        for r in range(col + 1, len(a)):
            if a[r][col]:
                q = a[r][col] / lead
                for j in range(col, len(a)):
                    a[r][j] -= q * a[col][j]
    return sign * out


def enc(q):
    return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"


def point_json(x):
    return [enc(q) for q in x]


def main():
    all_forms = [(Fraction(1), tuple(Fraction(s) for s in signs))
                 for signs in product((-1, 1), repeat=5)]
    section_vertices = [tuple(Fraction(1 if j < k else 0) for j in range(N))
                        for k in range(N + 1)]
    crossing = [f for f in all_forms
                if min(value(f, x) for x in section_vertices) < 0 < max(value(f, x) for x in section_vertices)]
    assert len(crossing) == 12
    chamber = (
        (Fraction(1), (Fraction(-1), 0, 0, 0, 0)),
        (Fraction(0), (Fraction(1), Fraction(-1), 0, 0, 0)),
        (Fraction(0), (0, Fraction(1), Fraction(-1), 0, 0)),
        (Fraction(0), (0, 0, Fraction(1), Fraction(-1), 0)),
        (Fraction(0), (0, 0, 0, Fraction(1), Fraction(-1))),
        (Fraction(0), (0, 0, 0, 0, Fraction(1))),
    )
    boundaries = chamber + tuple(crossing)
    assert len(boundaries) == 18

    vertices, singular = set(), 0
    for chosen in combinations(boundaries, N):
        x = solve([f[1] for f in chosen], [-f[0] for f in chosen])
        if x is None:
            singular += 1
        elif all(value(f, x) >= 0 for f in chamber):
            vertices.add(x)
    vertices = sorted(vertices)
    assert len(vertices) == 18
    assert singular == 3270

    compatible = [[] for _ in range(1 << len(crossing))]
    for index, x in enumerate(vertices):
        positive = free = 0
        for bit, f in enumerate(crossing):
            q = value(f, x)
            if q > 0:
                positive |= 1 << bit
            elif q == 0:
                free |= 1 << bit
        sub = free
        while True:
            compatible[positive | sub].append(index)
            if sub == 0:
                break
            sub = (sub - 1) & free

    independent = {}
    for word, ids in enumerate(compatible):
        if not ids:
            continue
        interior = tuple(sum((vertices[i][j] for i in ids), Fraction()) / len(ids)
                         for j in range(N))
        strict_chamber = all(value(f, interior) > 0 for f in chamber)
        strict_word = all((1 if word & (1 << bit) else -1) * value(f, interior) > 0
                          for bit, f in enumerate(crossing))
        if not (strict_chamber and strict_word):
            continue
        assert len(ids) == 6
        rays = [(Fraction(1),) + vertices[i] for i in ids]
        det = abs(determinant([[ray[row] for ray in rays] for row in range(6)]))
        assert det > 0
        # Bind the word to every compatible vertex, including all weak chamber
        # and crossing inequalities.  Then independently check every one of
        # the original 32 body forms has one weak sign throughout the simplex.
        for x in (vertices[i] for i in ids):
            assert all(value(f, x) >= 0 for f in chamber)
            assert all((1 if word & (1 << bit) else -1) * value(f, x) >= 0
                       for bit, f in enumerate(crossing))
        assert all(min(value(f, x) for x in (vertices[i] for i in ids)) >= 0
                   or max(value(f, x) for x in (vertices[i] for i in ids)) <= 0
                   for f in all_forms)
        rates = [sum(abs(value(f, ray[1:])) for f in all_forms) for ray in rays]
        independent[word] = {
            "ids": ids, "interior": interior, "rays": rays,
            "determinant": det, "rates": rates,
        }
    assert len(independent) == 21
    volume = sum((c["determinant"] / factorial(5) for c in independent.values()), Fraction())
    assert volume == Fraction(1, 120)
    base = sum((c["determinant"] / prod(c["rates"]) for c in independent.values()), Fraction())
    full = base * (2**6) * factorial(6)
    peak = f"({enc(full)})/(2*pi)^6"

    # Only now read and bind the producer's stored vertex/cell data.
    stored = json.loads(CERTIFICATE.read_text())
    stored_crossing = [(Fraction(a["constant"]), tuple(Fraction(z) for z in a["coefficients"]))
                       for a in stored["crossing_forms"]]
    assert stored_crossing == crossing
    assert {tuple(Fraction(z) for z in x) for x in stored["vertices"]} == set(vertices)
    assert len(stored["cells"]) == len(independent)
    for cell in stored["cells"]:
        word = cell["word"]
        own = independent[word]
        assert [Fraction(x) for x in cell["strict_interior"]] == list(own["interior"])
        stored_rays = [tuple(Fraction(z) for z in ray) for ray in cell["rays_a0_equals_one"]]
        assert stored_rays == own["rays"]
        assert Fraction(cell["determinant"]) == own["determinant"]
        assert [Fraction(q) for q in cell["rates"]] == own["rates"]
    assert Fraction(stored["section_volume_sum"]) == volume
    assert Fraction(stored["base_exp_minus_symbol_integral"]) == base
    assert Fraction(stored["full_exp_minus_symbol_integral"]) == full
    assert stored["heat_kernel_peak"] == peak

    output = {
        "status": "PASS",
        "raw_body_directions": 32,
        "crossing_forms_rebuilt": len(crossing),
        "boundary_hyperplanes": len(boundaries),
        "five_hyperplane_systems": comb(18, 5),
        "singular_systems": singular,
        "arrangement_vertices": len(vertices),
        "sign_words_reclassified": 1 << len(crossing),
        "full_dimensional_cells": len(independent),
        "each_cell_rays_bound_to_exact_compatible_vertices": True,
        "each_cell_has_six_affinely_independent_vertices": True,
        "all_chamber_and_crossing_weak_signs_checked": True,
        "all_32_body_form_weak_signs_checked": True,
        "section_volume": enc(volume),
        "base_exp_minus_symbol_integral": enc(base),
        "full_exp_minus_symbol_integral": enc(full),
        "heat_kernel_peak": peak,
        "producer_certificate_sha256": sha256(CERTIFICATE.read_bytes()).hexdigest(),
        "audit_source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "scope": "independent finite Fraction audit; neither a Green certificate nor a C6 upper bound"
    }
    (ROOT / "independent_fan_audit_results.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
