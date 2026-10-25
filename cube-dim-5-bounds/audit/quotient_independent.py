#!/usr/bin/env python3
"""Independent exact volume audit of the four-dimensional quotient ball.

This deliberately does not import any upper-theory code.  It enumerates the
15 subset forms with SymPy, gets all rank-three null rays, finds polar-support
facets, and triangulates every 3-face by a centroid-to-ridges decomposition.
All objects below are integers or ``fractions.Fraction``; no floating-point
convex-hull or coverage test is used.
"""
from __future__ import annotations

from fractions import Fraction
from itertools import combinations, product
from math import gcd
from functools import reduce
from pathlib import Path
import argparse
import hashlib
import json
import time

from sympy import Matrix, ilcm


EXPECTED = Fraction(14455, 191102976)


class BudgetExceeded(RuntimeError):
    pass


def as_fraction(value):
    return Fraction(int(value.p), int(value.q)) if hasattr(value, "p") else Fraction(value)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def mean(points):
    return tuple(sum(p[j] for p in points) / len(points) for j in range(len(points[0])))


def cross3(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def det4(columns):
    """Exact four-column determinant, independently delegated to SymPy."""
    return as_fraction(Matrix.hstack(*(Matrix(c) for c in columns)).det())


def affine_rank(points):
    if len(points) < 2:
        return 0
    return Matrix([list(sub(p, points[0])) for p in points[1:]]).rank()


def primitive_null_ray(rows):
    """A primitive integral generator of a SymPy one-dimensional nullspace."""
    vector = Matrix(rows).nullspace()[0]
    multiplier = 1
    for item in vector:
        multiplier = ilcm(multiplier, int(item.q))
    raw = [int(multiplier * item) for item in vector]
    divisor = reduce(gcd, (abs(x) for x in raw if x), 0)
    raw = [x // divisor for x in raw]
    if next(x for x in raw if x) < 0:
        raw = [-x for x in raw]
    return tuple(raw)


def local_coordinates(points, dimension):
    """Exact coordinates of points in their affine ``dimension``-flat.

    Choose an actual basis of difference vectors using SymPy rank, then a
    nonsingular coordinate minor.  This is a coordinate construction only,
    never a numerical projection.
    """
    origin = points[0]
    basis = []
    for point in points[1:]:
        candidate = basis + [sub(point, origin)]
        if Matrix.hstack(*(Matrix(v) for v in candidate)).rank() > len(basis):
            basis.append(sub(point, origin))
        if len(basis) == dimension:
            break
    if len(basis) != dimension:
        raise AssertionError("insufficient affine rank")
    columns = Matrix.hstack(*(Matrix(v) for v in basis))
    rows = next(rows for rows in combinations(range(len(origin)), dimension)
                if columns[list(rows), :].det() != 0)
    minor = columns[list(rows), :]
    answer = []
    for point in points:
        rhs = Matrix([sub(point, origin)[i] for i in rows])
        answer.append(tuple(as_fraction(v) for v in minor.inv() * rhs))
    return answer


def supporting_ridges(facet_vertices):
    """Return all exact polygonal 2-faces of a convex 3-polytope.

    The face test is purely sign-based in exact affine coordinates: a triple
    defines a ridge iff all vertices lie in one closed halfspace of its plane.
    Equalities are collected, so quadrilateral ridges are not split.
    """
    local = local_coordinates(facet_vertices, 3)
    ridges = set()
    for i, j, k in combinations(range(len(local)), 3):
        normal = cross3(sub(local[j], local[i]), sub(local[k], local[i]))
        if normal == (0, 0, 0):
            continue
        heights = [dot(normal, sub(point, local[i])) for point in local]
        if all(value >= 0 for value in heights) or all(value <= 0 for value in heights):
            equality = frozenset(index for index, value in enumerate(heights) if value == 0)
            if len(equality) >= 3 and affine_rank([local[q] for q in equality]) == 2:
                ridges.add(equality)
    return tuple(sorted(ridges, key=lambda face: tuple(sorted(face))))


def polygon_edges(points):
    """Exact boundary edges of a convex polygon in an affine two-flat."""
    local = local_coordinates(points, 2)
    edges = set()
    for i, j in combinations(range(len(local)), 2):
        signs = []
        for k in range(len(local)):
            area = ((local[j][0] - local[i][0]) * (local[k][1] - local[i][1])
                    - (local[j][1] - local[i][1]) * (local[k][0] - local[i][0]))
            signs.append(area)
        if (all(value >= 0 for value in signs) or all(value <= 0 for value in signs)):
            equality = [q for q, value in enumerate(signs) if value == 0]
            # A one-dimensional face has exactly its two endpoint vertices.
            if len(equality) == 2:
                edges.add((i, j))
    if len(edges) != len(points):
        raise AssertionError(("polygon edge count", len(points), len(edges)))
    return tuple(sorted(edges))


def checkpoint(label, start, limit):
    elapsed = time.monotonic() - start
    print(f"checkpoint {label}: {elapsed:.3f}s", flush=True)
    if elapsed > limit:
        raise BudgetExceeded(f"soft time limit exceeded at {label}: {elapsed:.3f}s")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--soft-seconds", type=float, default=115.0)
    args = parser.parse_args()
    started = time.monotonic()
    forms = tuple(tuple((mask >> j) & 1 for j in range(4)) for mask in range(1, 16))

    vertices = set()
    total_triples = 0
    rank_three = 0
    for triple in combinations(forms, 3):
        total_triples += 1
        matrix = Matrix(triple)
        if matrix.rank() != 3:
            continue
        rank_three += 1
        ray = primitive_null_ray(triple)
        psi = 2 * sum(abs(dot(form, ray)) for form in forms)
        for sign in (-1, 1):
            vertex = tuple(Fraction(sign * item, psi) for item in ray)
            zeros = [form for form in forms if dot(form, vertex) == 0]
            assert Matrix(zeros).rank() == 3
            vertices.add(vertex)
    vertices = tuple(sorted(vertices))
    assert (total_triples, rank_three, len(vertices)) == (455, 430, 90)
    checkpoint("rank-three rays and vertices", started, args.soft_seconds)

    supports = tuple(sorted({
        tuple(sum(signs[i] * form[j] for i, form in enumerate(forms)) for j in range(4))
        for signs in product((-1, 1), repeat=15)
    }))
    facets = []
    for support in supports:
        equal = tuple(i for i, vertex in enumerate(vertices)
                      if 2 * dot(support, vertex) == 1)
        if affine_rank([vertices[i] for i in equal]) == 3:
            facets.append((support, equal))
    assert (len(supports), len(facets)) == (1611, 370)
    checkpoint("exact supports and facets", started, args.soft_seconds)

    volume = Fraction(0)
    tetrahedra = 0
    ridge_total = 0
    edge_total = 0
    facet_vertex_histogram = {}
    ridge_vertex_histogram = {}
    for facet_number, (_, ids) in enumerate(facets, 1):
        facet = tuple(vertices[i] for i in ids)
        center3 = mean(facet)
        ridges = supporting_ridges(facet)
        if not ridges:
            raise AssertionError("facet has no ridges")
        ridge_total += len(ridges)
        facet_vertex_histogram[str(len(facet))] = facet_vertex_histogram.get(str(len(facet)), 0) + 1
        for ridge_ids in ridges:
            ridge = tuple(facet[i] for i in sorted(ridge_ids))
            center2 = mean(ridge)
            edges = polygon_edges(ridge)
            edge_total += len(edges)
            ridge_vertex_histogram[str(len(ridge))] = ridge_vertex_histogram.get(str(len(ridge)), 0) + 1
            for i, j in edges:
                determinant = det4((center3, center2, ridge[i], ridge[j]))
                if determinant == 0:
                    raise AssertionError("degenerate radial simplex")
                volume += abs(determinant) / 24
                tetrahedra += 1
        if facet_number in (1, 120, 240, 370):
            checkpoint(f"centroid decomposition {facet_number}/370", started, args.soft_seconds)

    assert volume == EXPECTED
    payload = {
        "status": "complete_independent_exact_centroid_ridge_audit",
        "quotient_ball": "{x in R^4: 2 sum_{empty != S subset [4]} |sum_{j in S} x_j| <= 1}",
        "form_count": len(forms),
        "triple_count": total_triples,
        "rank3_triple_count": rank_three,
        "vertex_count": len(vertices),
        "support_vector_count": len(supports),
        "facet_count": len(facets),
        "facet_vertex_count_histogram": facet_vertex_histogram,
        "ridge_count": ridge_total,
        "ridge_vertex_count_histogram": ridge_vertex_histogram,
        "boundary_edge_count": edge_total,
        "origin_coned_4simplex_count": tetrahedra,
        "exact_volume": str(volume),
        "expected_volume": str(EXPECTED),
        "coverage_certificate": {
            "vertex": "All 455 form triples were tested with SymPy rank. A rank-three zero set gives a one-dimensional cell whose normalized ray is an exact boundary vertex; conversely every vertex has three independent zero forms, so this list is complete.",
            "facet": "Every signed subset-form sum was tested exactly against all vertices. Equality sets of affine rank three are precisely the exposed 3-faces; their 370 cones cover the boundary of the convex ball up to lower-dimensional overlaps.",
            "ridge": "Within each affine 3-face, every vertex triple was sign-tested against all face vertices. The collected equality sets are exactly its supporting 2-faces. Each ridge is partitioned by its vertex centroid coned to its exact boundary edges.",
            "volume": "The 3-face centroid coned to those ridge triangles partitions the face. Since the origin is strictly interior, coning every boundary 3-face to the origin partitions the ball up to measure-zero intersections. Each summand is abs(det(c3,c2,v_i,v_j))/4!.",
        },
        "arithmetic": "SymPy exact rank/nullspace/Matrix determinant plus Python rational arithmetic only; no floating-point hull or coverage computation.",
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    raw = json.dumps(payload, indent=2, sort_keys=True).encode() + b"\n"
    path = Path(__file__).with_name("quotient_independent_results.json")
    path.write_bytes(raw)
    print(raw.decode(), end="")
    print("sha256", hashlib.sha256(raw).hexdigest())


if __name__ == "__main__":
    main()
