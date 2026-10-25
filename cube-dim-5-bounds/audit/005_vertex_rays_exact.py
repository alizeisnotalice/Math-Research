#!/usr/bin/env python3
"""Exact vertex test for normalized rays in the seven-cone n=5 body fan.

This is deliberately independent of the numerical heat-kernel implementation.
For V={(1,s_2,...,s_5): s_j in {-1,1}}, it uses

    psi(x) = sum_{v in V}|v.x| = max_{eps in {-1,1}^V} (sum eps_v v).x.

Thus P={psi<=1} is the intersection of the 2^16 half-spaces a.x<=1.
For a listed fan ray r, p=r/psi(r) is a true vertex exactly when its
active half-space normals span R^5.  Everything below is integer arithmetic.
"""
from __future__ import annotations

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


def directions() -> tuple[tuple[int, ...], ...]:
    return tuple((1,) + signs for signs in itertools.product((-1, 1), repeat=4))


def dot(a: tuple[int, ...], b: tuple[int, ...]) -> int:
    return sum(x * y for x, y in zip(a, b))


def integer_rank(rows: list[tuple[int, ...]]) -> int:
    """Fraction-free enough here: exact rational row reduction with Python ints."""
    matrix = [list(row) for row in rows if any(row)]
    rank = 0
    columns = 5
    for column in range(columns):
        pivot = next((i for i in range(rank, len(matrix)) if matrix[i][column]), None)
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        pivot_value = matrix[rank][column]
        for i in range(rank + 1, len(matrix)):
            entry = matrix[i][column]
            if entry == 0:
                continue
            # Bareiss-style cross elimination preserves exact integrality; gcd
            # normalization controls coefficient growth without changing rank.
            newrow = [pivot_value * matrix[i][j] - entry * matrix[rank][j] for j in range(columns)]
            nonzero = [abs(value) for value in newrow if value]
            if nonzero:
                divisor = nonzero[0]
                for value in nonzero[1:]:
                    import math
                    divisor = math.gcd(divisor, value)
                newrow = [value // divisor for value in newrow]
            matrix[i] = newrow
        rank += 1
        if rank == columns:
            return rank
    return rank


def primitive_positive(ray: tuple[int, ...]) -> tuple[int, ...]:
    import math
    divisor = 0
    for entry in ray:
        divisor = math.gcd(divisor, abs(entry))
    ray = tuple(entry // divisor for entry in ray)
    for entry in ray:
        if entry:
            return ray if entry > 0 else tuple(-value for value in ray)
    raise ValueError("zero ray")


def all_halfspace_normals(vectors: tuple[tuple[int, ...], ...]) -> tuple[tuple[int, ...], ...]:
    normals = set()
    for signs in itertools.product((-1, 1), repeat=len(vectors)):
        normals.add(tuple(sum(sign * vector[j] for sign, vector in zip(signs, vectors)) for j in range(5)))
    return tuple(sorted(normals))


def main() -> None:
    vectors = directions()
    normals = all_halfspace_normals(vectors)
    rays = tuple(sorted({primitive_positive(ray) for cone in BASE_CONES for ray in cone}))
    records = []
    for ray in rays:
        q = sum(abs(dot(vector, ray)) for vector in vectors)
        active = [normal for normal in normals if dot(normal, ray) == q]
        # A direct maximization check guards against an accidental incorrect
        # sign convention in the half-space representation.
        assert max(dot(normal, ray) for normal in normals) == q
        rank = integer_rank(active)
        records.append({
            "primitive_ray": list(ray),
            "q_equals_psi_ray": q,
            "active_distinct_halfspace_normals": len(active),
            "active_normal_rank": rank,
            "normalized_ray_is_true_vertex": rank == 5,
        })
    result = {
        "status": "PASS_exact_halfspace_vertex_test",
        "dimension": 5,
        "direction_count": len(vectors),
        "sign_choices_enumerated": 2 ** len(vectors),
        "distinct_halfspace_normals": len(normals),
        "ray_count_in_base_cones_after_projective_deduplication": len(rays),
        "records": records,
        "true_vertex_ray_count": sum(record["normalized_ray_is_true_vertex"] for record in records),
        "nonvertex_ray_count": sum(not record["normalized_ray_is_true_vertex"] for record in records),
        "interpretation": "rank<5 means the normalized fan ray lies in a positive-dimensional face of P={psi<=1}; it is an artificial star-triangulation ray, not a vertex of P.",
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    target = Path(__file__).with_name("005_vertex_rays_exact_results.json")
    target.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
