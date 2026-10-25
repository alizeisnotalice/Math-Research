#!/usr/bin/env python3
"""Exact finite audit for the n=5 body-frame face short-time interface."""
from __future__ import annotations
from itertools import combinations, product
import hashlib
import json
from pathlib import Path

def det(matrix):
    a = [list(row) for row in matrix]
    n = len(a); sign = 1; previous = 1
    for k in range(n - 1):
        pivot = next((i for i in range(k, n) if a[i][k]), None)
        if pivot is None: return 0
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]; sign = -sign
        pivot_value = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                value = a[i][j] * pivot_value - a[i][k] * a[k][j]
                assert value % previous == 0
                a[i][j] = value // previous
        previous = pivot_value
        for i in range(k + 1, n): a[i][k] = 0
    return sign * a[-1][-1]

def rank(rows):
    for r in range(min(len(rows), 5), 0, -1):
        for selected in combinations(rows, r):
            for columns in combinations(range(5), r):
                if det(tuple(tuple(row[j] for j in columns) for row in selected)):
                    return r
    return 0

def main():
    vectors = tuple((1,) + signs for signs in product((-1, 1), repeat=4))
    pair_ranks = [rank(x) for x in combinations(vectors, 2)]
    triple_ranks = [rank(x) for x in combinations(vectors, 3)]
    determinants = [abs(det(x)) for x in combinations(vectors, 5)]
    complements = []
    for i, v in enumerate(vectors):
        complements.append(max(abs(det((v,) + tuple(vectors[j] for j in c))) for c in combinations(tuple(j for j in range(16) if j != i), 4)))
    assert pair_ranks == [2] * 120
    assert triple_ranks == [3] * 560
    assert all(x == 48 for x in complements)
    payload = {
        "model": "n=5 body Cauchy frame V={(1,s2,s3,s4,s5): sj in {-1,1}}, c_v=1",
        "direction_count": 16,
        "pair_rank_histogram": {str(r): pair_ranks.count(r) for r in sorted(set(pair_ranks))},
        "triple_rank_histogram": {str(r): triple_ranks.count(r) for r in sorted(set(triple_ranks))},
        "five_subset_abs_det_histogram": {str(d): determinants.count(d) for d in sorted(set(determinants))},
        "per_direction_best_complement_det": complements,
        "k": 16, "A": 16, "C_face": "256/(3*pi^5)", "q": "5/4",
        "E_epsilon_positive_orthant": "16^(-4/5)*(5/6)*(256/(3*pi^5))^(1/5)*(16/pi)^(4/5)*epsilon^(6/5)"
    }
    encoded = json.dumps(payload, sort_keys=True, indent=2).encode() + b"\n"
    output = Path(__file__).with_name("round001_shorttime_results.json")
    output.write_bytes(encoded)
    print(encoded.decode(), end="")
    print("sha256", hashlib.sha256(encoded).hexdigest())

if __name__ == "__main__": main()
