#!/usr/bin/env python3
"""Exact finite diagnostics for the unoriented six-dimensional body frame.

The 32 directions are (1,s_2,...,s_6), s_j in {-1,1}.  We exactly certify one
displayed determinant-160 basis, enumerate every six-element candidate only as
a floating diagnostic screen, and exactly check every one-, two-, and
three-direction span rank. This is a frame diagnostic, not a Green-function
or maximal-operator certificate.
"""
from __future__ import annotations

from itertools import combinations, product
from hashlib import sha256
from pathlib import Path
import json
from time import perf_counter
import numpy as np
import sympy as sp


def main() -> None:
    started = perf_counter()
    vectors = [(1,) + signs for signs in product((-1, 1), repeat=5)]
    # This exact basis, not the floating global screen below, is the
    # certificate used in the analysis.
    certified_ids = (0, 3, 5, 14, 22, 25)
    certified_det = int(sp.Matrix([vectors[i] for i in certified_ids]).det())
    if abs(certified_det) != 160:
        raise RuntimeError("the advertised determinant-160 basis failed exact verification")
    float_screen_max = -1
    float_screen_count = 0
    tested = 0
    for ids in combinations(range(len(vectors)), 6):
        det = int(round(np.linalg.det(np.asarray([vectors[i] for i in ids], dtype=float))))
        absolute = abs(det)
        tested += 1
        if absolute > float_screen_max:
            float_screen_max, float_screen_count = absolute, 1
        elif absolute == float_screen_max:
            float_screen_count += 1
    ranks = {}
    for size in (1, 2, 3):
        values = [sp.Matrix([vectors[i] for i in ids]).rank()
                  for ids in combinations(range(len(vectors)), size)]
        ranks[str(size)] = {"count": len(values), "minimum": min(values), "maximum": max(values)}
    gram = sp.Matrix(vectors).T * sp.Matrix(vectors)
    payload = {
        "frame": "V6_body_unoriented={(1,s2,...,s6): sj in {-1,1}}",
        "vector_count": len(vectors),
        "basis_candidates": tested,
        "det160_exact_certificate": {"ids": list(certified_ids), "det": certified_det},
        "floating_screen_only": {
            "largest_rounded_abs_det": float_screen_max,
            "candidate_count": float_screen_count,
            "warning": "No exact determinant verification was performed for every screened basis; this is not a global maximum claim."
        },
        "subset_rank_summary": ranks,
        "gram": [list(map(int, row)) for row in gram.tolist()],
        "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "elapsed_seconds": perf_counter() - started,
        "interpretation": "All 1-, 2-, and 3-subsets have full subset rank. The exact determinant-160 basis supplies a C t^-6 heat bound with C=1/(160 pi^6).",
    }
    target = Path(__file__).with_name("body_frame_rank_audit.json")
    target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
