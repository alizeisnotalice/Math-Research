#!/usr/bin/env python3
"""Exact finite audit for the sharp body-frame/face-quotient constants."""
from __future__ import annotations

from hashlib import sha256
from itertools import product
from pathlib import Path
import json
import sympy as sp


def main() -> None:
    signs = list(product((-1, 1), repeat=5))
    one = (1, 1, 1, 1, 1)
    quotient = [tuple(s[j] - 1 for j in range(5)) for s in signs if s != one]
    # r=-s turns s-1 into -(1+r); the extra r=-1 vector is zero.
    transformed = [tuple(-one[j] - r[j] for j in range(5)) for r in signs]
    rank = sp.Matrix(quotient).rank()
    same_multiset = sorted(quotient + [(0, 0, 0, 0, 0)]) == sorted(transformed)
    full_equal = sum(abs(1 + sum(s)) for s in signs)
    free_equal = sum(abs(sum(s)) for s in signs)
    Cf = sp.Rational(32) ** 3 / sp.Rational(12) ** 5
    assert rank == 5 and same_multiset
    assert full_equal == 60 and free_equal == 60 and Cf == sp.Rational(32, 243)
    result = {
        "full_body_l1_constant": 10,
        "full_equal_sign_sum": full_equal,
        "free_l1_constant": 12,
        "free_equal_sign_sum": free_equal,
        "quotient_directions": len(quotient),
        "quotient_rank": rank,
        "quotient_reindexing_identity": same_multiset,
        "face_sup_prefactor_over_pi6": str(Cf),
        "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    Path(__file__).with_name("sharp_body_face_audit.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
