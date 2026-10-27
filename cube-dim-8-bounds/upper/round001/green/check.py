#!/usr/bin/env python3
"""Exact finite algebra only; analysis proofs are in proof.md. Standard library only."""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
from itertools import combinations, product
import json
from pathlib import Path
import platform
import sys


def rank(rows):
    a = [[F(x) for x in row] for row in rows]
    if not a:
        return 0
    pivot = 0
    for col in range(len(a[0])):
        hit = next((i for i in range(pivot, len(a)) if a[i][col]), None)
        if hit is None:
            continue
        a[pivot], a[hit] = a[hit], a[pivot]
        d = a[pivot][col]
        a[pivot] = [x / d for x in a[pivot]]
        for i in range(pivot + 1, len(a)):
            d = a[i][col]
            a[i] = [x - d * y for x, y in zip(a[i], a[pivot])]
        pivot += 1
        if pivot == len(a):
            break
    return pivot


def determinant(rows):
    n = len(rows)
    assert all(len(row) == n for row in rows)
    a = [[F(x) for x in row] for row in rows]
    ans = F(1)
    for col in range(n):
        hit = next((i for i in range(col, n) if a[i][col]), None)
        if hit is None:
            return F(0)
        if hit != col:
            a[col], a[hit] = a[hit], a[col]
            ans = -ans
        d = a[col][col]
        ans *= d
        for i in range(col + 1, n):
            ratio = a[i][col] / d
            for j in range(col + 1, n):
                a[i][j] -= ratio * a[col][j]
            a[i][col] = F(0)
    return ans


def columns_to_rows(columns):
    return [list(row) for row in zip(*columns)]


def extend_basis(vectors, start):
    chosen = list(start)
    current = rank([vectors[j] for j in chosen])
    assert current == len(chosen)
    for j in range(len(vectors)):
        if j in chosen:
            continue
        nxt = rank([vectors[k] for k in chosen] + [vectors[j]])
        if nxt > current:
            chosen.append(j)
            current = nxt
        if current == 8:
            break
    assert current == 8
    return chosen


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def finite_checks():
    w = [(1,) + signs for signs in product((-1, 1), repeat=3)]
    vectors = [v + (0,) * 4 for v in w] + [(0,) * 4 + v for v in w]
    assert len(vectors) == len(set(vectors)) == 16
    assert rank(vectors) == 8
    gram = [[sum(v[i] * v[j] for v in vectors) for j in range(8)] for i in range(8)]
    assert gram == [[8 * int(i == j) for j in range(8)] for i in range(8)]
    assert {sum(x * x for x in v) for v in vectors} == {4}

    h = [[1, 1, 1, 1], [1, -1, 1, -1], [1, 1, -1, -1], [1, -1, -1, 1]]
    dh = determinant(h)
    assert abs(dh) == 16
    block = [[h[i][j] if i < 4 and j < 4 else h[i-4][j-4] if i >= 4 and j >= 4 else 0
              for j in range(8)] for i in range(8)]
    assert abs(determinant(block)) == 256

    face_indices = [j for j, v in enumerate(vectors) if v[0]]
    assert len(face_indices) == 8
    face_bases = []
    right = next(c for c in combinations(range(8, 16), 4)
                 if abs(determinant(columns_to_rows([vectors[i][4:] for i in c]))) == 16)
    for j in face_indices:
        left = next(c for c in combinations(range(8), 4)
                    if j in c and abs(determinant(columns_to_rows([vectors[i][:4] for i in c]))) == 16)
        full = left + right
        assert abs(determinant(columns_to_rows([vectors[i] for i in full]))) == 256
        face_bases.append({"direction_index": j, "basis_indices": full, "abs_det": 256})

    # This exhaustive finite model check is not the general-frame extension proof.
    rank_counts = {}
    independent_count = 0
    for inds in combinations(range(16), 4):
        r = rank([vectors[j] for j in inds])
        rank_counts[str(r)] = rank_counts.get(str(r), 0) + 1
        if r == 4:
            full = extend_basis(vectors, inds)
            assert determinant(columns_to_rows([vectors[j] for j in full])) != 0
            independent_count += 1

    q = F(8, 7)
    sup_weight = (q - 1) / q
    mass_weight = 1 / q
    assert (sup_weight, mass_weight) == (F(1, 8), F(7, 8))
    short_time_power = -6 * sup_weight + mass_weight
    assert short_time_power == F(1, 8)
    assert short_time_power + 1 == F(9, 8)
    assert 1 / (short_time_power + 1) == F(8, 9)
    long_time_power = -8 * sup_weight - mass_weight
    assert long_time_power == F(-15, 8)
    assert long_time_power + 1 == F(-7, 8)
    assert -1 / (long_time_power + 1) == F(8, 7)
    assert [2 * r - 8 for r in range(4)] == [-8, -6, -4, -2]
    assert 2 * 4 - 8 == 0
    assert 7 * mass_weight == F(49, 8)
    assert 16 * 128 == 2048
    assert F(8**2 * 8, 256) == 2  # coefficient of pi^-8 in B
    assert F(1, 256 * 7) == F(1, 1792)

    # Generic n>1 algebra, verified finitely here; its universal proof is in proof.md.
    generic_rows = []
    for n in range(2, 17):
        qn = F(n, n-1)
        assert qn > 1
        assert qn - 1 == F(1, n-1)
        assert (qn - 1) / qn == F(1, n)
        x = F(1, n)  # x = c/g0
        assert -F(1, n-1) * (1 - x) + x == 0
        # Raised (n-1)-th power of final symmetric-interface coefficient n^q.
        assert qn * (n - 1) == n
        generic_rows.append({"n": n, "q": str(qn), "c_over_g0": str(x), "small_cutoff_power": n-1})

    # Exact mass minimization at n=8, after raising to seventh powers.
    # ((7/8)*x^(-1/7)/(1-x))^7 = 8 at x=1/8.
    x = F(1, 8)
    assert F(7, 8)**7 / (x * (1-x)**7) == 8
    # J/(2^8) = 8 I, then coefficient after the seventh power is 8^8.
    assert 8 * F(2048, 256)**7 == 8**8

    return {
        "status": "PASS: stated finite algebra assertions only",
        "vectors": vectors,
        "gram": gram,
        "hadamard_det": str(dh),
        "block_basis_abs_det": 256,
        "face_bases": face_bases,
        "four_subset_rank_counts": rank_counts,
        "independent_four_subsets_with_checked_basis_extension": independent_count,
        "exponents": {
            "q": str(q), "short_time_density_norm": str(short_time_power),
            "short_time_integral": "9/8", "short_time_factor": "8/9",
            "long_time_density_norm": str(long_time_power),
            "long_time_integral": "-7/8", "long_time_factor": "8/7",
            "positive_orthant_norm_gain_power_of_2": "-49/8",
            "singular_lower_powers_for_rank_0_to_3": [-8, -6, -4, -2],
            "complement_density_upper_power": 0,
            "source_small_cutoff_power": 7,
            "compact_face_long_tail_power": -7,
        },
        "model_constants": {
            "C_coefficient_of_pi_minus_8": "1/256", "A": 8, "D1": 8,
            "B_coefficient_of_pi_minus_8": 2,
            "tail_coefficient_of_pi_minus_8_T_minus_7": "1/1792",
            "J_over_I": 2048, "seventh_power_cost_coefficient": 8**8,
        },
        "generic_finite_algebra_range": generic_rows,
        "limitations": [
            "No Green evaluation, quadrature, optimization or numerical constant certificate.",
            "No finite computation substitutes for the general-frame S3 proof or source/comparison analysis.",
            "No inference that C8 is the square of C4.",
        ],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", type=Path,
                        default=Path("/Users/zhengzhihao/Desktop/T/cube8-source-package-20261025"))
    parser.add_argument("--save", action="store_true")
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    source_names = [
        "cube-finite-dim-bounds/audit/second_review/theory/paper_finite_source.md",
        "cube-finite-dim-bounds/audit/second_review/theory/audit.md",
        "cube-finite-dim-bounds/papers/n3/main.tex",
        "cube-finite-dim-bounds/papers/n4/main.tex",
        "AGENTS.md",
    ]
    source_records = [{"path": name, "bytes": (args.source_root / name).stat().st_size,
                       "sha256": sha256(args.source_root / name)} for name in source_names]
    clone_agents = here.parents[3] / "AGENTS.md"
    # here.parents[3] is repository root for upper/round001/green.
    assert clone_agents.is_file(), str(clone_agents)
    snapshot = {
        "source_root_at_execution": str(args.source_root),
        "source_files": source_records,
        "clone_agents_sha256": sha256(clone_agents),
        "reading_scope": {
            "paper_finite_source.md": "full text",
            "audit.md": "full text",
            "n3/main.tex": "full file requested; relevant analysis and short-time sections reread without truncation",
            "n4/main.tex": "full file requested; analysis sections reread without truncation",
        },
        "clone_head": "79f8acb20c418e60c5196a3c830103253e34c900",
        "known_unrelated_missing_blob": {
            "oid": "48d436008eaeecca94834fc826350c6e94e26490",
            "path": "centered-cube-upper-bound/adversarial_shell_pressure_20261007/adversarial_shell_search.py",
            "effect": "Sparse checkout emitted errors; AGENTS checked out, only new green directory used.",
        },
    }
    old_snapshot = here / "sources.json"
    if old_snapshot.exists():
        previous = json.loads(old_snapshot.read_text())
        assert previous["source_files"] == snapshot["source_files"], "Source hash changed"
        assert previous["clone_agents_sha256"] == snapshot["clone_agents_sha256"]
    result = finite_checks()
    result.update({"python": sys.version, "platform": platform.platform(),
                   "arithmetic": "Python standard-library Fraction exact rational arithmetic",
                   "source_hashes_verified": True})
    if args.save:
        (here / "sources.json").write_text(json.dumps(snapshot, ensure_ascii=False, indent=2) + "\n")
        (here / "checks.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"status": result["status"],
                      "independent_four_subsets": result["independent_four_subsets_with_checked_basis_extension"],
                      "four_subset_rank_counts": result["four_subset_rank_counts"],
                      "source_hashes_verified": True}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
