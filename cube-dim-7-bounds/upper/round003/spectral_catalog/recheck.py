#!/usr/bin/env python3
"""Independent stdout-only replay. Never writes a file or uses a rank library."""
import hashlib
import itertools
import json
import math
import time
from collections import Counter
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
V = [(1,) + s for s in itertools.product((-1, 1), repeat=6)]
MAX_BITS = 0

def exact_rank(rows):
    """Integer row elimination with gcd reduction; preserves rational row span."""
    global MAX_BITS
    a = [list(row) for row in rows]
    r = 0
    for c in range(7):
        p = next((i for i in range(r, len(a)) if a[i][c]), None)
        if p is None:
            continue
        a[r], a[p] = a[p], a[r]
        for i in range(r + 1, len(a)):
            if a[i][c] == 0:
                continue
            u, w = a[r][c], a[i][c]
            a[i] = [u * y - w * z for y, z in zip(a[i], a[r])]
            MAX_BITS = max(MAX_BITS, *(abs(y).bit_length() for y in a[i]))
            g = math.gcd(*a[i])
            if g:
                a[i] = [y // g for y in a[i]]
        r += 1
        if r == len(a):
            break
    assert MAX_BITS <= 768
    return r

def determinant(rows):
    """Small exact Leibniz determinant; independent of symbolic libraries."""
    n = len(rows)
    s = 0
    for p in itertools.permutations(range(n)):
        inversions = sum(p[i] > p[j] for i in range(n) for j in range(i + 1, n))
        s += (-1) ** inversions * math.prod(rows[i][p[i]] for i in range(n))
    return s

def main():
    tic, cpu_tic = time.perf_counter(), time.process_time()
    catalog = json.loads((HERE / "catalog.json").read_text())["orbits"]
    results = json.loads((HERE / "results.json").read_text())
    receipt = json.loads((HERE / "receipt.json").read_text())
    for section in ["source_sha256", "output_sha256"]:
        for name, expected in receipt[section].items():
            assert hashlib.sha256((HERE / name).read_bytes()).hexdigest() == expected
    accepted = {}
    hist = Counter()
    raw = primitive = nonprimitive = zero = 0
    for a in itertools.combinations_with_replacement(range(7), 7):
        raw += 1
        x = tuple(reversed(a))
        g = math.gcd(*x)
        if g == 0:
            zero += 1
            continue
        if g != 1:
            nonprimitive += 1
            continue
        primitive += 1
        dots = [sum(vj * xj for vj, xj in zip(v, x)) for v in V]
        z = [i for i, d in enumerate(dots) if d == 0]
        r = exact_rank([V[i] for i in z])
        hist[r] += 1
        if r == 6:
            accepted[x] = (z, sum(abs(d) for d in dots))
    assert raw == 1716 == math.comb(13, 7)
    assert raw == zero + nonprimitive + primitive
    assert len(catalog) == len(accepted)
    assert len({tuple(r["primitive_xi"]) for r in catalog}) == len(catalog)
    total = 0
    for record in catalog:
        x = tuple(record["primitive_xi"])
        z, rate = accepted[x]
        assert list(z) == record["zero_V_indices"]
        assert rate == record["psi_rate"]
        assert record["unit_ball_vertex"] == [str(Fraction(c, rate)) for c in x]
        assert sum(abs(sum(Fraction(vj * xj, rate) for vj, xj in zip(v, x))) for v in V) == 1
        w = record["rank6_witness_indices"]
        assert len(w) == 6 and len(set(w)) == 6 and set(w) <= set(z)
        B = [V[i] for i in w]
        assert [list(b) for b in B] == record["rank6_witness_directions"]
        assert exact_rank(B) == record["zero_rank"] == 6
        cols = record["rank6_minor_columns"]
        det = determinant([[row[j] for j in cols] for row in B])
        assert det != 0 and det == record["rank6_minor_determinant"]
        cofactors = [(-1) ** j * determinant([[row[k] for k in range(7) if k != j] for row in B]) for j in range(7)]
        assert cofactors == record["witness_signed_cofactors"]
        assert all(c % 32 == 0 and abs(c) <= 216 for c in cofactors)
        assert all(sum(c * v for c, v in zip(cofactors, row)) == 0 for row in B)
        reduced = [c // 32 for c in cofactors]
        gd = math.gcd(*reduced)
        reduced = [c // gd for c in reduced]
        assert reduced == list(x) or reduced == [-c for c in x]
        size = 2 ** sum(c != 0 for c in x) * math.factorial(7) // math.prod(math.factorial(m) for m in Counter(x).values())
        assert size == record["orbit_cardinality"]
        total += size
    assert results["raw_tuple_count"] == raw
    assert results["zero_tuple_count"] == zero
    assert results["nonprimitive_nonzero_count"] == nonprimitive
    assert results["primitive_tuple_count"] == primitive
    assert results["primitive_zero_rank_histogram"] == {str(k): v for k, v in sorted(hist.items())}
    assert results["vertex_orbit_count"] == len(catalog)
    assert results["rejected_primitive_count"] == primitive - len(catalog)
    assert results["total_vertex_count"] == total
    assert results["accepted_psi_rate_histogram"] == {str(k): v for k, v in sorted(Counter(r["psi_rate"] for r in catalog).items())}
    assert results["accepted_max_coordinate_histogram"] == {str(k): v for k, v in sorted(Counter(r["primitive_xi"][0] for r in catalog).items())}
    print(json.dumps({"status": "PASS", "raw_tuples": raw, "primitive_tuples": primitive, "orbits": len(catalog), "vertices": total, "independent_rank_max_integer_bits": MAX_BITS, "wall_seconds": time.perf_counter() - tic, "process_cpu_seconds": time.process_time() - cpu_tic}, sort_keys=True))

if __name__ == "__main__":
    main()
