#!/usr/bin/env python3
"""Complete exact catalog. Run without arguments to write only this directory."""
from __future__ import annotations
import hashlib
import itertools
import json
import math
import platform
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
import flint
import sympy as sp

HERE = Path(__file__).resolve().parent
V = [(1,) + s for s in itertools.product((-1, 1), repeat=6)]

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def rank(rows):
    return int(flint.fmpz_mat(rows).rank()) if rows else 0

def main():
    registration = json.loads((HERE / "registration.json").read_text())
    assert registration["raw_tuple_bound"] == 1716
    start_utc = datetime.now(timezone.utc).isoformat()
    tic = time.perf_counter()
    cpu_tic = time.process_time()
    records = []
    histogram = Counter()
    raw = zero = nonprimitive = primitive = 0
    for a in itertools.combinations_with_replacement(range(7), 7):
        raw += 1
        x = tuple(reversed(a))
        d = math.gcd(*x)
        if d == 0:
            zero += 1
            continue
        if d != 1:
            nonprimitive += 1
            continue
        primitive += 1
        dots = [sum(vj * xj for vj, xj in zip(v, x)) for v in V]
        zindices = [i for i, dot in enumerate(dots) if dot == 0]
        zr = rank([V[i] for i in zindices])
        assert 0 <= zr <= 6
        histogram[zr] += 1
        if zr != 6:
            continue
        zero_matrix = sp.Matrix([V[i] for i in zindices])
        independent_rows = zero_matrix.T.rref()[1]
        assert len(independent_rows) == 6
        witness_indices = [zindices[j] for j in independent_rows]
        B = sp.Matrix([V[i] for i in witness_indices])
        pivot_columns = list(B.rref()[1])
        minor_det = int(B[:, pivot_columns].det())
        assert minor_det != 0
        cofactors = [((-1) ** j) * int(B[:, [k for k in range(7) if k != j]].det())
                     for j in range(7)]
        assert all(c % 32 == 0 and abs(c) <= 216 for c in cofactors)
        scaled = [c // 32 for c in cofactors]
        gd = math.gcd(*scaled)
        kernel_primitive = [c // gd for c in scaled]
        assert kernel_primitive == list(x) or kernel_primitive == [-c for c in x]
        rate = sum(abs(dot) for dot in dots)
        mult = Counter(x)
        orbit_size = 2 ** sum(c != 0 for c in x) * math.factorial(7)
        orbit_size //= math.prod(math.factorial(m) for m in mult.values())
        records.append({
            "primitive_xi": list(x), "psi_rate": rate,
            "unit_ball_vertex": [str(Fraction(c, rate)) for c in x],
            "zero_V_indices": zindices, "zero_rank": zr,
            "rank6_witness_indices": witness_indices,
            "rank6_witness_directions": [list(V[i]) for i in witness_indices],
            "rank6_minor_columns": pivot_columns,
            "rank6_minor_determinant": minor_det,
            "witness_signed_cofactors": cofactors,
            "orbit_cardinality": orbit_size,
        })
    records.sort(key=lambda r: tuple(r["primitive_xi"]), reverse=True)
    for i, record in enumerate(records, 1):
        record["orbit_id"] = f"psi7-{i:03d}"
    assert raw == math.comb(13, 7) == 1716
    assert raw == zero + nonprimitive + primitive
    assert len({tuple(r["primitive_xi"]) for r in records}) == len(records)
    summary = {
        "dimension": 7, "body_direction_count": len(V),
        "raw_tuple_count": raw, "zero_tuple_count": zero,
        "nonprimitive_nonzero_count": nonprimitive,
        "primitive_tuple_count": primitive,
        "primitive_zero_rank_histogram": dict(sorted(histogram.items())),
        "rejected_primitive_count": primitive - len(records),
        "vertex_orbit_count": len(records),
        "total_vertex_count": sum(r["orbit_cardinality"] for r in records),
        "accepted_psi_rate_histogram": dict(sorted(Counter(r["psi_rate"] for r in records).items())),
        "accepted_max_coordinate_histogram": dict(sorted(Counter(r["primitive_xi"][0] for r in records).items())),
        "scope": registration["scope_exclusions"],
    }
    (HERE / "catalog.json").write_text(json.dumps({"index_convention": registration["body_directions"], "orbits": records}, indent=2) + "\n")
    (HERE / "results.json").write_text(json.dumps(summary, indent=2) + "\n")
    receipt = {
        "status": "success", "started_utc": start_utc,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": time.perf_counter() - tic,
        "process_cpu_seconds": time.process_time() - cpu_tic,
        "python_executable": sys.executable, "python": platform.python_version(),
        "python_flint": flint.__version__, "sympy": sp.__version__,
        "command": f"{sys.executable} -B {Path(__file__).resolve()}",
        "source_sha256": {name: sha(HERE / name) for name in ["registration.json", "catalog.py"]},
        "output_sha256": {name: sha(HERE / name) for name in ["catalog.json", "results.json"]},
        "exact_rank_calls": primitive, "body_dot_products": 64 * primitive,
        "old_baseline_evaluations": 0, "green_nodes": 0,
        "floating_rank_calls": 0, "integer_arithmetic": "unbounded exact fmpz/SymPy Integer; all input entries <=42 in absolute value; witness determinants <=216",
        "integer_budget_note": "All determinants used to certify rank are order <=7 of rows with +/-1 entries; Hadamard bounds them below 2^10, far below 768 bits. Library internal allocation width is not exposed or claimed.",
    }
    (HERE / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"status": "success", **{k: summary[k] for k in ["raw_tuple_count", "primitive_tuple_count", "vertex_orbit_count", "total_vertex_count"]}}, sort_keys=True))

if __name__ == "__main__":
    main()
