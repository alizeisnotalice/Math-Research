#!/usr/bin/env python3
"""One fixed Hermite candidate: rational Vandermonde and integer Bernstein audit."""
from collections import Counter
from fractions import Fraction as Q
from itertools import product
from math import comb, factorial, lcm
from pathlib import Path
import argparse
import hashlib
import json
import os
import sys
import time

HERE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def solve_two_rhs(matrix, rhs):
    a = [list(map(Q, row))+list(map(Q, value)) for row, value in zip(matrix, rhs)]
    for k in range(4):
        p = next(i for i in range(k, 4) if a[i][k])
        a[k], a[p] = a[p], a[k]
        divisor = a[k][k]
        a[k] = [x/divisor for x in a[k]]
        for i in range(4):
            if i != k:
                coefficient = a[i][k]
                a[i] = [x-coefficient*y for x, y in zip(a[i], a[k])]
    return [[a[i][j] for i in range(4)] for j in (4, 5)]


def root_floor(n, degree):
    lo, hi = 0, 1 << ((n.bit_length()+degree-1)//degree)
    while lo+1 < hi:
        mid = (lo+hi)//2
        if mid**degree <= n:
            lo = mid
        else:
            hi = mid
    return lo


def orbit_size(partition):
    answer = factorial(5)
    for count in Counter(partition).values():
        answer //= factorial(count)
    return answer


def unique_permutations(partition):
    counts = Counter(partition)
    keys = sorted(counts)
    def rec(prefix):
        if len(prefix) == 5:
            yield tuple(prefix)
            return
        for key in keys:
            if counts[key]:
                counts[key] -= 1
                yield from rec(prefix+[key])
                counts[key] += 1
    yield from rec([])


def decode(poly):
    return {tuple(map(int, key.split(","))): Q(value) for key, value in poly.items()}


def build_hermite(nodes, denominator):
    r1, r2 = nodes
    gamma = Q(6, 11)
    matrix = []
    for r in nodes:
        matrix.extend([[Q(1), r, r*r, r**3], [Q(0), Q(1), 2*r, 3*r*r]])
    rhs = [[1, 0], [-gamma/r1, 0], [0, 1], [0, -gamma/r2]]
    bases = solve_two_rhs(matrix, rhs)
    for row, expected in zip(matrix, rhs):
        assert [sum((x*y for x, y in zip(row, base)), Q(0)) for base in bases] == list(map(Q, expected))
    intervals = []
    for r in nodes:
        target_floor = denominator**11*r.denominator**6//r.numerator**6
        k = root_floor(target_floor, 11)
        lo, hi = Q(k, denominator), Q(k+1, denominator)
        assert lo > 0 and lo**11*r**6 < 1 < hi**11*r**6
        intervals.append((lo, hi))
    lower, upper = [Q(0)], [Q(0)]
    for i in range(4):
        lower.append(sum((base[i]*(bounds[0] if base[i] >= 0 else bounds[1])
                          for base, bounds in zip(bases, intervals)), Q(0)))
        upper.append(sum((base[i]*(bounds[1] if base[i] >= 0 else bounds[0])
                          for base, bounds in zip(bases, intervals)), Q(0)))
    assert all(a <= b for a, b in zip(lower, upper))
    fourth_derivative_multiplier = gamma*(gamma+1)*(gamma+2)*(gamma+3)
    assert fourth_derivative_multiplier > 0
    return bases, intervals, lower, upper, fourth_derivative_multiplier


def bernstein_integer(poly):
    face = {}
    for partition, value in poly.items():
        for exponent in unique_permutations(partition):
            key = exponent[1:]
            face[key] = face.get(key, Q(0))+value
    face = {key: value for key, value in face.items() if value}
    degrees = tuple(max(key[k] for key in face) for k in range(4))
    assert degrees == (6, 6, 6, 6)
    denominator = 1
    for value in face.values():
        denominator = lcm(denominator, value.denominator)
    indices = list(product(range(7), repeat=4))
    values = [int(face.get(index, Q(0))*denominator) for index in indices]
    binomial_lcm = lcm(*(comb(6, j) for j in range(7)))
    assert binomial_lcm == 60
    matrix = [[comb(k, j)*(60//comb(6, j)) if j <= k else 0 for j in range(7)] for k in range(7)]
    for axis in range(4):
        stride = 7**(3-axis)
        transformed = []
        for position in range(2401):
            k = position//stride % 7
            transformed.append(sum(matrix[k][j]*values[position+(j-k)*stride] for j in range(k+1)))
        values = transformed
    common = denominator*60**4
    assert len(values) == 2401 and min(values) > 0
    at = min(range(2401), key=values.__getitem__)
    tensor = {index: Q(value, common) for index, value in zip(indices, values)}
    serialized = json.dumps({",".join(map(str, k)): str(v) for k, v in sorted(tensor.items())}, sort_keys=True)
    integer_stream = json.dumps({"denominator": str(common), "numerators": list(map(str, values))}, separators=(",", ":"))
    return {"face_terms": len(face), "degrees": degrees, "count": len(values),
            "common_denominator": str(common), "binomial_matrix_scale_per_axis": 60,
            "minimum_numerator": str(values[at]), "minimum_index": indices[at],
            "minimum": str(tensor[indices[at]]), "all_integer_numerators_strictly_positive": True,
            "normalized_full_table_sha256": hashlib.sha256(serialized.encode()).hexdigest(),
            "common_integer_stream_sha256": hashlib.sha256(integer_stream.encode()).hexdigest()}


def calculate(reg, sources):
    candidate = json.loads(sources["candidate_result"].read_text())
    nodes = list(map(Q, reg["nodes"]))
    eta, a = Q(reg["eta"]), Q(reg["a"])
    assert 0 < nodes[0] < nodes[1] < 1 and 0 < eta < Q(1, 2) and a > 0
    assert eta == Q(candidate["eta"]) and a == Q(candidate["a"]) and nodes == list(map(Q, candidate["nodes"]))
    bases, intervals, ell, ell_upper, derivative_multiplier = build_hermite(nodes, 10**18)
    assert bases == [list(map(Q, row)) for row in candidate["Hermite_h_basis_coefficients"]]
    assert intervals == [tuple(map(Q, row)) for row in candidate["h_node_intervals"]]
    assert ell == list(map(Q, candidate["ell_d_coefficients"]))
    corner = json.loads(sources["corner_d_moments"].read_text())
    bcorner = json.loads(sources["corner_b_moments"].read_text())
    raw = json.loads(sources["space_moments"].read_text())
    families = {name: [decode(row) for row in rows] for name, rows in raw["families"].items()}
    assert set(families) == {"d1", "d2", "d3", "d4", "b"}
    checked_rows = 0
    for name, rows in families.items():
        expected = list(map(Q, bcorner["moments"]["E_b_corner_W_power"])) if name == "b" else list(map(Q, corner["moments"][name[1:]]))
        assert len(rows) == len(expected) == (12 if name == "b" else 8)
        for n, (row, target) in enumerate(zip(rows, expected)):
            assert all(len(key) == 5 and tuple(sorted(key, reverse=True)) == key
                       and sum(key) == n+(name == "b") for key in row)
            assert sum((value*orbit_size(key) for key, value in row.items()), Q(0)) == target
            checked_rows += 1
    mu = [sum((comb(7, n)*eta**n*Q(corner["moments"][str(k)][n]) for n in range(8)), Q(0)) for k in range(1, 5)]
    U = sum((ell[k]*mu[k-1] for k in range(1, 5)), Q(0))
    V = sum((comb(11, n)*eta**n*Q(bcorner["moments"]["E_b_corner_W_power"][n]) for n in range(12)), Q(0))
    assert U == Q(candidate["U"]) and V == Q(candidate["V"]) and U > 0 and V > 0
    power = {}
    for k in range(1, 5):
        for n, row in enumerate(families["d"+str(k)]):
            factor = Q(11, 4)*a**7*ell[k]*comb(7, n)*eta**n
            for key, value in row.items():
                power[key] = power.get(key, Q(0))+factor*value
    for n, row in enumerate(families["b"]):
        factor = -Q(7, 4)*a**11*comb(11, n)*eta**n
        for key, value in row.items():
            power[key] = power.get(key, Q(0))+factor*value
    power = {key: value for key, value in power.items() if value}
    assert power == decode(candidate["symmetric_power_coefficients"])
    corner_value = sum((value*orbit_size(key) for key, value in power.items()), Q(0))
    assert corner_value == Q(11, 4)*a**7*U-Q(7, 4)*a**11*V
    bernstein = bernstein_integer(power)
    minimum = Q(bernstein["minimum"])
    assert minimum > 0
    assert bernstein["normalized_full_table_sha256"] == candidate["Bernstein_sha256"]
    assert list(bernstein["minimum_index"]) == candidate["minimum_index"] and minimum == Q(candidate["minimum"])
    # Only after every global coefficient is positive, consume c_lower.
    cL = Q(json.loads(sources["constant_c"].read_text())["constants"]["c_lower"])
    import flint
    from flint import arb, ctx
    ctx.prec = 256
    pi = arb.pi()
    cball = 3/(128*arb(2).sqrt()*pi**2*pi.sqrt())
    assert cL > 0 and cL < Q(str(cball.lower().fmpq()))
    g = cL*minimum
    simple = Q((g*10**8).__floor__(), 10**8)
    previous = json.loads(sources["previous_g_result"].read_text())
    gref = max(Q(row["global_g_exact_rational"]) for row in previous["candidates"])
    difference = g-gref
    assert g == Q(candidate["exact_g"]) and difference == Q(candidate["new_minus_ref_exact"]) > 0
    assert simple == Q(1479, 25000000) == Q(candidate["simple_g"]) <= g
    return {"nodes": list(map(str, nodes)), "eta": str(eta), "a": str(a),
            "positive_f_lower": str(1-2*eta), "h_fourth_derivative_multiplier": str(derivative_multiplier),
            "Vandermonde_h_basis_coefficients": [list(map(str, row)) for row in bases],
            "integer_root_node_enclosures": [list(map(str, row)) for row in intervals],
            "ell_lower_coefficients": list(map(str, ell)), "ell_upper_coefficients": list(map(str, ell_upper)),
            "shared_moment_rows_checked": checked_rows, "original_graph_moments_rebuilt": False,
            "symmetric_polynomial_terms": len(power), "U": str(U), "V": str(V),
            "Bernstein": bernstein, "c_lower_consumed_after_global_positivity": True,
            "c_lower": str(cL), "c_actual_Arb_ball": cball.str(60), "flint_version": flint.__version__,
            "c_Arb_bits": ctx.prec, "exact_g": str(g), "simple_g": str(simple),
            "exact_previous_009_g": str(gref), "exact_improvement": str(difference),
            "g_display_only": float(g), "old_g_tables_recomputed": False,
            "Gauss_node_or_a_optimality_checked": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--run", action="store_true", help="execute the one new fixed candidate audit; default hashes only")
    args = parser.parse_args()
    root = args.project_root.resolve()
    if (root/"cube-dim-5-bounds").is_dir():
        root = root/"cube-dim-5-bounds"
    regpath = HERE/"registration.json"
    reg = json.loads(regpath.read_text())
    assert sha(Path(__file__)) == reg["script_sha256"]
    sources = {}
    for source in reg["sources"]:
        path = root/source["project_relative"]
        assert sha(path) == source["sha256"], str(path)
        sources[source["id"]] = path
    if not args.run:
        print(json.dumps({"status": "HASH_ONLY_PASS", "sources": len(sources), "run": False}))
        return
    started = time.monotonic()
    print(json.dumps({"status": "AUDIT_STARTED", "actual_pid": os.getpid()}), flush=True)
    record = calculate(reg, sources)
    result = {"status": "PASS_FIXED_HERMITE_GLOBAL_G_ON_LOCKED_MOMENT_INPUTS",
              "actual_pid": os.getpid(), "python": sys.executable,
              "registration_sha256": sha(regpath), "script_sha256": sha(Path(__file__)),
              "source_locks": reg["sources"], "record": record,
              "seconds": time.monotonic()-started,
              "author_program_or_helper_imported": False, "new_C5_bound": False}
    (HERE/"results.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": result["status"], "actual_pid": result["actual_pid"],
                      "seconds": result["seconds"], "simple_g": record["simple_g"],
                      "Bernstein_count": record["Bernstein"]["count"],
                      "results_sha256": sha(HERE/"results.json")}))


if __name__ == "__main__":
    main()
