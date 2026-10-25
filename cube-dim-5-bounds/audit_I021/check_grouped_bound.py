#!/usr/bin/env python3
"""Independent normal/Gram audit and two fixed Mellin formulas; never enumerates 4-subsets."""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations_with_replacement, permutations, product
from math import gcd
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


def determinant(matrix):
    a = [list(row) for row in matrix]
    previous, sign = 1, 1
    for k in range(len(a)-1):
        pivot = next((i for i in range(k, len(a)) if a[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            a[pivot], a[k] = a[k], a[pivot]
            sign *= -1
        p = a[k][k]
        for i in range(k+1, len(a)):
            for j in range(k+1, len(a)):
                numerator = p*a[i][j]-a[i][k]*a[k][j]
                assert numerator % previous == 0
                a[i][j] = numerator//previous
            a[i][k] = 0
        previous = p
    return sign*a[-1][-1]


def rank(rows):
    a = [list(map(Q, row)) for row in rows]
    r = 0
    for k in range(5):
        p = next((i for i in range(r, len(a)) if a[i][k]), None)
        if p is None:
            continue
        a[p], a[r] = a[r], a[p]
        divisor = a[r][k]
        a[r] = [x/divisor for x in a[r]]
        for i in range(r+1, len(a)):
            coefficient = a[i][k]
            a[i] = [x-coefficient*y for x, y in zip(a[i], a[r])]
        r += 1
    return r


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def covariance(rows):
    return [[sum(v[i]*v[j] for v in rows) for j in range(5)] for i in range(5)]


def sigma4(a):
    return sum(determinant([[a[i][j] for j in range(5) if j != k]
                            for i in range(5) if i != k]) for k in range(5))


def normal_orbit(normal):
    result = set()
    for p in set(permutations(normal)):
        for signs in product((-1, 1), repeat=5):
            n = tuple(a*s for a, s in zip(p, signs))
            if next(x for x in n if x) < 0:
                n = tuple(-x for x in n)
            result.add(n)
    return result


def geometry(saved):
    vectors = [(1,)+x for x in product((-1, 1), repeat=4)]
    frame = covariance(vectors)
    assert frame == [[16*int(i == j) for j in range(5)] for i in range(5)]
    candidates, valid = [], []
    for n in combinations_with_replacement(range(3), 5):
        divisor = 0
        for x in n:
            divisor = gcd(divisor, x)
        if divisor != 1:
            continue
        n = tuple(reversed(n))
        zeros = [v for v in vectors if dot(n, v) == 0]
        r = rank(zeros)
        candidates.append({"absolute_normal": n, "zero_set_size": len(zeros), "zero_set_rank": r})
        if r == 4:
            valid.append((n, zeros))
    assert len(candidates) == 15 and len(valid) == 3
    W, total_C = 16**5, 0
    rows = []
    for n, zeros in valid:
        orbit = normal_orbit(n)
        a = covariance(zeros)
        e4 = sigma4(a)
        C = 16*e4
        h = len(zeros)
        # Rank-one determinant differentiation gives the weighted inclusion sum,
        # independently of every choice of four columns.
        for v in zeros:
            updated = [[a[i][j]+v[i]*v[j] for j in range(5)] for i in range(5)]
            assert sigma4(updated)-e4 == Q(4, h)*e4
        probabilities = Counter(Q(dot(n, v)**2, 16*dot(n, n))
                                for v in vectors if dot(n, v))
        assert sum((p*count for p, count in probabilities.items()), Q(0)) == 1
        occupancy = [Q(0)]*16
        for direction in orbit:
            assert sum(dot(direction, v) == 0 for v in vectors) == h
            for i, v in enumerate(vectors):
                occupancy[i] += C if dot(direction, v) else C*Q(4, h)
        sum_C = len(orbit)*C
        total_C += sum_C
        profile = [x/sum_C for x in occupancy]
        assert len(set(profile)) == 1 and profile[0] == Q(20-h, 16)
        expected = next(t for t in saved["normal_types"] if tuple(t["normal_type_abs_sorted"]) == n)
        assert sum_C == expected["sum_completion_weight_C"]
        assert {str(p): count for p, count in probabilities.items()} == expected["probability_counts"]
        assert 20-h == expected["scope_size"]
        rows.append({"normal_type": n, "normal_directions": len(orbit), "zero_directions": h,
                     "sigma4_covariance": e4, "C_per_normal_direction": C,
                     "sum_C": sum_C, "scope_size": 20-h,
                     "conditional_scope_profile": list(map(str, profile)),
                     "probability_counts": {str(p): count for p, count in sorted(probabilities.items())}})
    rows.sort(key=lambda row: row["scope_size"])
    assert total_C == 5*W == saved["sum_C"]
    weights = [Q(row["sum_C"], total_C) for row in rows]
    assert weights == [Q(1, 2), Q(3, 8), Q(1, 8)]
    default = [sum((w*Q(row["conditional_scope_profile"][i]) for w, row in zip(weights, rows)), Q(0)) for i in range(16)]
    assert default == [Q(53, 64)]*16
    return {"normal_candidates": candidates, "types": rows, "sum_C": total_C,
            "default_type_weights": list(map(str, weights)),
            "default_coordinate_scope_occupancy": list(map(str, default)),
            "all_four_subsets_enumerated": False, "unweighted_four_subset_counts_rebuilt": False}


def parse_ball(s):
    c, r = map(Q, s.strip()[1:-1].split(" +/- "))
    return c-r, c+r


def root_floor(n, degree):
    low, high = 0, 1 << ((n.bit_length()+degree-1)//degree)
    while low+1 < high:
        mid = (low+high)//2
        if mid**degree <= n:
            low = mid
        else:
            high = mid
    return low


def arithmetic(geo, author, weighted, old, selected):
    import flint
    from flint import arb, ctx
    ctx.prec = 512
    def ar(q):
        q = Q(q)
        return arb(q.numerator)/q.denominator
    def interval(ball):
        return Q(str(ball.lower().fmpq())), Q(str(ball.upper().fmpq()))
    def contained(bounds, printed):
        a, b = parse_ball(printed)
        assert a <= bounds[0] <= bounds[1] <= b
    beta, b = Q(3, 4), Q(3, 20)
    assert 5*b == beta
    # Unshifted Mellin Gamma formulas, rather than the production shifted form.
    Z = ar(1).gamma()/(16*ar(beta)*ar(beta).gamma())
    R = [ar(8)**ar(Q(1, 4)), 1+1/ar(2).sqrt(),
         (7+8*ar(2).sqrt()+3*ar(3).sqrt())/ar(2)**ar(Q(15, 4))]
    types = sorted(author["geometry"]["normal_types"], key=lambda row: row["scope_size"])
    for value, row in zip(R, types):
        contained(interval(value), row["R_ball"])
    base = [Q(1, 2), Q(3, 8), Q(1, 8)]
    saved_weighted = next(c for c in weighted["candidates"] if list(map(Q, c["weights"])) == selected)
    output = []
    for name, weights, saved in [("default", base, None), ("fixed_weighted", selected, saved_weighted)]:
        assert sum(weights) == 1 and min(weights) > 0
        occupancy = [sum((w*Q(t["conditional_scope_profile"][i]) for w, t in zip(weights, geo["types"])), Q(0)) for i in range(16)]
        assert len(set(occupancy)) == 1
        rho = occupancy[0]
        assert 0 < rho <= 1
        if name == "fixed_weighted":
            assert rho == Q(66127, 80000) == Q(saved["rho"])
        else:
            assert rho == Q(53, 64)
        s = b*rho
        moment = ar(s/beta).gamma()/(ar(beta)*ar(s).gamma())
        KL = sum((ar(w)*(ar(w)/ar(t)).log() for w, t in zip(weights, base)), arb(0))
        D = arb(1)
        for value, w in zip(R, weights):
            D *= value**ar(-w/5)
        J = ar(16**5)**ar(-b)/Z*(ar(b)*KL).exp()*moment**ar(5/rho)*D
        # I^7 = 3^10 J^10 / (2^75 60^7 pi^11), avoiding a fractional I power.
        I7 = 3**10*J**10/(arb(2)**75*arb(60)**7*arb.pi()**11)
        L, U = interval(I7)
        den = 10**105
        lo = Q(root_floor(L.numerator*den**7//L.denominator, 7), den)
        hi = Q(root_floor(U.numerator*den**7//U.denominator, 7)+1, den)
        assert lo**7 <= L <= U <= hi**7
        if name == "default":
            contained(interval(J), author["arb_384bit"]["fractional_Holder_J_upper"])
            contained((lo, hi), author["arb_384bit"]["fractional_Holder_complete_I_upper"])
            assert hi < Q(author["strict_rational_I_upper"])
        else:
            contained(interval(J), saved["J_ball"])
            contained((lo, hi), saved["I_ball"])
        output.append({"name": name, "weights": list(map(str, weights)), "rho": str(rho),
                       "coordinate_scope_occupancy": list(map(str, occupancy)),
                       "qJ_over_original_thetaJ_by_type": list(map(str, [w/t for w, t in zip(weights, base)])),
                       "negative_moment_exponent": str(s), "mellin_gamma_arguments": [str(s/beta), str(s)],
                       "moment_power": str(5/rho), "KL_ball": KL.str(110), "J_ball": J.str(110),
                       "I_interval": [str(lo), str(hi)], "I_display_only": [float(lo), float(hi)]})
    default_interval = list(map(Q, output[0]["I_interval"]))
    new_interval = list(map(Q, output[1]["I_interval"]))
    assert new_interval[1] < default_interval[0]
    assert default_interval[1] < parse_ball(old["I_upper_ball"])[0]
    cap = Q((new_interval[1]*10**13).__ceil__(), 10**13)
    assert new_interval[1] < cap < default_interval[0]
    return {"flint_version": flint.__version__, "arb_bits": ctx.prec, "records": output,
            "fixed_weighted_strictly_below_default": True,
            "default_strictly_below_saved_old_bound_lower": True,
            "fixed_weighted_strict_rational_I_upper": str(cap),
            "old_Gamma_baseline_evaluated": False, "author_program_imported_or_executed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--run", action="store_true", help="run the one registered new audit; default hashes only")
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
    author = json.loads(sources["candidate_results"].read_text())
    geo = geometry(author["geometry"])
    numerical = arithmetic(geo, author, json.loads(sources["weighted_results"].read_text()),
                           json.loads(sources["saved_old_ball"].read_text()), list(map(Q, reg["fixed_weights"])))
    result = {"status": "PASS_GROUPED_CB_SCOPE_MELLIN_AND_ONE_FIXED_WEIGHT_EXTENSION",
              "actual_pid": os.getpid(), "python": sys.executable,
              "registration_sha256": sha(regpath), "script_sha256": sha(Path(__file__)),
              "geometry": geo, "arithmetic": numerical,
              "source_locks": reg["sources"], "seconds": time.monotonic()-started,
              "C5_upper_changed": False, "SLSQP_or_parameter_search_run": False,
              "spatial_or_probability_samples": 0}
    (HERE/"results.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": result["status"], "actual_pid": result["actual_pid"],
                      "seconds": result["seconds"], "weighted_I_upper": numerical["fixed_weighted_strict_rational_I_upper"],
                      "results_sha256": sha(HERE/"results.json")}))


if __name__ == "__main__":
    main()
