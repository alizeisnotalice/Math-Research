#!/usr/bin/env python3
"""One cutoff box: integer sqrt guards, Fraction Jbar and partial-fraction check."""
from fractions import Fraction as Q
from itertools import product
from math import comb, isqrt, lcm, prod
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import os
import sys
import time

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)
HERE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def add_term(row, key, value):
    if value:
        total = row.get(key, Q(0))+value
        if total:
            row[key] = total
        else:
            row.pop(key, None)


def atanh_upper(u):
    y = (u-1)/(u+1)
    assert 0 <= y <= Q(1, 3)
    return 2*sum((y**(2*j+1)/Q(2*j+1) for j in range(8)), Q(0))+2*y**17/(17*(1-y*y))


def log_upper(x):
    assert x >= 1
    e = max(0, x.numerator.bit_length()-x.denominator.bit_length())
    if Q(1 << e) > x:
        e -= 1
    u = x/Q(1 << e)
    assert 1 <= u < 2
    return e*atanh_upper(Q(2))+atanh_upper(u)


def time_upper(a, points):
    b = sorted(points)
    assert b[0] >= a > 0
    result = (b[0]-a)/prod(b)+log_upper(b[1]/b[0])/prod(b[1:])
    for m in range(2, 5):
        result += (b[m-1]**(1-m)-b[m]**(1-m))/((m-1)*prod(b[m:]))
    return result+Q(1)/(4*b[4]**4)


def direct_stencils(nodes, slopes, determinant, a, constant, derivative):
    weights, sigma = [], []
    for j in range(5):
        weights.append(Q(1, prod(nodes[j]-nodes[m] for m in range(5) if m != j)))
        sigma.append(sum((Q(slopes[j]-slopes[m], nodes[j]-nodes[m]) for m in range(5) if m != j), Q(0)))
    assert sum(weights) == 0 and sum((w*s for w, s in zip(weights, sigma)), Q(0)) == 0
    for k, slope, w, s in zip(nodes, slopes, weights, sigma):
        logkey = a.numerator**2+(a.denominator*k)**2
        add_term(constant, logkey, -determinant*w/2)
        add_term(derivative, logkey, determinant*w*s/2)
        add_term(derivative, 0, -determinant*w*slope*k/(a*a+k*k))


def parse_ball(s):
    c, r = map(Q, s[1:-1].split(" +/- "))
    return c-r, c+r


def calculate(reg, sources, cachepath, started):
    import flint
    from flint import arb, ctx, fmpq
    ctx.prec = 512
    def ar(q):
        q = Q(q)
        return arb(fmpq(q.numerator, q.denominator))
    def bounds(ball):
        return Q(str(ball.lower().fmpq())), Q(str(ball.upper().fmpq()))
    def contained(ball, printed):
        lo, hi = bounds(ball)
        a, b = parse_ball(printed)
        assert a <= lo <= hi <= b
    def power54(z):
        z = ar(z)
        return z*z.sqrt().sqrt()
    def chord(L, U, m, volume):
        assert 0 <= L < U and L <= m <= U
        return ar(volume)*(power54(L)+(power54(U)-power54(L))/ar(U-L)*ar(m-L))
    candidate = json.loads(sources["candidate_results"].read_text())
    prodreg = json.loads(sources["candidate_registration"].read_text())
    qualification = json.loads(sources["qualification"].read_text())
    cache = json.loads(cachepath.read_text())
    fingerprint = hashlib.sha256(json.dumps({"source": sha(sources["candidate_program"]),
              "registration": sha(sources["candidate_registration"]), "inputs": prodreg["sources"]}, sort_keys=True).encode()).hexdigest()
    assert cache["fingerprint"] == fingerprint == candidate["checkpoint"]["fingerprint"]
    assert cache["status"] == "complete_stencils" and cache["next_simplex"] == 3504
    saved = [{int(k): Q(v) for k, v in row.items()} for row in cache["stencil"]]
    assert len(saved) == 35 and sum(map(len, saved)) == 4934
    table = json.loads(gzip.decompress(sources["cone_table"].read_bytes()))
    vertices, cones = table["vertices"], table["simplices"]
    assert table["denominator"] == 144 and len(vertices) == 280 and len(cones) == 3504
    center = tuple(map(Q, reg["accepted_center"]))
    D, eps, h = lcm(*(x.denominator for x in center)), Q(reg["eps"]), Q(reg["halfwidth"])
    a = 144*D*eps
    assert D == 104729 and a == Q(942561, 32) == Q(candidate["integer_scaling"]["a_144D_eps"])
    assert center == tuple(map(Q, prodreg["center_candidates_in_order"][1])) == tuple(map(Q, candidate["box_center"]))
    assert all(0 < x-h < x+h < 1 for x in center[1:])
    n = [int(D*x) for x in center]
    k = [sum(x*y for x, y in zip(n, z)) for z in vertices]
    primary = tuple(map(Q, prodreg["center_candidates_in_order"][0]))
    primary_n = [int(D*x) for x in primary]
    primary_k = [sum(x*y for x, y in zip(primary_n, z)) for z in vertices]
    assert len(set(primary_k)) == 278 and len(set(k)) == 280
    assert qualification["candidate_index"] == 1 and qualification["rejected_candidates"][0]["distinct_count"] == 278
    den = 1 << 64
    b = []
    for node in k:
        squared = a*a+node*node
        root = isqrt(squared.numerator*den**2//squared.denominator)
        lo = Q(root, den)
        assert a <= lo and lo*lo <= squared < (lo+Q(1, den))**2
        b.append(lo)
    q = [h*D*sum(map(abs, z[1:]))/bv for z, bv in zip(vertices, b)]
    assert max(q) < Q(1, 32)
    witness = candidate["zero_and_repeated_node_witness"]
    i, j = witness["vertex_index"], witness["opposite_vertex_index"]
    point = tuple(map(Q, witness["zero_point"]))
    assert vertices[i] == [0, -9, 9, 0, 0] and vertices[j] == [0, 9, -9, 0, 0]
    assert point[1] == point[2] and all(abs(x-y) < h for x, y in zip(point[1:], center[1:]))
    assert sum(x*y for x, y in zip(point, vertices[i])) == sum(x*y for x, y in zip(point, vertices[j])) == 0
    crossings = sum(abs(node) < h*D*sum(map(abs, z[1:])) for z, node in zip(vertices, k))
    assert crossings == candidate["strict_zero_crossing_vertex_count"] == 4
    constant, derivative = {}, {}
    denominator_sum = 1 << 256
    rounded_sum = 0
    Jcache = {}
    for ci, cone in enumerate(cones):
        if ci % 64 == 0:
            assert time.monotonic()-started < 120, "fixed audit time cap"
        ids, det = cone["vertices"], int(cone["det_numerator"])
        assert len(ids) == len(set(ids)) == 5 and det > 0
        nodes = [k[i] for i in ids]
        slopes = [D*vertices[i][4] for i in ids]
        direct_stencils(nodes, slopes, det, a, constant, derivative)
        key = tuple(sorted(abs(v) for v in nodes))
        if key not in Jcache:
            Jcache[key] = time_upper(a, [b[i] for i in ids])
        qc = max(q[i] for i in ids)
        tail = (1-qc)**(-5)-sum((comb(n+4, 4)*qc**n for n in range(4)), Q(0))
        term = det*Jcache[key]*tail
        assert term >= 0
        rounded_sum += (term*denominator_sum).__ceil__()
    assert constant == saved[0] and derivative == saved[1]
    new_raw_upper = Q(rounded_sum, denominator_sum)
    alphas = sorted((a for a in product(range(4), repeat=4) if sum(a) <= 3), key=lambda a: (sum(a), a))
    assert [list(a) for a in alphas] == [row["alpha"] for row in candidate["Taylor_coefficients"]]
    logs = {key: arb(key).log() for row in saved for key in row if key}
    factor = arb(D)**4/(144*(2*arb.pi())**5)
    coefficients = []
    for row, old in zip(saved, candidate["Taylor_coefficients"]):
        value = sum((ar(c)*(logs[k] if k else 1) for k, c in row.items()), arb(0))
        ball = factor*value
        contained(ball, old["ball"])
        coefficients.append(ball)
    cached_R = factor*ar(Q(cache["remainder_sum_upper"]))
    contained(cached_R, candidate["cutoff_time_resolvent_remainder_upper_ball"])
    varU = sum((max(map(abs, bounds(c)))*h**sum(alpha) for alpha, c in zip(alphas[1:], coefficients[1:])), Q(0))
    meanP = sum((c*ar(h**sum(alpha)/prod(v+1 for v in alpha))
                 for alpha, c in zip(alphas, coefficients) if all(v % 2 == 0 for v in alpha)), arb(0))
    contained(meanP, candidate["mean_Taylor_polynomial_ball"])
    A0lo, A0hi = bounds(coefficients[0])
    cached_R_up = bounds(cached_R)[1]
    meanPlo, meanPhi = bounds(meanP)
    original_L, original_U = Q(candidate["H_eps_uniform_lower_dyadic"]), Q(candidate["H_eps_uniform_upper_dyadic"])
    original_m = Q(candidate["mean_H_eps_upper_dyadic"])
    assert original_L > Q(218951499, 10**12) and original_U < Q(235519403, 10**12)
    assert original_L <= max(Q(0), A0lo-varU-cached_R_up)
    assert A0hi+varU+cached_R_up <= original_U
    assert min(original_U, meanPhi+cached_R_up) <= original_m <= original_U
    volume = (2*h)**4
    assert volume == Q(1, 2**40) == Q(candidate["box_volume"])
    original_chord = chord(original_L, original_U, original_m, volume)
    contained(original_chord, candidate["chord_integral_upper_expression_ball"])
    assert bounds(original_chord)[1] <= Q(candidate["integral_H_eps_5over4_upper_dyadic"])
    assert Q(candidate["integral_H_eps_5over4_upper_dyadic"]) < Q(53, 2*10**18)
    new_R = factor*ar(new_raw_upper)
    R_up = bounds(new_R)[1]
    L, U = max(Q(0), A0lo-varU-R_up), A0hi+varU+R_up
    m = min(U, meanPhi+R_up)
    new_chord = chord(L, U, m, volume)
    integral_up = bounds(new_chord)[1]
    assert time.monotonic()-started < 120, "fixed audit time cap"
    return {"eps": str(eps), "a_144D_eps": str(a), "center": list(map(str, center)),
            "halfwidth": str(h), "box_volume": str(volume),
            "primary_distinct_projections": len(set(primary_k)), "accepted_distinct_projections": len(set(k)),
            "sqrt_lower_method": "integer isqrt at denominator2^64, all280 squared comparisons",
            "new_q_max": str(max(q)), "all_new_q_below_one_over32": True,
            "strict_zero_crossings": crossings, "zero_repeated_witness": witness,
            "incorrect_full_G_q_at_zero_ray_center": str(h*D*sum(map(abs, vertices[i][1:]))/abs(k[i])),
            "new_J_method": "8 positive atanh terms plus rational tail, complete integral froma",
            "pureFraction_remainder_rounding_denominator": str(denominator_sum),
            "raw_remainder_upper": str(new_raw_upper), "total_rounding_slack_upper": str(Q(3504, denominator_sum)),
            "new_cutoff_R_ball": new_R.str(105), "new_H_lower": str(L), "new_H_upper": str(U),
            "new_H_bounds_display_only": [float(L), float(U)],
            "new_mean_upper": str(m), "new_chord_integral_ball": new_chord.str(105),
            "new_box_integral_upper": str(integral_up), "new_integral_below_53_over2e18": integral_up < Q(53, 2*10**18),
            "original_chord_expression_reproduced": True, "original_saved_chord_cap_valid": True,
            "independent_global_partial_fraction_stencils": ["constant", "alpha0001"],
            "remaining_shared_stencil_count": 33, "all35_reassessed_balls_contained": True,
            "Arb_bits": ctx.prec, "flint_version": flint.__version__, "complete_cones": len(cones),
            "production_and_frozen_producers_executed": False,
            "object_is_H_eps_not_full_G_or_positivepart": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--stencil-cache", type=Path, required=True)
    parser.add_argument("--run", action="store_true", help="one fixed audit; default hashes only")
    args = parser.parse_args()
    root = args.project_root.resolve()
    if (root/"cube-dim-5-bounds").is_dir():
        root = root/"cube-dim-5-bounds"
    regpath = HERE/"registration.json"
    reg = json.loads(regpath.read_text())
    assert sha(Path(__file__)) == reg["script_sha256"] and sha(args.stencil_cache) == reg["cache_sha256"]
    sources = {}
    for source in reg["sources"]:
        path = root/source["project_relative"]
        assert sha(path) == source["sha256"], str(path)
        sources[source["id"]] = path
    if not args.run:
        print(json.dumps({"status": "HASH_ONLY_PASS", "sources": len(sources), "cache": True, "run": False}))
        return
    started = time.monotonic()
    print(json.dumps({"status": "AUDIT_STARTED", "actual_pid": os.getpid()}), flush=True)
    record = calculate(reg, sources, args.stencil_cache, started)
    result = {"status": "FIXED_CUTOFF_BOX_AUDIT_COMPLETE_WITH_SHARED_33_STENCILS", "actual_pid": os.getpid(),
              "python": sys.executable, "script_sha256": sha(Path(__file__)), "registration_sha256": sha(regpath),
              "source_locks": reg["sources"], "cache_sha256": reg["cache_sha256"], "record": record,
              "seconds": time.monotonic()-started, "global_I_g_or_C5_claim": False}
    (HERE/"results.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": result["status"], "actual_pid": result["actual_pid"],
                      "seconds": result["seconds"], "new_H_bounds": record["new_H_bounds_display_only"],
                      "new_integral_below_53_over2e18": record["new_integral_below_53_over2e18"],
                      "results_sha256": sha(HERE/"results.json")}))


if __name__ == "__main__":
    main()
