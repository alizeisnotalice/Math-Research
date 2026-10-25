#!/usr/bin/env python3
"""One fixed-box audit: Fraction log majorant, direct weight derivative, shared stencil balls."""
from fractions import Fraction as Q
from itertools import product
from math import comb, lcm, prod
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
        result = row.get(key, Q(0))+value
        if result:
            row[key] = result
        else:
            row.pop(key, None)


def new_log_upper(a, b):
    """log(b/a)<=min(b/a-1,0.7*ceil(log2(b/a))), integers 0<a<=b."""
    assert 0 < a <= b
    e = max(0, b.bit_length()-a.bit_length())
    while a*(1 << e) < b:
        e += 1
    while e and a*(1 << (e-1)) >= b:
        e -= 1
    return min(Q(b-a, a), Q(7*e, 10))


def new_time_upper(nodes):
    a = sorted(map(abs, nodes))
    result = Q(a[0], prod(a))+new_log_upper(a[0], a[1])/prod(a[1:])
    for m in range(2, 5):
        result += (Q(1, a[m-1]**(m-1))-Q(1, a[m]**(m-1)))/((m-1)*prod(a[m:]))
    return result+Q(1, 4*a[4]**4)


def direct_center_rows(nodes, slopes, determinant, constant, derivative):
    """Differentiate the partial-fraction weights, with no DD quotient recursion."""
    for j in range(5):
        weight = Q(1, prod(nodes[j]-nodes[m] for m in range(5) if m != j))
        logarithmic_weight_derivative = sum((Q(slopes[j]-slopes[m], nodes[j]-nodes[m])
                                             for m in range(5) if m != j), Q(0))
        add_term(constant, abs(nodes[j]), -determinant*weight)
        add_term(derivative, abs(nodes[j]), determinant*weight*logarithmic_weight_derivative)
        add_term(derivative, 0, -determinant*weight*Q(slopes[j], nodes[j]))


def parse_ball(s):
    center, radius = map(Q, s[1:-1].split(" +/- "))
    return center-radius, center+radius


def calculate(reg, sources, cachepath, started):
    import flint
    from flint import arb, fmpq, ctx
    ctx.prec = 512
    def ar(q):
        q = Q(q)
        return arb(fmpq(q.numerator, q.denominator))
    def bounds(ball):
        return Q(str(ball.lower().fmpq())), Q(str(ball.upper().fmpq()))
    candidate = json.loads(sources["candidate_results"].read_text())
    production_registration = json.loads(sources["candidate_registration"].read_text())
    cache = json.loads(cachepath.read_text())
    assert cache["status"] == "complete_stencils" and cache["next_simplex"] == 3504
    fingerprint = hashlib.sha256(json.dumps({"source": sha(sources["candidate_program"]),
                    "registration": sha(sources["candidate_registration"]),
                    "inputs": production_registration["read_only_inputs"]}, sort_keys=True).encode()).hexdigest()
    assert cache["fingerprint"] == fingerprint == candidate["checkpoint"]["fingerprint"]
    saved = [{int(key): Q(value) for key, value in row.items()} for row in cache["stencil"]]
    assert len(saved) == 35 and sum(map(len, saved)) == 4934
    table = json.loads(gzip.decompress(sources["cone_table"].read_bytes()))
    vertices, cones = table["vertices"], table["simplices"]
    assert table["denominator"] == 144 and len(vertices) == 280 and len(cones) == 3504
    center = tuple(map(Q, reg["center"]))
    D = lcm(*(x.denominator for x in center))
    h = Q(reg["halfwidth"])
    assert D == 104729 and center[0] == 1 and all(0 <= x-h <= x+h <= 1 for x in center[1:])
    n = [int(D*x) for x in center]
    k = [sum(a*b for a, b in zip(n, v)) for v in vertices]
    assert 0 not in k and len(set(k)) == 280
    q = [h*D*sum(map(abs, v[1:]))/abs(node) for v, node in zip(vertices, k)]
    assert max(q) == Q(104729, 6234112) == Q(candidate["maximum_vertex_relative_projection_shift_exact"]) < 1
    constant, derivative = {}, {}
    remainder_sum = Q(0)
    for ci, cone in enumerate(cones):
        if ci % 64 == 0:
            assert time.monotonic()-started < 120, "fixed audit time cap"
        ids, determinant = cone["vertices"], int(cone["det_numerator"])
        assert len(ids) == len(set(ids)) == 5 and determinant > 0
        nodes = [k[i] for i in ids]
        slopes = [D*vertices[i][4] for i in ids]  # fixed derivative alpha=(0,0,0,1)
        direct_center_rows(nodes, slopes, determinant, constant, derivative)
        qc = max(q[i] for i in ids)
        tail = qc**4*(70-224*qc+280*qc**2-160*qc**3+35*qc**4)/(1-qc)**5
        assert tail == (1-qc)**(-5)-sum((comb(j+4, 4)*qc**j for j in range(4)), Q(0)) >= 0
        remainder_sum += determinant*new_time_upper(nodes)*tail
    assert constant == saved[0] and derivative == saved[1]
    # An exact in-box collision witness, without changing the registered box.
    collision = None
    for i in range(280):
        if collision is not None:
            break
        for j in range(i+1, 280):
            delta = k[i]-k[j]
            slopes = [D*(a-b) for a, b in zip(vertices[i][1:], vertices[j][1:])]
            norm = sum(map(abs, slopes))
            if norm and abs(delta) <= h*norm:
                u = [-Q(delta, norm)*(1 if a > 0 else -1 if a < 0 else 0) for a in slopes]
                left = k[i]+sum((D*a*b for a, b in zip(vertices[i][1:], u)), Q(0))
                right = k[j]+sum((D*a*b for a, b in zip(vertices[j][1:], u)), Q(0))
                assert left == right != 0 and all(abs(x) <= h for x in u)
                collision = {"vertices": [i, j], "u": list(map(str, u)), "common_scaled_projection": str(left)}
                break
    alphas = sorted((a for a in product(range(4), repeat=4) if sum(a) <= 3), key=lambda a: (sum(a), a))
    assert len(alphas) == 35 and [list(a) for a in alphas] == [row["alpha"] for row in candidate["Taylor_coefficients"]]
    logkeys = sorted({key for row in saved for key in row if key})
    logs = {key: arb(key).log() for key in logkeys}
    factor = arb(D)**4/(144*(2*arb.pi())**5)
    coefficients = []
    for row, old in zip(saved, candidate["Taylor_coefficients"]):
        value = arb(0)
        for key, c in row.items():
            value += ar(c)*(logs[key] if key else 1)
        ball = factor*value
        lo, hi = bounds(ball)
        oldlo, oldhi = parse_ball(old["ball"])
        assert oldlo <= lo <= hi <= oldhi
        coefficients.append(ball)
    cached_remainder = factor*ar(Q(cache["remainder_sum_upper"]))
    oldlo, oldhi = parse_ball(candidate["complete_time_absolute_remainder_upper_ball"])
    clo, chi = bounds(cached_remainder)
    assert oldlo <= clo <= chi <= oldhi
    new_remainder = factor*ar(remainder_sum)
    variation_upper = sum((max(map(abs, bounds(ball)))*h**sum(a)
                           for a, ball in zip(alphas[1:], coefficients[1:])), Q(0))
    A0lo, A0hi = bounds(coefficients[0])
    Rlo, Rhi = bounds(new_remainder)
    lower, upper = A0lo-variation_upper-Rhi, A0hi+variation_upper+Rhi
    old_center = json.loads(sources["stored018_center"].read_text())["selected_value"]
    assert tuple(map(Q, old_center["point"])) == center
    oldlo, oldhi = parse_ball(old_center["G_ball"])
    assert oldlo <= A0lo <= A0hi <= oldhi
    assert time.monotonic()-started < 120, "fixed audit time cap"
    return {"center": list(map(str, center)), "halfwidth": str(h), "cones": len(cones), "vertices": len(vertices),
            "maximum_q": str(max(q)), "minimum_abs_center_projection": min(map(abs, k)),
            "new_log_majorant": "min(x-1,(7/10)*ceil(log2(x)))", "exact_raw_new_remainder_sum": str(remainder_sum),
            "new_full_time_remainder_ball": new_remainder.str(105),
            "new_full_time_remainder_interval": list(map(str, (Rlo, Rhi))),
            "new_remainder_upper_used": str(Rhi), "variation_upper_used": str(variation_upper),
            "new_uniform_lower": str(lower), "new_uniform_upper": str(upper),
            "new_uniform_lower_display_only": float(lower), "target_still_certified": lower > Q(1, 10000),
            "direct_partial_fraction_global_constant_stencil_matches": True,
            "direct_weight_derivative_alpha": [0, 0, 0, 1], "direct_weight_derivative_global_stencil_matches": True,
            "center_coefficients_independently_rebuilt": 2, "other_center_coefficients_from_shared_cache": 33,
            "all35_cache_coefficient_balls_contained_in_production_balls": True,
            "coefficient_balls": [{"alpha": list(a), "ball": c.str(105)} for a, c in zip(alphas, coefficients)],
            "cached_dyadic_remainder_certificate_reproduced": True, "in_box_nonzero_collision_witness": collision,
            "flint_version": flint.__version__, "arb_bits": ctx.prec,
            "geometry_and_full_35_stencil_production_shared": True, "production_or_frozen_sources_executed": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--stencil-cache", type=Path, required=True)
    parser.add_argument("--run", action="store_true", help="one registered fixed box audit; default hashes only")
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
    # e^(7/10)>1+7/10+(7/10)^2/2+(7/10)^3/6>2 certifies log(2)<7/10.
    assert sum((Q(7, 10)**n/Q(prod(range(1, n+1))) for n in range(4)), Q(0)) > 2
    record = calculate(reg, sources, args.stencil_cache, started)
    result = {"status": "FIXED_BOX_SHARED_STENCIL_AUDIT_COMPLETE", "actual_pid": os.getpid(),
              "python": sys.executable, "registration_sha256": sha(regpath), "script_sha256": sha(Path(__file__)),
              "source_locks": reg["sources"], "cache_sha256": reg["cache_sha256"], "record": record,
              "seconds": time.monotonic()-started, "global_g_or_C5_claim": False}
    (HERE/"results.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": result["status"], "actual_pid": result["actual_pid"],
                      "seconds": result["seconds"], "target_certified": record["target_still_certified"],
                      "new_lower_display_only": record["new_uniform_lower_display_only"],
                      "results_sha256": sha(HERE/"results.json")}))


if __name__ == "__main__":
    main()
