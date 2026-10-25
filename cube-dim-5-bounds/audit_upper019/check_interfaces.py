#!/usr/bin/env python3
"""Independent finite audit of shape I bounds; no basis enumeration or sampling."""
from fractions import Fraction as F
from itertools import product
from math import factorial
from pathlib import Path
import argparse
import hashlib
import json
import re
import sys

import flint
from flint import arb, ctx
import sympy as sp

HERE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def det(rows):
    a = [list(map(F, row)) for row in rows]
    ans = F(1)
    for i in range(len(a)):
        j = next((j for j in range(i, len(a)) if a[j][i]), None)
        if j is None:
            return F(0)
        if j != i:
            a[i], a[j] = a[j], a[i]
            ans = -ans
        pivot = a[i][i]
        ans *= pivot
        for j in range(i + 1, len(a)):
            factor = a[j][i] / pivot
            for k in range(i + 1, len(a)):
                a[j][k] -= factor * a[i][k]
            a[j][i] = F(0)
    return ans


def ar(x):
    x = F(x)
    return arb(x.numerator) / x.denominator


def exact_arb_endpoints(x):
    return F(str(x.lower().fmpq())), F(str(x.upper().fmpq()))


def iroot(n, k):
    lo, hi = 0, 1 << ((n.bit_length() + k - 1) // k)
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid**k <= n:
            lo = mid
        else:
            hi = mid
    return hi if hi**k == n else lo


def root_interval(x, k, digits=90):
    assert x >= 0
    scale = 10**digits
    j = iroot(x.numerator * scale**k // x.denominator, k)
    lo = F(j, scale)
    hi = lo if lo**k == x else F(j + 1, scale)
    assert lo**k <= x <= hi**k
    return lo, hi


def endpoint_decimal(x, upper=False, digits=44):
    scale = 10**digits
    j = (x * scale).__ceil__() if upper else (x * scale).__floor__()
    sign = "-" if j < 0 else ""
    whole, rem = divmod(abs(j), scale)
    return f"{sign}{whole}.{rem:0{digits}d}"


def encode(bounds):
    lo, hi = bounds
    return {"lower_fraction": str(lo), "upper_fraction": str(hi),
            "decimal_outward": [endpoint_decimal(lo), endpoint_decimal(hi, True)]}


def saved_ball(s):
    m = re.fullmatch(r"\[([^ ]+) \+/- ([^ ]+)\]", s)
    assert m, s
    c, r = map(F, m.groups())
    return c - r, c + r


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=HERE.parents[1])
    parser.add_argument("--hash-only", action="store_true")
    args = parser.parse_args()
    registration = json.loads((HERE / "registration.json").read_text())
    sources = {}
    for row in registration["inputs"]:
        # Prefer the archival project path, then the locked author source path.
        archived = args.source_root / row["path"]
        author = Path(row["author_source_root"]) / row["path"]
        selected = next((p for p in (archived, author) if p.is_file() and sha(p) == row["sha256"]), None)
        assert selected is not None, row["path"]
        sources[row["path"]] = selected
    if args.hash_only:
        print(json.dumps({"status": "pass", "locked_inputs": len(sources), "calculation_run": False}))
        return

    def inp(relative):
        return sources["cube-dim-5-bounds/" + relative]

    old = json.loads(inp("upper/round011_shape_integral_results.json").read_text())
    new = json.loads(inp("round006/orbit_minkowski.json").read_text())
    assert old["script_sha256"] == sha(inp("upper/round011_shape_integral.py"))
    assert new["source_sha256"] == sha(inp("round006/orbit_minkowski.py"))
    assert new["registration_sha256"] == sha(inp("round006/orbit_minkowski_registration.md"))

    n, beta, q = F(5), F(3, 4), F(10, 7)
    b, a, gamma = F(3, 14), F(15, 224), F(9, 28)
    bn, an = F(3, 20), F(3, 64)
    assert b == q / 2 - F(1, 2)
    assert a == 5 * b / 16 and gamma == 16 * a - beta and 0 < gamma < beta
    assert bn == b / q and an == 5 * bn / 16 and 16 * an - beta == 0
    assert -5 * b + 1 + gamma / beta == F(5, 14)
    assert -5 * bn + 1 == F(1, 4) and 2**4 == 16
    group, patches, stabilizer = factorial(5) * 2**5, 2 * 5 * 2**4, factorial(4)
    assert group == patches * stabilizer == 3840

    V = [(1,) + tail for tail in product((-1, 1), repeat=4)]
    frame = [[sum(v[i] * v[j] for v in V) for j in range(5)] for i in range(5)]
    assert frame == [[16 * int(i == j) for j in range(5)] for i in range(5)]
    W = det(frame)
    leverages = [sum(F(z * z, 16) for z in v) for v in V]
    assert W == 16**5 and leverages == [F(5, 16)] * 16
    # Matrix determinant lemma supplies exact Cauchy--Binet inclusion weights.
    for v in V:
        rank_one_update = [[frame[i][j] + v[i] * v[j] for j in range(5)] for i in range(5)]
        assert det(rank_one_update) == W * (1 + F(5, 16))

    cases = [
        ("uniform", [F(1, 16)] * 16),
        ("one_heavy", [F(241, 256)] + [F(1, 256)] * 15),
        ("integer_weights_1_to_16", [F(i, 136) for i in range(1, 17)]),
        ("rank_one_boundary", [F(1)] + [F(0)] * 15),
    ]
    shape_tests = []
    for name, probabilities in cases:
        assert sum(probabilities) == 1
        M = [[sum(p * v[i] * v[j] for p, v in zip(probabilities, V)) for j in range(5)] for i in range(5)]
        d = det(M)
        prod_p = F(1)
        for p in probabilities:
            prod_p *= p
        # Exact sixteenth powers remove all fractional exponents from AM--GM.
        residual = d**16 - 16**80 * prod_p**5
        assert residual >= 0
        if name == "uniform":
            assert d == 1 and residual == 0
        if name == "one_heavy":
            eta = F(1, 256)
            assert d == (16 * eta)**4 * (5 - 64 * eta)
        if name == "rank_one_boundary":
            assert d == residual == 0
        shape_tests.append({"case": name, "determinant": str(d), "AMGM_power16_residual_nonnegative": True})

    # Independent symbolic time/Gaussian and angular constants.
    t, z = sp.symbols("t z", positive=True)
    raw_time = sp.integrate(t**(-sp.Rational(10, 3)) * sp.exp(-z / (4 * t**sp.Rational(4, 3))), (t, 0, sp.oo))
    kappa = sp.simplify((4 * sp.pi)**(-sp.Rational(5, 2)) * raw_time * z**sp.Rational(7, 4))
    expected_kappa = 3 * sp.gamma(sp.Rational(7, 4)) / (2**sp.Rational(7, 2) * sp.pi**sp.Rational(5, 2))
    assert sp.simplify(kappa - expected_kappa) == 0
    area = sp.simplify(2 * sp.pi**sp.Rational(5, 2) / sp.gamma(sp.Rational(5, 2)))
    assert sp.simplify(area - 8 * sp.pi**2 / 3) == 0
    assert sp.simplify(area / patches - sp.pi**2 / 60) == 0
    rising = F(1)
    for j in range(5):
        rising *= beta + j
    first_det_moment = beta**5 * factorial(5) / rising
    assert first_det_moment == F(648, 1463)

    ctx.prec = 384
    pi = arb.pi()
    # Mellin formulas without shifted Gamma arguments, independent of source code.
    Z = (ar(1).gamma() / (ar(beta) * ar(beta).gamma())) / 16
    neg_old = ar(a / beta).gamma() / (ar(beta) * ar(a).gamma())
    neg_new = ar(an / beta).gamma() / (ar(beta) * ar(an).gamma())
    pos_old = ar(1 - gamma / beta).gamma() / ar(1 - gamma).gamma()
    Qold = ar(16)**ar(-5 * b) / Z * ar(16)**ar(gamma / beta) * pos_old * neg_old**16
    Jnew = ar(16)**ar(-5 * bn) / Z * neg_new**16
    c = (3 * ar(F(7, 4)).gamma() / (ar(2)**ar(F(7, 2)) * pi**ar(F(5, 2)))) * Z
    # I^7 closed formulas; all noninteger final roots are verified by Fraction powers.
    base = ar(3)**10 / (ar(2)**65 * ar(60)**7 * pi**11)
    Iold7 = base * ar(F(7, 4)).gamma()**7 * pos_old**7 * neg_old**112
    Inew7 = base * ar(F(7, 4)).gamma()**10 * neg_new**160
    old7, new7 = exact_arb_endpoints(Iold7), exact_arb_endpoints(Inew7)
    Iold = root_interval(old7[0], 7)[0], root_interval(old7[1], 7)[1]
    Inew = root_interval(new7[0], 7)[0], root_interval(new7[1], 7)[1]
    ratio = root_interval(new7[0] / old7[1], 7)[0], root_interval(new7[1] / old7[0], 7)[1]
    needed10 = (F(10, 3)**10 * new7[0] / 32**7, F(10, 3)**10 * new7[1] / 32**7)
    needed_g = root_interval(needed10[0], 10)[0], root_interval(needed10[1], 10)[1]
    comparisons = {}
    for label, own, source in [
        ("Z", exact_arb_endpoints(Z), old["arb_256bit"]["normalization_ET_minus_beta"]),
        ("c", exact_arb_endpoints(c), old["arb_256bit"]["shape_prefactor_c"]),
        ("Q_old", exact_arb_endpoints(Qold), old["arb_256bit"]["negative_determinant_moment_upper_Q"]),
        ("I_old", Iold, old["arb_256bit"]["complete_I_upper"]),
        ("J_new", exact_arb_endpoints(Jnew), new["J_upper_ball"]),
        ("I_new", Inew, new["I_upper_ball"]),
        ("new_over_old", ratio, new["new_over_old_I_upper_ball"]),
        ("needed_g", needed_g, new["g_needed_for_cost32_ball"]),
    ]:
        saved = saved_ball(source)
        assert saved[0] <= own[0] <= own[1] <= saved[1], label
        comparisons[label] = {"independent_interval": encode(own), "contained_in_source_printed_ball": True}
    assert Iold[1] < F(1, 50000)
    assert Inew[1] < F(1, 60000)
    assert Inew[1] < F(1579288, 10**11)
    assert Inew[1] < F(157929, 10**10)
    assert Inew[0] > F(11, 10**6)
    assert needed_g[0] > F(1, 10000)
    assert ratio[1] < 1
    costs = {}
    for g in [F(1, 65000), F(1, 10000)]:
        cost7 = F(10, 3)**10 * new7[0] / g**10, F(10, 3)**10 * new7[1] / g**10
        cost = root_interval(cost7[0], 7)[0], root_interval(cost7[1], 7)[1]
        assert cost[0] > 32
        costs[str(g)] = encode(cost)

    result = {
        "status": "pass_complete_I_interfaces_and_orbit_improvement_no_new_C5_bound",
        "script_sha256": sha(Path(__file__)), "registration_sha256": sha(HERE / "registration.json"),
        "locked_input_sha256": {row["path"]: row["sha256"] for row in registration["inputs"]},
        "runtime": {"python": sys.executable, "flint": flint.__version__, "sympy": sp.__version__, "arb_bits": ctx.prec},
        "angular_constants": {"sphere_area": str(area), "pieces": patches, "group": group,
                              "face_patch_stabilizer": stabilizer, "single_patch_factor": "pi^2/60",
                              "isotropic_calibration": "M=I: both the fixed-shape Jensen and orbit norm bounds are equalities"},
        "Cauchy_Binet_weights": {"W": str(W), "inclusion_weight_fraction": "5/16",
                                 "method": "tight-frame derivative / rank-one determinant identity; no all-bases enumeration"},
        "deterministic_shape_tests": shape_tests,
        "parameters": {"old_b": str(b), "old_a": str(a), "old_gamma": str(gamma),
                       "new_b": str(bn), "new_a": str(an), "new_T_exponent": "0",
                       "first_det_moment": str(first_det_moment)},
        "symbolic_time_kappa": str(kappa), "independent_comparisons": comparisons,
        "certified_simple_I_upper": "1/60000", "certified_finer_I_upper": "1579288/100000000000",
        "parent_requested_fine_I_upper_verified": "157929/10000000000",
        "theoretical_I_bound_exceeds_target_not_actual_I_lower_bound": True,
        "conditional_cost_of_new_formula_at_fixed_g": costs,
        "certified_C5_upper": 32,
        "scope_exclusions": ["Astra trace exclusion and polynomial minorants", "4368-basis enumeration",
                             "spatial/quadrature/QMC sampling", "target g or target I certification"],
    }
    (HERE / "check_results.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"status": result["status"], "I_new": encode(Inew), "needed_g": encode(needed_g),
                      "deterministic_cases": len(cases), "locked_sources": len(sources)}, indent=2))


if __name__ == "__main__":
    main()
