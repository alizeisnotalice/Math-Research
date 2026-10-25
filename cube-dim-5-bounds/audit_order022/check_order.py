#!/usr/bin/env python3
"""Small independent series/Bromwich and exact order-weight checks."""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
import time
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import flint
from flint import acb, arb, ctx


def show(x):
    return x.str(80, more=True)


def ar(q):
    q = F(q)
    return arb(q.numerator) / q.denominator


def direct_series(z):
    """Independent coefficient evaluation via reciprocal Gamma; no recurrence."""
    b = arb(3) / 4
    tail, first, inverse_tail = arb(0), arb(0), arb(0)
    factorial = 1
    first_k1_k5 = {}
    for k in range(1, 257):
        factorial *= k
        if k % 4 == 0:
            continue
        sign = 1 if k % 2 else -1
        tk = sign * (1 - b * k).rgamma() * z ** (-b * k) / factorial
        # Direct inverse-Laplace coefficient, not multiplier of the tail term.
        mk = b * sign * (2 - b * k).rgamma() * z ** (1 - b * k) / (factorial // k)
        tail += tk
        first += mk
        inverse_tail += tk * b * k / (1 + b * k) / z
        if k in (1, 5):
            first_k1_k5[k] = mk
    for k in (257, 258, 259):
        assert z ** (-3) / k < 1
        factorial = arb(k + 1).gamma()
        sign = 1 if k % 2 else -1
        tk = sign * (1 - b * k).rgamma() * z ** (-b * k) / factorial
        mk = b * sign * (2 - b * k).rgamma() * z ** (1 - b * k) / (factorial / k)
        nk = tk * b * k / (1 + b * k) / z
        tail += arb(0).union(tk)
        first += arb(0).union(mk)
        inverse_tail += arb(0).union(nk)
    assert all(term > 0 for term in first_k1_k5.values())
    cdf = 1 - tail
    n0 = (1 / b).gamma() / b
    inverse = n0 - inverse_tail
    assert 0 < cdf < 1 and first > 0 and inverse > 0
    return {"F": cdf, "M": first, "N": inverse, "m": first / cdf, "a": cdf / inverse}


def bromwich(z):
    """c=1 Fourier inversion of exp(-s^beta)/s and beta*s^(beta-2)exp(-s^beta)."""
    b, R, ii = arb(3) / 4, arb(256), acb(0, 1)

    def integrand(w, analytic, moment):
        s = 1 + ii * w
        sb = s.sqrt(analytic=analytic).sqrt(analytic=analytic) ** 3
        common = (ii * w * z - sb).exp()
        return common * (b * sb / s ** 2 if moment else 1 / s)

    values = {}
    kappa = (arb.pi() * 3 / 8).cos()
    common_tail = z.exp() / arb.pi() * (-kappa * R ** b).exp()
    tails = {"F": common_tail / (b * kappa * R ** b),
             "M": common_tail / (kappa * R)}
    for key, moment in (("F", False), ("M", True)):
        integral = sum((acb.integral(lambda w, flag: integrand(w, flag, moment),
                                     j * 16, (j + 1) * 16) for j in range(16)), acb(0))
        finite = (z.exp() / arb.pi() * integral).real
        assert finite.is_finite()
        values[key] = finite + arb(0, tails[key].upper())
    return values, tails


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--order001-dir", type=Path, required=True)
    parser.add_argument("--order002-dir", type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    inputs = {}
    for label, directory, script in (("001", args.order001_dir, "order_certificate.py"),
                                      ("002", args.order002_dir, "top2_certificate.py")):
        for name in ("proof.md", "report.md", script, "registration.json", "result.json"):
            inputs[label + "/" + name] = (directory / name).read_bytes()
    source1 = json.loads(inputs["001/result.json"])
    source2 = json.loads(inputs["002/result.json"])
    assert hashlib.sha256(inputs["001/order_certificate.py"]).hexdigest() == source1["source_sha256"]
    assert hashlib.sha256(inputs["002/top2_certificate.py"]).hexdigest() == source2["source_sha256"]
    ctx.prec = 384
    quarter = direct_series(ar(F(1, 4)))
    unit = direct_series(arb(1))
    for key in ("F", "m", "a"):
        assert (unit[key] - arb(source1["conditional_stats_at1"][key])).contains(0)
    inversion, inversion_tails = bromwich(ar(F(1, 4)))
    for key in ("F", "M"):
        assert (quarter[key] - inversion[key]).contains(0)

    # Three rational CDF bins, distinct from any source spatial/stable grid.
    knots = (F(1, 7), F(2, 5), F(3, 4), F(9, 10))
    pair_sum, maximum_sum, labelled_full, labelled_triangles = F(0), F(0), F(0), F(0)
    for i, (L, H) in enumerate(zip(knots, knots[1:])):
        maximum_sum += (H ** 16 - L ** 16) / 16
        for j, (l, h) in enumerate(zip(knots, knots[1:])):
            if i == j:
                weight = 2 * (H * (H ** 15 - L ** 15) / 15 - (H ** 16 - L ** 16) / 16)
            elif i > j:
                weight = (H - L) * (h ** 15 - l ** 15) / 15
            else:
                weight = (h - l) * (H ** 15 - L ** 15) / 15
            assert weight > 0
            pair_sum += weight
            labelled_full += weight * (10 * i + j + 1)
            if i == j:
                labelled_triangles += weight * (10 * i + i + 1)
            elif i > j:
                labelled_triangles += weight * ((10 * i + j + 1) + (10 * j + i + 1))
    L, H = knots[0], knots[-1]
    pair_closed = H ** 16 - 16 * H * L ** 15 + 15 * L ** 16
    assert 16 * maximum_sum == H ** 16 - L ** 16
    assert 120 * pair_sum == pair_closed
    assert labelled_full == labelled_triangles
    assert 2 * (F(1, 15) - F(1, 16)) == F(1, 120)

    directions = tuple((1,) + x for x in itertools.product((-1, 1), repeat=4))
    groups1 = Counter(F(sum(v[:4]) ** 2, 5) for v in directions)
    assert groups1 == Counter({F(16, 5): 2, F(4, 5): 8, F(0): 6})
    classes = Counter()
    for v, w in itertools.combinations(directions, 2):
        r = sum(x * y for x, y in zip(v, w))
        b, e = sum(v[:4]), sum(w[:4])
        B, E = sorted((b * b, e * e), reverse=True)
        key = (r * r, B, E, r * b * e)
        classes[key] += 1
        R, B, E, C = key
        coefficients = (484 - 4 * R + 11 * (B + E) + 2 * C,
                        220 + 4 * R - 11 * B + 5 * E - 2 * C,
                        220 + 4 * R + 5 * B - 11 * E - 2 * C,
                        100 - 4 * R - 5 * (B + E) + 2 * C)
        assert all(c >= 0 for c in coefficients)
        assert 16 - (5 + abs(r)) >= 8
    assert classes == Counter({tuple(row["invariants"]): row["count"] for row in source2["class_contributions"]})

    # Independent low-tail evaluations; high-tail formulas are audited analytically.
    beta = arb(3) / 4
    n0 = (1 / beta).gamma() / beta
    kappa = 3 * (arb(7) / 4).gamma() / (arb(2) ** (arb(7) / 2) * arb.pi() ** (arb(5) / 2))
    C = arb(4) ** (arb(7) / 4) / (256 * arb(5).sqrt())
    low1 = kappa * 16 * C * n0 ** 2 * unit["F"] ** 14 / 14
    C0 = arb(20) ** (arb(7) / 4) / (4 * arb(16) ** (arb(3) / 2))
    low2 = (120 * kappa * 2 * C0 * n0 ** (arb(3) / 2)
            * ar(F(1, 4)) ** beta * quarter["F"] ** (arb(27) / 2) / (arb(27) / 2))
    assert (low1 - arb(source1["Green_low_omission_upper"])).contains(0)
    assert (low2 - arb(source2["Green_low_omission_upper"])).contains(0)
    for source, key in ((source1, "projection_groups"), (source2, "class_contributions")):
        component_name = "Green_finite_upper_contribution" if key == "projection_groups" else "finite_upper"
        total = sum((arb(row[component_name]) for row in source[key]), arb(0))
        total += arb(source["Green_low_omission_upper"]) + arb(source["Green_high_omission_upper"])
        assert (total - arb(source["Green_formula_upper_Arb"])).contains(0)
    assert arb(source2["Green_lower_Arb"]) > arb(source1["Green_formula_upper_Arb"])
    assert arb(source2["Green_formula_upper_Arb"]) < arb(1) / 10000
    payload = {
        "status": "independent_small_formula_and_weight_audit_passed",
        "precision_bits": ctx.prec, "python_flint": flint.__version__,
        "input_sha256": {name: hashlib.sha256(raw).hexdigest() for name, raw in inputs.items()},
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "series_direct_reciprocal_gamma_at_quarter": {k: show(v) for k, v in quarter.items()},
        "series_direct_reciprocal_gamma_at_one": {k: show(v) for k, v in unit.items()},
        "unit_series_overlaps_source_stats": True,
        "quarter_Bromwich_balls": {k: show(v) for k, v in inversion.items()},
        "quarter_Bromwich_explicit_tail_balls": {k: show(v) for k, v in inversion_tails.items()},
        "Bromwich_and_series_F_M_overlap": True,
        "exact_three_bin_weights": {"knots": list(map(str, knots)),
                                   "maximum_residual": "0", "top2_closed_form_residual": "0",
                                   "labelled_cost_residual": "0", "pair_full_square_mass": "1/120"},
        "exact_projection_classes": [{"key": list(k), "count": v} for k, v in sorted(classes.items())],
        "residual_frame_minimum_eigenvalue_lower": 8,
        "recomputed_low_tail_balls": {"order001": show(low1), "order002": show(low2)},
        "source_formula_upper_ledgers_overlap": True,
        "top2_formula_lower_exceeds_top1_formula_upper": True,
        "top2_formula_upper_below_1_over_10000": True,
        "cost": {"direct_series_nodes": 2, "Bromwich_nodes": 1, "Bromwich_integrals": 2,
                 "Bromwich_subintegrals": 32, "exact_probability_bins": 3,
                 "source_grids_rerun": False, "random_samples": 0,
                 "elapsed_seconds": time.monotonic() - started},
        "scope": "Audits order-conditioned analytic lower formulas and their enclosures/method exclusions, not upper bounds on actual G or C5. The source full grids and grid-dependent high-tail A_U sum were not rerun."}
    output = Path(__file__).with_name("check_results.json")
    output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({key: payload[key] for key in ("status", "quarter_Bromwich_balls", "quarter_Bromwich_explicit_tail_balls", "recomputed_low_tail_balls", "cost")}, indent=2))


if __name__ == "__main__":
    main()
