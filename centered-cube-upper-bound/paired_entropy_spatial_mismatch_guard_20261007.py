#!/usr/bin/env python3
"""Exact new slow-entropy/paired-constant guard; not a flow simulation."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math

ROOT = Path(__file__).resolve().parent
REG = ROOT / "paired_entropy_spatial_mismatch_registration_20261007.json"
OUT = ROOT / "paired_entropy_spatial_mismatch_results_20261007.json"
registration = json.loads(REG.read_text())
assert registration["status"] == "registered_before_execution"
assert registration["rounds"] == [8, 32, 128]


def frac(x):
    return f"{x.numerator}/{x.denominator}"


def log2_interval(terms):
    q = F(1, 3)
    lo = 2 * sum((q ** (2 * j + 1) / (2 * j + 1)
                  for j in range(terms)), F(0))
    tail = 2 * q ** (2 * terms + 1) / ((2 * terms + 1) * (1 - q*q))
    return lo, lo + tail


def exp_interval(x, terms):
    lo = sum((x ** j / math.factorial(j) for j in range(terms + 1)), F(0))
    first = x ** (terms + 1) / math.factorial(terms + 1)
    tail = first / (1 - x / (terms + 2))
    return lo, lo + tail


lam, eta, offset = F(1, 2), F(1, 4), F(3, 2)
coarse_d = F(1, 3024)
coarse_expr = (lam * eta**2 * F(2, 3) * F(1, 3) * F(1, 4)
               / (2 * offset * (offset + eta)))
rounds = []
total_checks = 0
for n, log_terms, exp_terms in zip(
        registration["rounds"], registration["artanh_terms"],
        registration["exp_terms"]):
    count = 0
    check_names = []

    def check(name, predicate):
        assert predicate, (n, name)
        check_names.append(name)
        return 1

    gl, gu = log2_interval(log_terms)
    count += check("log2_lower_above_2_over_3", gl > F(2, 3))
    count += check("log2_upper_below_3_over_4", gu < F(3, 4))
    count += check("log2_interval_positive", F(0) < gl < gu < 1)
    e1lo, e1hi = exp_interval(F(1), exp_terms)
    count += check("e_upper_below_3", e1hi < 3)
    count += check("e_lower_above_2", e1lo > 2)
    count += check("coarse_slow_constant_exact", coarse_expr == coarse_d)

    node_records = []
    for zs in registration["time_nodes"]:
        z = F(zs)
        elo_exp, ehi_exp = exp_interval(z, exp_terms)
        el, eu = 1 / ehi_exp, 1 / elo_exp
        count += check(f"exp_inverse_interval_{zs}", 0 < el < eu <= 1)
        count += check(f"exp_inverse_above_1_over_3_{zs}", el > F(1, 3))
        bl = gl + el * (1 - gl)
        bu = gu + eu * (1 - gu)
        count += check(f"response_b_interval_{zs}", 0 < bl < bu < 1)
        lower_d = (lam * eta**2 * bl * el * (1 - gu)
                   / (2 * offset * (offset + eta)))
        count += check(f"slow_d_lower_above_1_over_3024_{zs}", lower_d > coarse_d)
        node_records.append({
            "z": zs,
            "exp_minus_z_interval": [frac(el), frac(eu)],
            "b_interval": [frac(bl), frac(bu)],
            "slow_d_certified_lower": frac(lower_d),
        })
        # The rational identity used before phase averaging. It is algebraic,
        # not a quadrature or a claim that cos phase has these finite values.
        for b in [bl, bu]:
            a = eta * b
            for j in range(-8, 9):
                q = F(j, 8)
                den = offset + a*q
                count += check(f"den_positive_{zs}_{frac(b)}_{j}", den > 0)
                residual = (-q/den + q/offset
                            - a*q*q/(offset*den))
                count += check(f"rational_phase_identity_{zs}_{frac(b)}_{j}",
                               residual == 0)

    # These are only new constant chains. The original fast a>=n/4 and
    # gap a>=n/12 are dependencies already reviewed in the cited paper.
    width = F(1, 4)
    fast_a, fast_gap_a = F(n, 4), F(n, 12)
    full_limit_ratio_squared = width**2 * coarse_d * fast_a / 4
    gap_limit_ratio_squared = width**2 * coarse_d * fast_gap_a / 4
    count += check("full_pair_limit_square_exact",
                   full_limit_ratio_squared == F(n, 256*3024))
    count += check("gap_pair_limit_square_exact",
                   gap_limit_ratio_squared == F(n, 64*36288))
    full_finite_ratio_squared = F(n, 2048**2)
    gap_finite_ratio_squared = F(n, 4096**2)
    count += check("full_finite_bound_below_half_limit",
                   full_finite_ratio_squared < full_limit_ratio_squared/4)
    count += check("gap_finite_bound_below_half_limit",
                   gap_finite_ratio_squared < gap_limit_ratio_squared/4)
    count += check("full_normalized_sqrt_n_coefficient",
                   full_finite_ratio_squared/n == F(1, 2048**2))
    count += check("gap_normalized_sqrt_n_coefficient",
                   gap_finite_ratio_squared/n == F(1, 4096**2))
    total_checks += count
    rounds.append({
        "n": n, "exact_checks": count, "status": "exact_PASS",
        "log2_interval": [frac(gl), frac(gu)],
        "e_interval": [frac(e1lo), frac(e1hi)],
        "node_records": node_records,
        "full_pair_limit_ratio_squared": frac(full_limit_ratio_squared),
        "gap_pair_limit_ratio_squared": frac(gap_limit_ratio_squared),
        "full_finite_claim_ratio_squared": frac(full_finite_ratio_squared),
        "gap_finite_claim_ratio_squared": frac(gap_finite_ratio_squared),
        "check_name_sha256": hashlib.sha256(
            "\n".join(check_names).encode()).hexdigest(),
    })

payload = {
    "status": "terminal_exact_PASS",
    "total_checks": total_checks,
    "random_seed": None,
    "scope": "new slow entropy and global paired constants; no fast-mode rerun",
    "analytic_dependency": "finite-L1 two-box lifting in companion md",
    "not_certified": ["numerical finite-L flow", "actual FIRST/CP/GP",
                      "signed commutator lower bound"],
    "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "registration_sha256": hashlib.sha256(REG.read_bytes()).hexdigest(),
    "rounds": rounds,
}
OUT.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n")
print(json.dumps({"status": payload["status"], "total_checks": total_checks,
                  "per_round": [r["exact_checks"] for r in rounds],
                  "output": str(OUT)}, ensure_ascii=False))
