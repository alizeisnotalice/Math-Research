#!/usr/bin/env python3
"""Certified Gaussian dilation components only; no actual FIRST simulation."""
from fractions import Fraction as F
from math import factorial, isqrt
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "physical_scale_frozen_transfer_exact_guard_20261007_results.json"


def exp_small(x):
    """Rational enclosure exp(x), 0<=x<=1; remainder <=3 x^(N+1)/(N+1)!."""
    assert 0 <= x <= 1
    count = 36
    total = sum((x**j / factorial(j) for j in range(count + 1)), F(0))
    return total, total + 3*x**(count+1)/factorial(count+1)


def exp_negative_small(x):
    lower, upper = exp_small(x)
    return 1/upper, 1/lower


def log_four_thirds():
    x = F(1, 3)
    lower = sum(((-1)**(j+1)*x**j/j for j in range(1, 41)), F(0))
    upper = lower + x**41/41
    return lower, upper


def log_density_interval(m, q, u):
    elo, ehi = exp_negative_small(2*u)
    return -m*u-q*ehi, -m*u-q*elo


def frac_hash(x):
    return hashlib.sha256(f"{x.numerator}/{x.denominator}".encode()).hexdigest()


lnlo, lnhi = log_four_thirds()
assert lnlo > F(1, 4) > F(3, 16)
rows = []
count_total = 1
for n in (576, 1024, 4096):
    rootn = isqrt(n)
    assert rootn*rootn == n
    delta = F(1, 2*rootn)
    checks = 0
    signatures = []
    margins = []
    # Uniform continuous certificate, including every face dimension m<=n.
    for m in range(n+1):
        assert m*delta**2 <= F(1, 4)
        assert m*delta**2*3/4 <= F(3, 16) < lnlo
        checks += 2
    # Non-asymptotic log-density enclosures at multiple actual Gaussian q values.
    for m in (1, n//4, n//2, n):
        qlist = [m*F(1,4), m*F(1,2), m*(F(1,2)+delta/4),
                 m*(F(1,2)+delta/2), m*(F(1,2)+delta), F(m), F(2*m)]
        for q in qlist:
            end0 = log_density_interval(m, q, F(0))
            end1 = log_density_interval(m, q, delta)
            endpoint_lower = max(end0[0], end1[0])
            for j in range(17):
                u = delta*j/16
                lo, hi = log_density_interval(m, q, u)
                margin = endpoint_lower + lnlo-hi
                assert margin > 0
                margins.append(margin)
                signatures.append(frac_hash(margin))
                checks += 1
            # The two endpoint values are included; maximum cannot be rounded down.
            assert end0[0] <= end0[1] and end1[0] <= end1[1]
            checks += 2
    # Holding component, m=0, is literally independent of scale.
    for j in range(17):
        assert log_density_interval(0, F(0), delta*j/16) == (F(0), F(0))
        checks += 1
    # Positive-mixture algebra independent of prior labels or receiver selection.
    weights = [F(1,19), F(2,19), F(5,19), F(11,19)]
    left = [F(1), F(7,3), F(1,5), F(11,2)]
    right = [F(5,2), F(1,7), F(4), F(3,2)]
    component_peak = [F(4,3)*max(a,b) for a,b in zip(left,right)]
    mix_peak = sum((w*p for w,p in zip(weights,component_peak)), F(0))
    mix_end_sum = F(4,3)*sum((w*(a+b) for w,a,b in zip(weights,left,right)), F(0))
    assert sum(weights) == 1 and mix_peak <= mix_end_sum
    checks += 2
    # Registered fixed-source endpoint union cost.  log(b/a)<=log2<1.
    assert 2*rootn+2 <= 4*rootn
    assert F(8,3)*(2*rootn+2) <= F(32,3)*rootn
    checks += 2
    count_total += checks
    rows.append({"n":n, "sqrt_n":rootn, "delta":str(delta), "checks":checks,
                 "log_guard_count":len(margins), "min_log_margin_approx":float(min(margins)),
                 "min_log_margin_sha256":frac_hash(min(margins)),
                 "all_margin_hashes_sha256":hashlib.sha256("".join(signatures).encode()).hexdigest(),
                 "mixture_left":str(mix_peak), "mixture_right":str(mix_end_sum),
                 "status":"PASS_EXACT_GAUSSIAN_SCALE_COMPONENTS_ONLY"})

result = {"scope":"Gaussian centered dilation, exact constants, positive mixture and fixed-source union algebra only; not Gc quadrature, nonhomogeneous transfer, cube hard average or actual FIRST/history",
          "random_seed":"none; deterministic", "arithmetic":"Fraction rational Taylor and alternating-log enclosures",
          "log_4_3_lower_sha256":frac_hash(lnlo), "log_4_3_upper_sha256":frac_hash(lnhi),
          "checks":count_total, "rounds":rows,
          "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
OUT.write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n")
print(json.dumps({"checks":count_total,"rounds":[{"n":r["n"],"checks":r["checks"],"status":r["status"]} for r in rows],"output":str(OUT)},ensure_ascii=False))
