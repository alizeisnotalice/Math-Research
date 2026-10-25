#!/usr/bin/env python3
"""One certified conditional-point pilot; no old order grid is executed."""
from __future__ import annotations
import argparse
import hashlib
import json
import time
from fractions import Fraction as F
from pathlib import Path
import flint
from flint import arb, ctx


def ar(q):
    q = F(q)
    return arb(q.numerator) / q.denominator


def show(x):
    return x.str(80, more=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--audited-stats-json", type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    input_raw = args.audited_stats_json.read_bytes()
    previous = json.loads(input_raw)
    assert previous["status"] == "independent_small_formula_and_weight_audit_passed"
    ctx.prec = 384
    m = arb(previous["series_direct_reciprocal_gamma_at_one"]["m"])
    a = arb(previous["series_direct_reciprocal_gamma_at_one"]["a"])
    beta = arb(3) / 4
    coefficients = [arb(0)] * 260
    factorial = 1
    for k in range(1, 260):
        factorial *= k
        if k % 4:
            coefficients[k] = (1 if k % 2 else -1) * (1 - beta * k).rgamma() / factorial

    def cdf(z):
        assert z >= arb(1) / 4
        r, partial = z ** (-beta), arb(0)
        for k in range(256, 0, -1):
            partial = partial * r + coefficients[k]
        tail = partial * r
        for k in (257, 258, 259):
            assert z ** (-3) / k < 1
            tail += arb(0).union(coefficients[k] * r ** k)
        result = 1 - tail
        assert 0 < result < 1
        return result

    moment_bins = 128
    ell = F(1, 4)
    step = (1 - ell) / moment_bins
    endpoints = [ell + j * step for j in range(moment_bins + 1)]
    Fs = []
    for endpoint in endpoints:
        if time.monotonic() - started > 120:
            raise TimeoutError("registered 120s checkpoint")
        Fs.append(cdf(ar(endpoint)))
    F1 = Fs[-1]
    assert (F1 - arb(previous["series_direct_reciprocal_gamma_at_one"]["F"])).contains(0)
    masses = [right - left for left, right in zip(Fs, Fs[1:])]
    assert all(mass > 0 for mass in masses)
    raw_moment_balls = {}
    omitted_low_moment_bounds = {}
    for power in (2, 3):
        low = sum((ar(left) ** power * mass for left, mass in zip(endpoints, masses)), arb(0))
        high = sum((ar(right) ** power * mass for right, mass in zip(endpoints[1:], masses)), arb(0))
        low_omission = ar(ell) ** power * Fs[0]
        raw_moment_balls[power] = low.union(high + low_omission)
        omitted_low_moment_bounds[power] = low_omission
    e2, e3 = raw_moment_balls[2] / F1, raw_moment_balls[3] / F1
    variance = e2 - m ** 2
    mu3 = e3 - 3 * m * e2 + 2 * m ** 3
    assert variance > 0

    # For the other 15 body directions, |v.w|^2=9 has count5 and =1 count10.
    def c_leverage(frame_c, r2):
        return (25 - r2) / (80 * frame_c) + r2 / (5 * (11 * frame_c + 5))

    c9, c1 = c_leverage(m, 9), c_leverage(m, 1)
    sum_c2, sum_c3 = 5 * c9 ** 2 + 10 * c1 ** 2, 5 * c9 ** 3 + 10 * c1 ** 3
    cubic_gap = variance / 2 * sum_c2 - mu3 / 3 * sum_c3

    hessian_lower, hessian_upper = arb(0), arb(0)
    hessian_bins = 64
    for j in range(hessian_bins):
        left, right = F(j, hessian_bins), F(j + 1, hessian_bins)
        weight = ar((right - left) * (1 - (left + right) / 2))
        # B_u has residual coefficient c(u)=m+u(1-m); leverage decreases in u.
        c_left, c_right = m + ar(left) * (1 - m), m + ar(right) * (1 - m)
        sum_low = 5 * c_leverage(c_right, 9) ** 2 + 10 * c_leverage(c_right, 1) ** 2
        sum_high = 5 * c_leverage(c_left, 9) ** 2 + 10 * c_leverage(c_left, 1) ** 2
        hessian_lower += weight * sum_low
        hessian_upper += weight * sum_high
    hessian_integral = hessian_lower.union(hessian_upper)
    hessian_gap = variance * hessian_integral
    assert hessian_gap > 0
    # Universal Amax=16I endpoint check for this one conditional event.
    assert hessian_lower > arb(375) / 512
    safe_gap_lower = arb(0)
    for candidate in (cubic_gap.lower(), hessian_gap.lower()):
        if candidate > safe_gap_lower:
            safe_gap_lower = candidate
    safe_multiplier = (safe_gap_lower / 2).exp()
    assert safe_multiplier > 1
    q = ar(F(4, 5)) / (16 * a) + ar(F(16, 5)) / (11 * a + 5)
    baseline_phi = (16 * m) ** (-2) * (11 * m + 5) ** (-arb(1) / 2) * q ** (-arb(7) / 4)
    conditional_phi_with_safe_multiplier = baseline_phi * safe_multiplier
    output = {
        "status": "one_conditional_point_detgap_certified_NOT_a_global_g_or_G_certificate",
        "precision_bits": ctx.prec, "python_flint": flint.__version__,
        "input_sha256": hashlib.sha256(input_raw).hexdigest(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "condition": {"point": [1, 1, 1, 1, 0], "winner_direction": [1, 1, 1, 1, 1],
                      "winner_value": 1, "remaining_variables": 15, "truncation_threshold": 1,
                      "Amax": "16I", "leverage_classes": [{"r_squared": 9, "count": 5}, {"r_squared": 1, "count": 10}]},
        "CDF_moment_partition": {"left": "1/4", "right": "1", "bins": moment_bins,
                                 "step": str(step), "all_stieltjes_weights_positive": True,
                                 "full_low_interval_included_by_bounds": True},
        "F1_ball": show(F1), "m_ball": show(m), "a_ball": show(a),
        "raw_M2_ball": show(raw_moment_balls[2]), "raw_M3_ball": show(raw_moment_balls[3]),
        "low_raw_moment_omission_bounds": {str(k): show(v) for k, v in omitted_low_moment_bounds.items()},
        "conditional_second_raw_moment_ball": show(e2), "conditional_third_raw_moment_ball": show(e3),
        "conditional_variance_ball": show(variance), "conditional_central_mu3_ball": show(mu3),
        "leverage_sum_c_squared_ball": show(sum_c2), "leverage_sum_c_cubed_ball": show(sum_c3),
        "cubic_gap_ball": show(cubic_gap), "hessian_positive_integral_ball": show(hessian_integral),
        "hessian_gap_ball": show(hessian_gap), "safe_gap_lower_ball": show(safe_gap_lower),
        "safe_conditional_Green_multiplier_ball": show(safe_multiplier),
        "baseline_conditional_phi_ball": show(baseline_phi),
        "conditional_phi_times_safe_multiplier_ball": show(conditional_phi_with_safe_multiplier),
        "cost": {"conditional_retained_value_nodes": 1, "CDF_nodes": len(Fs),
                 "hessian_u_bins": hessian_bins, "old_order_grids_rerun": False,
                 "random_samples": 0, "elapsed_seconds": time.monotonic() - started},
        "scope": "A single conditional event parameter, not a positive-probability retained-value box or an integrated order formula. No new g, actual G upper bound or C5 result."}
    path = Path(__file__).with_name("detgap_results.json")
    path.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({key: output[key] for key in ("status", "conditional_variance_ball", "conditional_central_mu3_ball", "cubic_gap_ball", "hessian_gap_ball", "safe_conditional_Green_multiplier_ball", "cost")}, indent=2))


if __name__ == "__main__":
    main()
