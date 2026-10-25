#!/usr/bin/env python3
"""Single exact Fraction/integer audit; no dual engine or Gamma evaluation."""
from __future__ import annotations
import argparse
import hashlib
import json
import platform
import re
import time
from fractions import Fraction as F
from pathlib import Path


def parse_decimal_ball(text):
    match = re.fullmatch(r"\[([^\s]+) \+/- ([^\s]+)\]", text)
    if match is None:
        raise ValueError("unsupported displayed-ball syntax")
    midpoint, radius = map(F, match.groups())
    assert radius >= 0
    return midpoint - radius, midpoint + radius


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--round011-results", type=Path, required=True)
    parser.add_argument("--orbit-results", type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    raw_inputs = {"Q_source": args.round011_results.read_bytes(),
                  "J_source": args.orbit_results.read_bytes()}
    Q_json = json.loads(raw_inputs["Q_source"])
    J_json = json.loads(raw_inputs["J_source"])
    delta, m, s = F(1, 4000), F(648, 1463), F(126404864, 474619145)
    t = (s - m * delta) / (m - delta)
    p = (m - delta) / (t - delta)
    assert 0 < delta < t < 1
    assert 0 < p < 1
    assert (1 - p) * delta + p * t == m
    assert (1 - p) * delta ** 2 + p * t ** 2 == s
    n, d = t.numerator, t.denominator
    # Each positive integer residual proves the displayed rational root bound.
    residuals = {
        "delta_f_upper_593_over_100": 593 ** 14 - 4000 ** 3 * 100 ** 14,
        "t_f_upper_9_over_8": 9 ** 14 * n ** 3 - 8 ** 14 * d ** 3,
        "delta_g_lower_69_over_20": 4000 ** 3 * 20 ** 20 - 69 ** 20,
        "t_g_lower_43_over_40": 40 ** 20 * d ** 3 - 43 ** 20 * n ** 3,
    }
    assert all(value > 0 for value in residuals.values())
    f_upper = (1 - p) * F(593, 100) + p * F(9, 8)
    g_lower = (1 - p) * F(69, 20) + p * F(43, 40)
    assert f_upper < F(12, 5)
    assert g_lower > F(17, 10)
    Q_lower, Q_upper = parse_decimal_ball(Q_json["arb_256bit"]["negative_determinant_moment_upper_Q"])
    J_lower, J_upper = parse_decimal_ball(J_json["J_upper_ball"])
    assert Q_lower > F(12, 5)
    assert J_upper < F(17, 10)
    payload = {
        "status": "exact_rational_scalar_family_obstruction_passed",
        "arithmetic": "Python integers and Fraction only; no floating point, root functions, Gamma, LP or dual engine",
        "python": platform.python_version(),
        "input_sha256": {name: hashlib.sha256(raw).hexdigest() for name, raw in raw_inputs.items()},
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "rational_data": {"delta": str(delta), "m": str(m), "s": str(s), "t": str(t),
                          "p_high_atom": str(p), "1_minus_p_low_atom": str(1 - p)},
        "exact_positive_support_gaps": {"delta": str(delta), "t_minus_delta": str(t - delta),
                                        "1_minus_t": str(1 - t), "p": str(p), "1_minus_p": str(1 - p)},
        "moments": {"first_moment_equals_m": True, "second_moment_equals_s": True,
                    "first_residual": "0", "second_residual": "0"},
        "positive_integer_power_residuals": {name: str(value) for name, value in residuals.items()},
        "negative_power_moment_witnesses": {
            "f_exponent": "3/14", "f_upper_rational": str(f_upper),
            "12_over_5_minus_f_upper": str(F(12, 5) - f_upper),
            "g_exponent": "3/20", "g_lower_rational": str(g_lower),
            "g_lower_minus_17_over_10": str(g_lower - F(17, 10))},
        "source_ball_comparisons": {
            "Q_lower_exact_decimal_rational": str(Q_lower), "Q_lower_minus_12_over_5": str(Q_lower - F(12, 5)),
            "J_upper_exact_decimal_rational": str(J_upper), "17_over_10_minus_J_upper": str(F(17, 10) - J_upper),
            "Q_greater_than_12_over_5": True, "J_less_than_17_over_10": True},
        "family": "For every alpha>=0 and real A,B,C whose pointwise envelope d^(-3/20)<=alpha*d^(-3/14)+A+B*d+C*d^2 holds on all d in (0,1], alpha*Qstar+A+B*m+C*s>17/10>Jstar.",
        "scope": "Abstract two-atom probability witness for this scalar envelope/moment family only; not the actual det(M) law, not an impossibility theorem for other bounds, and not a C5 result.",
        "cost": {"finite_batches": 1, "integer_power_comparisons": 4, "root_evaluations": 0,
                 "Gamma_evaluations": 0, "LP_calls": 0, "random_samples": 0,
                 "elapsed_seconds": time.monotonic() - started}}
    output = Path(__file__).with_name("exact_scalar_results.json")
    output.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({key: payload[key] for key in ("status", "rational_data", "negative_power_moment_witnesses", "cost")}, indent=2))


if __name__ == "__main__":
    main()
