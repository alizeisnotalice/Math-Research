#!/usr/bin/env python3
"""One pre-registered nested positive 1D quadrature diagnostic, not a certificate."""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import platform
import time
import warnings
from pathlib import Path
import scipy
from scipy.integrate import IntegrationWarning, quad
from flint import arb, ctx

PARAMETERS = {"beta": "3/4", "a": "15/224", "gamma": "9/28", "dimension": 16,
              "inner_logu": [-64, 6], "outer_loglambda": [-40, "log(16)"],
              "outer_breakpoints": [-40, -20, -10, -5, 0, "log(16)"],
              "inner_epsabs": 2e-14, "inner_epsrel": 2e-11, "inner_limit": 200,
              "outer_epsabs": 5e-10, "outer_epsrel": 2e-9, "outer_limit": 160,
              "max_distinct_inner_integrals": 2000}


def show(value):
    return value.str(70, more=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--round011-results", type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    source_raw = args.round011_results.read_bytes()
    source = json.loads(source_raw)
    beta, a, gamma = 3 / 4, 15 / 224, 9 / 28
    Gamma_a = math.gamma(a)
    m0 = math.gamma(1 + a / beta) / math.gamma(1 + a)
    prefactor = 16 ** (1 - 15 / 14) * math.gamma(1 + beta)
    positive_moment_factor = gamma / math.gamma(1 - gamma)
    q_factor = prefactor * positive_moment_factor * m0 ** 16
    cache = {}
    error_messages = []
    max_inner_estimated_error = 0.0
    inner_function_evaluations = 0
    clamps = 0

    def outer(v):
        nonlocal max_inner_estimated_error, inner_function_evaluations, clamps
        if v in cache:
            return cache[v][0]
        if len(cache) >= PARAMETERS["max_distinct_inner_integrals"]:
            raise RuntimeError("pre-registered inner-integral budget exhausted")
        lam = math.exp(v)

        def inner(w):
            nonlocal inner_function_evaluations
            inner_function_evaluations += 1
            u = math.exp(w)
            u_beta = math.exp(beta * w)
            difference = u_beta * math.expm1(beta * math.log1p(lam / u))
            return math.exp(a * w - u_beta) * (-math.expm1(-difference))

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", IntegrationWarning)
            integral, estimate = quad(inner, -64, 6, epsabs=2e-14, epsrel=2e-11, limit=200)
        error_messages.extend(str(item.message) for item in caught)
        left_proxy = (-math.expm1(-(lam ** beta))) * math.exp(-64 * a) / (a * Gamma_a)
        d = integral / Gamma_a + left_proxy
        estimated_error = estimate / Gamma_a
        max_inner_estimated_error = max(max_inner_estimated_error, estimated_error)
        r = d / m0
        if r >= 1:
            if r > 1 + 1e-10:
                raise ArithmeticError("normalized positive difference exceeds one materially")
            phi = 1.0
            clamps += 1
        elif r < 0:
            raise ArithmeticError("negative difference in positive representation")
        else:
            phi = -math.expm1(16 * math.log1p(-r))
        value = phi * math.exp(-gamma * v)
        cache[v] = (value, estimated_error, d)
        return value

    pieces = []
    breaks = [-40, -20, -10, -5, 0, math.log(16)]
    for left, right in zip(breaks, breaks[1:]):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", IntegrationWarning)
            value, error = quad(outer, left, right, epsabs=5e-10, epsrel=2e-9, limit=160)
        error_messages.extend(str(item.message) for item in caught)
        pieces.append({"loglambda_left": left, "loglambda_right": right,
                       "normalized_outer_integral": value, "quadpack_error_estimate": error})
    finite = sum(piece["normalized_outer_integral"] for piece in pieces)
    large_tail_main = 16 ** (-gamma) / gamma
    q_estimate = q_factor * (finite + large_tail_main)
    c = float(arb(source["arb_256bit"]["shape_prefactor_c"]))
    I_factor = c ** (10 / 7) * math.pi ** 2 / 60
    baseline_Q = float(arb(source["arb_256bit"]["negative_determinant_moment_upper_Q"]))
    target_Q = float(arb(source["arb_256bit"]["target_sufficient_negative_moment_upper_gate"]))

    # These analytic tail/transfer formula evaluations are certified Arb balls.
    # The finite nested integral above remains floating-point and uncertified.
    ctx.prec = 256
    bb, aa, gg = arb(3) / 4, arb(15) / 224, arb(9) / 28
    mm = (1 + aa / bb).gamma() / (1 + aa).gamma()
    pp = arb(16) ** (1 - arb(15) / 14) * (1 + bb).gamma()
    factor = pp * gg / (1 - gg).gamma() * mm ** 16
    delta, L = arb(-40).exp(), arb(16)
    u0, upper_u = arb(-64).exp(), arb(6).exp()
    # d - (left_proxy + finite inner integral) lies in [-EL, EU].
    EL = u0 ** (aa + bb) / ((aa + bb) * aa.gamma())
    t0, kappa = upper_u ** bb, aa / bb
    EU = t0 ** (kappa - 1) * (-t0).exp() / (bb * aa.gamma())
    lambda_weight = (delta ** (-gg) - L ** (-gg)) / gg
    inner_truncation_Q_error = factor * 16 / mm * (EL + EU) * lambda_weight
    ES_beta_minus_a = (aa / bb).gamma() / (1 - bb + aa).gamma()
    small_outer_Q_upper = (pp * gg / (1 - gg).gamma() * 16 * mm ** 15
                           * ES_beta_minus_a * delta ** (bb - gg) / (bb - gg))
    concavity_c = arb(2) ** (bb - 1)
    ratio_L_upper = concavity_c ** (-aa / bb) * (-concavity_c * L ** bb).exp()
    large_outer_Q_correction_upper = factor * ratio_L_upper ** 16 * L ** (-gg) / gg
    heuristic_inner_error_Q_proxy = q_factor * 16 / m0 * max_inner_estimated_error * float(lambda_weight)
    payload = {
        "status": "finite_nested_quadrature_diagnostic_NOT_a_rigorous_Q_or_I_certificate",
        "parameters": PARAMETERS, "python": platform.python_version(), "scipy": scipy.__version__,
        "source_json_sha256": hashlib.sha256(source_raw).hexdigest(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "diagnostic": {"m0": m0, "q_prefactor_including_m0_power": q_factor,
                       "normalized_finite_outer_integral": finite,
                       "normalized_large_lambda_main_tail": large_tail_main,
                       "coupled_Q_surrogate_estimate": q_estimate,
                       "coupled_I_bound_surrogate_estimate": I_factor * q_estimate,
                       "imported_separated_Q_baseline": baseline_Q,
                       "ratio_coupled_to_separated": q_estimate / baseline_Q,
                       "imported_target_sufficient_Q_gate": target_Q,
                       "coupled_Q_over_target_gate": q_estimate / target_Q,
                       "I_over_new_orbit_bound_as_reported_by_parent": I_factor * q_estimate / 1.5792875288e-5},
        "finite_outer_pieces": pieces,
        "analytic_error_interface_arb256": {
            "inner_left_proxy_overestimate_bound_EL": show(EL),
            "inner_right_omitted_bound_EU": show(EU),
            "induced_Q_error_from_inner_analytic_truncation": show(inner_truncation_Q_error),
            "small_lambda_omitted_Q_upper_bound": show(small_outer_Q_upper),
            "large_lambda_main_tail_overestimate_Q_bound": show(large_outer_Q_correction_upper),
            "interpretation": "If finite nested integral were evaluated exactly, add the small-lambda upper bound, allow the inner-truncation bound in both directions, and subtract at most the large-tail correction. These do not bound floating quadrature error."},
        "floating_error_diagnostics_NOT_certified": {
            "outer_quadpack_error_sum_in_Q_units": q_factor * sum(piece["quadpack_error_estimate"] for piece in pieces),
            "maximum_sampled_inner_quadpack_error_estimate": max_inner_estimated_error,
            "crude_Q_proxy_from_maximum_sampled_inner_error": heuristic_inner_error_Q_proxy,
            "warnings": error_messages, "clamps_to_phi_one": clamps,
            "interpretation": "QUADPACK estimates and sampled maxima are diagnostics, not rigorous enclosures or supremum bounds."},
        "cost": {"distinct_inner_integrals": len(cache), "inner_function_evaluations": inner_function_evaluations,
                 "outer_subintervals": len(pieces), "random_samples": 0,
                 "elapsed_seconds": time.monotonic() - started},
        "scope": "Two nested positive one-dimensional integrals retain T/S coupling in the AM-GM surrogate; actual determinant moment may be smaller. No new certified global I, g or C5."}
    output = Path(__file__).with_name("coupled_diagnostic_results.json")
    output.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({key: payload[key] for key in ("status", "diagnostic", "analytic_error_interface_arb256", "floating_error_diagnostics_NOT_certified", "cost")}, indent=2))


if __name__ == "__main__":
    main()
