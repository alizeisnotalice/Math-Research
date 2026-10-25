#!/usr/bin/env python3
"""One registered positive 1D integral pilot; no sampling-based claims."""
from __future__ import annotations
import argparse
import hashlib
import json
import platform
import time
from fractions import Fraction
from pathlib import Path
import flint
from flint import acb, arb, ctx


def ball(q):
    q = Fraction(q)
    return arb(q.numerator) / q.denominator


def show(x):
    return x.str(90, more=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--hermite-json", type=Path, required=True)
    parser.add_argument("--middle-budget-json", type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    inputs = {"hermite": args.hermite_json.read_bytes(), "middle_budget": args.middle_budget_json.read_bytes()}
    hermite = json.loads(inputs["hermite"])
    budget = json.loads(inputs["middle_budget"])
    assert hermite["status"] == "passed"
    ctx.prec = 384
    A = (arb(2) / 3).gamma() * 2 / (3 * arb.pi())
    K = arb(25) / 6 * (arb(10) / (3 * arb.pi())).sqrt()
    Cf = K * arb(16) ** (arb(5) / 2) * A ** 4 / 3
    B = arb(5) / 6 * Cf ** (arb(3) / 10) * K ** (arb(7) / 10)
    epsilon = ball(Fraction(1, 65536))
    T = arb(4)
    short = B * epsilon ** (arb(6) / 5)
    long = arb(3) / 7 * A ** 5 / 48 * T ** (-arb(7) / 3)
    middle_norm_gate = ball(Fraction(11, 1000000)) ** (arb(7) / 10) - short - long
    middle_integral_gate = middle_norm_gate ** (arb(10) / 7)
    source_gate = next(row for row in budget["rows"] if row["epsilon"] == "1/65536" and row["T"] == 4)
    assert middle_integral_gate in arb(source_gate["middle_integral_limit_ball"])
    downward_gate = "1.06143562300460297428169091520606831886653080773414323264992525041314912e-5"
    assert ball(downward_gate) < middle_integral_gate
    U, R, N = arb(4), Fraction(12), 32
    h = R / N
    scalar_tail = 2 / (3 * arb.pi() * U) * (-U ** 3).exp()
    reused, new = [], []
    cache = {}

    def density(x):
        if x in cache:
            return cache[x]
        for segment in hermite["intervals"]:
            a, b = Fraction(segment["a"]), Fraction(segment["b"])
            if a <= x <= b:
                t = ball(x - a) / ball(b - a)
                c0, c1, c2, c3 = map(arb, segment["coefficient_balls_c0_to_c3"])
                result = ((c3 * t + c2) * t + c1) * t + c0
                result += arb(0, arb(segment["hermite_remainder_absolute_bound"]).upper())
                cache[x] = result
                reused.append(str(x))
                return result
        xb = ball(x)
        integrand = lambda u, analytic: u * (-u ** 3).exp() * (xb * u ** 2).cos()
        truncated = (2 / acb.pi() * sum((acb.integral(integrand, k, k + 1) for k in range(4)), acb(0))).real
        assert truncated.is_finite()
        cache[x] = truncated + arb(0, scalar_tail.upper())
        new.append(str(x))
        return cache[x]

    values = []
    for k in range(N + 1):
        z = k * h
        values.append(density(z) * density(z / 3) ** 5)
    lower_rectangle = ball(h) / 24 * sum(values[1:], arb(0))
    upper_rectangle = ball(h) / 24 * sum(values[:-1], arb(0))
    tail = K ** 6 * arb(3) ** (arb(25) / 2) / (24 * 14 * ball(R) ** 14)
    c6_lower = lower_rectangle.lower()
    c6_upper = (upper_rectangle + tail).upper()
    c6_enclosure = c6_lower.union(c6_upper)
    width = c6_upper - c6_lower
    monotone_mesh_gap_bound = ball(h) * A ** 6 / 24
    baseline = A ** 5 / 48
    late_norm = c6_upper * arb(3) / 7 * (1 - arb(4) ** (-arb(7) / 3))
    late_integral_bound = late_norm ** (arb(10) / 7)
    early_norm_gate = middle_norm_gate - late_norm
    assert early_norm_gate > 0
    early_integral_gate = early_norm_gate ** (arb(10) / 7)
    direct_whole_norm = c6_upper * arb(3) / 7 * (epsilon ** (-arb(7) / 3) - T ** (-arb(7) / 3))
    payload = {
        "status": "finite_positive_integral_and_late_time_interface_certified_not_full_I",
        "python": platform.python_version(), "python_flint": flint.__version__, "precision_bits": ctx.prec,
        "input_sha256": {name: hashlib.sha256(raw).hexdigest() for name, raw in inputs.items()},
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "target": {"epsilon": "1/65536", "T": 4, "I_total": "11/1000000", "q": "10/7",
                   "middle_Lq_limit_ball": show(middle_norm_gate),
                   "middle_integral_limit_ball": show(middle_integral_gate),
                   "usable_downward_integral_limit_exact_decimal": downward_gate,
                   "independent_gate_contained_by_main_ball": True},
        "registered_trial": {"spatial_face_cell": "[0,1]^4", "time_slab": "[1,4]", "R": str(R), "N": N,
                             "h": str(h), "integrand": "s(z)*s(z/3)^5", "factor": "1/24"},
        "kernel_table": [{"x_exact": str(x), "density_ball": show(value),
                          "origin": "reused_continuous_Hermite_certificate" if str(x) in reused else "new_Arb_scalar_integral"}
                         for x, value in sorted(cache.items())],
        "rectangle_lower_ball_before_tail": show(lower_rectangle),
        "rectangle_upper_ball_before_tail": show(upper_rectangle),
        "analytic_tail_upper_ball": show(tail), "C6_enclosure_ball": show(c6_enclosure),
        "C6_enclosure_width_ball": show(width), "mesh_gap_analytic_upper_ball": show(monotone_mesh_gap_bound),
        "baseline_A5_over_48_ball": show(baseline), "C6_upper_over_baseline_ball": show(c6_upper / baseline),
        "late_time_Lq_upper_ball": show(late_norm), "late_time_integral_q_upper_ball": show(late_integral_bound),
        "remaining_early_time_interval": "[1/65536,1]",
        "required_early_Lq_limit_ball": show(early_norm_gate),
        "required_early_integral_limit_ball": show(early_integral_gate),
        "failure_diagnostic": {"direct_peak_bound_on_full_middle_Lq_ball": show(direct_whole_norm),
                               "ratio_to_middle_norm_budget_ball": show(direct_whole_norm / middle_norm_gate),
                               "uncovered": "No certified early-time spatial integral. A peak bound loses off-origin small-time decay."},
        "cost": {"unique_density_arguments": len(cache), "reused_Hermite_arguments": reused,
                 "new_density_integrals": len(new), "acb_integral_subcalls": 4 * len(new),
                 "elapsed_seconds": time.monotonic() - started},
        "scope": "One positive one-dimensional peak integral, one full-face cell and one late-middle time slab. No finite samples promoted to full I or C5."
    }
    output = Path(__file__).with_name("positive_peak_results.json")
    output.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({key: payload[key] for key in ("status", "C6_enclosure_ball", "C6_enclosure_width_ball", "analytic_tail_upper_ball",
                                                 "late_time_Lq_upper_ball", "required_early_integral_limit_ball", "failure_diagnostic", "cost")}, indent=2))


if __name__ == "__main__":
    main()
