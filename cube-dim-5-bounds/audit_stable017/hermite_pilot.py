#!/usr/bin/env python3
"""Certified Hermite pilot on exactly six registered one-dimensional intervals.

Requires python-flint; reads locked scalar JSON without importing source code.
python hermite_pilot.py --locked-kernel-json /path/010_arb_stable_density.json
"""
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

INTERVALS = (("0", "0.01"), ("1", "1.01"), ("2", "2.01"),
             ("8", "8.01"), ("32", "32.01"), ("127.99", "128"))


def as_ball(rational):
    rational = Fraction(rational)
    return arb(rational.numerator) / rational.denominator


def dump_ball(value):
    # Force ample digits to avoid needlessly inflating radii on text roundtrips.
    return value.str(80, more=True)


def evaluate_certificate(record, x_left, x_right=None):
    """Return a ball enclosing s(x) for every x in the exact rational input set.

    Coefficients include all endpoint quadrature/tail uncertainties. The one
    final symmetric remainder ball is the analytic Hermite interpolation error.
    """
    a, b = Fraction(record["a"]), Fraction(record["b"])
    left = Fraction(x_left)
    right = left if x_right is None else Fraction(x_right)
    if not a <= left <= right <= b:
        raise ValueError("requested input interval is outside the certified domain")
    x = as_ball(left).union(as_ball(right))
    t = (x - as_ball(a)) / as_ball(b - a)
    c0, c1, c2, c3 = map(arb, record["coefficient_balls_c0_to_c3"])
    polynomial = ((c3 * t + c2) * t + c1) * t + c0
    error = arb(record["hermite_remainder_absolute_bound"])
    return polynomial + arb(0, error.upper())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--locked-kernel-json", type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    locked_raw = args.locked_kernel_json.read_bytes()
    locked = json.loads(locked_raw)
    assert locked["status"] == "complete"
    assert [row["x"] for row in locked["rows"]] == [0, 1, 2, 8, 32, 128]
    source_rows = {Fraction(row["x"]): row for row in locked["rows"]}
    ctx.prec = 384
    U = arb(4)
    density_tail = 2 / (3 * arb.pi() * U) * (-U ** 3).exp()
    derivative_tail = (2 * U / (3 * arb.pi() * (1 - 1 / (3 * U ** 3)))
                       * (-U ** 3).exp())
    M4 = (arb(10) / 3).gamma() / (arb(3) / 2 * arb.pi())
    integral_counts = {"new_endpoint_densities": 0, "endpoint_derivatives": 0,
                       "midpoint_densities_for_check_only": 0}

    def integrate(x_fraction, derivative=False):
        x = as_ball(x_fraction)
        if derivative:
            f = lambda u, analytic: -(u ** 3) * (-u ** 3).exp() * (x * u ** 2).sin()
            tail = derivative_tail
        else:
            f = lambda u, analytic: u * (-u ** 3).exp() * (x * u ** 2).cos()
            tail = density_tail
        truncated = (2 / acb.pi() * sum((acb.integral(f, k, k + 1)
                                        for k in range(4)), acb(0))).real
        assert truncated.is_finite()
        return truncated + arb(0, tail.upper())

    nodes = {}
    for a, b in INTERVALS:
        for raw_x in (a, b):
            x_fraction = Fraction(raw_x)
            if x_fraction in nodes:
                continue
            if x_fraction in source_rows:
                density = arb(source_rows[x_fraction]["s1_arb_ball_with_tail"])
                origin = "locked_scalar_ball"
            else:
                density = integrate(x_fraction)
                integral_counts["new_endpoint_densities"] += 1
                origin = "new_384_bit_subdivided_Arb_integral"
            if x_fraction == 0:
                derivative = arb(0)
                derivative_origin = "exact_symmetric_derivative_zero"
            else:
                derivative = integrate(x_fraction, derivative=True)
                integral_counts["endpoint_derivatives"] += 1
                derivative_origin = "new_384_bit_subdivided_Arb_integral"
            nodes[x_fraction] = {"x": raw_x, "density_ball": dump_ball(density),
                                 "density_origin": origin,
                                 "derivative_ball": dump_ball(derivative),
                                 "derivative_origin": derivative_origin}
    assert len(nodes) == 12
    records = []
    for raw_a, raw_b in INTERVALS:
        a, b = Fraction(raw_a), Fraction(raw_b)
        h = as_ball(b - a)
        fa, da = map(arb, (nodes[a]["density_ball"], nodes[a]["derivative_ball"]))
        fb, db = map(arb, (nodes[b]["density_ball"], nodes[b]["derivative_ball"]))
        coefficients = (fa, h * da,
                        3 * (fb - fa) - h * (2 * da + db),
                        2 * (fa - fb) + h * (da + db))
        error = M4 * h ** 4 / 384
        rec = {"a": raw_a, "b": raw_b, "h_exact": str(b - a),
               "coordinate": "t=(x-a)/h; H(t)=c0+c1*t+c2*t^2+c3*t^3",
               "coefficient_balls_c0_to_c3": list(map(dump_ball, coefficients)),
               "hermite_remainder_absolute_bound": dump_ball(error)}
        # Convex-hull property of the cubic Bernstein basis gives a tighter
        # whole-domain enclosure than dependent interval Horner evaluation.
        controls = (fa, fa + h * da / 3, fb - h * db / 3, fb)
        hull = controls[0]
        for control in controls[1:]:
            hull = hull.union(control)
        rec["bernstein_control_balls"] = list(map(dump_ball, controls))
        rec["whole_domain_density_ball"] = dump_ball(hull + arb(0, error.upper()))
        check_a = evaluate_certificate(rec, a)
        check_b = evaluate_certificate(rec, b)
        assert fa in check_a and fb in check_b
        c0, c1, c2, c3 = map(arb, rec["coefficient_balls_c0_to_c3"])
        slope_a, slope_b = c1 / h, (c1 + 2 * c2 + 3 * c3) / h
        assert da in slope_a and db in slope_b
        mid = (a + b) / 2
        direct_mid = integrate(mid)
        integral_counts["midpoint_densities_for_check_only"] += 1
        interpolated_mid = evaluate_certificate(rec, mid)
        assert direct_mid in interpolated_mid
        polynomial_mid = c0 + c1 / 2 + c2 / 4 + c3 / 8
        midpoint_residual = direct_mid - polynomial_mid
        assert abs(midpoint_residual) < error
        rec["checks"] = {"endpoint_density_balls_contained": True,
                         "endpoint_derivative_balls_contained": True,
                         "midpoint_exact": str(mid),
                         "midpoint_direct_density_ball": dump_ball(direct_mid),
                         "midpoint_H_plus_error_ball": dump_ball(interpolated_mid),
                         "midpoint_residual_ball": dump_ball(midpoint_residual),
                         "midpoint_direct_ball_contained": True,
                         "midpoint_residual_within_M4_bound": True}
        records.append(rec)
    payload = {"status": "passed", "precision_bits": ctx.prec,
               "python": platform.python_version(), "python_flint": flint.__version__,
               "locked_input_sha256": hashlib.sha256(locked_raw).hexdigest(),
               "locked_kernel_source_sha256": locked["source_sha256"],
               "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               "integration_split": [0, 1, 2, 3, 4], "U": 4,
               "density_tail_absolute_bound": dump_ball(density_tail),
               "derivative_tail_absolute_bound": dump_ball(derivative_tail),
               "global_M4_absolute_bound": dump_ball(M4),
               "distinct_endpoint_nodes": len(nodes),
               "reused_locked_density_nodes": 6,
               "actual_integral_counts": integral_counts,
               "nodes": list(nodes.values()), "intervals": records,
               "scope": "Six continuous one-dimensional intervals only. Midpoint checks do not generate interpolation nodes. No global 5D g/I/C5 certificate.",
               "independence_limit": "Fresh Hermite construction and new derivative integrals, but same Arb backend and u-substitution as locked scalar calculation.",
               "elapsed_seconds": time.monotonic() - started}
    output = Path(__file__).with_name("hermite_results.json")
    output.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({key: payload[key] for key in ("status", "actual_integral_counts", "distinct_endpoint_nodes",
                                                  "global_M4_absolute_bound", "derivative_tail_absolute_bound", "elapsed_seconds")}, indent=2))


if __name__ == "__main__":
    main()
