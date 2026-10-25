#!/usr/bin/env python3
"""Arb cross-check of the C++ double heat-kernel batch at locked points."""
from __future__ import annotations

import json
import time
from fractions import Fraction
from pathlib import Path

from flint import acb, arb, ctx, fmpq

from body_fan_arb_anchor import exact_terms
from body_fan_green_arb import DENOMINATOR, TAILS


def ball(value: Fraction) -> arb:
    return arb(fmpq(value.numerator, value.denominator))


def heat(point, time_value, terms):
    total = acb(0)
    for rays, determinant, rates in terms:
        term = acb(determinant)
        for ray, rate in zip(rays, rates):
            dot = sum(Fraction(coordinate) * value for coordinate, value in zip(ray, point))
            term /= acb(ball(time_value * rate), -ball(dot))
        total += term
    return (total / (2 * acb.pi()) ** 5).real


def main() -> None:
    ctx.prec = 256
    cpp = json.loads(Path("004_cpp_results.json").read_text(encoding="utf-8"))
    times = tuple(Fraction(1, value) for value in (128, 32, 8, 2, 1))
    terms = tuple(exact_terms())
    started = time.monotonic()
    rows = []
    worst_abs = 0.0
    worst_rel = 0.0
    worst_cancel = 0.0
    for index, tail in enumerate(TAILS):
        point = (Fraction(1),) + tuple(Fraction(value, DENOMINATOR) for value in tail)
        values = []
        for position, time_value in enumerate(times):
            reference = heat(point, time_value, terms)
            reference_float = float(reference)
            double_value = cpp["locked_values"][index]["values"][position]["p"]
            abs_error = abs(double_value - reference_float)
            rel_error = abs_error / abs(reference_float)
            cancellation = cpp["locked_values"][index]["values"][position]["abs_sum"] / abs(double_value)
            worst_abs = max(worst_abs, abs_error)
            worst_rel = max(worst_rel, rel_error)
            worst_cancel = max(worst_cancel, cancellation)
            values.append({"time": str(time_value), "arb_ball": str(reference), "arb_midpoint_float": reference_float, "double": double_value, "absolute_error": abs_error, "relative_error": rel_error, "cancellation_factor": cancellation})
        rows.append({"tail_numerators": list(tail), "values": values})
    result = {"status": "double_batch_crosschecked_against_arb", "precision_bits": 256, "points": 8, "times": [str(value) for value in times], "worst_absolute_error": worst_abs, "worst_relative_error": worst_rel, "worst_cancellation_factor": worst_cancel, "elapsed_seconds": time.monotonic() - started, "rows": rows, "scope": "Agreement is a finite numerical cross-check only; double output is not a proof certificate."}
    Path("004_crosscheck.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in ("points", "times", "worst_absolute_error", "worst_relative_error", "worst_cancellation_factor", "elapsed_seconds")}, indent=2))


if __name__ == "__main__":
    main()
