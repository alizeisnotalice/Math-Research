#!/usr/bin/env python3
"""Certified finite-point alpha=3/2 stable-density kernel values via Arb."""
from __future__ import annotations

import hashlib
import json
import math
import time
from pathlib import Path

from flint import acb, arb, ctx

POINTS = (0, 1, 2, 8, 32, 128)
U = 4


def main() -> None:
    ctx.prec = 256
    output = {
        "status": "in_progress", "precision_bits": ctx.prec, "formula": "s1(x)=2/pi*integral_0^infty u*exp(-u^3)*cos(x*u^2) du",
        "substitution": "xi=u^2", "U": U,
        "tail_proof": "abs(tail)<=2/(3*pi*U)*exp(-U^3)",
        "points_requested": list(POINTS), "rows": [],
        "scope": "Finite certified 1D kernel values only; not a g, I, or C5 result.",
    }
    started = time.monotonic()
    tail = arb(2) / (3 * arb.pi() * U) * (-arb(U) ** 3).exp()
    for point in POINTS:
        began = time.monotonic()
        try:
            x = arb(point)
            integrand = lambda u, _: u * (-u ** 3).exp() * (x * u ** 2).cos()
            truncated = (2 * acb.integral(integrand, arb(0), arb(U)) / acb.pi()).real
            certified = truncated + arb(0, tail.upper())
            row = {"x": point, "truncated_integral_arb_ball": str(truncated),
                   "tail_absolute_upper_arb": str(tail), "s1_arb_ball_with_tail": str(certified),
                   "elapsed_seconds": time.monotonic() - began, "status": "success"}
        except Exception as error:
            row = {"x": point, "status": "failed", "error_type": type(error).__name__,
                   "error": str(error), "elapsed_seconds": time.monotonic() - began}
        output["rows"].append(row)
        output["elapsed_seconds_at_checkpoint"] = time.monotonic() - started
        output["source_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        Path("010_arb_stable_density.json").write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
        if row["status"] == "failed":
            break
    peak = 2 * (arb(2) / 3).gamma() / (3 * arb.pi())
    zero = next((row for row in output["rows"] if row["x"] == 0 and row["status"] == "success"), None)
    output["gamma_peak_formula_arb"] = str(peak)
    output["gamma_peak_contained_by_x0_tail_ball"] = None if zero is None else bool(peak in arb(zero["s1_arb_ball_with_tail"]))
    output["status"] = "complete" if len(output["rows"]) == len(POINTS) else "failed_at_first_unfinished_point"
    output["elapsed_seconds"] = time.monotonic() - started
    output["source_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    Path("010_arb_stable_density.json").write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
