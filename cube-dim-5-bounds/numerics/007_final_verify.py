#!/usr/bin/env python3
"""Add Arb analytic error/tail to the 007 bit-pattern heat quadrature."""
from __future__ import annotations

import hashlib
import json
import math
import struct
from fractions import Fraction
from pathlib import Path

from flint import arb, ctx

from body_fan_arb_anchor import exact_terms
from body_fan_green_arb import TAILS, point_from_tail
from body_fan_h_epsilon_arb import truncated_green


def bits_to_float(text: str) -> float:
    return struct.unpack(">d", int(text, 16).to_bytes(8, "big"))[0]


def as_bits(value: float) -> str:
    return f"0x{struct.unpack('>Q', struct.pack('>d', value))[0]:016x}"


def lower_binary64(value: arb) -> float:
    lo = float(value)
    while arb(lo) > value:
        lo = math.nextafter(lo, -math.inf)
    while arb(math.nextafter(lo, math.inf)) <= value:
        lo = math.nextafter(lo, math.inf)
    return lo


def upper_binary64(value: arb) -> float:
    hi = float(value)
    while arb(hi) < value:
        hi = math.nextafter(hi, math.inf)
    while arb(math.nextafter(hi, -math.inf)) >= value:
        hi = math.nextafter(hi, -math.inf)
    return hi


def enclosure(value: arb) -> tuple[float, float]:
    lower, upper = value.lower(), value.upper()
    lo = float(lower)
    while arb(lo) > lower:
        lo = math.nextafter(lo, -math.inf)
    while arb(math.nextafter(lo, math.inf)) <= lower:
        lo = math.nextafter(lo, math.inf)
    hi = float(upper)
    while arb(hi) < upper:
        hi = math.nextafter(hi, math.inf)
    while arb(math.nextafter(hi, -math.inf)) >= upper:
        hi = math.nextafter(hi, -math.inf)
    return lo, hi


def main() -> None:
    ctx.prec = 256
    heat = json.loads(Path("007_heat_quadrature.json").read_text(encoding="utf-8"))
    inputs = json.loads(Path("007_fejer_inputs.json").read_text(encoding="utf-8"))
    assert heat["input_nodes_sha256"] == hashlib.sha256(Path("007_fejer_nodes.bits").read_bytes()).hexdigest()
    assert inputs["nodes_bits_sha256"] == heat["input_nodes_sha256"]
    rho = arb(49) / 10
    cpeak = arb(35509) / (arb(764411904) * arb.pi() ** 5)
    interpolation = arb(0)
    for k in range(7):
        a = arb(1) / (2 ** (7 - k))
        sigma = a * (arb(3) / 2 - (rho + 1 / rho) / 4)
        interpolation += a * 4 * cpeak * sigma ** -5 * rho ** (-24) / (rho - 1)
    tail = cpeak / 4
    q = heat["quadrature"]
    qlo, qhi = bits_to_float(q["lo_bits"]), bits_to_float(q["hi_bits"])
    # Interpolation is two-sided; p_t >= 0 yields the one-sided [0,tail] tail.
    lower_expression = arb(qlo) - interpolation.upper()
    upper_expression = arb(qhi) + interpolation.upper() + tail.upper()
    h_lo = lower_binary64(lower_expression.lower())
    h_hi = upper_binary64(upper_expression.upper())
    direct_complex, _ = truncated_green(point_from_tail(TAILS[0]), Fraction(1, 128), tuple(exact_terms()))
    direct = direct_complex.real
    contains = bool(direct.lower() >= arb(h_lo) and direct.upper() <= arb(h_hi))
    result = {
        "status": "complete_onepoint_H_epsilon_interval_certificate_conditional_on_Cpeak",
        "precision_bits": 256,
        "point_tail_numerators": list(TAILS[0]),
        "epsilon": "1/128", "T": 1, "degree": 24, "rho": "49/10", "nodes": 175,
        "quadrature_bits": q,
        "Cpeak_arb": str(cpeak),
        "interpolation_error_arb": str(interpolation),
        "tail_interval_arb": f"[0, {tail.upper()}]",
        "H_interval_bits": {"lo_bits": as_bits(h_lo), "hi_bits": as_bits(h_hi)},
        "H_interval_decimal_display": [h_lo, h_hi],
        "direct_H_1_128_arb_ball": str(direct),
        "direct_arb_ball_contained": contains,
        "input_nodes_sha256": heat["input_nodes_sha256"],
        "heat_output_sha256": hashlib.sha256(Path("007_heat_quadrature.json").read_bytes()).hexdigest(),
        "node_generator_source_sha256": hashlib.sha256(Path("007_generate_fejer_inputs.py").read_bytes()).hexdigest(),
        "final_verifier_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "scope": "A single-point interval certificate conditional on the stated complex-time Cpeak inequality and the audited 005b binary64 kernel assumptions. It proves no spatial uniformity, g, I, or C5 bound.",
    }
    Path("007_final_verify.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
