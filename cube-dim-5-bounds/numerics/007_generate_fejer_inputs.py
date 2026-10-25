#!/usr/bin/env python3
"""Generate Arb-enclosed Fejer-I time nodes and positive weights for 007."""
from __future__ import annotations

import hashlib
import json
import math
import struct
from pathlib import Path

from flint import arb, ctx

EPSILON_POWER = 7
SEGMENTS = 7
D = 24
N = D + 1
RHO_NUMERATOR, RHO_DENOMINATOR = 49, 10


def bits(value: float) -> str:
    return f"0x{struct.unpack('>Q', struct.pack('>d', value))[0]:016x}"


def binary64_enclosure(value: arb) -> tuple[float, float]:
    """Find binary64 endpoints by Arb comparisons, without decimal endpoints."""
    lower, upper = value.lower(), value.upper()
    lo = float(value)
    while arb(lo) > lower:
        lo = math.nextafter(lo, -math.inf)
    while arb(math.nextafter(lo, math.inf)) <= lower:
        lo = math.nextafter(lo, math.inf)
    hi = float(value)
    while arb(hi) < upper:
        hi = math.nextafter(hi, math.inf)
    while arb(math.nextafter(hi, -math.inf)) >= upper:
        hi = math.nextafter(hi, -math.inf)
    assert arb(lo) <= lower and arb(hi) >= upper and lo <= hi
    return lo, hi


def main() -> None:
    ctx.prec = 256
    pi = arb.pi()
    rho = arb(RHO_NUMERATOR) / RHO_DENOMINATOR
    cpeak = arb(35509) / (arb(764411904) * pi ** 5)
    error = arb(0)
    lines, records = [], []
    min_weight = None
    for segment in range(SEGMENTS):
        a = arb(1) / (2 ** (EPSILON_POWER - segment))
        sigma = a * (arb(3) / 2 - (rho + 1 / rho) / 4)
        error += a * 4 * cpeak * sigma ** -5 * rho ** (-D) / (rho - 1)
        for j in range(N):
            theta = arb(2 * j + 1) * pi / (2 * N)
            cosine = theta.cos()
            time_value = a * (arb(3) / 2 + cosine / 2)
            weight = a / N
            for k in range(1, D // 2 + 1):
                weight -= (a / N) * 2 * (2 * k * theta).cos() / (4 * k * k - 1)
            assert weight.lower() > 0
            min_weight = weight if min_weight is None or weight.lower() < min_weight.lower() else min_weight
            tlo, thi = binary64_enclosure(time_value)
            wlo, whi = binary64_enclosure(weight)
            lines.append(f"{segment} {j} {bits(tlo)} {bits(thi)} {bits(wlo)} {bits(whi)}")
            records.append({"segment": segment, "j": j, "time_bits": [bits(tlo), bits(thi)], "weight_bits": [bits(wlo), bits(whi)]})
    assert len(lines) == SEGMENTS * N
    bits_path = Path("007_fejer_nodes.bits")
    bits_path.write_text("\n".join(lines) + "\n", encoding="ascii")
    tail = cpeak / 4
    result = {
        "status": "arb_fejer_inputs_complete",
        "precision_bits": 256,
        "epsilon": "1/128", "T": 1, "segments": SEGMENTS, "degree": D, "N": N,
        "rho": "49/10", "Cpeak": str(cpeak),
        "interpolation_error_arb": str(error), "tail_error_arb": str(tail),
        "minimum_weight_arb": str(min_weight),
        "all_weights_certified_positive": True,
        "nodes_bits_sha256": hashlib.sha256(bits_path.read_bytes()).hexdigest(),
        "generator_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "records": records,
        "scope": "Arb-generated input enclosures only; heat values and final H enclosure are computed separately.",
    }
    Path("007_fejer_inputs.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in ("segments", "degree", "N", "rho", "interpolation_error_arb", "tail_error_arb", "minimum_weight_arb", "nodes_bits_sha256")}, indent=2))


if __name__ == "__main__":
    main()
