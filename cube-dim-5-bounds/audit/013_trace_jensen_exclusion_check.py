#!/usr/bin/env python3
"""Exact rational exclusion check for the trace-Jensen Gaussian interface."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json


def atan_partial(x: Q, last: int) -> Q:
    return sum(((-1) ** k) * x ** (2*k + 1) / (2*k + 1) for k in range(last + 1))


def main():
    # Machin: pi/4=4 atan(1/5)-atan(1/239).  Alternating-series remainders
    # give atan(1/5)>S_3 and atan(1/239)<S_0, hence this strict pi lower bound.
    pi_lower = 4 * (4 * atan_partial(Q(1, 5), 3) - atan_partial(Q(1, 239), 0))
    assert pi_lower > Q(333, 106)
    peak_file = Path(__file__).with_name("004_body_fan_independent_results.json")
    peak = json.loads(peak_file.read_text())
    assert peak["peak_coefficient"] == "35509/764411904"
    P_rational, target = Q(35509, 764411904), Q(11, 100000)
    # If P=p_1(0)=P_rational*pi^-5, then the trace-Jensen upper bound U obeys
    # U^5=(2/(25(2pi)^2))^5 P.  Replace pi by a certified lower bound.
    upper_fifth = Q(2, 25) ** 5 * P_rational / ((2 * pi_lower) ** 10 * pi_lower ** 5)
    assert upper_fifth < target ** 5
    result = {
        "status": "PASS_exact_trace_Jensen_exclusion",
        "peak": "35509/(764411904*pi^5)",
        "machin_pi_lower": str(pi_lower),
        "upper_bound_fifth_power": str(upper_fifth),
        "target_fifth_power": str(target ** 5),
        "strictly_below_11_over_100000": True,
        "numeric_upper_from_rational_pi_lower_display_only": float(upper_fifth) ** 0.2,
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "peak_source_sha256": hashlib.sha256(peak_file.read_bytes()).hexdigest(),
    }
    output = Path(__file__).with_name("013_trace_jensen_exclusion_results.json")
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
