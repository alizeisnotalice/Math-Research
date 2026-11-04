#!/usr/bin/env python3
"""Pure-Fraction sign certificate for two critical points at t=11/50.

For each rational base b>0 this computes rational lo< b^(11/50) <hi
by integer 50th-root inequalities.  It then propagates those intervals through
G=(1+d)N'-N using the sign of every rational coefficient.
"""
import json
from fractions import Fraction as Q
from pathlib import Path


BITS = 320
TAU = Q(11, 50)
POINTS = [Q(17, 100), Q(18, 100), Q(32, 100), Q(33, 100)]


def integer_root(number, degree):
    """floor(number^(1/degree)) with integer arithmetic only."""
    assert number >= 0 and degree >= 1
    if number == 0 or degree == 1:
        return number
    low, high = 0, 1 << ((number.bit_length() + degree - 1) // degree)
    while low + 1 < high:
        middle = (low + high) // 2
        if middle**degree <= number:
            low = middle
        else:
            high = middle
    return low


def power_bounds(base):
    """Strict rational bounds on base^(11/50), verified by 50th powers."""
    assert base > 0
    scale = 1 << BITS
    # (base^(11/50)*scale)^50 = base^11 * scale^50.
    numerator = base.numerator**11 * scale**50
    denominator = base.denominator**11
    floor_value = integer_root(numerator // denominator, 50)
    while (floor_value + 1)**50 * denominator <= numerator:
        floor_value += 1
    while floor_value**50 * denominator > numerator:
        floor_value -= 1
    assert floor_value**50 * denominator <= numerator < (floor_value + 1)**50 * denominator
    lo, hi = Q(floor_value, scale), Q(floor_value + 1, scale)
    # Explicitly certify the displayed strict inequality after clearing denominators.
    assert lo**50 < base**11 < hi**50
    return lo, hi


def multiply_interval(coefficient, lo, hi):
    if coefficient > 0:
        return coefficient * lo, coefficient * hi
    if coefficient < 0:
        return coefficient * hi, coefficient * lo
    return Q(0), Q(0)


def certificate_at(delta):
    bases = [1 + 2*delta, 1 - 2*delta, delta, 2 - delta]
    coefficients = [
        (1 + delta) * 2*TAU / (1 + 2*delta) - 1,
        (1 + delta) * 2*TAU / (1 - 2*delta) + 1,
        (1 + delta) * TAU / delta - 1,
        -(1 + delta) * TAU / (2 - delta) - 1,
    ]
    assert all(base > 0 for base in bases)
    pieces = []
    lower = upper = Q(0)
    for base, coefficient in zip(bases, coefficients):
        lo, hi = power_bounds(base)
        part_lo, part_hi = multiply_interval(coefficient, lo, hi)
        lower += part_lo
        upper += part_hi
        pieces.append({"base": str(base), "coefficient": str(coefficient),
                       "power_lower": str(lo), "power_upper": str(hi),
                       "contribution_lower": str(part_lo),
                       "contribution_upper": str(part_hi)})
    assert lower <= upper
    return {"delta": str(delta), "G_lower": str(lower), "G_upper": str(upper),
            "sign": "+" if lower > 0 else "-" if upper < 0 else "indeterminate",
            "terms": pieces}


def main():
    rows = [certificate_at(delta) for delta in POINTS]
    signs = [row["sign"] for row in rows]
    assert signs == ["+", "-", "-", "+"], signs
    result = {
        "status": "PASS", "method": "Fraction plus integer 50th-root inequalities",
        "tau": str(TAU), "bits": BITS, "rows": rows,
        "consequence": "G changes sign on (17/100,18/100) and (32/100,33/100); continuity gives at least two distinct roots in (0,1/3).",
    }
    destination = Path(__file__).with_name("critical_fraction_certificate.json")
    destination.write_text(json.dumps(result, indent=2))
    print("PASS: signs +,-,-,+; two disjoint root intervals certified with Fraction arithmetic.")


if __name__ == "__main__":
    main()
