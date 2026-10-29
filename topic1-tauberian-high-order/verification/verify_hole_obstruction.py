#!/usr/bin/env python3
"""Exact rational smoke checks for lower_round002.md; no external packages.

These finite checks support transcription and implementation consistency.
They do not replace the analytic proofs and do not certify a global search.
"""

from decimal import Decimal, localcontext
from fractions import Fraction as F
from itertools import product


def decimal(q, places=22):
    with localcontext() as ctx:
        ctx.prec = places + 8
        return f"{Decimal(q.numerator) / Decimal(q.denominator):.{places}g}"


def one_dimensional_mass(x, radius, delta):
    """Exact mass of Z+[-delta,delta] in [x-radius,x+radius]."""
    left, right = x - radius, x + radius
    k_min = (left - delta).__floor__() - 1
    k_max = (right + delta).__ceil__() + 1
    result = F(0)
    for k in range(k_min, k_max + 1):
        result += max(F(0), min(right, k + delta) - max(left, k - delta))
    return result


def defect_fraction(x, radius, delta):
    result = F(1)
    for coordinate in x:
        result *= one_dimensional_mass(coordinate, radius, delta) / (2 * radius)
    return result


def check_scalar_inequality():
    count = 0
    for n in range(2, 9):
        for v in (F(1, 10), F(1, 50), F(1, 200)):
            c = 1 / (2 * v)
            for k in range(1, 101):
                q = F(k, 100)
                # q^n / ((q+c)^n-1) <= 1 / ((1+c)^n-1).
                assert q**n * ((1 + c) ** n - 1) <= (q + c) ** n - 1
                count += 1
    print(f"scalar monotonicity inequality: {count} exact samples PASS")


def check_all_radius_formula_samples():
    count = 0
    for n in (2, 3, 4):
        for delta in (F(1, 4 * n), F(1, 20 * n)):
            for x in product((F(0), delta / 2, delta), repeat=n):
                r0 = 1 - delta - max(x)
                expected = (delta / r0) ** n
                assert defect_fraction(x, r0, delta) == expected
                radii = (
                    delta / 2,
                    delta,
                    2 * delta,
                    r0 / 2,
                    r0,
                    r0 + delta / 2,
                    F(1),
                    1 + delta,
                    F(3, 2),
                    F(9, 4),
                    F(33, 8),
                    F(10),
                )
                for radius in radii:
                    actual = defect_fraction(x, radius, delta)
                    assert actual >= expected, (n, delta, x, radius, actual, expected)
                    count += 1
    print(f"all-radius formula: {count} exact finite samples PASS")


def check_finite_construction():
    count = 0
    for n in (2, 3, 4):
        for u in (F(1, 20), F(1, 100), F(1, 500)):
            t = u**n
            rho = u * F(999999, 1000000)
            N = 10000000
            delta = rho / (N * (1 + 2 * rho))
            radius = 1 / (N * (1 + 2 * rho))
            assert radius + 2 * delta == F(1, N)
            assert radius + delta < F(1, N)
            assert 2 * delta < radius
            assert rho**n < t
            removed = (1 - F(1, N)) ** n * (2 * rho / (1 + 2 * rho)) ** n
            ratio = 1 / (1 - removed)
            normalized_remainder = (ratio - 1 - 2**n * t) / (t * u)
            print(
                f"finite n={n} u={u}: ratio={decimal(ratio)} "
                f"R_lower={decimal(normalized_remainder)} "
                f"baseline_limit={-n * 2 ** (n + 1)}"
            )
            count += 1
    print(f"finite strict-witness constructions: {count} PASS")


def main():
    check_scalar_inequality()
    check_all_radius_formula_samples()
    check_finite_construction()


if __name__ == "__main__":
    main()
