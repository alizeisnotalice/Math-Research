#!/usr/bin/env python3
"""Numerical checks for the dimension-only second-order cutoff window.

This verifies algebra used by round17/norm_geometry.tex; it is not a proof of
the weak-type hypothesis or of the norm-ball geometric theorem.
"""
from decimal import Decimal, getcontext

getcontext().prec = 70


def sqrt_decimal(x: Decimal) -> Decimal:
    return x.sqrt()


def threshold(n: int) -> Decimal:
    nn = Decimal(n)
    return min(Decimal(1), (nn + 4 - sqrt_decimal(nn * nn + 16)) / 2)


def main() -> None:
    print("NORM-GEOMETRY CUTOFF CERTIFICATE")
    print("n,T_n,tau,tested_interval_nonempty,D(n,tau)")
    for n in range(1, 13):
        tn = threshold(n)
        tau = tn * Decimal("0.8")
        low = tau / (2 * n)
        high = (Decimal(n) - 2 * tau) / (Decimal(n) * (Decimal(n) - tau))
        dconst = tau * (Decimal(n) - tau) * (Decimal(2) ** (Decimal(n) + 2 - tau))
        assert Decimal(0) < tau < tn
        assert low < high
        print(f"{n},{tn},{tau},({low},{high}),{dconst}")

    # Exact algebraic form of the cutoff inequality, evaluated numerically.
    for n in (1, 2, 3, 8, 32):
        tau = threshold(n) * Decimal("0.9")
        left = tau * (Decimal(n) - tau)
        right = 2 * Decimal(n) - 4 * tau
        assert left < right
    print("PASS: all cutoff windows and quadratic-error inequalities verified")


if __name__ == "__main__":
    main()
