#!/usr/bin/env python3
"""Deterministic adaptive-quadrature checks for the layer-cake expansion."""
from math import sqrt
from scipy.integrate import quad


def corner_mass(n: int, tau: float, delta: float) -> float:
    """Mass of Q((delta,...,delta),1) for ||x||_infty^(tau-n) dx.

    Only the finite boundary shell is numerically integrated.  The inner and
    outer layer-cake pieces are evaluated in closed form.
    """
    lo, hi = 1 - delta, 1 + delta
    first = 2.0**n * lo**tau / tau
    middle, error = quad(
        lambda s: (s + 1 - delta)**n * s**(tau - n - 1),
        lo, hi, epsabs=1e-13, epsrel=1e-13, limit=200,
    )
    assert error < 1e-11
    tail = 2.0**n * hi**(tau - n) / (n - tau)
    return (n - tau) * (first + middle + tail)


def main() -> None:
    print("LAYER-CAKE CORNER EXPANSION CHECK")
    print("n,tau,delta,normalized quadratic quotient,expected,residual/delta^3")
    cases = ((2, 0.5), (2, 1.25), (3, 0.75), (4, 0.5))
    for n, tau in cases:
        kappa = 2.0**n * n / tau
        expected = -tau * (n - tau) / 2
        for power in (6, 8, 10):
            delta = 2.0**(-power)
            value = corner_mass(n, tau, delta) / kappa
            quotient = (value - 1) / delta**2
            residual = (quotient - expected) / delta
            print(
                f"{n},{tau:.8g},2^-{power},{quotient:.15g},"
                f"{expected:.15g},{residual:.8g}"
            )

    print("\nCUTOFF COMPATIBILITY")
    print("n,b_crit,b(tau=1/2),admissible a interval")
    for n in range(2, 9):
        bcrit = (n + 4 - sqrt(n * n + 16)) / (2 * n)
        tau = 0.5
        b = tau / n
        lower = b / 2
        upper = (1 - 2 * b) / (n - tau)
        assert b < bcrit and lower < upper
        print(
            f"{n},{bcrit:.14g},{b:.10g},({lower:.10g},{upper:.10g})"
        )


if __name__ == "__main__":
    main()
