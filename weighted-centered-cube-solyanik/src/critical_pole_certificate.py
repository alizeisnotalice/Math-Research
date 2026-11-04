#!/usr/bin/env python3
"""Exact scalar certificates for the supercritical critical-pole expansion.

All arithmetic is rational.  The checks include a switch at s=1, so no
claim below relies on the critical water level avoiding an exponent.
"""
from fractions import Fraction
from math import comb


def critical_data(weights: tuple[int, ...]):
    """Return the left active set and exact data at total defect s=1."""
    order = sorted(range(len(weights)), key=lambda i: weights[i])
    for q in range(2, len(weights) + 1):
        active = order[:q]
        inactive = order[q:]
        total = sum(weights[i] for i in active)
        rho = Fraction(total, q - 1)
        if all(weights[i] < rho for i in active) and all(weights[i] >= rho for i in inactive):
            defects = tuple(1 - Fraction(a, rho) if a < rho else Fraction(0)
                            for a in weights)
            lam = Fraction(1)
            for a, d in zip(weights, defects):
                lam *= (1 - d) ** a
            return active, total, rho, defects, lam
    raise RuntimeError("no finite left water-filling regime at s=1")


def water_fill(weights: tuple[int, ...], s: Fraction):
    """Exact KKT optimizer at an admissible rational total defect s."""
    order = sorted(range(len(weights)), key=lambda i: weights[i])
    for q in range(1, len(weights) + 1):
        if s == q:
            continue
        active = order[:q]
        inactive = order[q:]
        rho = Fraction(sum(weights[i] for i in active), q) / (1 - s / q)
        if all(weights[i] < rho for i in active) and all(weights[i] >= rho for i in inactive):
            defects = tuple(1 - Fraction(a, rho) if a < rho else Fraction(0)
                            for a in weights)
            return rho, defects
    raise RuntimeError("no KKT regime")


def product(weights: tuple[int, ...], defects: tuple[Fraction, ...]) -> Fraction:
    ans = Fraction(1)
    for a, d in zip(weights, defects):
        ans *= (1 - d) ** a
    return ans


def main() -> None:
    cases = ((1, 1), (1, 1, 2), (1, 2, 3), (1, 1, 3, 5),
             (2, 2, 2, 7), (1, 1, 1))
    print("CRITICAL LEFT-WATER-LEVEL CERTIFICATE")
    print("weights,q,A,rho_c,lambda_c,defects,linear beta coefficient,pass")
    for weights in cases:
        active, total, rho, defects, lam = critical_data(weights)
        q = len(active)
        assert sum(defects) == 1
        assert product(weights, defects) == lam
        # Choose a rational beta inside the final active stratum.
        strict_gap = min(Fraction(total, a) - (q - 1) for a in (weights[i] for i in active))
        beta = min(Fraction(1, 10), strict_gap / 2)
        rho_beta, defects_beta = water_fill(weights, 1 - beta)
        expected_rho = Fraction(total, q - 1 + beta)
        expected = tuple(1 - Fraction(a, expected_rho) if a < rho else Fraction(0)
                         for a in weights)
        exact = product(weights, defects_beta)
        formula = lam * (1 + beta / (q - 1)) ** total
        linear = lam * Fraction(total, q - 1)
        assert rho_beta == expected_rho
        assert defects_beta == expected
        assert exact == formula
        assert linear == rho * lam
        # The binomial coefficient is the exact coefficient of beta.
        assert lam * Fraction(comb(total, 1), q - 1) == rho * lam
        print(f"{weights},{q},{total},{rho},{lam},{defects},{linear},PASS")

    # A visibly switching case: a_3=rho_c, yet it remains inactive from the left.
    active, total, rho, defects, lam = critical_data((1, 1, 2))
    assert rho == 2 and defects == (Fraction(1, 2), Fraction(1, 2), Fraction(0))
    print("switch-at-s=1 case (1,1,2): a_3=rho_c=2 and d_3=0, PASS")


if __name__ == "__main__":
    main()
