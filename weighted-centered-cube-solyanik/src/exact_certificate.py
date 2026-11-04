#!/usr/bin/env python3
"""Exact rational certificate for the scalar parameters in Topic 2.

This file deliberately uses only Python's standard library.  It proves exact
identities at the algebraic threshold without treating a floating point scan as
a proof of any maximal-operator statement.
"""
from fractions import Fraction
from itertools import product


MS = (2, 3, 4, 8)
NS = (1, 2, 3)
ALPHAS = (Fraction(1, 2), Fraction(3, 4), Fraction(9, 10))


def lambda_star(m: int) -> Fraction:
    return Fraction(m - 1, m) ** m


def s_at_rational_root(m: int, root: Fraction) -> Fraction:
    """s_m(lambda) when lambda^(1/m) has been supplied exactly as root."""
    return m * (1 - root)


def halo_target(n: int, alpha: Fraction) -> Fraction:
    return 1 + (2**n) * (1 - alpha) / alpha


def star_ratio(n: int, alpha: Fraction, eta: Fraction) -> tuple[Fraction, Fraction]:
    """Return q and the exact finite star ratio 1+2^n q.

    Central mass is one and every signed-corner satellite has mass q.  Taking
    q=(1-eta)(1-alpha)/alpha makes every satellite witness strict, since
    1/(1+q)>alpha.  The construction is a finite model of the lower-bound
    geometry; a separate smoothing lemma is required for its use with a
    density w dx.
    """
    q = (1 - eta) * (1 - alpha) / alpha
    assert Fraction(1, 1) / (1 + q) > alpha
    return q, 1 + (2**n) * q


def signed_corners(n: int) -> list[tuple[int, ...]]:
    return list(product((-1, 1), repeat=n))


def main() -> None:
    print("EXACT SCALAR CERTIFICATE")
    print("m,lambda_star,s(lambda_star),Phi(lambda_star),checks")
    for m in MS:
        root = Fraction(m - 1, m)
        lam = lambda_star(m)
        s = s_at_rational_root(m, root)
        phi = 1 - s
        assert s == 1 and phi == 0
        print(f"{m},{lam},{s},{phi},PASS")

    print("\nSTAR LOWER-BOUND FINITE MODELS")
    print("n,alpha,eta,q,ratio,target,gap,strict_witness")
    eta = Fraction(1, 1000)
    for n in NS:
        assert len(signed_corners(n)) == 2**n
        for alpha in ALPHAS:
            q, ratio = star_ratio(n, alpha, eta)
            target = halo_target(n, alpha)
            assert ratio < target and target - ratio == (2**n) * eta * (1-alpha)/alpha
            print(f"{n},{alpha},{eta},{q},{ratio},{target},{target-ratio},PASS")

    print("\nINTERPRETATION")
    print("All displayed equalities are rational scalar identities or finite atomic")
    print("geometry.  They do not prove a covering theorem or a smoothing transfer.")


if __name__ == "__main__":
    main()
