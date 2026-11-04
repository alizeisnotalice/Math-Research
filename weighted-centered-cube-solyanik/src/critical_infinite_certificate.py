#!/usr/bin/env python3
"""Exact scalar certificate for the one-dimensional critical construction.

The proof chooses a centered interval with integer total mass N.  Periodicity
then gives exactly p_1=c+f/N and p_i=c (i>1), c=(m-1)/m.  Fractions certify
the identity prod(p_i)-c**m=c**(m-1)*f/N>0 for sample parameters.
This is only an algebraic reproduction check; the measure argument is in TeX.
"""
from fractions import Fraction


def check(m: int, f: Fraction, total_mass: int) -> Fraction:
    c = Fraction(m - 1, m)
    p = [c + f / total_mass] + [c] * (m - 1)
    value = Fraction(1)
    for x in p:
        value *= x
    excess = value - c**m
    assert excess == c**(m - 1) * f / total_mass
    assert excess > 0
    return excess


def main() -> None:
    print("ONE-DIMENSIONAL CRITICAL EXCESS CERTIFICATE")
    print("m,mu(F),integer_total_mass,exact_excess,status")
    for m, f, total_mass in ((2, Fraction(1, 7), 17),
                             (3, Fraction(2, 9), 29),
                             (4, Fraction(1, 3), 41),
                             (8, Fraction(1, 5), 97)):
        print(f"{m},{f},{total_mass},{check(m, f, total_mass)},PASS")


if __name__ == "__main__":
    main()
