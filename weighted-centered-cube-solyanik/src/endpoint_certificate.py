#!/usr/bin/env python3
"""Exact checks for the first water-filling interval of Topic 2.

The program checks only scalar identities.  The common-centre transmission
theorem used in the report is proved separately in extension.tex.
"""
from fractions import Fraction


def first_interval(weights: tuple[int, ...]):
    amin = min(weights)
    k = weights.count(amin)
    larger = sorted(a for a in weights if a > amin)
    if not larger:
        return amin, k, None, None
    next_a = larger[0]
    return amin, k, next_a, Fraction(k * (next_a - amin), next_a)


def water_fill(weights: tuple[int, ...], total: Fraction):
    active = sorted(range(len(weights)), key=lambda i: weights[i])
    for length in range(1, len(weights) + 1):
        group = active[:length]
        rest = active[length:]
        # At total=length this candidate would force a zero denominator;
        # it cannot be a finite KKT water level and must be skipped.
        if total == length:
            continue
        rho = Fraction(sum(weights[i] for i in group), length) / (1 - total / length)
        if all(weights[i] < rho for i in group) and all(weights[i] >= rho for i in rest):
            ds = tuple(Fraction(0) if i not in group else 1 - Fraction(weights[i], 1) / rho
                       for i in range(len(weights)))
            return rho, ds
    raise RuntimeError("no KKT regime")


def product(weights: tuple[int, ...], ds: tuple[Fraction, ...]) -> Fraction:
    ans = Fraction(1)
    for a, d in zip(weights, ds):
        ans *= (1 - d) ** a
    return ans


def main() -> None:
    cases = ((1, 1, 3, 5), (2, 2, 2, 7), (1, 4, 4, 9),
             (1, 1, 2), (3, 3, 3))
    print("FIRST WATER-FILLING INTERVAL CERTIFICATE")
    print("weights,a_min,k,a_next,s_break,s_test,rho,deficits,P(s),formula,pass")
    for weights in cases:
        amin, k, next_a, s_break = first_interval(weights)
        s_test = Fraction(1, 2) if s_break is None else min(s_break, Fraction(1)) / 2
        rho, ds = water_fill(weights, s_test)
        expected = tuple(Fraction(s_test, k) if a == amin else Fraction(0) for a in weights)
        exact = product(weights, ds)
        formula = (1 - s_test / k) ** (k * amin)
        assert ds == expected
        assert exact == formula
        if s_break is not None and s_break <= 1:
            rho_break, ds_break = water_fill(weights, s_break)
            assert rho_break == next_a
            assert all(ds_break[i] == 0 for i, a in enumerate(weights) if a >= next_a)
        print(f"{weights},{amin},{k},{next_a},{s_break},{s_test},{rho},{ds},{exact},{formula},PASS")

    # Exact near-endpoint instances: lambda=(1-epsilon) and the claimed
    # inverse s=k(1-lambda^(1/(k a_min))) when the root is rational.
    print("\nINVERSE CHECKS")
    print("weights,root,lambda,s,inverse,pass")
    for weights, root in (((1, 1, 3), Fraction(9, 10)), ((2, 2, 7), Fraction(19, 20))):
        amin, k, _, s_break = first_interval(weights)
        lam = root ** (k * amin)
        s = k * (1 - root)
        assert s_break is None or s < min(s_break, Fraction(1))
        rho, ds = water_fill(weights, s)
        assert product(weights, ds) == lam
        assert rho == Fraction(amin, 1) / root
        print(f"{weights},{root},{lam},{s},{s},PASS")


if __name__ == "__main__":
    main()
