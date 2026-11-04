#!/usr/bin/env python3
"""Exact KKT water-filling certificate for weighted deficit optimization.

For a_i>0 and 0<=d_i<1, maximize
    prod_i (1-d_i)^a_i
subject to sum_i d_i=s.
The logarithm is strictly concave on the active affine section.  This program
enumerates the possible active sets in exact rational arithmetic and then
checks the KKT inequalities.  This is a scalar optimisation result; it is not
a replacement for the geometric common-witness lemma needed by the paper.
"""
from fractions import Fraction
from itertools import product


def water_fill(weights: tuple[Fraction, ...], total: Fraction):
    m = len(weights)
    assert all(a > 0 for a in weights)
    assert 0 <= total < m
    order = sorted(range(m), key=lambda i: weights[i])
    for k in range(1, m + 1):
        active = order[:k]
        inactive = order[k:]
        lam = sum((weights[i] for i in active), Fraction(0)) / (Fraction(k, 1) - total)
        if all(weights[i] < lam for i in active) and all(weights[i] >= lam for i in inactive):
            deficits = [Fraction(0) for _ in weights]
            for i in active:
                deficits[i] = 1 - weights[i] / lam
            assert sum(deficits) == total
            assert all(0 <= d < 1 for d in deficits)
            return lam, tuple(deficits)
    raise RuntimeError("no KKT active set found")


def objective(weights, deficits):
    value = Fraction(1)
    for a, d in zip(weights, deficits):
        # Only integer weights can remain entirely rational.
        assert a.denominator == 1
        value *= (1-d) ** a.numerator
    return value


def grid_sanity(weights, total, denominator=24):
    """Coarse independent rational grid; equality need not be attained on it."""
    m = len(weights)
    numer = total * denominator
    if numer.denominator != 1:
        return None
    target = numer.numerator
    best = None
    for ds in product(range(target + 1), repeat=m-1):
        last = target - sum(ds)
        if last < 0:
            continue
        candidate = tuple(Fraction(x, denominator) for x in (*ds, last))
        if max(candidate) >= 1:
            continue
        val = objective(weights, candidate)
        if best is None or val > best[0]:
            best = (val, candidate)
    return best


def main():
    cases = [
        ((1, 1), Fraction(1, 2)),
        ((1, 2, 3), Fraction(3, 4)),
        ((1, 3, 5, 7), Fraction(5, 4)),
        ((2, 2, 5, 9), Fraction(7, 4)),
    ]
    print("WEIGHTED DEFICIT KKT CERTIFICATE")
    print("weights,total,lambda,deficits,objective,grid_objective,grid_leq_kkt")
    for raw, total in cases:
        weights = tuple(Fraction(x) for x in raw)
        lam, ds = water_fill(weights, total)
        exact = objective(weights, ds)
        grid = grid_sanity(weights, total)
        assert grid is not None and grid[0] <= exact
        print(f"{raw},{total},{lam},{ds},{exact},{grid[0]},{grid[0] <= exact}")


if __name__ == "__main__":
    main()
