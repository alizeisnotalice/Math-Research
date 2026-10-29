#!/usr/bin/env python3
"""Exact rational checks of weighted and common-radius model identities.

This script checks algebra and test values. It is not the proof of the
geometric or convolution inequalities; those proofs are in the report.
"""
from fractions import Fraction as F
rows = []
delta = F(1, 4)
eps = F(1, 100)
R = 1 + delta + eps
assert 0 < eps < 1 - 3 * delta
assert 2 * delta < R < 2 - 2 * delta

for n in range(2, 11):
    m = 2**n
    b_delta = ((1 + delta) / delta)**n
    b_eps = (R / delta)**n
    b_star = 4**n
    for c in (1, 10, 1000, 10**6):
        B_delta = m * (1 + c) / (F(c) + b_delta)
        B_eps = m * (1 + c) / (F(c) + b_eps)
        B_star = m * F(1 + c, c + b_star)
        K_delta = m * (b_delta - 1)
        K_star = m * (b_star - 1)
        assert F(m) - B_delta == K_delta / (F(c) + b_delta)
        assert F(m) - B_delta <= K_delta / c
        assert F(m) - B_star == F(K_star, c + b_star)
        assert F(m) - B_star <= F(K_star, c)
        assert B_eps < B_delta < B_star < m
        for t in (F(1, 2 * b_star), F(1, b_star)):
            transferred = B_star * (1 + c*t) / ((1+c) * (1-t))
            baseline = m*t / (1-t)
            gap = m * (1-b_star*t) / ((c+b_star)*(1-t))
            assert transferred - baseline == gap >= 0
        rows.append({"n": n, "c": c, "b_delta": str(b_delta),
                     "W_lower_delta": float(B_delta),
                     "m_minus_W_lower_delta": float(m-B_delta),
                     "K_delta_over_c": float(K_delta/c),
                     "W_lower_star": float(B_star)})
    a = F(n, 2*(n+1))
    kappa = a * 2**(n+1)
    assert kappa == F(n*2**n, n+1)
    for u in (F(1, 10), F(1, 100), F(1, 1000)):
        t = u**n
        q_root = 2*u/(1+2*u)
        q = q_root**n
        # Existing cubic-hole family is consistent with the new restricted bound.
        assert q + a*q*q_root <= m*t

print("PASS: exact rational geometry, 36 weighted examples, 72 transfer checks,")
print("      27 fixed-radius consistency checks, and coefficient identities for n=2,...,10.")
print("delta=1/4, epsilon=1/100, R=", R)
print("n=2, c=1000000: W >=", rows[3]["W_lower_delta"])
print("n=3, c=1000000: W >=", rows[7]["W_lower_delta"])
