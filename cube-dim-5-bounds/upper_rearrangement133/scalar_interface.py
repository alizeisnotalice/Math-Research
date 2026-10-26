#!/usr/bin/env python3
"""Exact rational interface only; no physical budget is run on import.

The caller must supply a qualified rational Cbar>=C3+C1*epsilon^2, a
qualified pi interval, and (for mixed use) a qualified same-epsilon H sup M.
There is no H/R/kernel/log/Gamma/Arb evaluation or H-sup optimization here.
All lower endpoints below concern the abstract comparator, not actual R/I.
"""
from fractions import Fraction as F
from math import isqrt

def fourth_root_interval(value, bits=160):
    """Dyadic enclosing fourth root, verified by exact integer powers."""
    value = F(value)
    assert value >= 0 and type(bits) is int and bits >= 1
    den = 1 << bits
    integer = value.numerator * den**4 // value.denominator
    k = isqrt(isqrt(integer))
    lower = F(k, den)
    assert lower**4 <= value < F(k+1, den)**4
    upper = lower if lower**4 == value else F(k+1, den)
    return lower, upper

def norm_fifth_power(epsilon, Cbar):
    """Upper on ||R_epsilon||_(5/4)^5 when true B=t/pi is used.

    Exact cancellation of the analytic beta integral removes pi. This is
    a new certificate-family budget, not a measured norm or its lower bound.
    """
    epsilon, Cbar = F(epsilon), F(Cbar)
    assert epsilon > 0 and Cbar > 0
    return F(625, 524288) * Cbar * epsilon**6

def mixed_increment_interval(epsilon, Cbar, pi_interval, M, intervals=4096, bits=160):
    """Enclose abstract J_M=integral Phi_M(Rstar), q=5/4.

    Positive nonsingular formula:
      J_M=5*m*D*integral_0^1 z^2[D+(M-D)z^4]^(1/4) dz,
      D=Cbar/(2epsilon^2), m=epsilon^4/(pi*Cbar).
    On each rational z interval, enclose the two nonnegative factors
    separately. This is conservative even when M<D and the product is not
    monotone. Caller may use only the upper endpoint for an actual I bound.
    """
    epsilon, Cbar, M = map(F, (epsilon, Cbar, M))
    pi_lower, pi_upper = map(F, pi_interval)
    assert epsilon > 0 and Cbar > 0 and 0 < pi_lower <= pi_upper and M >= 0
    assert type(intervals) is int and intervals >= 1
    D = Cbar/(2*epsilon**2)
    ml, mu = epsilon**4/(pi_upper*Cbar), epsilon**4/(pi_lower*Cbar)
    lower = upper = F(0)
    for index in range(intervals):
        a, b = F(index, intervals), F(index+1, intervals)
        wa, wb = D+(M-D)*a**4, D+(M-D)*b**4
        wlower, wupper = min(wa, wb), max(wa, wb)
        assert wlower >= 0
        root_lower, _ = fourth_root_interval(wlower, bits)
        _, root_upper = fourth_root_interval(wupper, bits)
        lower += (b-a)*a*a*root_lower
        upper += (b-a)*b*b*root_upper
    result = (5*ml*D*lower, 5*mu*D*upper)
    assert 0 <= result[0] <= result[1]
    return result

def mixed_infinite_tail_upper(epsilon, Cbar, pi_lower, M, T, bits=160):
    """Upper on the positive survival-form tail r in [T,infinity).

    Tail <= (5/3)*m*D^2*(M+T)^(1/4)/T.
    This is an alternate tail interface, not used by the finite z rule.
    """
    epsilon, Cbar, pi_lower, M, T = map(F, (epsilon, Cbar, pi_lower, M, T))
    assert epsilon > 0 and Cbar > 0 and pi_lower > 0 and M >= 0 and T > 0
    D = Cbar/(2*epsilon**2)
    mupper = epsilon**4/(pi_lower*Cbar)
    _, root_upper = fourth_root_interval(M+T, bits)
    return F(5,3)*mupper*D*D*root_upper/T

def complete_I_upper(IH_upper, epsilon, Cbar, pi_interval, M, intervals=4096, bits=160):
    """Conditional interface: needs independent complete IH and global H<=M.

    No qualification of these caller-supplied quantities is inferred here.
    """
    IH_upper = F(IH_upper)
    assert IH_upper >= 0
    _, increment_upper = mixed_increment_interval(epsilon, Cbar, pi_interval, M, intervals, bits)
    return IH_upper+increment_upper
