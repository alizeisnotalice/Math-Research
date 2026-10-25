#!/usr/bin/env python3
"""Exact rational enclosure for Aldaz (2000), Proposition 1.4, at n=6.

No floating-point arithmetic is used.  ``floor_nth_root`` supplies lower
bounds for 2^(1/6) and (3/2)^(1/6), and the next rational grid point supplies
upper bounds.  The interval propagation is deliberately dependency-blind,
so its result is wider than necessary but remains rigorous.
"""

from fractions import Fraction


def floor_nth_root(num: int, den: int, n: int, scale: int) -> int:
    """floor(scale * (num / den) ** (1/n)) by integer bisection."""
    lo, hi = 0, 2 * scale
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if den * mid**n <= num * scale**n:
            lo = mid
        else:
            hi = mid
    return lo


def interval_pow(lo: Fraction, hi: Fraction, n: int) -> tuple[Fraction, Fraction]:
    assert 0 <= lo <= hi
    return lo**n, hi**n


def aldaz_interval(n: int, scale: int) -> tuple[Fraction, Fraction]:
    """Dependency-blind rigorous enclosure for the A_n formula."""
    uraw = floor_nth_root(2, 1, n, scale)
    vraw = floor_nth_root(3, 2, n, scale)
    u0, u1 = Fraction(uraw, scale), Fraction(uraw + 1, scale)
    v0, v1 = Fraction(vraw, scale), Fraction(vraw + 1, scale)
    assert u0**n <= 2 < u1**n
    assert v0**n <= Fraction(3, 2) < v1**n

    # a=(1+u)/2; h=3a-2v; b=(3-u)/2.
    a0, a1 = (1 + u0) / 2, (1 + u1) / 2
    h0, h1 = 3 * a0 - 2 * v1, 3 * a1 - 2 * v0
    b0, b1 = (3 - u1) / 2, (3 - u0) / 2
    assert 0 < h0 < h1 and 0 < b0 < b1
    power = n - 1
    ap = interval_pow(a0, a1, power)
    hp = interval_pow(h0, h1, power)
    bp = interval_pow(b0, b1, power)

    # A6=(3/2)a^5-(1/2)(h^5-b^5).  Interval arithmetic with correlations
    # omitted is still valid: h^5-b^5 lies in [h5.lo-b5.hi,h5.hi-b5.lo].
    lo = Fraction(3, 2) * ap[0] - Fraction(1, 2) * (hp[1] - bp[0])
    hi = Fraction(3, 2) * ap[1] - Fraction(1, 2) * (hp[0] - bp[1])
    return lo, hi


def main() -> None:
    # This precision is intentionally far beyond the displayed endpoints.
    scale = 10**80
    lo, hi = aldaz_interval(6, scale)
    target_lo = Fraction(1763857472539509489679047135, 10**27)
    target_hi = Fraction(1763857472539509489679047136, 10**27)
    assert target_lo < lo <= hi < target_hi
    lo5, hi5 = aldaz_interval(5, scale)
    # A5<1.749 while A6>1.763; this is a strict, exact comparison, so the
    # Aldaz n=5 formula cannot improve the direct n=6 formula by monotonicity.
    assert hi5 < Fraction(1749, 1000) < Fraction(1763, 1000) < lo
    print("PASS exact enclosure")
    print("A6 in (1763857472539509489679047135/10^27,")
    print("       1763857472539509489679047136/10^27)")
    print("The two endpoints are strict rational bounds.")
    print("PASS exact comparison: Aldaz A5 < A6.")


if __name__ == "__main__":
    main()
