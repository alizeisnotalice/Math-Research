"""Exact outward interval evaluation of Aldaz (2000), Proposition 1.4, at n=8."""
from fractions import Fraction as F

SCALE = 10**50

def root_interval(numerator: int, denominator: int, degree: int) -> tuple[F, F]:
    """Return adjacent SCALE-grid rationals enclosing (numerator/denominator)^(1/degree)."""
    lo, hi = 0, 2 * SCALE
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if denominator * mid**degree <= numerator * SCALE**degree:
            lo = mid
        else:
            hi = mid
    assert denominator * lo**degree <= numerator * SCALE**degree
    assert denominator * hi**degree > numerator * SCALE**degree
    return F(lo, SCALE), F(hi, SCALE)

def add(x: tuple[F, F], y: tuple[F, F]) -> tuple[F, F]:
    return x[0] + y[0], x[1] + y[1]

def sub(x: tuple[F, F], y: tuple[F, F]) -> tuple[F, F]:
    return x[0] - y[1], x[1] - y[0]

def mul(x: tuple[F, F], y: tuple[F, F]) -> tuple[F, F]:
    products = (x[0]*y[0], x[0]*y[1], x[1]*y[0], x[1]*y[1])
    return min(products), max(products)

def scale(c: F, x: tuple[F, F]) -> tuple[F, F]:
    return mul((c, c), x)

def power(x: tuple[F, F], exponent: int) -> tuple[F, F]:
    assert x[0] >= 0 and exponent >= 0
    return x[0]**exponent, x[1]**exponent

def root_floor(numerator: int, denominator: int, degree: int) -> int:
    return root_interval(numerator, denominator, degree)[0].numerator * SCALE // root_interval(numerator, denominator, degree)[0].denominator

u = root_interval(2, 1, 8)                 # 2^(1/8)
v = root_interval(3, 2, 8)                 # (3/2)^(1/8)
a = scale(F(1, 2), add((F(1), F(1)), u))
b = sub(scale(F(3), a), scale(F(2), v))
c = scale(F(1, 2), sub((F(3), F(3)), u))
# Prop. 1.4: A_n=(3/2)a^(n-1) - (1/2)[b^(n-1)-c^(n-1)].
A8 = sub(scale(F(3, 2), power(a, 7)), scale(F(1, 2), sub(power(b, 7), power(c, 7))))
lower, upper = F(178390524162816, 10**14), F(178390524162817, 10**14)
assert lower < A8[0] <= A8[1] < upper
print(f"root_2_floor={u[0].numerator * SCALE // u[0].denominator}")
print(f"root_3_over_2_floor={v[0].numerator * SCALE // v[0].denominator}")
print(f"A8_outward_lower={A8[0]}")
print(f"A8_outward_upper={A8[1]}")
print("PASS: 178390524162816/10^14 < A8 < 178390524162817/10^14")
