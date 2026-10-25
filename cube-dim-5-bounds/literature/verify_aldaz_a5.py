#!/usr/bin/env python3
"""Exact rational enclosure for Aldaz (2000), Proposition 1.4 at n=5."""
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import json


def floor_nth_root(value: int, n: int) -> int:
    lo, hi = 0, 1
    while hi**n <= value:
        hi *= 2
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid**n <= value:
            lo = mid
        else:
            hi = mid
    return lo


def add(a, b):
    return (a[0] + b[0], a[1] + b[1])


def neg(a):
    return (-a[1], -a[0])


def sub(a, b):
    return add(a, neg(b))


def mul(a, b):
    values = [x * y for x in a for y in b]
    return (min(values), max(values))


def scale(a, q):
    return mul(a, (q, q))


def power(a, n):
    out = (Fraction(1), Fraction(1))
    for _ in range(n):
        out = mul(out, a)
    return out


def root_interval(numer: int, denom: int, degree: int, digits: int):
    d = 10**digits
    k = floor_nth_root(numer * d**degree // denom, degree)
    assert denom * k**degree <= numer * d**degree < denom * (k + 1) ** degree
    return (Fraction(k, d), Fraction(k + 1, d))


def certificate():
    n = 5
    digits = 70
    u = root_interval(2, 1, n, digits)
    w = root_interval(3, 2, n, digits)
    one = (Fraction(1), Fraction(1))
    a = scale(add(one, u), Fraction(1, 2))
    x = sub(scale(a, 3), scale(w, 2))
    y = scale(sub((Fraction(3), Fraction(3)), u), Fraction(1, 2))
    # Proposition 1.4: (3/2)a^(n-1) - (1/2)(x^(n-1)-y^(n-1)).
    value = sub(scale(power(a, 4), Fraction(3, 2)),
                scale(sub(power(x, 4), power(y, 4)), Fraction(1, 2)))
    return {
        "status": "PASS",
        "source": "Aldaz 2000, Proposition 1.4, p. 110",
        "dimension": 5,
        "formula": "(3/2)*a^4-(1/2)*((3*a-2*(3/2)^(1/5))^4-((3-2^(1/5))/2)^4), a=(1+2^(1/5))/2",
        "root_enclosure_decimal_denominator": str(10**digits),
        "interval": [str(value[0]), str(value[1])],
        "display_decimal": "1.748124634348480040151671905588294427328468848965173734372867563426809...",
    }


if __name__ == "__main__":
    result = certificate()
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    target = Path(__file__).with_name("aldaz_a5_certificate.json")
    target.write_text(encoded, encoding="utf-8")
    print(encoded, end="")
    print("sha256", sha256(target.read_bytes()).hexdigest())
