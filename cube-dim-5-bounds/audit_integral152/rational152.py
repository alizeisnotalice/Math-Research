"""Integer/gcd rational path; no Fraction, scalar R, or power roots."""
from math import gcd


class Q:
    def __init__(self, n=0, d=1):
        if isinstance(n, Q):
            n, d = n.n, n.d
        elif isinstance(n, str):
            parts = n.split('/')
            assert len(parts) in (1, 2)
            n, d = int(parts[0]), int(parts[1]) if len(parts) == 2 else 1
        assert type(n) is int and type(d) is int and d > 0
        g = gcd(n, d)
        self.n, self.d = n // g, d // g

    def __add__(self, other):
        other = Q(other)
        g = gcd(self.d, other.d)
        return Q(self.n * (other.d // g) + other.n * (self.d // g),
                 self.d * (other.d // g))

    __radd__ = __add__

    def __neg__(self):
        return Q(-self.n, self.d)

    def __sub__(self, other):
        return self + -Q(other)

    def __rsub__(self, other):
        return Q(other) + -self

    def __mul__(self, other):
        other = Q(other)
        g, h = gcd(self.n, other.d), gcd(other.n, self.d)
        return Q((self.n // g) * (other.n // h),
                 (self.d // h) * (other.d // g))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = Q(other)
        assert other.n != 0
        sign = 1 if other.n > 0 else -1
        return self * Q(sign * other.d, sign * other.n)

    def __rtruediv__(self, other):
        return Q(other) / self

    def __pow__(self, power):
        assert type(power) is int and power >= 0
        return Q(self.n ** power, self.d ** power)

    def __eq__(self, other):
        other = Q(other)
        return self.n == other.n and self.d == other.d

    def __lt__(self, other):
        other = Q(other)
        return self.n * other.d < other.n * self.d

    def __le__(self, other):
        other = Q(other)
        return self.n * other.d <= other.n * self.d

    def __str__(self):
        return str(self.n) if self.d == 1 else f'{self.n}/{self.d}'


def ceil_div(n, d):
    assert type(n) is int and type(d) is int and d > 0
    quotient, remainder = divmod(n, d)
    return quotient + bool(remainder)


def dyadic_bits(q):
    q = Q(q)
    assert q.d & (q.d - 1) == 0
    return q.d.bit_length() - 1


def on_grid(q, bits):
    q = Q(q)
    n, r = divmod(q.n << bits, q.d)
    assert r == 0
    return n
