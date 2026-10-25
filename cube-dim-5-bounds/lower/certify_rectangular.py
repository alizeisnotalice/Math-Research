#!/usr/bin/env python3
"""Exact five-dimensional rectangular-witness certificate.

The lattice has one unit atom at every point of a rectangular lattice.  For
each l in {1,2,3,4}^5 the centre set is the contiguous l-block whose common
cube radius is (prod(l))^(1/5)/2.  Reflection reduces the receiver cell to
prod [0,a_j/2], so the certified weak ratio is 32 times the union volume.

All decisions are in Q(u,v), u^5=2, v^5=3.  Rational isolating intervals are
used only to order endpoints and compare final algebraic quantities.
"""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from pathlib import Path
import hashlib
import json


def floor_root(n, k):
    lo, hi = 0, 1 << ((n.bit_length() + k - 1) // k)
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid ** k <= n:
            lo = mid
        else:
            hi = mid
    return lo


def field():
    d, den = 5, 10 ** 60
    roots = [floor_root(n * den ** d, d) for n in (2, 3)]
    ul, vl = [Q(x, den) for x in roots]
    uh, vh = [Q(x + 1, den) for x in roots]
    size = d * d

    @lru_cache(None)
    def add(x, y):
        return tuple(a + b for a, b in zip(x, y))

    @lru_cache(None)
    def neg(x):
        return tuple(-a for a in x)

    @lru_cache(None)
    def mul(x, y):
        out = [Q(0)] * size
        for i, a in enumerate(x):
            if not a:
                continue
            ei, ej = divmod(i, d)
            for j, b in enumerate(y):
                if not b:
                    continue
                fi, fj = divmod(j, d)
                ii, jj = ei + fi, ej + fj
                out[(ii % d) * d + jj % d] += a * b * 2 ** (ii // d) * 3 ** (jj // d)
        return tuple(out)

    @lru_cache(None)
    def bounds(x):
        lo = hi = Q(0)
        for k, c in enumerate(x):
            if not c:
                continue
            i, j = divmod(k, d)
            low, high = ul ** i * vl ** j, uh ** i * vh ** j
            if c > 0:
                lo += c * low
                hi += c * high
            else:
                lo += c * high
                hi += c * low
        return lo, hi

    @lru_cache(None)
    def positive(x):
        if not any(x):
            return False
        lo, hi = bounds(x)
        if lo > 0:
            return True
        if hi < 0:
            return False
        raise ArithmeticError("isolating interval did not separate an endpoint")

    class F:
        __slots__ = ("a",)
        def __init__(self, n=0, a=None):
            self.a = tuple(a) if a is not None else (Q(n),) + (Q(0),) * (size - 1)
        def __hash__(self): return hash(self.a)
        def __eq__(self, other):
            return self.a == (other.a if isinstance(other, F) else F(other).a)
        def __add__(self, other):
            return F(a=add(self.a, (other.a if isinstance(other, F) else F(other).a)))
        __radd__ = __add__
        def __neg__(self): return F(a=neg(self.a))
        def __sub__(self, other): return self + -other
        def __rsub__(self, other): return F(other) + -self
        def __mul__(self, other):
            return F(a=mul(self.a, (other.a if isinstance(other, F) else F(other).a)))
        __rmul__ = __mul__
        def __truediv__(self, other): return F(a=tuple(x / Q(other) for x in self.a))
        def __pow__(self, n):
            out = F(1)
            for _ in range(n): out = out * self
            return out
        def __lt__(self, other):
            other = other if isinstance(other, F) else F(other)
            return positive((other - self).a)
        def __le__(self, other): return self == other or self < other

    ua = [Q(0)] * size; ua[d] = Q(1)
    va = [Q(0)] * size; va[1] = Q(1)
    u, v = F(a=ua), F(a=va)
    assert u ** d == F(2) and v ** d == F(3)
    return F, u, v, bounds, roots, den


def calculate(delta1=Q(0), delta2=Q(0)):
    F, u, v, bounds, roots, den = field()
    a = [1 + u ** 3 / 2 + F(delta1), (2 + u ** 4 * v) / 4 + F(delta2)] + [(1 + u) / 2] * 3
    def root(n):
        e2 = e3 = 0
        while n % 2 == 0: n //= 2; e2 += 1
        while n % 3 == 0: n //= 3; e3 += 1
        assert n == 1
        return u ** e2 * v ** e3
    boxes = []
    for ls in product(range(1, 5), repeat=5):
        r, box = root(__import__('math').prod(ls)) / 2, []
        for aa, ell in zip(a, ls):
            if ell % 2:
                lo, hi = F(0), min(aa / 2, r - (ell - 1) * aa / 2)
            else:
                lo, hi = max(F(0), ell * aa / 2 - r), aa / 2
            if hi <= lo:
                break
            box.extend((lo, hi))
        else:
            boxes.append((ls, tuple(box)))
    unique = {box: ls for ls, box in boxes}
    kept = []
    for box, ls in unique.items():
        if not any(other != box and all(other[2*j] <= box[2*j] and box[2*j+1] <= other[2*j+1] for j in range(5)) for other in unique):
            kept.append((ls, box))
    @lru_cache(None)
    def volume(boxes):
        if not boxes:
            return F(0)
        if len(boxes[0]) == 2:
            left, right, out = None, None, F(0)
            for lo, hi in sorted(boxes):
                if left is None:
                    left, right = lo, hi
                elif right < lo:
                    out += right - left; left, right = lo, hi
                elif right < hi:
                    right = hi
            return out + right - left
        cuts = sorted({x for box in boxes for x in box[:2]})
        out = F(0)
        for lo, hi in zip(cuts, cuts[1:]):
            sub = tuple(sorted({box[2:] for box in boxes if box[0] <= lo and hi <= box[1]}))
            out += (hi - lo) * volume(sub)
        return out
    retained = tuple(sorted(box for _, box in kept))
    recursive_volume = volume(retained)
    # Separate coordinate-cell enumeration: a cell is covered exactly when the
    # intersection of its five one-dimensional covering masks is nonempty.
    segments = []
    for j in range(5):
        cuts = sorted({box[2*j+k] for box in retained for k in (0, 1)})
        one_dim = []
        for lo, hi in zip(cuts, cuts[1:]):
            mask = 0
            for i, box in enumerate(retained):
                if box[2*j] <= lo and hi <= box[2*j+1]:
                    mask |= 1 << i
            if mask:
                one_dim.append((hi - lo, mask))
        segments.append(one_dim)
    def grid_volume(j=0, mask=None):
        if j == 5:
            return F(1)
        total = F(0)
        for length, next_mask in segments[j]:
            overlap = next_mask if mask is None else mask & next_mask
            if overlap:
                total += length * grid_volume(j + 1, overlap)
        return total
    enumerated_volume = grid_volume()
    assert enumerated_volume == recursive_volume
    ratio = 32 * recursive_volume
    lo, hi = bounds(ratio.a)
    return dict(dimension=5, delta1=str(delta1), delta2=str(delta2),
                mechanism="one-orbit rectangular lattice; common-radius contiguous blocks", roots_denominator=str(den),
                root_lower_numerators=roots, retained_witness_counts=[ls for ls, _ in kept],
                union_engines_agree=True, coordinate_segment_counts=[len(x) for x in segments],
                coefficient_basis="u^i v^j, u^5=2, v^5=3",
                coefficients={f"{i // 5},{i % 5}": str(c) for i, c in enumerate(ratio.a) if c},
                ratio_interval=[str(lo), str(hi)], decimal_display_only=[float(lo), float(hi)])


def main():
    candidates = [calculate(Q(i, 64), Q(j, 64)) for i in (-1, 0, 1) for j in (-1, 0, 1)]
    base = next(x for x in candidates if x['delta1'] == '0' and x['delta2'] == '0')
    base_upper = Q(base['ratio_interval'][1])
    for row in candidates:
        row['strictly_above_base'] = Q(row['ratio_interval'][0]) > base_upper
    out = dict(status="exact finite witness certificates", candidates=candidates,
               source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    Path(__file__).with_name('rectangular_certificate.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
