#!/usr/bin/env python3
"""Exact rational verification; standard library only. No floating point required."""
from fractions import Fraction as F
from pathlib import Path
import json

class I:
    def __init__(self, lo, hi=None): self.lo, self.hi = F(lo), F(lo if hi is None else hi)
    def __add__(self, x):
        if not isinstance(x, I): x = I(x)
        return I(self.lo+x.lo, self.hi+x.hi)
    __radd__ = __add__
    def __neg__(self): return I(-self.hi, -self.lo)
    def __sub__(self, x): return self + (-x if isinstance(x,I) else -F(x))
    def __rsub__(self, x): return -self+x
    def __mul__(self, x):
        if not isinstance(x, I): x = I(x)
        z = [self.lo*x.lo,self.lo*x.hi,self.hi*x.lo,self.hi*x.hi]
        return I(min(z), max(z))
    __rmul__ = __mul__
    def __pow__(self, n):
        assert n >= 0 and self.lo >= 0
        return I(self.lo**n, self.hi**n)

def verify(data):
    checked = []
    for row in data['aldaz_intervals']:
        n = row['n']
        t = I(*row['root_2']); s = I(*row['root_3_over_2'])
        assert 0 < t.lo < t.hi and t.lo**n < 2 < t.hi**n
        assert 0 < s.lo < s.hi and s.lo**n < F(3,2) < s.hi**n
        # Transcribed from Aldaz (2000), Proposition 1.4, printed p. 110.
        a = F(3,2)*((1+t)*F(1,2))**(n-1) - F(1,2)*((F(3,2)*(1+t)-2*s)**(n-1) - ((3-t)*F(1,2))**(n-1))
        lo, hi = map(F,row['A_interval'])
        assert lo < a.lo <= a.hi < hi
        assert hi-lo <= F(1,10**49)
        checked.append((n,lo,hi))
    assert [n for n,lo,hi in checked] == list(range(2,8))
    assert all(hi < checked[-1][1] for n,lo,hi in checked[:-1])
    # Lin's own Challenge/README bounds Phi < 1.686; Melas c1 < 1.568.
    assert F(1686,1000) < checked[-1][1]
    assert F(1568,1000) < checked[-1][1]
    assert data['upper_bound'] == 2**7 == 128
    return checked

if __name__ == '__main__':
    data=json.loads(Path(__file__).with_name('a7-certificate.json').read_text())
    for n,lo,hi in verify(data): print(f'A{n}: {lo} < A{n} < {hi}')
    print('PASS: exact rational root brackets, A2..A7, comparisons, and 2^7=128')
