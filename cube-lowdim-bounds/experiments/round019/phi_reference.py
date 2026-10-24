"""Exact rational enclosure of the fixed Lin algebraic reference; no search."""
from fractions import Fraction as Q
from math import isqrt
from pathlib import Path
import json
S=10**50
class I:
    def __init__(self,a,b=None): self.lo=Q(a);self.hi=Q(a if b is None else b)
    def __add__(self,o):
        o=o if isinstance(o,I) else I(o);return I(self.lo+o.lo,self.hi+o.hi)
    __radd__=__add__
    def __neg__(self):return I(-self.hi,-self.lo)
    def __sub__(self,o):return self+-I(o) if not isinstance(o,I) else self+-o
    def __rsub__(self,o):return -self+o
    def __mul__(self,o):
        o=o if isinstance(o,I) else I(o);v=[a*b for a in (self.lo,self.hi) for b in (o.lo,o.hi)];return I(min(v),max(v))
    __rmul__=__mul__
    def __truediv__(self,o):
        o=o if isinstance(o,I) else I(o);assert o.lo>0;return self*I(1/o.hi,1/o.lo)
    def sqrt(self):
        assert self.lo>=0
        l=isqrt((self.lo.numerator*S*S)//self.lo.denominator)
        h=isqrt((self.hi.numerator*S*S)//self.hi.denominator)+1
        assert Q(l,S)**2<=self.lo and Q(h,S)**2>=self.hi
        return I(Q(l,S),Q(h,S))
def reference():
    u=(2+I(22).sqrt())/3;w=u*u-1;h=(1+u)/2;v=h+1
    a=2*h-(2*(2+w)).sqrt()/2-u/2
    b=v-I(2).sqrt()*u/2-w.sqrt()/2
    assert a.lo>0 and b.lo>0
    return (2*h*v-4*a*b)/(1+w)
if __name__=='__main__':
    p=reference();out={'source':'Lin fixed 72c022ba09032b522c878cca2f9b6c68217d77ea docs/PROOF.md section 3 phi_eq','formula':'(2*h*V-4*a*b)/(1+w), u=(2+sqrt(22))/3, w=u^2-1, h=(1+u)/2, V=h+1, a=2h-sqrt(2(2+w))/2-u/2, b=V-sqrt(2)*u/2-sqrt(w)/2','lower':str(p.lo),'upper':str(p.hi),'width':str(p.hi-p.lo),'sqrt_scale':str(S),'scope':'Reference enclosure only; zero candidate evaluations.'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print('Phi display',float(p.lo),'width',float(p.hi-p.lo))
