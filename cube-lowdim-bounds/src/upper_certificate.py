"""Exact Fraction certificate for the scalar Cauchy-kernel conditions.
The operator transfer is a separate analytic argument; no Lean build is asserted.
"""
from fractions import Fraction as F
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[1]

class I:
    def __init__(self,a,b=None):self.a=F(a);self.b=F(a if b is None else b)
    def __add__(self,o):
        if not isinstance(o,I):o=I(o)
        return I(self.a+o.a,self.b+o.b)
    __radd__=__add__
    def __neg__(self):return I(-self.b,-self.a)
    def __sub__(self,o):return self+-asI(o)
    def __rsub__(self,o):return asI(o)+-self
    def __mul__(self,o):
        o=asI(o);v=[x*y for x in (self.a,self.b) for y in (o.a,o.b)]
        return I(min(v),max(v))
    __rmul__=__mul__
    def __truediv__(self,o):
        o=asI(o);assert not o.a<=0<=o.b
        return self*I(1/o.b,1/o.a)
    def __rtruediv__(self,o):return asI(o)/self
    def pair(self):return [str(self.a),str(self.b)]
def asI(o):return o if isinstance(o,I) else I(o)

def logq(q,n=100):
    q=F(q);assert q>0
    z=(q-1)/(q+1);s=sum((z**(2*k+1)/F(2*k+1) for k in range(n)),F(0))*2
    e=2*abs(z)**(2*n+1)/(F(2*n+1)*(1-z*z))
    # Outward rational compression keeps the exported certificate readable.
    D=10**35;lo=s-e;hi=s+e
    return I(F(lo.numerator*D//lo.denominator,D),F(-((-hi.numerator*D)//hi.denominator),D))

def certify(R,e,r0):
    R,e,r0=map(F,(R,e,r0));s=2*r0/R;g=e*R*(R-1)
    L=logq(1-s/2);lr=logq(r0);lc=logq(1-r0);l2=logq(2)
    J=I(F(1,2)-2/s)-4/s**2*L
    T=2*(2-r0)*(1-l2+lr)-2*(1-r0)*lc
    Jp=I(2/s**2+2/(s**2*(1-s/2)))+8/s**3*L
    Tp=2*l2-2*lr+2*lc+4/r0-2
    Q=J-g*T;Qp=Jp*(2/R)-g*Tp
    # J convex and T''<=-27/2 give Q''>=18 when gamma>4/3.
    assert g>F(4,3)
    b=max(abs(Qp.a),abs(Qp.b));margin=Q.a-b*b/36
    ext=1-g*(2-2*l2)
    assert 1<R and 0<e<1 and F(1,2)<r0<1
    assert margin>0 and ext.a>0
    # K>=1: R - e(R-1) r(1-r) >= R-e(R-1)/4 > 0.
    dominance=R-e*(R-1)/4;assert dominance>0
    cost=R*R/(R-1)-e/6;old=F(2787954611,718800000)
    return {'R':str(R),'epsilon':str(e),'r0':str(r0),'gamma':str(g),'cost':str(cost),
      'cost_decimal_display':float(cost),'improvement_vs_public_exact':str(old-cost),
      'Q_at_r0':Q.pair(),'Qprime_at_r0':Qp.pair(),'global_Q_lower':str(margin),
      'global_Q_lower_display':float(margin),'exterior_slack':ext.pair(),
      'dominance_margin':str(dominance),'log_terms':100,
      'status':'scalar_conditions_certified; see separate proof for operator transfer'}
