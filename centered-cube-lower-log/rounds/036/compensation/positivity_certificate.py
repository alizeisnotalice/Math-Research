"""Outward Decimal enclosure for source positivity only (not Gram errors).

Inputs U are the exact dyadics saved by the NumPy diagnostic. Pi is enclosed
by Machin's identity and rational alternating-series bounds. Sin/cos use
degree 49 Taylor polynomials and a proven 1e-35 uniform remainder on [-pi,pi].
"""
from pathlib import Path
from decimal import Decimal as D,localcontext,ROUND_FLOOR,ROUND_CEILING
from fractions import Fraction as F
import json,math,hashlib

HERE=Path(__file__).resolve().parent
PREC=70
def down(f):
    with localcontext()as c:c.prec=PREC;c.rounding=ROUND_FLOOR;return f()
def up(f):
    with localcontext()as c:c.prec=PREC;c.rounding=ROUND_CEILING;return f()
class I:
    def __init__(self,a,b=None):self.a=D(a);self.b=D(a if b is None else b)
    def __add__(self,o):
        o=asiv(o);return I(down(lambda:self.a+o.a),up(lambda:self.b+o.b))
    __radd__=__add__
    def __neg__(self):return I(self.b.copy_negate(),self.a.copy_negate())
    def __sub__(self,o):return self+-asiv(o)
    def __rsub__(self,o):return asiv(o)+-self
    def __mul__(self,o):
        o=asiv(o)
        lows=[down(lambda:x*y)for x in [self.a,self.b]for y in [o.a,o.b]]
        highs=[up(lambda:x*y)for x in [self.a,self.b]for y in [o.a,o.b]]
        return I(min(lows),max(highs))
    __rmul__=__mul__
    def __truediv__(self,o):
        o=asiv(o);assert o.a>0 or o.b<0
        return self*I(down(lambda:1/o.b),up(lambda:1/o.a))
    def square(self):
        low=D(0)if self.a<=0<=self.b else min(down(lambda:self.a*self.a),down(lambda:self.b*self.b))
        return I(low,max(up(lambda:self.a*self.a),up(lambda:self.b*self.b)))
    def sqrt(self):
        assert self.a>=0
        with localcontext()as c:
            c.prec=PREC
            # Decimal.sqrt is correctly rounded to nearest; widen by one ulp.
            a=self.a.sqrt().next_minus();b=self.b.sqrt().next_plus()
        return I(max(D(0),a),b)
def asiv(x):return x if isinstance(x,I)else I(x)
def frac(x):return I(down(lambda:D(x.numerator)/D(x.denominator)),up(lambda:D(x.numerator)/D(x.denominator)))
def atan_bounds(q,n):
    s=sum((F((-1)**k,(2*k+1)*q**(2*k+1))for k in range(n)),F(0))
    nxt=F((-1)**n,(2*n+1)*q**(2*n+1))
    return min(s,s+nxt),max(s,s+nxt)
a,b=atan_bounds(5,50);c,d=atan_bounds(239,20)
pi=I(frac(16*a-4*d).a,frac(16*b-4*c).b)
# Machin angle identity and a generous rigorous remainder cap.
assert (F(120,119)-F(1,239))/(1+F(120,119)*F(1,239))==1
assert F(22,7)**50/math.factorial(50)<F(1,10**35)
assert F(22,7)**51/math.factorial(51)<F(1,10**35)
def trig(turn):
    turn=turn%1
    if turn>F(1,2):turn-=1
    x=2*pi*frac(turn);q=-x.square()
    ct=I(1);st=x;co=ct;si=st
    for k in range(1,25):
        ct=ct*q/((2*k-1)*(2*k));co=co+ct
        st=st*q/((2*k)*(2*k+1));si=si+st
    e=I('-1e-35','1e-35')
    return co+e,si+e
def magnitude(re,im):return (re.square()+im.square()).sqrt()

sources=json.loads((HERE/'source_jitters.json').read_text())
rows=[]
for source in sources:
    K=source['K'];U=[F(x)for x in source['U']];atoms=[(j+x)/K for j,x in enumerate(U)]
    mean=F(1,2)-sum(U,F(0))/K
    R0=sum((x*x-x+F(1,3)for x in U),F(0))/K
    gaps=[(atoms[(j+1)%K]-atoms[j])%1 for j in range(K)]
    cache={}
    for m in range(1,13):
        re=I(0);im=I(0)
        for atom in atoms:
            co,si=trig(m*atom);re=re+co;im=im-si
        cache[m]=(re,im)
    for v in [8,32,128]:
        M=math.ceil(math.sqrt(v));delta=F(1,4*K*v)
        EV=F(K)*delta/6-F(K*K)*delta*delta/12+sum(((delta-z)**3/(3*delta*delta)for z in gaps if z<delta),F(0))
        C0=frac(R0-EV-mean*mean)
        modes=[]
        for m in range(1,M+1):
            _,ss=trig(F(m)*delta/2)
            sinc=ss/(pi*frac(m*delta))
            zr,zi=(w*sinc for w in cache[m])
            phnorm=(zr.square()+zi.square())/(4*pi.square()*m*m)
            C0=C0-2*phnorm
            _,s1=trig(F(5*m,27));_,s2=trig(F(m,18))
            bcoef=4*frac(F(729,29))*((s1/(pi*m)).square().square()+3*(s2/(pi*m)).square().square())
            cm=(2*bcoef/3).sqrt();gamma=pi*m*cm
            modes.append((gamma,zr,zi))
        bs=(6*C0).sqrt();hs=I(1)/bs
        sums={name:I(0)for name in ['original','projected','B_calibrated','H_rescaled']}
        for ga,zr,zi in modes:
            sums['original']+=ga
            sums['projected']+=magnitude(ga-2*zr,-2*zi)
            sums['B_calibrated']+=magnitude(bs*ga-2*zr,-2*zi)
            sums['H_rescaled']+=magnitude(ga-2*hs*zr,-2*hs*zi)
        p=frac(F(1,2))
        lower={name:((1-p*hs if name=='H_rescaled'else 1-p)-p*s/K).a for name,s in sums.items()}
        assert all(z>0 for z in lower.values())
        rows.append({'K':K,'v':v,'positive_density_lower_enclosures':{n:str(z)for n,z in lower.items()},
                     'C0_high_enclosure':[str(C0.a),str(C0.b)]})
out={'scope':'Analytic modified-source positivity only; exact dyadic U; directed Decimal arithmetic. Does not certify NumPy Gram errors.',
     'decimal_precision':PREC,'pi_enclosure':[str(pi.a),str(pi.b)],
     'source_jitters_SHA256':hashlib.sha256((HERE/'source_jitters.json').read_bytes()).hexdigest(),
     'certificate_script_SHA256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
     'sin_cos_uniform_remainder':'1e-35','rows':rows,'all_36_source_bounds_positive':True}
(HERE/'positivity_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'all_36_source_bounds_positive':True,'minimum_positive_lower':str(min(D(z)for r in rows for z in r['positive_density_lower_enclosures'].values()))}))
