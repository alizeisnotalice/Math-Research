"""Exact interval sign certificate for the sixth-order all-input gap.
All assertions use rational arithmetic; numerical geometry is separate.
"""
from fractions import Fraction as F
from pathlib import Path
import json, platform, time
import mpmath as mp

POLYS=[(F(704,32),F(415,32),F(576,32),F(648,32)),
       (F(23,16),F(-330,16),F(-157,16),F(-252,16),F(-324,16)),
       (F(24,16),F(64,16),F(-232,16),F(135,16),F(-66,16),F(-13,16),F(396,16),F(-324,16)),
       (F(304,64),F(640,64),F(240,64),F(801,64),F(1080,64),F(1296,64))]

def evaluate(cs,t): return sum(c*t**j for j,c in enumerate(cs))
def certify_interval(cs,T):
    # Every negative monomial c*t^k >= c*T^k; positive monomials are >=0.
    lower=cs[0]+sum(min(c,F(0))*T**j for j,c in enumerate(cs[1:],1))
    assert lower>0
    return str(lower)

def margins(t):
    a=1-t;r=t/2+t*t/2+F(15,16)*t**3+F(9,4)*t**4
    q=a+t*t/4;ep=t**3/16;M=2*q-a
    g=lambda z:(z-2*r)**2-2*a*z+1
    variation=lambda z:4*ep*(z-2*r)+4*ep**2
    actual=[g(a)-variation(a),-g(q)-variation(q),-g(1/a)-variation(1/a),(a+2*r)**2-4*M*(M-a)-1-16*a*ep]
    expected=[t**5*evaluate(POLYS[0],t),t**4*evaluate(POLYS[1],t),t**3/a**2*evaluate(POLYS[2],t),t**3*evaluate(POLYS[3],t)]
    assert actual==expected and min(actual)>0
    assert 0<2*r+2*ep<a and a<q<(a+1/a)/2
    return {'t':str(t),'alpha':str(a),'r':str(r),'q':str(q),'epsilon':str(ep),'robust_margins':[str(v) for v in actual],'area_gap_lower':str(t**6/768)}

if __name__=='__main__':
    start=time.monotonic();T=F(1,100)
    lower=[certify_interval(cs,T) for cs in POLYS]
    # Rational regression of the independent original formulas and coefficient tables.
    ts=sorted(set([F(j,100000) for j in range(1,1001)]+[F(1,10**j) for j in range(2,13)]))
    rows=[margins(t) for t in ts]
    keep=[row for row in rows if F(row['t']) in {T,F(1,1000),F(1,10**6),F(1,10**12)}]
    out={'python':platform.python_version(),'arithmetic':'fractions.Fraction','interval':['0','1/100'],'positive_interval_lower_bounds':lower,'coefficient_tables':[[str(c) for c in row] for row in POLYS],'regression_count':len(rows),'selected_exact_certificates':keep,'seconds':time.monotonic()-start,'scope':'polynomial lower bounds certify the entire interval; the global geometric transfer is proved in the paper'}
    mp.mp.dps=80
    moment=[]
    for pv in [64,128,256,512,1024,4096,16384]:
        pp=mp.mpf(pv)
        scaled=mp.quad(lambda u: mp.exp((pp-1)*mp.log1p(-u/pp))*u**6,[0,pp/100])/768
        beta_scaled=mp.mpf(15)/16*mp.fprod(pp/(pp+j) for j in range(1,7))
        assert 0<scaled<beta_scaled
        moment.append({'p':pv,'p6_D':mp.nstr(scaled,70),'p6_full_beta':mp.nstr(beta_scaled,70)})
    out['moment_precision']=mp.mp.dps;out['moment_rows']=moment
    dest=Path(__file__).resolve().parents[1]/'results/quantitative_gap_python.json'
    dest.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
