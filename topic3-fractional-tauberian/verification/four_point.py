"""Exact rational four-point obstruction; floating exploration is separate."""
from fractions import Fraction as F
from math import isqrt
import json, platform, time
from pathlib import Path

def check(a,r,q,eps=F(0)):
    u=1/a; m=2*q-a
    def g(z): return (z-2*r)**2-2*a*z+1
    def variation(z): return 4*eps*(z-2*r)+4*eps**2
    vals={'alpha':a,'r':r,'q':q,'M':m,'epsilon':eps,
          'pair_margin':g(a)-variation(a),
          'q_negative_margin':-g(q)-variation(q),
          'upper_negative_margin':-g(u)-variation(u),
          'mass_excess':(a+2*r)**2-4*m*(m-a)-1-16*a*eps}
    assert F(0)<2*r<a and a<q<(a+u)/2
    assert 2*r+2*eps<a
    assert all(vals[k]>0 for k in ['pair_margin','q_negative_margin','upper_negative_margin','mass_excess'])
    # Diagonal is the worst of the six point pairs. A rational common pair scale
    # z>a with G(z)>0 proves each pair separately strictly feasible.
    step=(q-a)/2
    while g(a+step)<=0: step/=2
    z=a+step
    assert a<z<u and a*z<min(1,z*z) and 2*a*z<1+(z-2*r)**2
    vals['pair_witness_scale']=z
    return {k:str(v) for k,v in vals.items()}

if __name__=='__main__':
    start=time.monotonic()
    rows=[check(F(9,10),F(9,160),F(901,1000),F(1,100000))]
    # Exhaustive exact verification of 295 fixed thresholds, not a proof over alpha.
    # r is rational immediately inside r*, using integer square-root brackets.
    D=10**60
    alphas=sorted(set([F(j,1000) for j in range(708,1000)]+[1-F(1,10**j) for j in range(2,7)]))
    for a in alphas:
        c=2*a*a-1
        root_upper=F(isqrt(c.numerator*D*D//c.denominator)+1,D)
        r=(a-root_upper)/2
        step=(1/a-a)/8
        for _ in range(300):
            try:
                row=check(a,r,a+step)
                rows.append(row); break
            except AssertionError: step/=2
        else: raise RuntimeError(a)
    out={'environment':{'python':platform.python_version(),'arithmetic':'fractions.Fraction + integer isqrt'},'rows':rows,'count':len(rows),'seconds':time.monotonic()-start,'scope':'exact certificate checks at listed parameters; all-alpha theorem is analytic'}
    dest=Path(__file__).resolve().parents[1]/'results/four_point_python.json'
    dest.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'count':len(rows),'seconds':out['seconds'],'explicit_robust_certificate':rows[0]},indent=2))
