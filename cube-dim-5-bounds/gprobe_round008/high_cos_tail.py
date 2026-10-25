"""Exact cone moments for the m=3,5 cosine tail diagnostics."""
from __future__ import annotations
from fractions import Fraction as F
from itertools import permutations,product
import hashlib,json,time,sys
from pathlib import Path
from flint import arb,ctx
ctx.prec=256
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT.parent/'numerics'))
from body_fan_geometry import base_terms

def h(n,values):
    # complete homogeneous symmetric polynomial h_n by one-variable updates.
    a=[F(0)]*(n+1);a[0]=F(1)
    for value in values:
        for degree in range(1,n+1): a[degree]+=value*a[degree-1]
    return a[n]

def raw_moments(point,maximum):
    answer=[F(0)]*(maximum+1);count=0
    for perm in permutations(range(5)):
      for signs in product((-1,1),repeat=5):
       for cone,det,rates in base_terms():
        a=[]
        for ray,q in zip(cone,rates):
            transformed=tuple(signs[row]*ray[perm[row]] for row in range(5))
            a.append(sum(F(coordinate)*value for coordinate,value in zip(transformed,point))/q)
        factor=F(det)
        for q in rates:factor/=q
        for k in range(maximum+1):answer[k]+=factor*h(2*k,a)
        count+=1
    assert count==26880
    return answer

def minorant_value(m,moments,T):
    return sum(F((-1)**k,4+2*k)*moments[k]/T**(4+2*k) for k in range(m+1))

def over_pi5(value):
    return arb(value.numerator)/value.denominator/arb.pi()**5

def main():
    began=time.monotonic();points={'corner':(F(1),)*5,'low':(F(1),F(223,4096),F(1021,1024),F(2937,4096),F(3639,4096))};Ts=(F(1,32),F(1,16),F(1,8),F(1,4));rows=[]
    for label,x in points.items():
        moments=raw_moments(x,5)
        entry={'point':[str(v) for v in x],'raw_h_even_0_to_10':[str(v) for v in moments],'values':{}}
        for m in (3,5):
            entry['values'][str(m)]={str(T):{'exact_over_pi5':str(minorant_value(m,moments,T))+' / pi^5','arb_ball':str(over_pi5(minorant_value(m,moments,T)))} for T in Ts}
        rows.append({label:entry})
    # Exact negative/positive signs are decisive because pi^5>0.
    signs={}
    for item in rows:
        label,entry=next(iter(item.items()));signs[label]={str(m):{str(T):('positive' if minorant_value(m,[F(v) for v in entry['raw_h_even_0_to_10']],T)>0 else 'nonpositive') for T in Ts} for m in (3,5)}
    out={'status':'exact_finite_high_cos_tail_diagnostic','cones_per_point':26880,'m_values':[3,5],'T_values':[str(T) for T in Ts],'rows':rows,'signs':signs,'cosine_minorant_scope':'For odd m, cos(u) is bounded below on R by its degree-2m Taylor polynomial; proof recorded in report before any use.','seconds':time.monotonic()-began,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (ROOT/'high_cos_tail.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
