"""Exact pressure of a proposed effective-group bound; not actual geom traffic.
Uniform L1 source on [0,1]^n, legal midpoint balanced tree, query side h.
Three fixed rounds; no statistical fitting, no asymptotic claim from numbers.
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json, hashlib

def run():
    rows=[]
    for stage,(n,k) in enumerate(((4,2),(16,8),(64,32)),1):
        h=F(1,2**k); z=F(1,2); A=h**n
        for name in ('center','boundary_alternating','interior_alternating'):
            off=[F(0) if name=='center' else (F(1,4) if name=='boundary_alternating' else F(3,16))*(-1)**i for i in range(n)]
            pieces=[]
            for d in off:
                lo=z+d*h-h/2;hi=z+d*h+h/2
                left=max(F(0),min(hi,z)-max(lo,z-h))
                right=max(F(0),min(hi,z+h)-max(lo,z))
                assert left+right==h and 0<=lo<hi<=1
                pieces.append((left/h,right/h))
            j=F(1);pmax=F(1);posterior_mass=F(1)
            for l,r in pieces:
                j*=l**3+r**3;pmax*=max(l,r);posterior_mass*=l+r
            assert posterior_mass==1 and j<=pmax**2
            assert pmax<=F(3,4)**n and j<F(1,n)
            explicit=None
            if n==4:
                vals=[]
                for bits in product((0,1),repeat=n):
                    mass=F(1)
                    for i,b in enumerate(bits):mass*=pieces[i][b]*h
                    vals.append(mass/A)
                explicit=sum(p**3 for p in vals)
                assert explicit==j and sum(vals)==1
            rows.append(dict(stage=stage,n=n,k=k,pattern=name,query_side=str(h),capture=str(A),parent_mass=str(A),lambda_='1/2',actual_hard_density='1',J2=str(j),pmax=str(pmax),Keff_squared=str(1/j),Keff_gt_sqrt_n=True,explicit_enumeration_J2=None if explicit is None else str(explicit),scope='legal midpoint tree; no FIRST/CP/GP/early-exit gate certification'))
    return rows
if __name__=='__main__':
    out=Path(__file__).with_name('frontier_pressure_results.json')
    if out.exists():raise SystemExit('Refuse to overwrite an existing first result')
    rows=run()
    result={'status':'PASS_EXACT_PRESSURE','records':len(rows),'seed':None,'arithmetic':'exact rational','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'rows':rows,'conclusion':'Counters universal small effective-group premise for this legal tree, not cube weak bound or full geom remainder.'}
    out.write_text(json.dumps(result,ensure_ascii=False,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k!='rows'},ensure_ascii=False))
