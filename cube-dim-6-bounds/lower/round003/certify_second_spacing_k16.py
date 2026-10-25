#!/usr/bin/env python3
"""Exact promotion check for the pre-registered line winner k=16."""
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import hashlib,json,sys,time
HERE=Path(__file__).resolve();sys.path.insert(0,str(HERE.parents[1]))
import round001_rectangular_certificate as A

def main():
    t=time.monotonic();u,v,F=A.u,A.v,A.F
    a=(F(1)+u**4/2,((F(2)+u**5*v)/4)*Q(10016,10000),*((F(1)+u)/2 for _ in range(4)))
    raw=[]
    for ell in product(range(1,5),repeat=6):
        m=1
        for z in ell:m*=z
        r=A.root(m)/2;b=[]
        for x,k in zip(a,ell):
            lo,hi=(F(0),min(x/2,r-(k-1)*x/2)) if k%2 else (max(F(0),k*x/2-r),x/2)
            if hi<=lo:break
            b.extend((lo,hi))
        else:
            s=[-(k//2) if k%2 else 1-k//2 for k in ell]
            for j,(x,k,q) in enumerate(zip(a,ell,s)):
                for c in b[2*j:2*j+2]:
                    for z in (q*x,(q+k-1)*x):assert -r<=c-z<=r
            raw.append((ell,tuple(b)))
    U={b:e for e,b in raw};keep=[(e,b) for b,e in U.items() if not any(c!=b and all(c[2*j]<=b[2*j] and b[2*j+1]<=c[2*j+1] for j in range(6)) for c in U)]
    R=64*A.box_union(tuple(sorted(b for _,b in keep)));lo,hi=A.bounds(R.a)
    ledger=HERE.with_name('second_spacing_certificate_results.json');L={'lower':json.loads(ledger.read_text())}
    old=F(a=(Q(L['lower']['exact_interval'][1]),)+(Q(0),)*35);assert old<R
    scale=10**14
    out={'status':'primary_exact_pass_pending_independent_audit','candidate':'a2*(1+16/10000)','raw_nonempty_boxes':len(raw),'maximal_boxes':len(keep),'exact_interval':[str(lo),str(hi)],'decimal_enclosure':[f'{lo.numerator*scale//lo.denominator}/10^14',f'{(hi.numerator*scale+hi.denominator-1)//hi.denominator}/10^14'],'strictly_exceeds_frozen_current_lower_upper':True,'frozen_upper':L['lower']['exact_interval'][1],'witness_counts':[list(e) for e,_ in keep],'coefficients':{f'{i},{j}':str(R.a[i*6+j]) for i in range(6) for j in range(6) if R.a[i*6+j]},'ledger_sha256':hashlib.sha256(ledger.read_bytes()).hexdigest(),'code_sha256':hashlib.sha256(HERE.read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-t}
    HERE.with_name('second_spacing_k16_certificate_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
