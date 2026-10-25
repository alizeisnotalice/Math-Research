#!/usr/bin/env python3
"""Exact promotion check for a_2 -> (10001/10000) a_2.

This initial certificate reuses the already audited Q(u,v) arithmetic layer
but rebuilds all witness boxes and the six-dimensional union.  A separate
auditor must still use a non-importing implementation before a new bound is
accepted.
"""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from pathlib import Path
import hashlib,json,sys,time

HERE=Path(__file__).resolve();sys.path.insert(0,str(HERE.parents[1]))
import round001_rectangular_certificate as A

def main():
    started=time.monotonic();u,v,F=A.u,A.v,A.F
    spacing=(F(1)+u**4/2, ((F(2)+u**5*v)/4)*Q(10001,10000), *((F(1)+u)/2 for _ in range(4)))
    records=[]
    for ell in product(range(1,5),repeat=6):
        mass=1
        for z in ell:mass*=z
        r=A.root(mass)/2; box=[]
        for a,k in zip(spacing,ell):
            lo,hi=(F(0),min(a/2,r-(k-1)*a/2)) if k%2 else (max(F(0),k*a/2-r),a/2)
            if hi<=lo:break
            box.extend((lo,hi))
        else:
            starts=[-(k//2) if k%2 else 1-k//2 for k in ell]
            for j,(a,k,s) in enumerate(zip(spacing,ell,starts)):
                for x in box[2*j:2*j+2]:
                    for atom in (s*a,(s+k-1)*a):assert -r<=x-atom<=r
            records.append((ell,tuple(box)))
    unique={box:ell for ell,box in records}
    kept=[(ell,box) for box,ell in unique.items() if not any(other!=box and all(other[2*j]<=box[2*j] and box[2*j+1]<=other[2*j+1] for j in range(6)) for other in unique)]
    R=64*A.box_union(tuple(sorted(box for _,box in kept)));lo,hi=A.bounds(R.a)
    ledger_path=HERE.parents[2]/'rounds/003/frozen_lower.json'
    ledger=json.loads(ledger_path.read_text());old=F(a=tuple([Q(ledger['lower']['exact_interval'][1])]+[Q(0)]*35))
    assert old<R
    scale=10**16
    out={'status':'primary_exact_pass_pending_independent_audit','mechanism':'R6-001 rectangular lattice with second spacing multiplied by 10001/10000',
         'dimension':6,'raw_nonempty_boxes':len(records),'maximal_boxes':len(kept),
         'exact_interval':[str(lo),str(hi)],
         'decimal_enclosure':[f'{lo.numerator*scale//lo.denominator}/10^16',f'{(hi.numerator*scale+hi.denominator-1)//hi.denominator}/10^16'],
         'strictly_exceeds_frozen_R6_upper':True,'frozen_R6_upper':ledger['lower']['exact_interval'][1],
         'witness_counts':[list(x) for x,_ in kept],
         'polynomial_coefficients':{f'{i},{j}':str(R.a[i*6+j]) for i in range(6) for j in range(6) if R.a[i*6+j]},
         'ledger_sha256':hashlib.sha256(ledger_path.read_bytes()).hexdigest(),
         'imported_arithmetic_sha256':hashlib.sha256((HERE.parents[1]/'round001_rectangular_certificate.py').read_bytes()).hexdigest(),
         'code_sha256':hashlib.sha256(HERE.read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-started}
    HERE.with_name('second_spacing_certificate_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
