"""Fourcoset old-block embedding; new4cell check reuses original exact R4."""
from pathlib import Path
from itertools import product
from functools import cmp_to_key
from fractions import Fraction as F
import json,time,numpy as np
from round005_exposure import rvec
from round005_pure_satellite_check import AA,ZERO,add,scale,compare
from round008_embedding import union_scaled
P=Path(__file__).resolve().parent
COSETS=[(0,0),(0,1),(1,0),(1,1)];SIGNS=[1,-1,-1,1]

def split_block(lo,hi):
    lows=[];highs=[]
    for p,q in COSETS:
        l=[-((-(lo[0]-p))//2),-((-(lo[1]-q))//2),lo[2],lo[3]];h=[(hi[0]-p)//2,(hi[1]-q)//2,hi[2],hi[3]]
        if any(b<a for a,b in zip(l,h)):l=[0]*4;h=[-1]*4
        lows+=l;highs+=h
    return tuple(lows+highs)

def blocks():
    for ell in product(range(1,5),repeat=4):
        for lo in product(range(-2,3),range(-2,3),range(-2,2),range(-2,2)):
            yield lo,tuple(lo[j]+ell[j]-1 for j in range(4)),ell

def baseline_keys():return {split_block(lo,hi) for lo,hi,ell in blocks()}

def main():
    st=time.monotonic();checked=0;boxes=set();Q=[scale(AA[0],2),scale(AA[1],2),AA[2],AA[3]]
    for lo,hi,ell in blocks():
        key=split_block(lo,hi);ll=[key[4*i:4*i+4] for i in range(4)];hh=[key[16+4*i:20+4*i] for i in range(4)];counts=[int(np.prod([max(0,hh[i][j]-ll[i][j]+1) for j in range(4)])) for i in range(4)];n=int(np.prod(ell));assert sum(counts)==n
        for axis in [0,1]:
            ext=[]
            for i,coset in enumerate(COSETS):
                if counts[i]:ext.extend([2*ll[i][axis]+coset[axis],2*hh[i][axis]+coset[axis]])
            assert min(ext)==lo[axis] and max(ext)==hi[axis]
        checked+=1
        # The exact parity bijection identifies these old physical extrema, radius and mass.
        rv=tuple(map(int,rvec(n)));b=[]
        for j in range(4):
            l=add(scale(AA[j],hi[j]),scale(rv,-1));h=add(scale(AA[j],lo[j]),rv)
            if compare(l,ZERO)<0:l=ZERO
            if compare(h,Q[j])>0:h=Q[j]
            if compare(h,l)<=0:break
            b.extend([l,h])
        else:boxes.add(tuple(b))
    assert checked==102400
    kept=[b for b in boxes if not any(c!=b and all(compare(c[2*j],b[2*j])<=0 and compare(b[2*j+1],c[2*j+1])<=0 for j in range(4)) for c in boxes)]
    coefficient=union_scaled(tuple(sorted(kept)));old=json.loads((P.parents[1]/'baseline/original_certificates.json').read_text())[1];expected=[F(0)]*16
    for label,value in old['coefficients'].items():i,j=map(int,label.split(','));expected[4*i+j]=F(value)
    normalized=[F(c,4*4**4) for c in coefficient];assert normalized==expected
    out=dict(status='new fourcoset embedding exact; old R4 proof reused, no old baseline calculator rerun',original_blocks=checked,coset_order=COSETS,signs=SIGNS,all_integer_counts_extrema_passed=True,four_old_residues_retained=True,nonempty_boxes=len(boxes),kept_boxes=len(kept),cell_mass=4,relation='V_Q=4*R4; V_Q/4=R4 exactly',old16_coefficients_exact_match=True,normalized_coefficients={f'{i//4},{i%4}':str(c) for i,c in enumerate(normalized) if c},shared_arithmetic='old root brackets/field union helper reused; not independent recertification of R4',seconds=time.monotonic()-st)
    (P/'round011_embedding_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
