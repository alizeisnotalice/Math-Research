"""Exact old-block parity embedding and full-supercell Q(u,v) union."""
from pathlib import Path
from fractions import Fraction as F
from functools import lru_cache,cmp_to_key
from itertools import product
import json,time
import numpy as np
from round005_exposure import AV,rvec
from round005_pure_satellite_check import AA,ZERO,interval,compare,add,scale
P=Path(__file__).resolve().parent

def split_block(lo,hi):
    lower=[];upper=[]
    for parity in [0,1]:
        l=[-((-(lo[0]-parity))//2)]+list(lo[1:]);h=[(hi[0]-parity)//2]+list(hi[1:])
        if h[0]<l[0]:l=[0]*4;h=[-1]*4
        lower.extend(l);upper.extend(h)
    return tuple(lower+upper)

def original_blocks():
    for ell in product(range(1,5),repeat=4):
        for lo in product(range(-2,3),range(-2,2),range(-2,2),range(-2,2)):
            hi=tuple(lo[j]+ell[j]-1 for j in range(4));yield lo,hi,ell

def baseline_keys():return set(split_block(lo,hi) for lo,hi,ell in original_blocks())

@lru_cache(None)
def mul(a,b):
    out=[0]*16
    for k,x in enumerate(a):
        if not x:continue
        i,j=divmod(k,4)
        for l,y in enumerate(b):
            if not y:continue
            u,v=divmod(l,4);out[4*((i+u)%4)+(j+v)%4]+=x*y*2**((i+u)//4)*3**((j+v)//4)
    return tuple(out)

@lru_cache(None)
def union_scaled(boxes):
    if not boxes:return ZERO
    if len(boxes[0])==2:
        ss=sorted(boxes,key=cmp_to_key(lambda a,b:compare(a[0],b[0]) or compare(a[1],b[1])));lo,hi=ss[0];value=ZERO
        for a,b in ss[1:]:
            if compare(hi,a)<0:value=add(value,add(hi,scale(lo,-1)));lo,hi=a,b
            elif compare(b,hi)>0:hi=b
        return add(value,add(hi,scale(lo,-1)))
    cuts=sorted({x for b in boxes for x in b[:2]},key=cmp_to_key(compare));value=ZERO
    for lo,hi in zip(cuts,cuts[1:]):
        tail=tuple(sorted({b[2:] for b in boxes if compare(b[0],lo)<=0 and compare(hi,b[1])<=0}))
        value=add(value,mul(add(hi,scale(lo,-1)),union_scaled(tail)))
    return value

def main():
    st=time.monotonic();checked=0;key_to_count={};Q=[scale(AA[0],2)]+AA[1:]
    for lo,hi,ell in original_blocks():
        key=split_block(lo,hi);ll=[key[:4],key[4:8]];hh=[key[8:12],key[12:]];counts=[];extrema=[]
        for parity in [0,1]:
            n=int(np.prod([max(0,hh[parity][j]-ll[parity][j]+1) for j in range(4)]));counts.append(n)
            if n:extrema.extend([2*ll[parity][0]+parity,2*hh[parity][0]+parity])
        n=int(np.prod(ell));assert sum(counts)==n and min(extrema)==lo[0] and max(extrema)==hi[0]
        key_to_count[key]=n;checked+=1
    assert checked==81920
    boxes=set()
    for key,n in key_to_count.items():
        rv=tuple(map(int,rvec(n)));ll=[key[:4],key[4:8]];hh=[key[8:12],key[12:]];active=[o for o in [0,1] if all(hh[o][j]>=ll[o][j] for j in range(4))];b=[]
        for j in range(4):
            low_points=[add(scale(Q[j],ll[o][j]),AA[0] if j==0 and o==1 else ZERO) for o in active]
            high_points=[add(scale(Q[j],hh[o][j]),AA[0] if j==0 and o==1 else ZERO) for o in active]
            mn=min(low_points,key=cmp_to_key(compare));mx=max(high_points,key=cmp_to_key(compare));lo=add(mx,scale(rv,-1));hi=add(mn,rv)
            if compare(lo,ZERO)<0:lo=ZERO
            if compare(hi,Q[j])>0:hi=Q[j]
            if compare(hi,lo)<=0:break
            b.extend([lo,hi])
        else:boxes.add(tuple(b))
    # Safe exact containment pruning; full supercell, no coordinate-reflection shortcut.
    kept=[]
    for b in boxes:
        if not any(c!=b and all(compare(c[2*j],b[2*j])<=0 and compare(b[2*j+1],c[2*j+1])<=0 for j in range(4)) for c in boxes):kept.append(b)
    polynomial=union_scaled(tuple(sorted(kept)));old=json.loads((P.parents[1]/'baseline/original_certificates.json').read_text())[1];expected=[F(0)]*16
    for k,v in old['coefficients'].items():i,j=map(int,k.split(','));expected[4*i+j]=F(v)
    normalized=[F(x,2*4**4) for x in polynomial];assert normalized==expected
    out=dict(status='exact parity/full-supercell embedding passed; not a new bound',old_integer_blocks_checked=checked,distinct_parity_keys=len(key_to_count),raw_nonempty_supercell_boxes=len(boxes),kept_boxes=len(kept),cell_mass=2,supercell_volume_relation='V_super=2*R4 exactly',normalized_coefficients={f'{i//4},{i%4}':str(x) for i,x in enumerate(normalized) if x},old_certificate_coefficients_exact_match=True,root_comparisons='inherited exact rational root brackets, integer Q(u,v) equality; shared arithmetic disclosed',seconds=time.monotonic()-st)
    (P/'round008_embedding_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
