#!/usr/bin/env python3
"""Finite diagnostic for off-lattice mass insertion. No certification claim."""
from pathlib import Path
from itertools import product
from functools import lru_cache
from math import prod
import numpy as np
import json,time,hashlib,sys
HERE=Path(__file__).resolve().parent
A=np.array([1+2**(3/5)/2]+[(1+2**(1/5))/2]*4)
HALF=A/2
SECONDS=120;NODES=200000

def indices(shift,l):
    # nearest l integers or half-integers for x in[0,a/2].
    if not shift:
        return tuple(range(-(l-1)//2,(l-1)//2+1)) if l%2 else tuple(range(1-l//2,1+l//2))
    return tuple(range(-(l-1)//2,(l-1)//2+1)) if l%2 else tuple(range(-l//2,l//2))

def boxes_for(first,k,eps):
    shift=(first,)+(1,)*k+(0,)*(4-k)
    opts=[]
    for s in shift:
        opts.append([(p,q) for p in range(1,5) for q in range(1,5) if (p==q if not s else abs(p-q)<=1)])
    vectors=np.array(list(product(*opts)),dtype=int)
    p=vectors[:,:,0];q=vectors[:,:,1]
    n0=np.prod(p,axis=1);n1=np.prod(q,axis=1)
    lows=np.empty((len(p),5));highs=np.empty_like(lows)
    for j,s in enumerate(shift):
        e0=np.array([(min(indices(0,l)),max(indices(0,l))) for l in range(1,5)])
        e1=np.array([(min(indices(s,l))+s/2,max(indices(s,l))+s/2) for l in range(1,5)])
        lows[:,j]=np.minimum(e0[p[:,j]-1,0],e1[q[:,j]-1,0])*A[j]
        highs[:,j]=np.maximum(e0[p[:,j]-1,1],e1[q[:,j]-1,1])*A[j]
    mass=n0+eps*n1;r=mass**.2/2
    lo=np.maximum(0,highs-r[:,None]);hi=np.minimum(HALF,lows+r[:,None])
    good=np.all(hi>lo+1e-14,axis=1)
    rows=[tuple(np.stack((l,h),axis=1).flat) for l,h in zip(lo[good],hi[good])]
    metadata=[(tuple(pp),tuple(qq),float(mm),float(rr)) for pp,qq,mm,rr in zip(p[good],q[good],mass[good],r[good])]
    for orb,weight in ((0,1.),(1,eps)):
        for lengths in product(range(1,5),repeat=5):
            mass=weight*prod(lengths);r=mass**.2/2;b=[]
            for j,l in enumerate(lengths):
                s=shift[j] if orb else 0
                inds=indices(s,l);low=(min(inds)+s/2)*A[j];high=(max(inds)+s/2)*A[j]
                left=max(0,high-r);right=min(HALF[j],low+r)
                if right<=left+1e-14:break
                b.extend((left,right))
            else:rows.append(tuple(b));metadata.append((tuple(lengths) if orb==0 else None,tuple(lengths) if orb else None,mass,r))
    unique=dict(zip(rows,metadata));B=np.array(list(unique))
    kept=[]
    for i,b in enumerate(B):
        covers=np.all(B[:,::2]<=b[::2]+1e-14,axis=1)&np.all(B[:,1::2]>=b[1::2]-1e-14,axis=1)
        covers[i]=False
        # Avoid deleting all numerically identical boxes from float coincidence.
        for q in np.flatnonzero(covers):
            if np.max(np.abs(B[q]-b))<1e-14 and q>i:covers[q]=False
        if not covers.any():kept.append(tuple(b))
    return shift,unique,tuple(sorted(kept)),len(vectors)

def run(first,k,eps):
    start=time.monotonic();nodes=0
    shift,raw,kept,configs=boxes_for(first,k,eps)
    guards=0;max_indices=[0]*5
    for b in kept:
        l0,l1,mass,r=raw[b];atoms={}
        for orb,ls,weight in ((0,l0,1.),(1,l1,eps)):
            if ls is None:continue
            for inds in product(*(indices(shift[j] if orb else 0,l) for j,l in enumerate(ls))):
                pt=tuple((2*z+(shift[j] if orb else 0)) for j,z in enumerate(inds))
                assert pt not in atoms
                atoms[pt]=weight
                max_indices=[max(old,abs(z)) for old,z in zip(max_indices,inds)]
        assert abs(sum(atoms.values())-mass)<1e-10
        assert abs((2*r)**5-mass)<1e-9
        for j in range(5):
            coords=[p[j]*A[j]/2 for p in atoms]
            assert max(coords)-r<=b[2*j]+1e-12
            assert b[2*j+1]<=min(coords)+r+1e-12
            guards+=2
    @lru_cache(None)
    def sweep(boxes):
        nonlocal nodes
        nodes+=1
        if nodes>NODES or time.monotonic()-start>SECONDS:raise TimeoutError('registered cap')
        if not boxes:return 0.
        if len(boxes[0])==2:
            out=0.;left=right=None
            for lo,hi in sorted(boxes):
                if left is None:left,right=lo,hi
                elif right<lo:out+=right-left;left,right=lo,hi
                elif hi>right:right=hi
            return out+right-left
        cuts=sorted({x for b in boxes for x in b[:2]});out=0.
        for lo,hi in zip(cuts,cuts[1:]):
            sub=tuple(sorted({b[2:] for b in boxes if b[0]<=lo and hi<=b[1]}))
            out+=(hi-lo)*sweep(sub)
        return out
    vol=sweep(kept);ratio=32*vol/(1+eps)
    return dict(first_shift_half=first,tail_shift_count=k,epsilon=eps,ratio_float=ratio,raw_boxes=len(raw),kept_boxes=len(kept),mixed_configurations=configs,union_nodes=nodes,endpoint_guards=guards,max_period_indices=max_indices,distinct_atoms_checked=True,wall_seconds=time.monotonic()-start)

def main():
    baseline=json.loads((HERE.parent/'lower_round003/a2_equals_a3_results.json').read_text())
    C=baseline['display_only'][0];out=[]
    for first,k in product(range(2),range(5)):
        if first==k==0:continue
        for eps in (1/1024,1/512):
            try:
                row=run(first,k,eps)
                row['difference_float']=row['ratio_float']-C
                row['finite_difference_insertion_response']=(row['ratio_float']-C)/eps
            except TimeoutError as exc:row=dict(first_shift_half=first,tail_shift_count=k,epsilon=eps,status='COMPUTE_LIMIT',reason=str(exc))
            out.append(row);print(json.dumps(row),flush=True)
            (HERE/'screen_results.json').write_text(json.dumps(dict(status='floating_diagnostic_only',rows=out,baseline_source_sha256=hashlib.sha256((HERE.parent/'lower_round003/a2_equals_a3_results.json').read_bytes()).hexdigest(),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2)+'\n')
if __name__=='__main__':main()
