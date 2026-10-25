"""Batch 003: full-cell joint union, staggered periodic layers; numerical exploration."""
from pathlib import Path
from itertools import product
from functools import lru_cache
import json, time, hashlib
import numpy as np
import mpmath as mp
ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
mp.mp.dps=80
u=mp.root(2,4);v=mp.root(3,4)
A=np.array([float(1+u*u/2),float((2+u**3*v)/4),float((1+u)/2),float((1+u)/2)])

def boxes(q,w,cap=4):
    assert 0<=w<=.7 and np.all(q>=0) and np.all(q<=A[1:]/2)
    L=A.copy();L[0]*=2
    rows=[]
    for ell in product(range(1,cap+1),repeat=4):
        for k1 in range(-4,5):
            parity=(-1.)**k1
            ps=parity if ell[0]%2 else 0.
            mass=np.prod(ell[1:])*(ell[0]+w*ps)
            r=mass**.25/2
            lo1=A[0]*(k1+ell[0]-1)-r;hi1=A[0]*k1+r
            if hi1<=0 or lo1>=L[0] or hi1<=lo1:continue
            if ell[0]==1:sm=sM=parity*q/2
            else:sm=-abs(q)/2;sM=abs(q)/2
            choices=[]
            for j in range(1,4):
                opts=[]
                for k in range(-4,5):
                    lo=A[j]*(k+ell[j]-1)+sM[j-1]-r
                    hi=A[j]*k+sm[j-1]+r
                    lo=max(0.,lo);hi=min(L[j],hi)
                    if hi>lo:opts.append((lo,hi))
                choices.append(opts)
            for zz in product(*choices):rows.append((max(0.,lo1),min(L[0],hi1))+sum(zz,()))
    b=np.unique(np.array(rows),axis=0)
    # exact floating comparisons; all certified conclusions kept separate.
    lo=b[:,::2];hi=b[:,1::2];keep=[]
    for i in range(len(b)):
        included=np.all(lo<=lo[i],axis=1)&np.all(hi>=hi[i],axis=1)
        included[i]=False
        if not included.any():keep.append(i)
    return b[keep],len(b)

def volume(box):
    @lru_cache(None)
    def rec(B):
        if not B:return 0.
        if len(B[0])==2:
            a,b=B[0];s=0.
            for c,e in B[1:]:
                if c>b:s+=b-a;a,b=c,e
                else:b=max(b,e)
            return s+b-a
        cuts=sorted({x for b in B for x in b[:2]})
        return sum((v-u)*rec(tuple(sorted({b[2:] for b in B if b[0]<=u and b[1]>=v}))) for u,v in zip(cuts,cuts[1:]))
    return rec(tuple(sorted(set(map(tuple,box))))),rec.cache_info()._asdict()

def evaluate(f,w,pattern):
    q=A[1:]*f*np.array([1,0,0] if pattern=='single' else [1,1,1])
    t=time.monotonic();b,raw=boxes(q,w);vol,cache=volume(b)
    return dict(f=f,w=w,pattern=pattern,q=q.tolist(),volume=vol,ratio=vol/2,raw_boxes=raw,kept_boxes=len(b),seconds=time.monotonic()-t,cache=cache)

def main():
    t=time.monotonic();rows=[evaluate(0,0,'single')];print(rows[-1],flush=True)
    fs=[0,.01,.03,.08,.15,.3,.5];ws=[0,.02,.08,.2,.4,.7]
    for pattern in ['single','diagonal']:
        for f,w in product(fs,ws):
            if f==w==0:continue
            row=evaluate(f,w,pattern);rows.append(row);print(row,flush=True)
            result=dict(status='numerical_only',seed_formula=['1+u^2/2','(2+u^3*v)/4','(1+u)/2','(1+u)/2'],u='positive root u^4=2',v='positive root v^4=3',source_sha256=hashlib.sha256((ROOT/'baseline/original_certificates.json').read_bytes()).hexdigest(),elapsed=time.monotonic()-t,rows=rows)
            (OUT/'staggered_results.json').write_text(json.dumps(result,indent=2)+'\n')
            if row['ratio']>rows[0]['ratio']+1e-7:return
if __name__=='__main__':main()
