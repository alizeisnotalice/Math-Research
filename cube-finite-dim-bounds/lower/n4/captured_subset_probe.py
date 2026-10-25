"""Finite Sobol probe for nonrectangular cube-captured atom subsets."""
from pathlib import Path
from itertools import product
import json,time,hashlib
import numpy as np
from scipy.stats import qmc
from staggered_search import A, boxes, volume
P=Path(__file__).resolve().parent

def run(f,w,pattern,power=12):
    st=time.monotonic();q=A[1:]*f*np.array([1,0,0] if pattern=='single' else [1,1,1]);L=A.copy();L[0]*=2
    b,_=boxes(q,w);bv=volume(b)[0]/2
    centers=qmc.Sobol(4,scramble=True,seed=2026102504).random_base2(power)*L
    covered=np.zeros(len(centers),dtype=bool)
    for box in b:covered|=np.all(centers>=box[::2],axis=1)&np.all(centers<=box[1::2],axis=1)
    ids=np.array(list(product(range(-3,5),repeat=4)))
    parity=1-2*(ids[:,0]%2);points=ids*A;points[:,1:]+=parity[:,None]*q/2;weights=1+w*parity
    added={};discovered=0
    for c in centers[~covered]:
        dist=np.max(abs(points-c),axis=1);idx=np.argsort(dist)[:256];mass=np.cumsum(weights[idx]);score=mass/(2*dist[idx])**4
        best=np.argmax(score)
        if score[best]<=1+1e-10:continue
        subset=idx[:best+1];key=tuple(sorted(map(int,subset)))
        if key in added:continue
        r=mass[best]**.25/2
        lo=np.maximum(0,np.max(points[subset],axis=0)-r);hi=np.minimum(L,np.min(points[subset],axis=0)+r)
        if np.all(hi>lo):added[key]=np.stack([lo,hi],axis=1).ravel();discovered+=1
    extras=np.array(list(added.values()))
    if len(extras):joint=np.concatenate([b,extras]);vv,cache=volume(joint)
    else:vv=2*bv;cache={}
    nonrect=sum(np.prod([int(ids[list(k),j].max()-ids[list(k),j].min()+1) for j in range(4)])>len(k) for k in added)
    return dict(verified_nonrectangular_subsets=int(nonrect),f=f,w=w,pattern=pattern,sobol_power=power,seed=2026102504,point_index_range=[-3,4],prefix_cap=256,uncovered_centers=int(sum(~covered)),new_subsets=discovered,baseline_finite_ratio=bv,enriched_ratio=vv/2,seconds=time.monotonic()-st,cache=cache,subsets=[list(k) for k in added])

def main():
    rows=[];st=time.monotonic()
    for pattern,f,w in product(['single','diagonal'],[.03,.15],[0,.2]):
        row=run(f,w,pattern);rows.append(row);print({k:v for k,v in row.items() if k!='subsets'},flush=True)
        (P/'captured_subset_results.json').write_text(json.dumps(dict(status='finite numerical probe, not family exclusion',elapsed=time.monotonic()-st,rows=rows),indent=2)+'\n')
if __name__=='__main__':main()
