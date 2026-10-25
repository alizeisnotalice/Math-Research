"""Periodic sheared lattice, real axis-aligned shared-radius receiver unions."""
from pathlib import Path
from itertools import product
from fractions import Fraction as F
import sys,time,numpy as np
from scipy.stats import qmc
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent))
from staggered_search import A
from round007_union_fast import volume
IDS=np.array(list(product(range(-3,5),repeat=4)),int)

def coordinates(ids,t):
    xyz=np.asarray(ids,dtype=float)*A;xyz[...,0]+=float(t)*A[0]*np.asarray(ids)[...,1];return xyz

def canonical(ids):
    z=np.asarray(ids,int);z=z-z.min(axis=0);ell3=int(z[:,2].max())+1;ell4=int(z[:,3].max())+1;key=[ell3,ell4];n=0
    for j in np.unique(z[:,1]):
        row=z[z[:,1]==j,0];lo=int(row.min());hi=int(row.max());key.extend([int(j),lo,hi]);n+=hi-lo+1
    assert len(z)==n*ell3*ell4
    return tuple(key)

def old_generators():
    return {tuple([ell3,ell4]+[v for k2 in range(ell2) for v in [k2,0,ell1-1]]) for ell1,ell2,ell3,ell4 in product(range(1,5),repeat=4)}

def capture(t):
    xyz=coordinates(IDS,t);base=qmc.Sobol(4,scramble=True,seed=20261025114).random_base2(8);points=np.r_[base,base*np.array([-1,1,1,1])+np.array([1,0,0,0])]*A;keys=set();passed=0
    for point in points:
        dist=np.max(np.abs(xyz-point),axis=1);order=np.argsort(dist);dd=dist[order];ends=np.r_[np.flatnonzero(dd[1:]!=dd[:-1]),len(dd)-1]
        for end in ends:
            r=float(dd[end]);n=int(end)+1
            if r<=3 and n>(2*r)**4*(1+1e-10):
                keys.add(canonical(IDS[order[:n]]));passed+=1
    return keys,dict(centers=len(points),passing_event_batches=passed,unique_generators=len(keys))

def raw_box(key,t):
    ell3,ell4=key[:2];rows=np.array(key[2:],int).reshape(-1,3);n=int(sum(rows[:,2]-rows[:,1]+1)*ell3*ell4);r=n**.25/2
    low=np.array([np.min(rows[:,1]+float(t)*rows[:,0])*A[0],rows[:,0].min()*A[1],0,0]);high=np.array([np.max(rows[:,2]+float(t)*rows[:,0])*A[0],rows[:,0].max()*A[1],(ell3-1)*A[2],(ell4-1)*A[3]])
    return high-r,low+r,n

def periodic_boxes(keys,t,row=0):
    """Close each generator under ALL lattice shifts, then clip to row D."""
    origin=np.array([0,row*A[1],0,0]);end=origin+A;parts=[];visible=0
    for key in sorted(keys):
        lo,hi,n=raw_box(key,t)
        if np.any(hi<=lo):continue
        ranges=[]
        for j in [1,2,3]:ranges.append(range(int(np.floor((origin[j]-hi[j])/A[j]))-1,int(np.ceil((end[j]-lo[j])/A[j]))+2))
        for k2 in ranges[0]:
            first=range(int(np.floor((origin[0]-hi[0])/A[0]-float(t)*k2))-1,int(np.ceil((end[0]-lo[0])/A[0]-float(t)*k2))+2)
            shifts=np.array(list(product(first,ranges[1],ranges[2])),float)
            delta=np.c_[(shifts[:,0]+float(t)*k2)*A[0],np.full(len(shifts),k2*A[1]),shifts[:,1]*A[2],shifts[:,2]*A[3]]
            l=np.maximum(origin,lo+delta);h=np.minimum(end,hi+delta);good=np.all(h>l+1e-14,axis=1)
            if good.any():parts.append(np.stack([l[good]-origin,h[good]-origin],axis=2).reshape(-1,8));visible+=int(sum(good))
    if not parts:return np.empty((0,8)),dict(visible_boxes=0,unique_boxes=0,kept_boxes=0)
    boxes=np.unique(np.concatenate(parts),axis=0);l=boxes[:,::2];h=boxes[:,1::2];keep=[]
    for i in range(len(boxes)):
        inside=np.all(l<=l[i],axis=1)&np.all(h>=h[i],axis=1);inside[i]=False
        if not inside.any():keep.append(i)
    return boxes[keep],dict(visible_boxes=visible,unique_boxes=len(boxes),kept_boxes=len(keep))

def evaluate(t):
    start=time.monotonic();old=old_generators();captured,meta=capture(t);b0,oldmeta=periodic_boxes(old,t);b,fullmeta=periodic_boxes(old|captured,t);v0=volume(b0)[0];v=volume(b)[0]
    return dict(t=str(F(t)),physical_shear_coefficient=float(t)*A[0]/A[1],old_generators=len(old),new_generators=len(captured-old),capture=meta,old_boxes=oldmeta,full_boxes=fullmeta,old_finite_ratio=v0,full_finite_ratio=v,capture_gain=v-v0,complete_fundamental_mass=1,supercell_mass=F(t).denominator,seconds=time.monotonic()-start),captured
