"""Three-orbit actual witness engine: full masses and common-radius4D union."""
from pathlib import Path
from itertools import product
import time
import numpy as np
from scipy.stats import qmc
from round004_satellite_search import A,PRIMARY as PRIMARY2
from round007_union_fast import volume
PRIMARY={tuple(list(k[:8])+[0]*4+list(k[8:])+[-1]*4) for k in PRIMARY2}
PURE=set();PAIRS=set()
for k in product([-1,0,1],repeat=4):
    for orb in [1,2]:
        ll=np.zeros((3,4),int);hh=np.full((3,4),-1,int);ll[orb]=hh[orb]=k;PURE.add(tuple(map(int,np.r_[ll.ravel(),hh.ravel()])))
    for rel in product([-1,0,1],repeat=4):
        ll=np.zeros((3,4),int);hh=np.full((3,4),-1,int);ll[1]=hh[1]=k;ll[2]=hh[2]=np.array(k)+rel;PAIRS.add(tuple(map(int,np.r_[ll.ravel(),hh.ravel()])))

def capture(a,phases,weights,t):
    events=[]
    for orb in range(3):
        for j in range(4):events.extend(np.abs((t[:,j,None]-phases[orb,j]-np.arange(-2,4))*a[j]).T)
    rr=np.array(events).T
    ll=np.ceil(t[:,None,None,:]-phases[None,None,:,:]-rr[:,:,None,None]/a-1e-11).astype(np.int16)
    hh=np.floor(t[:,None,None,:]-phases[None,None,:,:]+rr[:,:,None,None]/a+1e-11).astype(np.int16)
    counts=np.prod(np.maximum(0,hh-ll+1),axis=3);mass=counts@weights;score=mass/np.maximum((2*rr)**4,1e-30);good=(score>1+1e-10)&(rr<=2.5)&(mass>0);score[~good]=-np.inf
    largest=np.where(good,rr,-np.inf);chosen=np.stack([np.argmax(score,axis=1),np.argmax(largest,axis=1)],axis=1);keys=set()
    for i in range(len(t)):
        for j in chosen[i]:
            if not good[i,j]:continue
            lo=ll[i,j].copy();hi=hh[i,j].copy()
            for orb in range(3):
                if np.any(hi[orb]<lo[orb]):lo[orb]=0;hi[orb]=-1
            keys.add(tuple(map(int,np.r_[lo.ravel(),hi.ravel()])))
    return keys

def boxes(a,phases,weights,keys):
    ordered=sorted(keys);data=np.array(ordered,int);ll=data[:,:12].reshape(-1,3,4);hh=data[:,12:].reshape(-1,3,4)
    counts=np.prod(np.maximum(0,hh-ll+1),axis=2);mass=counts@weights;r=mass**.25/2
    low=(ll+phases)*a;high=(hh+phases)*a;low[counts==0]=np.inf;high[counts==0]=-np.inf
    lo=np.maximum(0,np.max(high,axis=1)-r[:,None]);hi=np.minimum(a,np.min(low,axis=1)+r[:,None]);ok=np.all(hi>lo+1e-14,axis=1)
    active=np.flatnonzero(ok);b=np.stack([lo[ok],hi[ok]],axis=2).reshape(-1,8);b,ii=np.unique(b,axis=0,return_index=True);active=active[ii];lows=b[:,::2];highs=b[:,1::2];keep=[]
    for i in range(len(b)):
        inside=np.all(lows<=lows[i],axis=1)&np.all(highs>=highs[i],axis=1);inside[i]=False
        if not inside.any():keep.append(i)
    return b[keep],[ordered[int(active[i])] for i in keep]

def evaluate(a,phi1,phi2,sigma,share,pool=(),power=8,seed=20261025707):
    st=time.monotonic();phases=np.array([np.zeros(4),phi1,phi2]);weights=np.array([1,sigma*share,sigma*(1-share)])
    t=qmc.Sobol(4,scramble=True,seed=seed).random_base2(power);keys=PRIMARY|PURE|PAIRS|set(pool)|capture(np.array(a),phases,weights,t)
    b,selected=boxes(np.array(a),phases,weights,keys);v=float(volume(b)[0]/sum(weights))
    row=dict(a=list(map(float,a)),phi1=list(map(float,phi1)),phi2=list(map(float,phi2)),sigma=float(sigma),share1=float(share),epsilon1=float(weights[1]),epsilon2=float(weights[2]),cell_mass=float(sum(weights)),ratio=v,increment_quotient=(v-1.9017426111009554)/sigma,key_count=len(keys),selected_boxes=len(b),centers=2**power,pure_singletons=162,pure_pair_clusters=6561,seconds=time.monotonic()-st)
    return row,keys,selected

def face_keys(phases,patches):
    keys=set();mismatch=0;total=0
    for p in patches:
        j=p['axis'];other=[k for k in range(4) if k!=j];ns=np.array(p['counts']);rr=ns**.25/2;RR=rr[:,None]/A[other];cuts=[]
        for k,axis in enumerate(other):
            cc=[p['lo'][k],p['hi'][k]]
            for phi in phases:
                for R in RR[:,k]:
                    for shift in range(-3,4):
                        for x in [phi[axis]+shift-R,phi[axis]+shift+R]:
                            if p['lo'][k]+1e-12<x<p['hi'][k]-1e-12:cc.append(x)
            cuts.append(np.unique(cc))
        for ii in product(*(range(len(c)-1) for c in cuts)):
            c=np.zeros(4);c[j]=p['c'];c[other]=[(cuts[k][ii[k]]+cuts[k][ii[k]+1])/2 for k in range(3)];total+=1
            for n,r in zip(ns,rr):
                lo0=np.ceil(c-r/A-1e-10).astype(int);hi0=np.floor(c+r/A+1e-10).astype(int)
                mismatch+=int(np.prod(np.maximum(0,hi0-lo0+1)))!=n
                for side in [None,1,-1]:
                    ll=[lo0];hh=[hi0]
                    for phi in phases:
                        l=np.ceil(c-phi-r/A-1e-10).astype(int);h=np.floor(c-phi+r/A+1e-10).astype(int)
                        if side==1:h[j]=int(np.floor(c[j]-phi[j]+r/A[j]-1e-10))
                        if side==-1:l[j]=int(np.ceil(c[j]-phi[j]-r/A[j]+1e-10))
                        if np.any(h<l):l[:]=0;h[:]=-1
                        ll.append(l);hh.append(h)
                    keys.add(tuple(map(int,np.r_[np.array(ll).ravel(),np.array(hh).ravel()])))
    return keys,dict(face_keys=len(keys),face_partition_centers=total,primary_count_mismatches=int(mismatch))
