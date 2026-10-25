"""Free small satellite, two-orbit independent rectangular captures, full 4D union."""
from pathlib import Path
from itertools import product
import json,time,hashlib
import numpy as np
from scipy.stats import qmc
from scipy.optimize import minimize
from staggered_search import A,volume,boxes as old_half_boxes
P=Path(__file__).resolve().parent
OLD=1.9017426111009554

# All old rectangular subsets, translated to one full cell; no reflection shortcut.
def primary_keys():
    keys=[]
    for ell in product(range(1,5),repeat=4):
        for start in product(range(-2,2),repeat=4):
            lo=np.array([start,[0]*4]);hi=np.array([np.array(start)+ell-1,[-1]*4])
            keys.append(tuple(np.r_[lo.ravel(),hi.ravel()]))
    return keys
PRIMARY=primary_keys()

def capture(a,phi,eps,t):
    phases=np.array([np.zeros(4),phi]);weights=np.array([1.,eps]);events=[]
    for orbit in range(2):
        for axis in range(4):
            dist=np.abs((t[:,axis,None]-phases[orbit,axis]-np.arange(-2,4))*a[axis])
            events.extend(dist.T)
    rr=np.array(events).T
    # Vectorized over centers and entry radii. Cubes collect independent rectangles on both orbits.
    kl=np.ceil(t[:,None,None,:]-phases[None,None,:,:]-rr[:,:,None,None]/a-1e-11).astype(np.int16)
    kh=np.floor(t[:,None,None,:]-phases[None,None,:,:]+rr[:,:,None,None]/a+1e-11).astype(np.int16)
    count=np.prod(np.maximum(0,kh-kl+1),axis=3);mass=count@weights
    score=mass/np.maximum((2*rr)**4,1e-30);good=(score>1+1e-10)&(rr<=2.5)&(mass>0)
    score[~good]=-np.inf;large=np.where(good,rr,-np.inf)
    chosen=np.stack([np.argmax(score,axis=1),np.argmax(large,axis=1)],axis=1)
    keys=set()
    for i in range(len(t)):
        for j in chosen[i]:
            if not good[i,j]:continue
            lo=kl[i,j].copy();hi=kh[i,j].copy()
            for orb in range(2):
                if np.any(hi[orb]<lo[orb]):lo[orb]=0;hi[orb]=-1
            keys.add(tuple(map(int,np.r_[lo.ravel(),hi.ravel()])))
    return keys

def witness_boxes(a,phi,eps,keys):
    data=np.array(list(keys),dtype=int);loidx=data[:,:8].reshape(-1,2,4);hiidx=data[:,8:].reshape(-1,2,4)
    count=np.prod(np.maximum(0,hiidx-loidx+1),axis=2);mass=count@np.array([1,eps]);r=mass**.25/2
    phases=np.array([np.zeros(4),phi]);pmin=(loidx+phases)*a;pmax=(hiidx+phases)*a
    pmin[count==0]=np.inf;pmax[count==0]=-np.inf
    lo=np.maximum(0,np.max(pmax,axis=1)-r[:,None]);hi=np.minimum(a,np.min(pmin,axis=1)+r[:,None]);ok=np.all(hi>lo+1e-13,axis=1)
    b=np.stack([lo[ok],hi[ok]],axis=2).reshape(-1,8);b=np.unique(b,axis=0)
    lows=b[:,::2];highs=b[:,1::2];keep=[]
    for i in range(len(b)):
        inside=np.all(lows<=lows[i],axis=1)&np.all(highs>=highs[i],axis=1);inside[i]=False
        if not inside.any():keep.append(i)
    return b[keep]

def evaluate(a,phi,eps,power=9,seed=20261025404,pool=()):
    st=time.monotonic();t=qmc.Sobol(4,scramble=True,seed=seed).random_base2(power)
    keys=set(PRIMARY)|set(pool)|capture(a,phi,eps,t);b=witness_boxes(a,phi,eps,keys);v=volume(b)[0]/(1+eps)
    return dict(a=list(map(float,a)),phi=list(map(float,phi)),epsilon=float(eps),ratio=float(v),centers=2**power,key_count=len(keys),box_count=len(b),seconds=time.monotonic()-st),keys

def main():
    st=time.monotonic();rows=[];runs=[];candidates=[]
    def ev(a,phi,eps,power=9,pool=()):
        r,keys=evaluate(a,phi,eps,power=power,pool=pool);rows.append(r)
        print(len(rows),r,flush=True)
        if r['ratio']>OLD+1e-7:
            (P/'round004_candidate.json').write_text(json.dumps(dict(row=r,keys=[list(k) for k in keys]),indent=2)+'\n');raise RuntimeError('candidate > old+1e-7')
        if r['seconds']>5:raise RuntimeError('evaluation >5seconds; finite computational stop')
        return r,keys
    stop='finite batch completed'
    try:
        ev(A,np.zeros(4),0)
        phases=qmc.Sobol(4,scramble=True,seed=20261025404).random_base2(4)
        for eps in [.01,.05]:
            for phi in phases:
                r,k=ev(A,phi,eps);candidates.append((r,k))
        for r,k in sorted(candidates,key=lambda x:-x[0]['ratio'])[:4]:ev(A,r['phi'],.001)
        starts=sorted(candidates,key=lambda x:-x[0]['ratio'])[:2]
        for initial,pool in starts:
            eps=initial['epsilon'];x=np.r_[initial['a'],initial['phi']]
            def obj(y):return -ev(y[:4],y[4:],eps,power=7,pool=pool)[0]['ratio']
            res=minimize(obj,x,method='Nelder-Mead',bounds=[(.94*a,1.08*a) for a in A]+[(0,1)]*4,options={'maxfev':100,'xatol':1e-7,'fatol':1e-9})
            runs.append(dict(epsilon=eps,nfev=res.nfev,stop=str(res.message),ratio=float(-res.fun),x=res.x.tolist()))
        selected=sorted([r for r in rows if r['epsilon']>0],key=lambda x:-x['ratio'])[:2]
        for r in selected:
            old_capture=capture(np.array(r['a']),np.array(r['phi']),r['epsilon'],qmc.Sobol(4,scramble=True,seed=20261025404).random_base2(9))
            final,k=evaluate(np.array(r['a']),np.array(r['phi']),r['epsilon'],power=12,seed=20261025405,pool=old_capture)
            final['union_with_original_capture']=True;rows.append(final)
            if final['ratio']>OLD+1e-7:
                (P/'round004_candidate.json').write_text(json.dumps(dict(row=final,keys=[list(x) for x in k]),indent=2)+'\n');raise RuntimeError('final candidate > old+1e-7')
    except RuntimeError as e:stop=str(e)
    out=dict(status='numerical finite witness satellite search only',stop=stop,elapsed=time.monotonic()-st,rows=rows,runs=runs,best=max(rows,key=lambda r:r['ratio']))
    (P/'round004_results.json').write_text(json.dumps(out,indent=2)+'\n');print('FINAL',stop,out['best'],flush=True)
if __name__=='__main__':main()
