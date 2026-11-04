import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
import json,itertools,time,sys,multiprocessing as mp
import numpy as np
from scipy.stats import qmc
from numpy.polynomial import polynomial as P
from pathlib import Path
OUT=Path(__file__).resolve().parent / 'amplitude-batch'; OUT.mkdir(exist_ok=True); NPOINTS=int(os.environ.get('SOBOL_POWER','14'));CHUNK=1024

def fit(x,boxes,rmax=1.5):
    n=2; events=[1e-8,rmax]
    for center,d in boxes:
        for z,di in zip(np.abs(x-center),d):
            for p in (abs(z-di),z+di):
                if 1e-8<p<rmax:events.append(float(p))
    ev=np.unique(events); best=1.
    for lo,hi in zip(ev[:-1],ev[1:]):
        mid=(lo+hi)/2;poly=np.zeros(3)
        for center,d in boxes:
            length=np.maximum(0.,np.minimum(x+mid,center+d)-np.maximum(x-mid,center-d))
            if np.any(length<=0):continue
            slope=(x+mid<center+d).astype(float)+(x-mid>center-d).astype(float)
            f=length-slope*mid;poly[0]+=f[0]*f[1];poly[1]+=f[0]*slope[1]+f[1]*slope[0];poly[2]+=slope[0]*slope[1]
        if not np.any(poly):continue
        # For n=2 the stationary equation r F'(r)-2 F(r)=0 is linear.
        cand=[lo,hi]
        if abs(poly[1])>1e-25:
            rt=-2*poly[0]/poly[1]
            if lo<rt<hi:cand.append(rt)
        for rad in cand:
            vol=sum(float(np.prod(np.maximum(0.,np.minimum(x+rad,center+d)-np.maximum(x-rad,center-d)))) for center,d in boxes)
            best=min(best,vol/(2*rad)**2)
    return best

def job(spec):
    idx,k,u,seed,amplitude=spec;start=time.time();rng=np.random.default_rng(seed);opt=u/(1+2*u)
    # Site-dependent independent offsets and rectangular half-widths; no triangular layering assumption.
    c0=np.array(list(itertools.product(range(k),repeat=2)),float); offsets=rng.uniform(-1.2,1.2,size=c0.shape)*u*u*amplitude
    centers=c0+offsets; dims=opt*(1+rng.uniform(-2.,2.,size=c0.shape)*u*u*amplitude)
    q=float(np.sum(np.prod(2*dims,axis=1))/k**2);counts=[];progress=[]
    for site,(center,d) in enumerate(zip(centers,dims)):
        boxes=[]
        for ctr,w in zip(centers,dims):
            for ij in itertools.product((-1,0,1),repeat=2):
                shifted=ctr+k*np.array(ij)
                if np.max(np.abs(shifted-center))<1.5+2*np.max(dims):boxes.append((shifted,w))
        pts=center+(2*qmc.Sobol(2,scramble=True,seed=seed*100+site).random_base2(NPOINTS)-1)*d;count=0
        for b in range(0,len(pts),CHUNK):
            for x in pts[b:b+CHUNK]: count+=int(fit(x,boxes)<u*u*(1-1e-11))
            progress.append(dict(site=site,points=min(b+CHUNK,len(pts)),captured=count,elapsed=time.time()-start))
            partial=dict(amplitude=amplitude,id=idx,k=k,u=u,seed=seed,pid=os.getpid(),centers=centers.tolist(),halfwidths=dims.tolist(),counts_completed=counts,progress=progress[-1],note='Floating QMC, truncated radius 1.5; not certified')
            (OUT/f'job_{idx:02d}_progress.json').write_text(json.dumps(partial))
        counts.append(count)
    cap=sum(float(np.prod(2*dims[i]))*counts[i]/len(pts) for i in range(len(dims)))/k**2
    delta=opt;q0=(2*delta)**2;baseline=q0/(1-q0);gain=cap/(1-q)
    result=dict(amplitude=amplitude,id=idx,k=k,u=u,seed=seed,pid=os.getpid(),centers=centers.tolist(),halfwidths=dims.tolist(),q=q,counts=counts,samples_per_hole=len(pts),estimated_gain=gain,baseline=baseline,estimated_difference=gain-baseline,elapsed=time.time()-start,progress=progress,note='Floating QMC; radius limited to <=1.5; estimates are not certified')
    (OUT/f'job_{idx:02d}_result.json').write_text(json.dumps(result,indent=2));return result

if __name__=='__main__':
    OUT.mkdir(exist_ok=True)
    jobs=[(i,2,[.04,.06,.08][i%3],45101+i,[1/64,1/16,1/4,1/2,1][i//3]) for i in range(15)]
    print(json.dumps({'master':os.getpid(),'jobs':jobs,'samples_per_hole':2**NPOINTS,'purpose':'resolve amplitude dependence of loss in nonidentical multi-pore perturbations, exploratory only'}),flush=True)
    with mp.Pool(15) as pool:
        for result in pool.imap_unordered(job,jobs): print(json.dumps({k:result[k] for k in ('id','pid','amplitude','estimated_difference','elapsed')}),flush=True)
