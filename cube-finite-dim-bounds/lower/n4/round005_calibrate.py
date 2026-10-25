"""True finite-epsilon 4D joint unions for selected critical directions."""
from pathlib import Path
from itertools import product
import json,time
import numpy as np
from round005_exposure import A,geometry
from round004_satellite_search import evaluate
P=Path(__file__).resolve().parent

def face_keys(phi,patches):
    keys=set();mismatch=[];centers=0
    for p in patches:
        j=p['axis'];others=[k for k in range(4) if k!=j];ns=np.array(p['counts']);rr=ns**.25/2;RR=rr[:,None]/A[others];cuts=[]
        for k,axis in enumerate(others):
            cc=[p['lo'][k],p['hi'][k]]
            for R in RR[:,k]:
                for shift in range(-3,4):
                    for x in [phi[axis]+shift-R,phi[axis]+shift+R]:
                        if p['lo'][k]+1e-12<x<p['hi'][k]-1e-12:cc.append(x)
            cuts.append(np.unique(cc))
        for ii in product(*(range(len(c)-1) for c in cuts)):
            c=np.zeros(4);c[j]=p['c'];c[others]=[(cuts[k][ii[k]]+cuts[k][ii[k]+1])/2 for k in range(3)];centers+=1
            for n,r in zip(ns,rr):
                ll0=np.ceil(c-r/A-1e-10).astype(int);hh0=np.floor(c+r/A+1e-10).astype(int)
                count=int(np.prod(np.maximum(0,hh0-ll0+1)))
                if count!=n:mismatch.append([int(n),count])
                ll1=np.ceil(c-phi-r/A-1e-10).astype(int);hh1=np.floor(c-phi+r/A+1e-10).astype(int)
                for side in [None,1,-1]:
                    ll=ll1.copy();hh=hh1.copy()
                    if side==1:hh[j]=int(np.floor(c[j]-phi[j]+r/A[j]-1e-10))
                    elif side==-1:ll[j]=int(np.ceil(c[j]-phi[j]-r/A[j]+1e-10))
                    if np.any(hh<ll):ll[:]=0;hh[:]=-1
                    keys.add(tuple(map(int,np.r_[ll0,ll,hh0,hh])))
    return keys,dict(face_partition_centers=centers,face_keys=len(keys),primary_count_mismatches=len(mismatch),mismatch_examples=mismatch[:10])

def main():
    st=time.monotonic();data=json.loads((P/'round005_exposure_results.json').read_text());patches,_,geom=geometry();rows=[];models=[]
    for run in data['runs']:
        phi=np.array(run['phi']);d=np.array(run['direction']);pool,meta=face_keys(phi,patches);models.append(dict(phi=phi.tolist(),direction=d.tolist(),prediction=run['prediction'],**meta));print('KEYS',models[-1],flush=True)
        for eps in [.01,.003,.001,.0003]:
            a=A*(1+eps*d);r,keys=evaluate(a,phi,eps,power=10,seed=20261025505,pool=pool)
            r.update(direction=d.tolist(),prediction=run['prediction'],difference_quotient=(r['ratio']-geom['baseline_joint_ratio'])/eps,face_pool_size=len(pool));rows.append(r);print('ROW',r,flush=True)
            if r['ratio']>1.9017426111009554+1e-7:
                rr,kk=evaluate(a,phi,eps,power=12,seed=20261025506,pool=keys)
                (P/'round005_positive_candidate.json').write_text(json.dumps(dict(row=rr,keys=[list(k) for k in kk]),separators=(',',':'))+'\n')
                raise RuntimeError('positive true finite-union candidate: prioritize exact certification')
    out=dict(status='actual joint unions; single-sided difference quotients are numerical evidence only',models=models,rows=rows,elapsed=time.monotonic()-st)
    (P/'round005_calibration_results.json').write_text(json.dumps(out,indent=2)+'\n');print('FINAL',out['elapsed'],flush=True)
if __name__=='__main__':main()
