"""Nonrectangular two-layer stairs with phase/mass/period optimization."""
from pathlib import Path
from itertools import product
import time,json,sys,hashlib
import numpy as np
from scipy.optimize import minimize
from scipy.stats import qmc
from staggered_search import A,volume
P=Path(__file__).resolve().parent
OLD=1.9017426111009554
COUNTS=np.array(list(product(range(1,5),repeat=4)))
STAIRS=np.array(list(product(range(1,6),range(1,6),range(-2,3),range(1,5),range(1,5))))

def templates(z):
    a=np.array(z[:4]);f,h,w=z[4:];q=f*a[1]
    rows=[]
    for parity in [0,1]:
        sign=1-2*parity
        ell=COUNTS
        span1=2*a[0]*((parity+ell[:,0]-1)//2)+(parity+ell[:,0]-1)%2*2*a[0]*h-parity*2*a[0]*h
        mass=np.prod(ell[:,1:],axis=1)*(ell[:,0]+sign*w*(ell[:,0]%2))
        r=mass**.25/2
        pfirst=parity*2*a[0]*h
        min2=np.where(ell[:,0]==1,sign*q/2,-q/2)
        max2=np.where(ell[:,0]==1,sign*q/2,q/2)+(ell[:,1]-1)*a[1]
        lo=np.stack([pfirst+span1-r,max2-r,(ell[:,2]-1)*a[2]-r,(ell[:,3]-1)*a[3]-r],axis=1)
        hi=np.stack([pfirst+r,min2+r,r,r],axis=1)
        rows.append(np.stack([lo,hi],axis=2).reshape(-1,8))
        m,n,d,l3,l4=STAIRS.T
        mass=(m*(1+sign*w)+n*(1-sign*w))*l3*l4
        r=mass**.25/2
        min2=np.minimum(sign*q/2,d*a[1]-sign*q/2)
        max2=np.maximum((m-1)*a[1]+sign*q/2,(d+n-1)*a[1]-sign*q/2)
        span1=2*a[0]*(h if parity==0 else 1-h)
        lo=np.stack([pfirst+span1-r,max2-r,(l3-1)*a[2]-r,(l4-1)*a[3]-r],axis=1)
        hi=np.stack([pfirst+r,min2+r,r,r],axis=1)
        rows.append(np.stack([lo,hi],axis=2).reshape(-1,8))
    b=np.concatenate(rows);return b[np.all(b[:,1::2]>b[:,::2],axis=1)]

def full_boxes(z):
    a=np.array(z[:4]);L=a.copy();L[0]*=2
    b=templates(z);out=[]
    # each viable template has widths at most4; translated receiving boxes suffice in [-4,4].
    for j in range(4):
        new=[]
        for t in range(-4,5):
            x=b.copy();x[:,2*j:2*j+2]+=t*L[j]
            ok=(x[:,2*j+1]>0)&(x[:,2*j]<L[j]);x=x[ok]
            x[:,2*j]=np.maximum(0,x[:,2*j]);x[:,2*j+1]=np.minimum(L[j],x[:,2*j+1]);new.append(x)
        b=np.concatenate(new)
    b=np.unique(b,axis=0);lo=b[:,::2];hi=b[:,1::2];keep=[]
    for i in range(len(b)):
        inside=np.all(lo<=lo[i],axis=1)&np.all(hi>=hi[i],axis=1);inside[i]=False
        if not inside.any():keep.append(i)
    return b[keep],len(b)

def evaluate(z):
    st=time.monotonic();b,raw=full_boxes(z);v,cache=volume(b)
    return dict(z=list(map(float,z)),ratio=float(v/2),raw_boxes=raw,kept_boxes=len(b),seconds=time.monotonic()-st)

def main():
    start=time.monotonic();seed=np.r_[A,0,.5,0];rows=[];optimizations=[]
    bounds=[(.9*x,1.12*x) for x in A]+[(0,.5),(.3,.7),(-.6,.6)]
    def ev(z):
        r=evaluate(z);rows.append(r)
        if len(rows)%20==1:print(len(rows),r,flush=True)
        if r['seconds']>10:raise RuntimeError('evaluation exceeded finite time condition')
        if r['ratio']>OLD+1e-7:
            (P/'round003_candidate.json').write_text(json.dumps(r,indent=2)+'\n')
            raise RuntimeError('candidate > old +1e-7: stop for certification')
        return -r['ratio']
    reason='finite batch completed'
    try:
        ev(seed)
        sample=qmc.Sobol(3,scramble=True,seed=2026102504).random(12)
        for f,h,w in sample:ev(np.r_[A,.5*f,.3+.4*h,-.6+1.2*w])
        starts=sorted(rows,key=lambda x:-x['ratio'])[:2]
        for x in starts:
            rr=minimize(ev,x['z'],method='Powell',bounds=bounds,options={'maxfev':200,'xtol':1e-6,'ftol':1e-7})
            optimizations.append(dict(nfev=rr.nfev,ratio=float(-rr.fun),success=bool(rr.success),message=str(rr.message)))
            print('optimization',rr.message,rr.fun,rr.nfev,flush=True)
    except RuntimeError as e:reason=str(e)
    result=dict(status='numerical finite witness exploration only',stop_reason=reason,elapsed=time.monotonic()-start,exact_seed_formula=['1+u^2/2','(2+u^3*v)/4','(1+u)/2','(1+u)/2'],u='positive u^4=2',v='positive v^4=3',old_certificate_sha256=hashlib.sha256((P.parents[1]/'baseline/original_certificates.json').read_bytes()).hexdigest(),optimizations=optimizations,rows=rows,best=max(rows,key=lambda x:x['ratio']))
    (P/'round003_search_results.json').write_text(json.dumps(result,indent=2)+'\n');print('FINAL',reason,result['best'],flush=True)
if __name__=='__main__':main()
