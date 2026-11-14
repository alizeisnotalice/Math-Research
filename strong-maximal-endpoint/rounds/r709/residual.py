"""Local one-step comparisons for the median-wall potential.
Floating point tests only; no starting dimension is certified by these samples.
"""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[key]='1'
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import json,time
import numpy as np
from scipy import special as sp
from scipy.optimize import minimize_scalar
HERE=Path(__file__).resolve().parent
MU=float(1+sp.lambertw(1/np.e).real);SIG=np.sqrt(MU);CQ=(9-MU)/24;G=(MU-1)/8;LAM=1/32
RULES={j:sp.roots_laguerre(j) for j in [64,128]}

def run(job):
    t,u,offset,wall=job;start=time.time();K=G+offset;L=np.log(t);n0=t/2;C=6 if wall else 0;M=1 if wall else 0
    b=SIG*np.sqrt(t)*u
    def potential(time,x):
        z=x/(SIG*np.sqrt(time));f=np.exp(-z*z/2)/(SIG*np.sqrt(2*np.pi*time))
        w=2-z-np.exp(z)
        return sp.ndtr(z)-(CQ*z*z+C+K*np.log(time/n0))*f-M*L*f*(n0/time)**LAM*w
    terminal=potential(t,b);density=np.exp(-u*u/2)/(SIG*np.sqrt(2*np.pi*t))
    def value(c,order=128):
        nodes,weights=RULES[order];p=sp.expit(-c)
        avg=(1-p)*(weights@potential(t-1,b+MU-nodes))+p*(weights@potential(t-1,b+MU-nodes-c))
        return float(t*(avg-terminal)/density)
    opt=minimize_scalar(value,bounds=(.05,4),method='bounded',options={'xatol':1e-9,'maxiter':100})
    candidates=[(float(opt.x),float(opt.fun)),(0,value(0)),(10,value(10))]
    c,res=min(candidates,key=lambda a:a[1]);v64=value(c,64)
    zp=b/(SIG*np.sqrt(t-1));Ap=CQ*zp*zp+C+K*np.log((t-1)/n0);Bp=M*L*(n0/(t-1))**LAM
    wp=2-zp-np.exp(zp);wprime=-1-np.exp(zp)
    derivative_ratio=1+(zp*Ap-2*CQ*zp+Bp*(zp*wp-wprime))/(SIG*np.sqrt(t-1))
    return dict(t=t,u=u,K=K,fee_offset=offset,wall=wall,n0=n0,M=M,C2=C,c=c,scaled_residual=res,
       derivative_ratio_at_previous_state=float(derivative_ratio),uncorrected_asymptotic=offset,quadrature_64_128_difference=abs(v64-res),optimizer_success=bool(opt.success),
       optimizer_message=str(opt.message),evaluations=int(opt.nfev),seconds=time.time()-start)
if __name__=='__main__':
    jobs=[(t,u,k,wall) for t in [4096,16384,65536,262144,1048576,4194304,16777216]
          for u in [-4,-2,0] for k in [-.02,0,.02] for wall in [False,True]]
    out=[]
    with ProcessPoolExecutor(12) as pool:
        for fut in as_completed([pool.submit(run,j) for j in jobs]):
            row=fut.result();out.append(row)
            if row['t']==16777216:print(json.dumps(row),flush=True)
    (HERE/'residual.json').write_text(json.dumps(sorted(out,key=lambda r:(r['t'],r['u'],r['fee_offset'],r['wall'])),indent=2)+'\n')
    print('completed',len(out),'deterministic local comparisons',flush=True)
