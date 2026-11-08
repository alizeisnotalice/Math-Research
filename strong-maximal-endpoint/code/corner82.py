"""Cap-dependent optimized radial Fenchel charge for a physical corner Hardy operator."""
import os
for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[name]='1'
import json,sys,math,time,multiprocessing as mp
from pathlib import Path
import numpy as np
from scipy.special import gammainc,gammaincc,gammaln,ndtr
from scipy.integrate import quad
from scipy.optimize import brentq
DIMS=[2,4,8,16,32,64,128,256,512,1024,4096,16384,65536,262144,1048576]

def mean(n,A,q,eps=2e-13):
    if q<=0:return 1.0,0.0
    d=(q-n)/n
    if n<32:logfirst=n*math.log(q)-q-gammaln(n+1)
    else:
        corr=1/(12*n)-1/(360*n**3)+1/(1260*n**5)-1/(1680*n**7)+1/(1188*n**9)
        logfirst=-.5*math.log(2*math.pi*n)+n*(math.log1p(d)-d)-corr
    first=math.exp(logfirst);rt=math.sqrt(n)
    wmax=min((A-1)*rt,max(24.0,(n-q)/rt+24))
    def arg(w):
        l=math.log1p(w/rt)
        return q+l+(n-1)*math.log1p(l)
    val,err=quad(lambda w:gammaincc(n,arg(w))/rt,0,wmax,epsabs=eps,epsrel=eps,limit=300,points=[v for v in [1,2,4,8,16] if v<wmax])
    tail=(A-1)*gammaincc(n,arg(wmax)) if wmax<(A-1)*rt else 0.0
    return first-gammainc(n,q)+val,err+tail

def job(z):
    n,A=z;t=time.monotonic();q=brentq(lambda x:mean(n,A,x)[0],1e-8,2*n,xtol=1e-11,rtol=1e-14)
    residual,error=mean(n,A,q);zz=brentq(lambda z:ndtr(-z)-z/math.sqrt(n),0,math.sqrt(math.log(n)+10),xtol=1e-14)
    q0=n-math.sqrt(n)*zz;ze=math.log1p(math.log(A))/math.log(A);p=1+ze*(n-1);old=n*p*math.log1p(1/(p-1))
    return dict(N=n,A=A,q=q,mean_residual=residual,quadrature_and_tail_error_estimate=error,gaussian_z=zz,gaussian_reference=q0,
      logN_coefficient_diagnostic=(q-q0)/math.log(n),old_Lp_log_charge=old,log_charge_saving=old-q,pid=os.getpid(),seconds=time.monotonic()-t)

def main():
    out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True);t=time.monotonic()
    with mp.get_context('spawn').Pool(12) as p:rows=p.map(job,[(n,A) for A in [2,3] for n in DIMS])
    (out/'corner82.json').write_text(json.dumps(dict(workers=12,seconds=time.monotonic()-t,precision='float64 diagnostic; interval certificates are separate',rows=rows),indent=2)+'\n')
    print(json.dumps(dict(rows=len(rows),workers=len({r['pid'] for r in rows}),seconds=time.monotonic()-t,small=rows[0],last=rows[14])))
if __name__=='__main__':main()
