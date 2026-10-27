#!/usr/bin/env python3
"""Bounded stationary-point diagnostics; positive quadrature, no global-minimum claim."""
import os,sys,math,time,json,hashlib,resource,concurrent.futures
from pathlib import Path
os.environ.update({k:'1' for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS')})
import scipy
from scipy.integrate import quad
from scipy.optimize import brentq
import mpmath as mp
DIMS=(4,7,16,64,256,1024)

def run(n):
    start=time.monotonic(); cpu=time.process_time(); calls=0
    def evaluate(t):
        nonlocal calls
        calls+=1
        if calls>128: raise RuntimeError('128 evaluation limit')
        y=math.log(t); h=(1+1/t)/2
        I=quad(lambda s: 1-.5*(1+1/(2*t)-.5*math.exp(-s))**(n-1),0,y,epsabs=1e-10)[0]
        J=quad(lambda s: (1+1/(2*t)-.5*math.exp(-s))**(n-2),0,y,epsabs=1e-10)[0]
        B=1-h**n+n*I; A=(1+1/t)**n
        bp=n/(2*t)+n*h**(n-1)/(2*t*t)+n*(n-1)*J/(4*t*t)
        G=bp/math.sqrt(B)-n/t**2*(1+1/t)**(n/2-1)
        return (math.sqrt(A)+math.sqrt(B))**2,G
    lo=max(1.001,n/(3*math.log(n))) if n>=64 else 1.001
    hi=4*n+2
    grid=[math.exp(math.log(lo)+(math.log(hi)-math.log(lo))*i/31) for i in range(32)]
    vals=[evaluate(t)[1] for t in grid]; roots=[]
    for a,b,fa,fb in zip(grid,grid[1:],vals,vals[1:]):
        if fa*fb<0:
            t=brentq(lambda z:evaluate(z)[1],a,b,xtol=1e-12,rtol=1e-14,maxiter=50)
            roots.append({'t':repr(t),'F':repr(evaluate(t)[0]),'sign_change':('-+' if fa<0 else '+-')})
    mins=[r for r in roots if r['sign_change']=='-+']
    best=mins[0] if mins else None
    denominator=math.log(n)-3*math.log(math.log(n)); t0=n/denominator if denominator>0 else None
    out={'n':n,'pid':os.getpid(),'calls':calls,'scan_interval':[lo,hi],'roots':roots,'asymptotic_t0':t0,'asymptotic_t0_F':evaluate(t0)[0] if t0 and t0>1 else None}
    if best:
        mp.mp.dps=60; t=mp.mpf(best['t']); N=mp.mpf(n); y=mp.log(t); h=(1+1/t)/2
        ell=lambda s: 1+1/(2*t)-mp.exp(-s)/2
        B=1-h**n+N*mp.quad(lambda s: 1-ell(s)**(n-1)/2,[0,y]); A=(1+1/t)**n
        J=mp.quad(lambda s:ell(s)**(n-2),[0,y]); bp=N/(2*t)+N*h**(n-1)/(2*t*t)+N*(N-1)*J/(4*t*t)
        f=(mp.sqrt(A)+mp.sqrt(B))**2; g=bp/mp.sqrt(B)-N/t**2*(1+1/t)**(N/2-1)
        main=N*mp.log(N)-N*mp.log(mp.log(N))
        out.update({'high_precision_F':mp.nstr(f,45),'high_precision_derivative_bracket':mp.nstr(g,20),'remainder_F_minus_main_over_n':mp.nstr((f-main)/N,30),'x_n_over_t':mp.nstr(N/t,25),'difference_x_minus_log_n_plus_3loglog_n':mp.nstr(N/t-mp.log(N)+3*mp.log(mp.log(N)),25)})
    out.update({'calls':calls,'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'stop':'32-point bounded scan, brentq sign changes; <=128 calls; no global claim'})
    return out
if __name__=='__main__':
    t=time.monotonic()
    with concurrent.futures.ProcessPoolExecutor(max_workers=6) as pool: rows=list(pool.map(run,DIMS))
    out={'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'python':sys.version,'scipy':scipy.__version__,'mpmath':mp.__version__,'wall_seconds':time.monotonic()-t,'precision_note':'60-digit positive quadrature, roots found in float64; non-certified','workers':6,'results':rows}
    Path(__file__).with_name('results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
