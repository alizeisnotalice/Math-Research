"""True capped physical-corner weak norm via r028 r=0 exact scalar reduction."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[key]='1'
import json,math,sys,time,multiprocessing as mp
from pathlib import Path
from scipy.special import gammainc,gammaincc,ndtr
from scipy.integrate import quad
from scipy.optimize import brentq
DIMS=[2,4,8,16,32,64,128,256,512,1024,4096,16384,65536,262144,1048576]

def integrals(n,A,b,need_cost=True):
    rt=math.sqrt(n);wmax=(A-1)*rt
    def terms(w):
        l=math.log1p(w/rt);v=(n-2)*math.log1p(l)+math.log(n+l);d=b*math.expm1(v/(n-1));return b+d,v-d
    points=[v for v in [1,2,4,8,16,24] if v<wmax]
    ii,err=quad(lambda w:gammaincc(n,terms(w)[0])/rt,0,wmax,epsabs=2e-13,epsrel=2e-13,points=points,limit=250)
    cc,cerr=quad(lambda w:math.exp(terms(w)[1])/rt,0,wmax,epsabs=2e-13,epsrel=2e-13,points=points,limit=250) if need_cost else (0.,0.)
    tail=(A-1)*gammaincc(n,terms(wmax)[0]) if wmax<(A-1)*rt else 0.0
    return ii-gammainc(n,b),1+cc,err,cerr,tail

def job(z):
    n,A=z;t=time.monotonic();b=brentq(lambda x:integrals(n,A,x,False)[0],1e-10,2*n,xtol=1e-11,rtol=1e-14)
    residual,scaled,err,cerr,tail=integrals(n,A,b);lc=b-math.log(scaled)
    zg=brentq(lambda z:ndtr(-z)-z/math.sqrt(n),0,math.sqrt(math.log(n)+10),xtol=1e-14);ref=n-math.sqrt(n)*zg
    return dict(N=n,A=A,activation=b,log_true_corner_norm=lc,scaled_source_modular=scaled,mean_residual=residual,
      mean_quadrature_error_estimate=err,scaled_cost_quadrature_error_estimate=cerr,omitted_mean_tail_upper=tail,
      gaussian_z=zg,gaussian_reference=ref,asymptotic_constant_diagnostic=lc-ref-zg*zg/3,
      pid=os.getpid(),seconds=time.monotonic()-t)

def main():
    out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True);start=time.monotonic()
    with mp.get_context('spawn').Pool(12) as pool:rows=pool.map(job,[(n,A) for A in [2,3] for n in DIMS])
    (out/'exact83.json').write_text(json.dumps(dict(workers=12,seconds=time.monotonic()-start,precision='float64 diagnostics; true norm equality is proved in r028; directed certificates separate',rows=rows),indent=2)+'\n')
    print(json.dumps(dict(rows=len(rows),workers=len({r['pid'] for r in rows}),seconds=time.monotonic()-start,small=rows[0],last=rows[14])))
if __name__=='__main__':main()
