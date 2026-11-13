"""Positive Gamma series with interval Stirling remainder; no special-function CDF."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
from fractions import Fraction
import json,time,sys
from mpmath import iv
import scipy
from scipy.special import gammainc
HERE=Path(__file__).resolve().parent

def exact_endpoint(t):
    sign,man,exponent,_=t
    return Fraction((-1 if sign else 1)*man*(2**max(exponent,0)),2**max(-exponent,0))
def bounds(value,places=36):
    scale=10**places
    a=exact_endpoint(value._mpi_[0])*scale;b=exact_endpoint(value._mpi_[1])*scale
    lower=a.numerator//a.denominator;upper=-((-b.numerator)//b.denominator)
    def decimal(k):
        sign='-' if k<0 else '';digits=str(abs(k)).zfill(places+1)
        return sign+digits[:-places]+'.'+digits[-places:]
    return [decimal(lower),decimal(upper)]
def run(job):
    n,x,shift=job;start=time.time();iv.dps=60
    num,den=x.as_integer_ratio();xx=iv.mpf(num)/den;nn=iv.mpf(n)
    assert 0<x<n
    # Alternating Stirling: 0 < remainder after -1/(360n^3) < 1/(1260n^5).
    logpref=nn*iv.log(xx/nn)+nn-xx-iv.log(2*iv.pi*nn)/2-1/(12*nn)+1/(360*nn**3)
    remainder=1/(1260*nn**5)
    logpref+=iv.mpf([-remainder.b,0]);pref=iv.exp(logpref)
    term=iv.mpf(1);total=iv.mpf(1);tol=iv.mpf('1e-27');count=0
    while True:
        count+=1;term*=xx/(n+count);total+=term
        nextterm=term*xx/(n+count+1)
        tail=nextterm/(1-xx/(n+count+2))
        if pref.b*tail.b<tol:break
        if count>1000000:raise RuntimeError('series limit')
    certified=pref*(total+iv.mpf([0,tail.b]));sc=float(gammainc(n,x));sn,sd=sc.as_integer_ratio()
    error=iv.mpf(sn)/sd-certified
    return dict(n=n,x=x,x_numerator=num,x_denominator=den,shift=shift,
                standardized_x=(x-n)/n**.5,cdf_interval=bounds(certified),
                scipy_value=sc,scipy_error_interval=bounds(error),
                series_terms=count,tail_upper=bounds(pref*tail)[1],interval_precision=60,seconds=time.time()-start)
if __name__=='__main__':
    original=json.loads((HERE/'targeted_cdf.json').read_text())
    jobs=[(r['n'],r['worst_weighted_x']+delta,delta) for r in original for delta in [-2.,0.,2.]]
    out=[]
    with ProcessPoolExecutor(12) as pool:
        for future in as_completed([pool.submit(run,j) for j in jobs]):
            row=future.result();out.append(row);print(json.dumps(row),flush=True)
            (HERE/'gamma_certificates.json').write_text(json.dumps(dict(scipy_version=scipy.__version__,mpmath_interval_precision=60,
                deterministic=True,rows=sorted(out,key=lambda x:(x['n'],x['shift']))),indent=2)+'\n')
