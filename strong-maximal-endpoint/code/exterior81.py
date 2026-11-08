"""Independent positive-series evaluation of the N-dimensional product exterior halo.
The integral is numerical evidence; rational normalized-volume intervals are certificates.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
import json,sys,time,multiprocessing as mp
from fractions import Fraction as Q
from pathlib import Path
import mpmath as m
from soft81 import floorq,ceilq
DIMS=[2,4,8,16,32,64,128,256,512,1024,4096,16384,65536,262144,1048576]

def ln2_interval():
    lo=sum((Q(2,(2*k+1)*3**(2*k+1)) for k in range(110)),Q(0))
    hi=lo+Q(2,221*3**221)/(1-Q(1,9))
    return floorq(lo),ceilq(hi)

def positive_series_interval(n):
    tl,th=ln2_interval();bl=bh=Q(1);sl=sh=Q(1)
    for k in range(40):
        bl=floorq(bl*Q(n+k,k+1)*tl**2/((2*n+2*k+1)*(2*n+2*k+2)))
        bh=ceilq(bh*Q(n+k,k+1)*th**2/((2*n+2*k+1)*(2*n+2*k+2)))
        sl=floorq(sl+bl);sh=ceilq(sh+bh)
    rb=th**2/(2*(2*n+1))
    first=ceilq(bh*Q(n+40,41)*th**2/((2*n+81)*(2*n+82)))
    return sl,ceilq(sh+first/(1-rb))

def job(n):
    started=time.monotonic();iv=positive_series_interval(n);rows=[]
    for wp in (60,120):
        with m.workdps(wp):
            t=m.log(2);eta=m.mpf(1)/n**2
            b=m.mpf(1);coeff=[b]
            for k in range(59):
                b*=m.mpf(n+k)/(k+1)*t*t/((2*n+2*k+1)*(2*n+2*k+2));coeff.append(b)
            r=sum(coeff); ymax=-2*n*m.log1p(-m.log1p(eta)/t)
            def integrand(v):
                if v==0:return ymax*sum(coeff[k]*(1+m.mpf(k)/n) for k in range(len(coeff)))/r
                if v==1:return m.mpf(0)
                y=ymax*v;s=t*m.exp(-y/(2*n));z=m.expm1(t-s)/eta
                if z<=0:cut=0
                elif z>=1:cut=1
                else:cut=1/(1+m.exp(1/z-1/(1-z)))
                density=sum(coeff[k]*(1+m.mpf(k)/n)*m.exp(-(1+m.mpf(k)/n)*y) for k in range(len(coeff)))/r
                return ymax*(1-cut)*density
            deficit=m.quad(integrand,[0,m.mpf(1)/4,m.mpf(1)/2,m.mpf(3)/4,1])
            logv=2*n*m.log(t)-m.loggamma(2*n+1)+m.log(r)
            rows.append(dict(precision=wp,normalized_volume_series=m.nstr(r,55),log_exterior_volume=m.nstr(logv,55),
              log_exterior_ratio_upper=m.nstr(logv+(n-1)*t,55),smooth_to_hard_volume=m.nstr(1-deficit,55),
              scaled_smoothing_deficit=m.nstr(n*deficit,55),log_modular_relative_error_upper=m.nstr(n*m.log(2*t)-m.loggamma(n),55)))
    return dict(N=n,normalized_volume_interval=[str(x) for x in iv],series_terms_certificate=41,series_terms_numeric=60,
                rows=rows,pid=os.getpid(),seconds=time.monotonic()-started)

def main():
    out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True);start=time.monotonic()
    with mp.get_context('spawn').Pool(12) as p:rows=p.map(job,DIMS)
    data=dict(workers=12,worker_pids=sorted({r['pid'] for r in rows}),seconds=time.monotonic()-start,
      seed=None,precision=[60,120],quadrature='mpmath tanh-sinh on four fixed subintervals; positive 60-term density series',
      certificate='41-term exact rational positive series and geometric tail with outward rounding',rows=rows)
    (out/'exterior81.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(dict(dimensions=len(rows),workers=len(data['worker_pids']),seconds=data['seconds'],last=rows[-1]['rows'][-1])))
if __name__=='__main__':main()
