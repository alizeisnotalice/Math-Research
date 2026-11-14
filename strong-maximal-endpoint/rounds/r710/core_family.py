"""Exact rational common-core directories, evaluated through factorized formulas.
No exponential-size directory or high-dimensional grid is materialized.
"""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[key]='1'
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import json,time
from sympy import integer_nthroot
from mpmath import iv
from fractions import Fraction
HERE=Path(__file__).resolve().parent

def bounds(value,places=45):
    def exact(z):
        s,m,e,_=z
        return Fraction((-1 if s else 1)*m*2**max(e,0),2**max(-e,0))
    a=exact(value._mpi_[0])*10**places;b=exact(value._mpi_[1])*10**places
    lo=a.numerator//a.denominator;hi=-((-b.numerator)//b.denominator)
    def dec(v):
        sign='-' if v<0 else '';s=str(abs(v)).zfill(places+1)
        return sign+s[:-places]+'.'+s[-places:]
    return [dec(lo),dec(hi)]

def run(h):
    start=time.time();den=2**h;k=2**h;n=h*k
    root,_=integer_nthroot(2*den**h,h);num=int(root)-den
    assert (den+num)**h<=2*den**h<(den+num+1)**h
    assert 0<num<=den
    iv.dps=max(100,int(h*.31)+120)
    a=iv.mpf(num)/den;one=iv.mpf(1);two=iv.mpf(2);log2=iv.log(two);lnn=iv.log(iv.mpf(n));lnk=iv.log(iv.mpf(k))
    logunion=h*iv.log(1+iv.mpf(k)*a);logrow=h*iv.log(1+a);logsource=logrow-log2
    entropy=logrow+h*a/(1+a)*lnk
    minKL=logunion-entropy;testgain=logunion-logsource;gap=log2+h*a/(1+a)*lnk
    assert minKL.a>=0
    q=log2/n;Hb=-q*iv.log(q)-(1-q)*iv.log(1-q)
    universal=n*Hb+log2*log2+log2
    assert universal.a>gap.b
    exact_ratio=(Fraction(1)+Fraction(h*h*num,den+num))/Fraction(h+h.bit_length()-1)
    return dict(exact_gap_over_log_n=str(exact_ratio.numerator)+"/"+str(exact_ratio.denominator),h=h,k=str(k),n=str(n),row_count_log2=h*h,a_numerator=str(num),a_denominator=str(den),
        integer_row_volume_at_most_two=True,interval_precision=iv.dps,
        log_union_interval=bounds(logunion),log_source_interval=bounds(logsource),
        row_volume_interval=bounds(iv.exp(logrow)),min_mixture_KL_interval=bounds(minKL),
        true_median_log_score_interval=bounds(testgain),testing_minus_minKL_interval=bounds(gap),
        gap_over_log_n_interval=bounds(gap/lnn),universal_gap_upper_interval=bounds(universal),
        proof_status='exact factorized true LP and entropy optima by symmetry; intervals evaluate formula; no full directory enumerated',
        seconds=time.time()-start)
if __name__=='__main__':
    out=[]
    with ProcessPoolExecutor(12) as pool:
        for future in as_completed([pool.submit(run,h) for h in [2,4,8,16,32,64,128,256,512,1024]]):
            r=future.result();out.append(r);print(r['h'],r['gap_over_log_n_interval'],r['seconds'],flush=True)
            (HERE/'core_family.json').write_text(json.dumps(sorted(out,key=lambda r:r['h']),indent=2)+'\n')
