"""r081: exact flat-cutoff cover certificates and independent numerical checks.
Run: python soft81.py OUTPUT_DIRECTORY [--workers 12]
No random sampling. Spawned workers use one BLAS thread each.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
import sys, json, math, time, argparse, multiprocessing as mp
from fractions import Fraction as Q
from pathlib import Path
import mpmath as m
SCALE=10**100

def floorq(x): return Q((x.numerator*SCALE)//x.denominator,SCALE)
def ceilq(x): return -floorq(-x)
def interval_exp(x):
    assert x>=0
    b=0;t=x
    while t>Q(1,2):t/=2;b+=1
    term=Q(1);lo=term
    for k in range(1,81):term*=t/k;lo+=term
    nxt=term*t/81
    hi=lo+nxt/(1-t/82)
    lo,hi=floorq(lo),ceilq(hi)
    for _ in range(b):lo,hi=floorq(lo*lo),ceilq(hi*hi)
    return lo,hi

def cutoff_interval(z):
    if z<=0:return Q(0),Q(0)
    if z>=1:return Q(1),Q(1)
    if z>Q(1,2):
        lo,hi=cutoff_interval(1-z);return 1-hi,1-lo
    r=1/z-1/(1-z);lo,hi=interval_exp(r)
    return floorq(1/(1+hi)),ceilq(1/(1+lo))

def ceil_sqrt_q(x):
    a=math.isqrt(x.numerator//x.denominator)
    return a if a*a*x.denominator==x.numerator else a+1

def nodes(j):return ceil_sqrt_q(Q(2**(j+2),3))+1

def chi(z):
    if z<=0:return m.mpf(0)
    if z>=1:return m.mpf(1)
    return 1/(1+m.exp(1/z-1/(1-z)))

def source_interval():
    lo=sum((Q(2, (2*k+1)*3**(2*k+1)) for k in range(110)),Q(0))
    hi=lo+Q(2,221*3**221)/(1-Q(1,9))
    return floorq((1+5*lo)/8),ceilq((1+5*hi)/8)

def cover_job(k):
    started=time.monotonic(); eta=Q(1,2**k);jstop=k+4
    slo=Q(0);shi=Q(0)
    for j in range(jstop):
        z=Q(2**k,16*2**j);lo,hi=cutoff_interval(z);n=nodes(j)
        slo+=n*lo;shi+=n*hi
    # jstop gives m=jstop-k+4=8. For all j>=jstop:
    # chi(z_j)<=8*2^(-2^m), n_{j+1}<=2 n_j, hence tail<=2 first majorant.
    tail=Q(16*nodes(jstop),2**256)
    jflo,jfhi=source_interval();budget=[floorq(jflo*slo),ceilq(jfhi*(shi+tail))]
    values=[]
    for precision in (60,120):
        with m.workdps(precision):
            val=(1+5*m.log(2))/8*sum((nodes(j)*chi(m.mpf(2)**(k-4-j)) for j in range(jstop)),m.mpf(0))
            values.append(dict(precision=precision,prefix_budget=m.nstr(val,55),normalized_by_eta_sqrt=m.nstr(val*m.power(2,-m.mpf(k)/2),45)))
    return dict(k=k,eta=str(eta),prefix_families=jstop,tail_start=jstop,tail_nodes=nodes(jstop),
                tail_weighted_node_sum_upper=str(tail),source_modular_interval=[str(jflo),str(jfhi)],
                budget_interval=[str(v) for v in budget],precision_rows=values,
                pid=os.getpid(),seconds=time.monotonic()-started)

DIMS=[2,4,8,16,32,64,128,256,512,1024,4096,16384,65536,262144,1048576]

def dilation_job(n):
    rows=[]
    for wp in (60,120):
        with m.workdps(wp):
            eta=m.mpf(1)/n**2;tau=-m.log1p(-eta);zeta=m.log1p(m.log(2))/m.log(2);p=1+zeta*(n-1);q=1/(p-1)
            plus=m.log1p(eta)+(n-1)*m.log1p(m.log1p(eta))
            minus=tau+(n-1)*m.log1p(tau)
            crit=n*(q-m.log1p(q))
            rows.append(dict(precision=wp,log_dilation_plus=m.nstr(plus,52),log_dilation_minus=m.nstr(minus,52),
              cap_branch_increment=m.nstr(p*tau,52),tau_critical=m.nstr(crit,52),cap_branch=bool(tau<=crit),
              scaled_plus=m.nstr(n*plus,45),scaled_minus=m.nstr(n*minus,45),scaled_Hardy_increment=m.nstr(n*p*tau,45)))
    return dict(N=n,eta=str(Q(1,n*n)),rows=rows,pid=os.getpid())

def main():
    ap=argparse.ArgumentParser();ap.add_argument('output');ap.add_argument('--workers',type=int,default=12);a=ap.parse_args();out=Path(a.output);out.mkdir(parents=True,exist_ok=True)
    t=time.monotonic()
    with mp.get_context('spawn').Pool(a.workers) as pool:
        covers=pool.map(cover_job,list(range(2,25)))
        dilations=pool.map(dilation_job,DIMS)
    # Rational e bounds supporting the tail inequality 2<e<11/4 and e^2<8.
    elo=sum((Q(1,math.factorial(i)) for i in range(11)),Q(0))
    ehi=elo+Q(1,math.factorial(11))/(1-Q(1,12))
    assert 2<elo<=ehi<Q(11,4) and Q(11,4)**2<8
    data=dict(workers=a.workers,worker_pids=sorted({r['pid'] for r in covers+dilations}),seconds=time.monotonic()-t,
      arithmetic='exact rational outward rounding at 10^-100; numerical mpmath 60/120 dps',
      e_interval=[str(elo),str(ehi)],flat_weight_intervals=[dict(m=j,z=str(Q(1,2**j)),interval=[str(v) for v in cutoff_interval(Q(1,2**j))]) for j in range(1,8)],cover=covers,dilation=dilations,
      scope='Flat-cutoff true-column covers only on the explicit r080 two-dimensional curved receiver strip; no general C_N upper bound claimed.')
    (out/'soft81.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(dict(pass_all=True,workers=a.workers,worker_processes=len(data['worker_pids']),cover_scales=len(covers),dimensions=len(dilations),seconds=data['seconds'],last_normalized=covers[-1]['precision_rows'][-1]['normalized_by_eta_sqrt'])))
if __name__=='__main__':main()
