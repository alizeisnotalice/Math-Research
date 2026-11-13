"""Explicit finite-depth physical witnesses from certified rational monotone steps."""
import os
os.environ['SM_CERT_PANELS']='4096'
for name in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS']:os.environ[name]='1'
import sys,json,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from mpmath import mp,iv
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'r702'))
import certify as old
N=64;K=N-1;PANELS=4096

def outward(v,places=40):
    def exact(t):
        sign,m,e,_=t
        return Fraction((-1 if sign else 1)*m*2**max(e,0),2**max(-e,0))
    lo=exact(v._mpi_[0])*10**places;hi=exact(v._mpi_[1])*10**places
    a=lo.numerator//lo.denominator;b=-((-hi.numerator)//hi.denominator)
    def dec(i):
        s='-'if i<0 else '';d=str(abs(i)).zfill(places+1)
        return s+d[:-places]+'.'+d[-places:]
    return [dec(a),dec(b)]

def AH(H):
    value=iv.mpf(1);term=iv.mpf(1)
    for j in range(1,N):term*=iv.mpf(N-j)/H;value+=term
    return value

def main():
    start=time.monotonic()
    with ProcessPoolExecutor(12) as pool:nodes=list(pool.map(old.node,range(PANELS+1),chunksize=16))
    mp.dps=90;iv.dps=55
    c,weights,t0=old.setup();dt=(iv.mpf(old.END)-t0)/PANELS
    heights=[];mean=iv.mpf(0);budget=iv.mpf(0)
    for i in range(PANELS):
        left,right=nodes[i:i+2]
        h=(iv.mpf(left['height_mid'])+iv.mpf(right['height_mid']))/2
        hs=mp.nstr(old.midpoint(h),55);h=iv.mpf(hs)
        assert h.a>=old.decode(left['height']).b and h.b<=old.decode(right['height']).a
        mass=old.decode(left['survival'])-old.decode(right['survival']);assert mass.a>0
        mean+=h*mass
        budget+=old.phi(h)*iv.exp(-dt*i)*(-iv.expm1(-dt))
        heights.append(hs)
    last_s=nodes[-1]['height_mid'];last=iv.mpf(last_s)
    mean+=last*old.decode(nodes[-1]['survival']);budget+=old.phi(last)*iv.exp(t0-old.END)
    assert mean.a>1
    stepdef=dict(n=N,c_interval=outward(c,65),t0=old.START,T=old.END,panels=PANELS,heights=heights,tail_height=last_s,
       definition='g(t)=0 below t0; heights[i] on [t0+(T-t0)i/panels,t0+(T-t0)(i+1)/panels); tail_height on [T,infinity). Exact decimal rational heights.')
    (HERE/'steps.json').write_text(json.dumps(stepdef,indent=2)+'\n')
    depths=[2**k for k in [16,18,20,22,24,26,28,30]];results=[]
    flat= t0-iv.log(budget)
    for depth in depths:
        H=iv.mpf(depth);numerator=AH(H);physical_budget=iv.mpf(0)
        def tail_integral(t):
            return iv.exp(t0-t)*(1+t/H)**K*AH(H+t)
        left=tail_integral(t0)
        for i,hs in enumerate(heights):
            t=t0+dt*(i+1);right=tail_integral(t)
            mass=left-right;assert mass.a>0
            physical_budget+=old.phi(iv.mpf(hs))*mass;left=right
        physical_budget+=old.phi(last)*left
        score=t0+iv.log(numerator)-iv.log(physical_budget)
        crude_lower=flat-K*iv.log(1+iv.mpf(old.END)/H)
        assert score.a>crude_lower.b
        # Physical halo/input logs avoid underflow for exponentially small volumes.
        common=N*iv.log(2)-H+K*iv.log(H)-iv.log(iv.mpf(__import__('math').factorial(K)))
        results.append(dict(depth=depth,log_score_interval=outward(score),
            deep_score_minus_physical_interval=outward(flat-score),
            simple_log_score_lower_interval=outward(crude_lower),
            log_halo_volume_interval=outward(common+iv.log(numerator)),
            log_actual_orlicz_interval=outward(common-t0+iv.log(physical_budget))))
    out=dict(n=N,panels=PANELS,interval_precision=55,strict_average_interval=outward(mean),
       deep_profile_log_score_interval=outward(flat),tail_height_interval=outward(last),
       depths=results,all_heights_monotone=True,all_mean_conditions_strict=True,
       runtime_seconds=time.monotonic()-start)
    (HERE/'physical.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2),flush=True)
if __name__=='__main__':main()
