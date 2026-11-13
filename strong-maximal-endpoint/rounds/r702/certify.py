"""Interval primal-dual enclosure at n=64 using monotone likelihood endpoints."""
import json,time,math,os
from functools import lru_cache
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from mpmath import mp,iv
N=64
START='72.2564147337939'
END=250
PANELS=int(os.environ.get("SM_CERT_PANELS","8192"))

def setup():
    mp.dps=90;iv.dps=85
    c0=1+mp.lambertw(1/mp.e)
    lo=mp.floor(c0*10**65)/10**65;hi=lo+mp.mpf(10)**-65
    ls,hs=mp.nstr(lo,85),mp.nstr(hi,85)
    left,right=iv.mpf(ls),iv.mpf(hs)
    assert left-1-iv.exp(-left)<0 and right-1-iv.exp(-right)>0
    iv.dps=55
    c=iv.mpf([ls,hs])
    p=1/(1+iv.exp(c))
    weights=[(1-p)**N]
    for k in range(N):weights.append(weights[-1]*(N-k)/(k+1)*p/(1-p))
    return c,weights,iv.mpf(START)

def distribution(t,c,weights):
    density=iv.mpf(0);survival=iv.mpf(0)
    for k,w in enumerate(weights):
        x=t-c*k
        if x.b<=0:
            survival+=w
            continue
        assert x.a>0
        ex=iv.exp(-x)
        poly=iv.mpf(1)
        for j in range(N-1,0,-1):poly=1+x*poly/j
        survival+=w*ex*poly
        density+=w*ex*x**(N-1)/math.factorial(N-1)
    return density,survival

def midpoint(interval):
    return mp.mpf(interval.mid._mpi_[0])

def solve_height(target):
    if target.b<=N:return iv.mpf(1)
    logtarget=iv.log(target)
    y0=(midpoint(logtarget)-mp.log(N))/(N-1)
    y=mp.findroot(lambda y:(N-2)*mp.log1p(y)+mp.log(N+y)-midpoint(logtarget),y0)
    rad=mp.mpf('1e-35')
    yl=max(mp.mpf(0),y-rad);yh=y+rad
    left=iv.mpf(mp.nstr(yl,65));right=iv.mpf(mp.nstr(yh,65))
    value=lambda z:(N-2)*iv.log(1+z)+iv.log(N+z)
    assert value(left).b<=logtarget.a and value(right).a>=logtarget.b
    return iv.exp(iv.mpf([mp.nstr(yl,65),mp.nstr(yh,65)]))

@lru_cache(maxsize=1)
def reference():
    c,weights,t0=setup()
    f0,_=distribution(t0,c,weights)
    return c,weights,t0,f0

def node(index):
    c,weights,t0,f0=reference()
    t=t0+(iv.mpf(END)-t0)*index/PANELS
    f,Q=distribution(t,c,weights)
    target=iv.mpf(1) if index==0 else iv.exp(t-t0)*f/f0
    height=iv.mpf(1) if index==0 else solve_height(target)
    return dict(index=index,survival=str(Q),height=str(height),target=str(target),
                height_mid=mp.nstr(midpoint(height),60),f0=str(f0))

def decode(s):
    a,b=s[1:-1].split(',')
    return iv.mpf([a.strip(),b.strip()])

def phi(h):
    return h*(1+iv.log(h))**(N-1)

def main():
    started=time.monotonic()
    with ProcessPoolExecutor(12) as pool: nodes=list(pool.map(node,range(PANELS+1),chunksize=32))
    mp.dps=70;iv.dps=55
    c,weights,t0=setup();f0=decode(nodes[0]['f0'])
    assert iv.mpf(END)>N*c+18*iv.sqrt(N*c)
    dt=(iv.mpf(END)-t0)/PANELS
    mean=iv.mpf(0);budget=iv.mpf(0);gap=iv.mpf(0);pieces=[]
    for i in range(PANELS):
        left,right=nodes[i],nodes[i+1]
        h=(iv.mpf(left['height_mid'])+iv.mpf(right['height_mid']))/2
        # One constant exact decimal height per panel.
        h=iv.mpf(mp.nstr(midpoint(h),55))
        loh=decode(left['height']);hih=decode(right['height'])
        assert h.a>=loh.b and h.b<=hih.a
        mass=decode(left['survival'])-decode(right['survival'])
        assert mass.a>=0
        nu=iv.exp(-dt*i)*(-iv.expm1(-dt))
        mean+=h*mass;budget+=phi(h)*nu
        gap+=(decode(right['target']).b-decode(left['target']).a)*(hih.b-loh.a)*nu
        pieces.append((h,nu))
    dual=budget+(1-mean)/f0-gap-iv.mpf('4e-29')
    last=iv.mpf(nodes[-1]['height_mid']);tailQ=decode(nodes[-1]['survival'])
    tailnu=iv.exp(t0-END)
    totalmean=mean+last*tailQ
    lowermean=mp.mpf(totalmean.a._mpi_[0])
    scale=max(mp.mpf(1),mp.ceil((1+mp.mpf('1e-25'))/lowermean*10**40)/10**40)
    scale=iv.mpf(mp.nstr(scale,60))
    primal=sum((phi(scale*h)*nu for h,nu in pieces),iv.mpf(0))+phi(scale*last)*tailnu
    assert scale*totalmean>1 and dual>0
    score_lower=t0-iv.log(primal)
    score_upper=t0-iv.log(dual)
    out=dict(n=N,panels=PANELS,interval_precision=55,t0=START,T=END,
             normalized_mean=str(totalmean),scale=str(scale),
             normalized_primal=str(primal),normalized_dual=str(dual),gap_bound=str(gap),
             optimal_score_lower=str(score_lower),optimal_score_upper=str(score_upper),
             score_bracket_width=str(score_upper-score_lower),
             runtime_seconds=time.monotonic()-started)
    Path('certificate.json').write_text(json.dumps(out,indent=2)+'\n')
    Path('certificate_'+str(PANELS)+'.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2),flush=True)
if __name__=='__main__':main()
