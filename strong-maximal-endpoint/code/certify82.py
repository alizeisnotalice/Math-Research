"""Directed Decimal certificates: Gamma tails, logPhi and midpoint derivative ranges.
All transcendental bounds use positive series plus explicit remainders; no Decimal.exp/ln.
"""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[k]='1'
from decimal import Decimal as D,Context,ROUND_FLOOR,ROUND_CEILING
from functools import lru_cache
from pathlib import Path
import multiprocessing as mp,json,sys,time,math
F=Context(prec=70,rounding=ROUND_FLOOR);C=Context(prec=70,rounding=ROUND_CEILING)
Z=D(0);O=D(1);T=D(2)

def add(a,b):return F.add(a[0],b[0]),C.add(a[1],b[1])
def mul(a,b):return F.multiply(a[0],b[0]),C.multiply(a[1],b[1])
def div(a,b):return F.divide(a[0],b[1]),C.divide(a[1],b[0])
def exact(x):return D(x),D(x)
@lru_cache(None)
def logpoint(x):
    if not O<=x<=D(3):raise ValueError("log domain")
    if x>T:
        a=logpoint(F.divide(x,T))[0];b=logpoint(C.divide(x,T))[1];v=logpoint(T);return F.add(a,v[0]),C.add(b,v[1])
    wl=F.divide(F.subtract(x,O),C.add(x,O));wh=C.divide(C.subtract(x,O),F.add(x,O))
    low=high=Z;tl=wl;th=wh;w2l=F.multiply(wl,wl);w2h=C.multiply(wh,wh)
    for k in range(90):
        low=F.add(low,F.divide(F.multiply(T,tl),D(2*k+1)));high=C.add(high,C.divide(C.multiply(T,th),D(2*k+1)))
        tl=F.multiply(tl,w2l);th=C.multiply(th,w2h)
    rem=C.divide(C.multiply(T,th),F.multiply(D(181),F.subtract(O,w2h)))
    return low,C.add(high,rem)

def exp_point(x):
    if x<0:raise ValueError("exponential domain")
    b=0;tl=th=x
    while th>D('.5'):tl=F.divide(tl,T);th=C.divide(th,T);b+=1
    terml=termh=O;lo=hi=O
    for k in range(1,81):
        terml=F.divide(F.multiply(terml,tl),D(k));termh=C.divide(C.multiply(termh,th),D(k));lo=F.add(lo,terml);hi=C.add(hi,termh)
    nxt=C.divide(C.multiply(termh,th),D(81));rem=C.divide(nxt,F.subtract(O,C.divide(th,D(82))));hi=C.add(hi,rem)
    for _ in range(b):lo=F.multiply(lo,lo);hi=C.multiply(hi,hi)
    return lo,hi

def exp_iv(s):return exp_point(s[0])[0],exp_point(s[1])[1]

def gamma(n,s):
    # polynomial sum and gamma density evaluated with positive interval operations
    poly=exact(1)
    for k in range(n-1,0,-1):poly=add(exact(1),div(mul(s,poly),exact(k)))
    term=exact(1)
    for k in range(1,n):term=div(mul(term,s),exact(k))
    ee=exp_iv(s)
    return div(poly,ee),div(term,ee)

def run(z):
    n=z['N'];A=z['A'];q=D(z['candidate']);side=z['side'];started=time.monotonic()
    if n<2 or A not in [2,3] or q<=0 or side not in ['upper','lower']:raise ValueError('unsupported certificate parameters')
    @lru_cache(None)
    def node(t):
        ll=logpoint(t);uu=logpoint(F.add(O,ll[0]))[0],logpoint(C.add(O,ll[1]))[1]
        s=add(exact(q),add(ll,mul(exact(n-1),uu)))
        qq,ff=gamma(n,s)
        gd=div(add(exact(n),ll),mul(exact(t),add(exact(1),ll)))
        return dict(s=s,Q=qq,pdf=ff,gprime=gd)
    lowQ,lowpdf=gamma(n,exact(q));first=div(mul(lowpdf,exact(q)),exact(n))
    modepdf=gamma(n,exact(n-1))[1]
    for power in range(5,14):
        M=2**power;step=F.divide(D(A-1),D(M));il=ih=Z
        for i in range(M):
            a=F.add(O,F.multiply(D(i),step));b=F.add(a,step);mid=F.divide(F.add(a,b),T)
            na,nb,nm=node(a),node(b),node(mid)
            fl=min(na['pdf'][0],nb['pdf'][0]);fh=max(na['pdf'][1],nb['pdf'][1])
            if na['s'][0]<=n-1<=nb['s'][1]:fh=max(fh,modepdf[1])
            dm=F.multiply(fl,nb['gprime'][0]);dM=C.multiply(fh,na['gprime'][1])
            correction=C.divide(C.multiply(C.multiply(step,step),C.subtract(dM,dm)),D(8))
            il=F.add(il,F.subtract(F.multiply(step,nm['Q'][0]),correction));ih=C.add(ih,C.add(C.multiply(step,nm['Q'][1]),correction))
        ml=F.add(F.add(lowQ[0],first[0]),il);mh=C.add(C.add(lowQ[1],first[1]),ih)
        if (side=='upper' and mh<O) or (side=='lower' and ml>O):break
    else:raise ArithmeticError(('No certificate',z,str(ml),str(mh)))
    return dict(z,partition_cells=M,evaluated_nodes=node.cache_info().currsize,mean_interval=[str(ml),str(mh)],
      mean_minus_one_interval=[str(F.subtract(ml,O)),str(C.subtract(mh,O))],positive_gamma_tail_at_q=[str(v) for v in lowQ],
      proof='monotone logPhi, positive series, directed rounding; midpoint error <= (h^2/8)(sup derivative-inf derivative)',
      pid=os.getpid(),seconds=time.monotonic()-started)

def main():
    src=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=True)
    rows=json.loads(src.read_text())['rows'];jobs=[]
    for r in rows:
        if r['N']>256:continue
        q=r['q']
        for side,v in [('lower',(math.floor(100*q)-1)/100),('upper',(math.ceil(100*q)+1)/100)]:
            jobs.append(dict(N=r['N'],A=r['A'],candidate=str(D(str(v))),side=side))
    start=time.monotonic()
    with mp.get_context('spawn').Pool(12) as pool:cert=list(pool.map(run,jobs))
    data=dict(pass_all=True,workers=12,actual_processes=len({r['pid'] for r in cert}),precision=70,log_positive_series_terms=90,exp_positive_series_terms=81,
      seed=None,seconds=time.monotonic()-start,rows=cert)
    (out/'certificates82.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(dict(pass_all=True,brackets=len(cert)//2,actual_processes=data['actual_processes'],cells=sum(r['partition_cells'] for r in cert),max_cells=max(r['partition_cells'] for r in cert),seconds=data['seconds'])))
if __name__=='__main__':main()
