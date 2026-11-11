"""r085: exact entropy-margin and Hamming source-cost certificates, local 12-process run."""
import os
for z in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[z]='1'
import json,math,sys,time,hashlib,multiprocessing as mp
from pathlib import Path
from fractions import Fraction as Q
import mpmath as m
from verify84 import exp_interval,add,sub,scale
from sliver79 import logi

def interval_log(iv):return logi(iv[0])[0],logi(iv[1])[1]
def sqrt_interval(iv):
    d=10**70
    lo=math.isqrt((iv[0]*d*d).__floor__());hi=math.isqrt((iv[1]*d*d).__floor__())+1
    return Q(lo,d),Q(hi,d)
def negexp_interval(iv):
    a=exp_interval(iv[0]);b=exp_interval(iv[1]);return 1/b[1],1/a[0]
def log_integer(z):
    bits=z.bit_length();k=bits-1;d=2**128
    lead=z>>(k-128) if k>=128 else z
    iv=(Q(lead,d),Q(lead+1,d)) if k>=128 else (Q(z,2**k),Q(z,2**k))
    out=add(interval_log(iv),scale(logi(Q(2)),k))
    return out,dict(bits=bits,binary_log_scale=k,leading128_lower_integer=str(lead),SHA256=hashlib.sha256(z.to_bytes((bits+7)//8,'big')).hexdigest())
def constants85():
    with m.workdps(160):
        mu=m.findroot(lambda x:x-1-m.exp(-x),m.mpf('1.28'));d=10**60;a=int(m.floor(mu*d))
    iv=(Q(a,d),Q(a+1,d));elo=exp_interval(iv[0]);ehi=exp_interval(iv[1])
    fl=(iv[0]-1-1/elo[0],iv[0]-1-1/elo[1]);fh=(iv[1]-1-1/ehi[0],iv[1]-1-1/ehi[1])
    assert fl[1]<0<fh[0]
    lm=interval_log(iv);critical=add(lm,sub(iv,(Q(1),Q(1))));eps=Q(1,100);delta=Q(1,20);theta=Q(1,100)
    sr=theta*delta-2*theta**2/(1-theta)**3;margin=sub(sub(lm,(delta,delta)),scale(iv,eps))
    gates=dict(root_lower_sign=fl[1]<0,root_upper_sign=fh[0]>0,entropy_gap=lm[0]>Q(24,100),critical_exponent=critical[0]>Q(52,100),
      uniform_S_tail_rate=sr>Q(1,4000),pattern_Hoeffding_rate=2*eps**2==Q(1,5000),budget_entropy_margin=margin[0]>Q(9,50))
    assert all(gates.values())
    return dict(mu_interval=list(map(str,iv)),log_mu_interval=list(map(str,lm)),critical_entropy_interval=list(map(str,critical)),
      epsilon=str(eps),delta=str(delta),theta=str(theta),uniform_S_tail_rate=str(sr),budget_entropy_margin_interval=list(map(str,margin)),gates=gates)
def job(arg):
    n,mu=arg;mu=tuple(map(Q,mu));start=time.monotonic();ln=logi(Q(n));root=sqrt_interval((mu[0]*ln[0]/n,mu[1]*ln[1]/n));c=sub(mu,scale(root,Q(1,2)));r=negexp_interval(c);p=(r[0]/(1+r[0]),r[1]/(1+r[1]));ss=scale(sqrt_interval(scale(ln,n)),2)
    kth=add(scale(p,n),ss);k=kth[1].__ceil__();assert 0<k<n//2
    b=math.comb(n,k);lb,bdata=log_integer(b);corr=logi(Q(n-k+1,n-2*k+1));upper=add(lb,corr);nr=[]
    for wp in [80,160]:
        with m.workdps(wp):
            mm=m.findroot(lambda x:x-1-m.exp(-x),m.mpf('1.28'));nn=m.mpf(n);cc=mm-m.sqrt(mm*m.log(nn)/nn)/2;rr=m.exp(-cc);pp=rr/(1+rr);hh=-pp*m.log(pp)-(1-pp)*m.log1p(-pp)
            ll=m.loggamma(nn+1)-m.loggamma(k+1)-m.loggamma(nn-k+1);up=ll+m.log(m.mpf(n-k+1)/(n-2*k+1));nr.append(dict(working_precision=wp,c=m.nstr(cc,65),r=m.nstr(rr,65),p=m.nstr(pp,65),entropy=m.nstr(hh,65),log_binomial=m.nstr(ll,65),log_sphere_upper=m.nstr(up,65),entropy_limit=m.nstr(m.log(mm)+mm-1,65),canonical_lost_mean_upper=m.nstr(2/nn**8,65)))
    return dict(N=n,k=k,c_interval=list(map(str,c)),r_interval=list(map(str,r)),p_interval=list(map(str,p)),k_threshold_interval=list(map(str,kth)),
      log_binomial_interval=list(map(str,lb)),log_sphere_upper_interval=list(map(str,upper)),canonical_lost_mean_upper=str(Q(2,n**8)),binomial_integer=bdata,
      precision_rows=nr,pid=os.getpid(),seconds=time.monotonic()-start)
def main():
    out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True);cs=constants85();start=time.monotonic();dims=[2**k for k in range(10,22)]
    with mp.get_context('spawn').Pool(12) as pool:rr=pool.map(job,[(n,cs['mu_interval']) for n in dims])
    z=dict(pass_all=True,constants=cs,rows=rr,configured_workers=12,actual_processes=len({x['pid'] for x in rr}),seconds=time.monotonic()-start,
      scope='Entropy budget for retaining original canonical rectangles; no general upper/lower improvement. Sphere bounds use exact binomial integers and a rigorous geometric ratio bound.')
    (out/'certificates85.json').write_text(json.dumps(z,indent=2)+'\n');print(json.dumps({k:v for k,v in z.items() if k not in ['constants','rows']},indent=2))
if __name__=='__main__':main()
