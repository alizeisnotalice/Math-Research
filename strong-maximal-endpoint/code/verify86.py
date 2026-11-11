"""r086: actual clipped-log moment tests, directed action tolerance and mean-loss certificates."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[k]='1'
import sys,json,time,multiprocessing as mp,math
from pathlib import Path
from fractions import Fraction as Q
import mpmath as m
from scipy.integrate import quad
from verify85 import constants85,interval_log
from verify84 import add,sub,scale

def mul(a,b):
    z=[x*y for x in a for y in b];return min(z),max(z)
def divide(a,b):
    assert b[0]>0;return mul(a,(1/b[1],1/b[0]))
def mean_interval(r):return sub((Q(1),Q(1)),mul(divide(r,add((Q(1),Q(1)),r)),interval_log(r)))
def constants86():
    c=constants85();mu=tuple(map(Q,c['mu_interval']));r=sub(mu,(Q(1),Q(1)));aa=[Q(1,100),Q(1,1000),Q(1,10000),Q(1,1000000)];rows=[]
    for a in aa:
        left=mean_interval(sub(r,(a,a)));right=mean_interval(add(r,(a,a)));iv=[sub(mu,left),sub(mu,right),r];gap=(min(z[0] for z in iv),min(z[1] for z in iv));da=interval_log(divide(add(mu,(a,a)),sub(mu,(a,a))))
        assert gap[0]>0
        rows.append(dict(a=str(a),mean_loss_gap_interval=list(map(str,gap)),log_volume_inflation_per_dimension_interval=list(map(str,da))))
    assert 2*Q(rows[-1]['log_volume_inflation_per_dimension_interval'][1])<Q(1,40000)
    return dict(mu_interval=c['mu_interval'],rows=rows,epsilon='1/100',delta='1/20',tail_theta='1/100',uniform_S_tail_rate=c['uniform_S_tail_rate'],
      asymptotic_example_squared_mean_exponent='1/20000',asymptotic_example_final_success_fraction_exponent='1/100000',
      finite_N_onset_certified=False,scope='Positive mean-loss gaps and action tolerance only; the arbitrary-rectangle theorem is analytic and asymptotic.')
def job(arg):
    ds,bs=arg;b=float(Q(bs))
    with m.workdps(170):
        mu=m.findroot(lambda x:x-1-m.exp(-x),m.mpf('1.28'));dd=-(mu-1) if ds=='critical' else m.mpf(Q(ds).numerator)/Q(ds).denominator
    d=float(dd);nn=[]
    for wp in [80,160]:
        with m.workdps(wp):
            mu=m.findroot(lambda x:x-1-m.exp(-x),m.mpf('1.28'));dd=-(mu-1) if ds=='critical' else m.mpf(Q(ds).numerator)/Q(ds).denominator;bb=m.mpf(Q(bs).numerator)/Q(bs).denominator
            if dd<0:
                rr=-dd
                if rr<=1:mean=1-rr*m.log(rr)/(1+rr);mg=(1+rr**(1-bb))/((1-bb)*(1+rr))
                else:mean=2/(1+rr);mg=(rr-1+2/(1-bb))/(1+rr)
            elif dd<1:
                mean=1 if dd==0 else 1+dd*m.log(dd)/(1-dd);mg=1/(1-bb) if dd==0 else -m.expm1((1-bb)*m.log(dd))/((1-bb)*(1-dd))
            else:mean=m.mpf(0);mg=m.mpf(1)
            margin=bb*mean+2*bb**2/(1-bb)**3-m.log(mg)
            assert mean<=mu+m.mpf('1e-70') and margin>0
            nn.append(dict(working_precision=wp,mean=m.nstr(mean,65),log_mgf=m.nstr(m.log(mg),65),moment_slack=m.nstr(margin,65)))
    if d<0:
        rr=-d
        if rr<=1:
            c=-math.log(rr);qm=quad(lambda z:z*math.exp(-z),0,math.inf,epsabs=1e-12)[0]+quad(lambda z:z*math.exp(-z),c,math.inf,epsabs=1e-12)[0]
            qg=quad(lambda z:math.exp(-(1-b)*z),0,math.inf,epsabs=1e-12)[0]+quad(lambda z:math.exp(-(1-b)*z),c,math.inf,epsabs=1e-12)[0];qm/=1+rr;qg/=1+rr
        else:qm=2*quad(lambda z:z*math.exp(-z),0,math.inf,epsabs=1e-12)[0]/(1+rr);qg=(rr-1+2*quad(lambda z:math.exp(-(1-b)*z),0,math.inf,epsabs=1e-12)[0])/(1+rr)
    elif d<1:
        c=math.inf if d==0 else -math.log(d);qm=quad(lambda z:z*math.exp(-z)/(1-d),0,c,epsabs=1e-12)[0];qg=quad(lambda z:math.exp(-(1-b)*z)/(1-d),0,c,epsabs=1e-12)[0]
    else:qm=0.;qg=1.
    with m.workdps(80):errs=[float(abs(m.mpf(qm)-m.mpf(nn[-1]['mean']))),float(abs(m.mpf(math.log(qg))-m.mpf(nn[-1]['log_mgf'])))]
    assert max(errs)<1e-10
    return dict(endpoint=ds,b=bs,precision_rows=nn,scipy_mean=qm,scipy_log_mgf=math.log(qg),scipy_absolute_errors=errs,pid=os.getpid())
def main():
    out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True);cs=constants86();ds=['-100','-9','-2','-1','-3/4','critical','-1/10','-1/1000000','0','1/1000000','1/10','9/20','999/1000','999999/1000000','1000001/1000000','101/100','3','100'];jobs=[(d,b) for d in ds for b in ['1/100','1/25','1/10','1/4']];start=time.monotonic()
    with mp.get_context('spawn').Pool(12) as pool:rows=pool.map(job,jobs)
    z=dict(pass_all=True,constants=cs,rows=rows,configured_workers=12,actual_processes=len({r['pid'] for r in rows}),seconds=time.monotonic()-start,scope='Clipped-log interval means/MGFs and exact action-tolerance gates; no finite-N onset or new C_N bound.')
    (out/'moments86.json').write_text(json.dumps(z,indent=2)+'\n');print(json.dumps(dict(pass_all=True,cases=len(rows),actual_processes=z['actual_processes'],seconds=z['seconds'],max_scipy_error=max(max(r['scipy_absolute_errors']) for r in rows))))
if __name__=='__main__':main()
