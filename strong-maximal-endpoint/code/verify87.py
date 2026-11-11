"""r087: fixed-corner geometric constants and conditional finite-N efficiency bounds."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[k]='1'
import sys,json,math,time,multiprocessing as mp
from pathlib import Path
from fractions import Fraction as Q
import mpmath as m
from verify85 import constants85,interval_log
from verify86 import mul,divide
from verify84 import add,sub,scale
from sliver79 import logi

def constants87():
    mu=tuple(map(Q,constants85()['mu_interval']));r=sub(mu,(Q(1),Q(1)));ee=Q(1,50);rm=sub(r,(ee,ee));rp=add(r,(ee,ee));one=(Q(1),Q(1));rr=divide(mul(rp,add(one,rm)),mul(rm,add(one,rp)));pp=divide((2*ee,2*ee),mul(rp,add(one,rm)));cc=divide(mul(rr,rr),mul(rm,add(one,rm)))
    den=10**60
    clean=lambda x:(Q((x[0]*den).__floor__(),den),Q((x[1]*den).__ceil__(),den))
    rr,pp,cc=map(clean,[rr,pp,cc])
    gates=dict(critical_root=Q(27,100)<r[0]<r[1]<Q(28,100),regular_range=Q(1,4)<rm[0]<rp[1]<Q(3,10),curvature_uniform=Q(133000,65910)>2,
      irregular_gap=r[0]>ee**2,ball_constant=2*Q(22,7)*Q(11,4)==Q(121,7),outer_inflation=rr[0]>1,cell_ratio=0<pp[0]<pp[1]<1)
    assert all(gates.values())
    return dict(mu_interval=list(map(str,mu)),e=str(ee),mean_gap_lower=str(ee**2),R_interval=list(map(str,rr)),P_interval=list(map(str,pp)),C_interval=list(map(str,cc)),gates=gates)
def job(arg):
    n,gamma,cs=arg;gamma=Q(gamma);ll=logi(Q(n));dhi=gamma*ll[1];kk=(dhi/Q(cs['mean_gap_lower'])).__ceil__();nn=n-kk;cu=Q(cs['C_interval'][1]);aa=Q(121,7)*dhi
    base=dict(N=n,gamma=str(gamma),qualification_D_upper=str(dhi),exceptional_count_upper=kk,C_upper=str(cu),ball_squared_radius_upper=str(dhi),pid=os.getpid())
    if not (0<kk and 2*kk<n and nn>=cu**2*aa):
        return dict(**base,certified_formula=False,status='Excluded from simplified finite formula: subset/monotonicity gate fails; not counted as a validated efficiency bound.')
    q=Q(kk,n);ent=sub((Q(0),Q(0)),add(scale(logi(q),q),scale(logi(1-q),1-q)))
    iv=add(logi(Q(2)),scale(ent,n));iv=add(iv,scale(logi(Q(4)),kk));iv=add(iv,scale(logi(cu),nn));iv=add(iv,scale(sub(logi(aa),logi(Q(nn))),Q(nn,2)))
    pr=[]
    for wp in [80,160]:
        with m.workdps(wp):
            qq=m.mpf(kk)/n;hh=-qq*m.log(qq)-(1-qq)*m.log1p(-qq);c=m.mpf(cu.numerator)/cu.denominator;a=m.mpf(aa.numerator)/aa.denominator
            val=m.log(2)+n*hh+kk*m.log(4)+nn*m.log(c)+nn*(m.log(a)-m.log(nn))/2
            pr.append(dict(working_precision=wp,log_conditional_efficiency_upper=m.nstr(val,65)))
    return dict(**base,certified_formula=True,formula='2 exp(N H(K/N)) 4^K C^(N-K) ((121/7) D/(N-K))^((N-K)/2)',
      log_efficiency_upper_interval=list(map(str,iv)),precision_rows=pr)
def main():
    out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True);cs=constants87();args=[(2**k,str(g),cs) for k in range(12,24) for g in [Q(1,100),Q(1,10),Q(1),Q(10)]];start=time.monotonic()
    with mp.get_context('spawn').Pool(12) as pool:rows=pool.map(job,args)
    z=dict(pass_all=True,constants=cs,rows=rows,configured_workers=12,actual_processes=len({r['pid'] for r in rows}),seconds=time.monotonic()-start,
      scope='Conditional on a uniform scalar action deficit D<=gamma logN; the r022 implied constant is unknown, so gamma is not an unconditional finite-N certificate for the canonical source.')
    (out/'conditional87.json').write_text(json.dumps(z,indent=2)+'\n');print(json.dumps(dict(pass_all=True,inputs=len(rows),certified_cases=sum(r['certified_formula'] for r in rows),excluded_cases=sum(not r['certified_formula'] for r in rows),actual_processes=z['actual_processes'],seconds=z['seconds'])))
if __name__=='__main__':main()
