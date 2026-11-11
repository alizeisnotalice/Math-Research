"""r088 exact rational/log certificates, stable dyadic cover and physical-center volumes."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[k]='1'
import sys,json,math,time,multiprocessing as mp
from pathlib import Path
from fractions import Fraction as Q
import mpmath as m
from verify84 import exp_interval,logarithm_factorial,add,sub,scale,constants
from sliver79 import logi
EPS=Q(1,10**100)

def exp_at(x):
    if x>=0:return exp_interval(x)
    z=exp_interval(-x);return 1/z[1],1/z[0]
def logchi_at(x):
    if x<=-300:return x-EPS,x
    if x>=300:return -EPS,Q(0)
    z=exp_at(x);return x-logi(1+z[1])[1],x-logi(1+z[0])[0]
def logchi(delta,eta):
    ex=exp_interval(delta);sl=(ex[0]-1)/eta;sh=(ex[1]-1)/eta
    assert sl>0
    if sl>=1:return Q(0),Q(0)
    xl=-1/sl+1/(1-sl);lo=logchi_at(xl)[0]
    if sh>=1:return lo,Q(0)
    xh=-1/sh+1/(1-sh);return lo,logchi_at(xh)[1]
def truncated_exp(x):
    if x<=-300:return Q(0),EPS
    assert x<300
    return exp_at(x)
def logsum_interval(terms,tail_hi):
    peak=max(t[0] for t in terms);assert tail_hi<=peak-300
    lo=Q(0);hi=EPS
    for a,b in terms:
        lo+=truncated_exp(a-peak)[0];hi+=truncated_exp(b-peak)[1]
    assert lo>0
    return peak+logi(lo)[0],peak+logi(hi)[1]
def grid_subdivisions(j):
    if j%2==0:return 256*2**(j//2)
    a=2*(256*2**(j//2))**2;z=math.isqrt(a);assert z*z<a
    return z+1

def dyadic_job(n):
    start=time.monotonic();k=n.bit_length()-1;assert n==2**k;eta=Q(1,n*n);jj=3*k+20;ln=logi(Q(n));terms=[]
    for j in range(jj+1):
        d=Q(1,2**(j+2));mm=grid_subdivisions(j);assert Q(mm*mm)>=128**2/d and Q((mm-1)**2)<128**2/d if j%2 else Q(mm*mm)==128**2/d
        terms.append(add(scale(logi(Q(mm+1)),n),logchi(d,eta)))
    v=Q(n*2**23);assert v==eta/Q(1,2**(jj+3));tail=add(logi(Q(2)),(Q(2),Q(2)));tail=add(tail,scale(logi(Q(130)),n));tail=sub(tail,scale(logi(eta),Q(n,2)));tail=add(tail,scale(logi(v),Q(n,2)));tail=sub(tail,(v/2,v/2))
    assert v>(n+2)*logi(Q(2))[1]
    full=logsum_interval(terms,tail[1]);lf,fdata=logarithm_factorial(n//2);lower=[]
    for b in [2,3]:
        z=add(lf,scale(ln,Q(n*(b-1),2)-b));z=sub(z,scale(logi(Q(96)),n));z=sub(z,scale(logi(Q(132,7)),Q(n,2)));z=sub(z,logi(Q(3)));z=sub(z,scale(logi(Q(21,20)),n))
        if b==3:z=sub(z,add(logi(Q(2)),(Q(2*n,3),Q(2*n,3))))
        lower.append(dict(depth_power=b,log_budget_lower_interval=list(map(str,z))))
    us=scale(add(logi(Q(184)),ln),n);uf=scale(add(logi(Q(132)),scale(ln,Q(3,2))),n);assert full[1]<uf[0]
    nums=[]
    for wp in [80,160]:
        with m.workdps(wp):
            vals=[]
            for j in range(jj+1):
                d=m.mpf(2)**(-j-2);s=m.expm1(d)*n*n
                if s>=1:lq=m.mpf(0)
                else:
                    x=-1/s+1/(1-s);lq=min(x,m.mpf(0))-m.log1p(m.exp(-abs(x)))
                vals.append(n*m.log(grid_subdivisions(j)+1)+lq)
            mx=max(vals);ls=mx+m.log(m.fsum(m.exp(z-mx) for z in vals));peakj=vals.index(mx)
            nums.append(dict(working_precision=wp,finite_grid_log_sum=m.nstr(ls,65),peak_j=peakj,log_tail_upper=m.nstr(m.log(2)+2+n*m.log(130)-n*m.log(m.mpf(1)/(n*n))/2+n*m.log(n*2**23)/2-m.mpf(n*2**23)/2,65)))
    return dict(N=n,j_last=jj,eta=str(eta),factorial=fdata,lower_rows=lower,log_constructive_full_cost_interval=list(map(str,full)),
      log_tail_upper_interval=list(map(str,tail)),log_saturation_upper_interval=list(map(str,us)),log_full_upper_interval=list(map(str,uf)),precision_rows=nums,pid=os.getpid(),seconds=time.monotonic()-start)

def log_large_integer(z):
    k=z.bit_length()-1;d=2**192
    if k>=192:
        lead=z>>(k-192);iv=(Q(lead,d),Q(lead+1,d))
    else:iv=(Q(z,2**k),Q(z,2**k))
    return add((logi(iv[0])[0],logi(iv[1])[1]),scale(logi(Q(2)),k))
def log_large_fraction(z):
    assert z>0
    return sub(log_large_integer(z.numerator),log_large_integer(z.denominator))

def center_job(arg):
    n,ds=arg;delta=Q(ds);siv=logi(1+delta);lo=siv[0]**n/Q(math.factorial(n));hi=siv[1]**n/Q(math.factorial(n));tl=lo;th=hi
    for k in range(1,81):
        tl*=siv[0]*Q(n+k-1,k*(n+k));th*=siv[1]*Q(n+k-1,k*(n+k));lo+=tl;hi+=th
    nxt=th*siv[1]*Q(n+80,81*(n+81));tail=nxt/(1-siv[1]/82);hi+=tail;li=(log_large_fraction(lo)[0],log_large_fraction(hi)[1]);assert tail/lo<EPS;rr=[]
    for r in [Q(0),Q(1,4),Q(1)]:
        logv=sub(li,scale(logi(1+r),n));vals=[]
        for wp in [80,160]:
            with m.workdps(wp):
                ss=m.log1p(m.mpf(delta.numerator)/delta.denominator);scaled=m.quad(lambda u:m.exp(ss*u)*u**(n-1),[0,m.mpf(1)/2,1]);lv=n*m.log(ss)-m.loggamma(n)+m.log(scaled)-n*m.log1p(m.mpf(r.numerator)/r.denominator)
                vals.append(dict(working_precision=wp,log_center_volume=m.nstr(lv,65)))
        rr.append(dict(r=str(r),log_unit_rectangle_center_volume_interval=list(map(str,logv)),precision_rows=vals))
    return dict(N=n,delta=str(delta),positive_terms=81,relative_series_tail_upper=str(EPS),tail_relative_bound_certified=True,rows=rr,pid=os.getpid())
def main():
    out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True);start=time.monotonic();cs=constants();assert Q(2*130**2)<Q(184**2) and Q(8)*(Q(3,4)**8)<1 and Q(512,27)>10
    with mp.get_context('spawn').Pool(12) as pool:
        dr=pool.map(dyadic_job,[2**k for k in range(10,22)])
        (out/'dyadic_phase88.json').write_text(json.dumps(dr,indent=2)+'\n');print('dyadic phase complete',flush=True)
        cr=pool.map(center_job,[(2**k,str(d)) for k in range(1,9) for d in [Q(1,(2**k)**2),Q(1,10),Q(1,2)]])
    z=dict(pass_all=True,configured_workers=12,actual_dyadic_processes=len({r['pid'] for r in dr}),actual_center_processes=len({r['pid'] for r in cr}),seconds=time.monotonic()-start,
      constants=dict(prior_r084=cs,small_exp_cutoff=-300,small_exp_upper=str(EPS),grid_constant=130,saturation_upper_constant=184,full_upper_constant=132),dyadic_rows=dr,center_rows=cr,
      scope='Exact affine-patch reflected-center budget bounds and exact outward-certificate volume; no global input allocation theorem or new C_N bound.')
    (out/'certificates88.json').write_text(json.dumps(z,indent=2)+'\n');print(json.dumps({k:v for k,v in z.items() if k not in ['dyadic_rows','center_rows','constants']},indent=2))
if __name__=='__main__':main()
