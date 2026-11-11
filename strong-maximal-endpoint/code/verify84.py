"""r084: exact constants, factorial-volume certificates and independent high-precision diagnostics.
No N-dimensional grid or artificial CPU/memory load is required for the universal proof.
"""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[k]='1'
import sys,json,math,time,hashlib,multiprocessing as mp
from pathlib import Path
from fractions import Fraction as Q
import mpmath as m
from sliver79 import logi
DEN=10**80

def outward(lo,hi):return Q((lo*DEN).__floor__(),DEN),Q((hi*DEN).__ceil__(),DEN)
def exp_interval(x):
    if x<0:raise ValueError('positive exponential only')
    b=0;t=x
    while t>Q(1,2):t/=2;b+=1
    term=Q(1);lo=term
    for k in range(1,81):term*=t/k;lo+=term
    nxt=term*t/81;hi=lo+nxt/(1-t/82);lo,hi=outward(lo,hi)
    for _ in range(b):lo,hi=outward(lo*lo,hi*hi)
    return lo,hi

def constants():
    n0=1024
    exponents=[Q(1,2)+Q(2,n0**2),1+Q(2,n0-1),Q(1),Q(2,3)]
    intervals=[exp_interval(x) for x in exponents]
    gates=dict(v_upper=Q(3,4)+Q(9,32*n0)<Q(19,25),radius_under_tenth=Q(6,n0)<Q(1,100),
      actual_anchor_upper=Q(19,25)+Q(1,10)<Q(7,8),loss_constant=5*(1+Q(1,n0))<6,
      exterior_coordinate_upper=Q(9*n0,32*(n0-1))<1,
      source_height=intervals[0][1]<Q(33,20),moment_power_lower=intervals[1][1]<Q(200,73),
      e_upper=intervals[2][1]<Q(11,4),flat_weight=intervals[3][1]<2,
      ball_constant=12*Q(22,7)*Q(11,4)<121,moment_product=Q(33,20)*Q(127,200)==Q(4191,4000)<Q(21,20))
    if not all(gates.values()):raise ArithmeticError('constant gate failed')
    return dict(N_min=n0,gates=gates,positive_exp_terms=81,exponential_certificates=[dict(exponent=str(x),interval=[str(y) for y in v]) for x,v in zip(exponents,intervals)])

def logarithm_factorial(n):
    f=math.factorial(n);bits=f.bit_length();k=bits-1;d=2**128
    if k>=128:
        lead=f>>(k-128);tl=Q(lead,d);th=Q(lead+1,d)
    else:lead=f;tl=th=Q(f,2**k)
    l2=logi(Q(2));lo=logi(tl)[0]+k*l2[0];hi=logi(th)[1]+k*l2[1]
    return (lo,hi),dict(argument=n,bits=bits,binary_log_scale=k,leading128_lower_integer=str(lead),
                       SHA256=hashlib.sha256(f.to_bytes((bits+7)//8,'big')).hexdigest())

def add(a,b):return a[0]+b[0],a[1]+b[1]
def sub(a,b):return a[0]-b[1],a[1]-b[0]
def scale(a,x):return (x*a[0],x*a[1]) if x>=0 else (x*a[1],x*a[0])
def job(n):
    started=time.monotonic();ln=logi(Q(n));lf,fdata=logarithm_factorial(n//2);l2=logi(Q(2));rows=[]
    for p in [2,3]:
        elementary=sub(scale(sub(add(logi(Q(5)),scale(ln,Q(p,2))),logi(Q(924 if p==2 else 1848))),n),add(logi(Q(9 if p==2 else 18)),scale(ln,p)))
        exact=lf
        exact=add(exact,scale(ln,Q(n*(p-1),2)-p))
        exact=sub(exact,scale(logi(Q(132,7)),Q(n,2)))
        exact=sub(exact,scale(logi(Q(16)),n));exact=sub(exact,logi(Q(3)));exact=sub(exact,scale(logi(Q(21,20)),n));exact=sub(exact,(Q(3,4),Q(3,4)))
        if p==3:exact=sub(exact,add(l2,(Q(2*n,3),Q(2*n,3))))
        vals=[]
        for wp in [80,160]:
            with m.workdps(wp):
                nn=m.mpf(n);aa=m.mpf(1)/2+2/nn**p
                lm1=aa+nn*m.log1p(-1/(2*nn))
                ii=-m.expm1((nn+1)*m.log1p(-1/nn))*nn/(nn+1)
                lmn=nn*(aa+m.log(ii));z=max(lm1,lmn);lj=z+m.log1p(m.exp(min(lm1,lmn)-z))
                caplog= m.mpf(3)/4+nn*m.log(m.pi)/2+nn*m.log(6)/2-nn*(p-1)*m.log(nn)/2-m.loggamma(nn/2+1)
                lw=0 if p==2 else -m.log(2)-2*nn/3
                sharp=lw-nn*m.log(2)-p*m.log(nn)-caplog-nn*m.log(8)-lj
                cert=m.loggamma(nn/2+1)+nn*(p-1)*m.log(nn)/2-p*m.log(nn)-nn*m.log(m.mpf(132)/7)/2-nn*m.log(16)-m.log(3)-nn*m.log(m.mpf(21)/20)-m.mpf(3)/4+lw
                eco=nn*(m.log(5)+m.mpf(p)/2*m.log(nn)-m.log(924 if p==2 else 1848))-m.log(9 if p==2 else 18)-p*m.log(nn)
                vals.append(dict(working_precision=wp,log_elementary_budget_lower=m.nstr(eco,65),log_factorial_budget_lower=m.nstr(cert,65),
                  log_sharpened_budget_lower_numerical=m.nstr(sharp,65),log_actual_modular_upper=m.nstr(lj,65),
                  source_log_height=m.nstr(aa,65),log_physical_residual_ratio_upper=m.nstr(m.log(2)-nn*m.log(nn),65)))
        rows.append(dict(p=p,eta=str(Q(1,n*n)),parameter_layer_volume_lower=str(Q(1,2**n*n**p)) if n<=4096 else None,
          log_elementary_lower_interval=[str(x) for x in elementary],log_factorial_lower_interval=[str(x) for x in exact],precision_rows=vals))
    return dict(N=n,rows=rows,factorial=fdata,pid=os.getpid(),seconds=time.monotonic()-started)

def main():
    out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True);g=constants();start=time.monotonic();dims=[2**k for k in range(10,22)]
    with mp.get_context('spawn').Pool(12) as pool:rows=pool.map(job,dims)
    r=dict(pass_all=True,workers=12,actual_processes=len({z['pid'] for z in rows}),seconds=time.monotonic()-start,constants=g,rows=rows,
      scope='All-N argument in accompanying proof; exact factorial/log interval certificates for twelve finite dimensions. No new general C_N upper/lower.')
    (out/'certificates84.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(dict(pass_all=True,dimensions=len(rows),actual_processes=r['actual_processes'],seconds=r['seconds'],last_factorial_bits=rows[-1]['factorial']['bits'])))
if __name__=='__main__':main()
