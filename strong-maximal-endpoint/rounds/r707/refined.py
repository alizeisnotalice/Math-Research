"""Mills-resummed diagnostics and an independent Gamma Edgeworth cross-check.
The Gamma approximation is a floating-point cross-check, not an interval certificate.
"""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
import json,sys,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from scipy import special as sp
from scipy.stats import binom
from scipy.optimize import brentq
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'r706'))
import absolute
MU=absolute.MU;A0=(3-MU)/6;D0=(MU-1-np.log(MU))/2;AQ=(MU+3)/8-np.log(2*np.pi*MU)

def gamma_edgeworth(n,x,order=4):
    t=(x-n)/np.sqrt(n);H=[np.ones_like(t),t]
    for k in range(2,12):H.append(t*H[-1]-(k-1)*H[-2])
    poly=H[2]/(3*np.sqrt(n))
    if order>=2:poly+=(H[3]/4+H[5]/18)/n
    if order>=3:poly+=(H[4]/5+H[6]/12+H[8]/162)/n**1.5
    if order>=4:poly+=(H[5]/6+47*H[7]/480+H[9]/72+H[11]/1944)/float(n)**2
    return np.clip(sp.ndtr(t)-np.exp(-t*t/2)/np.sqrt(2*np.pi)*poly,0,1)

def full_soft(n,c,mode):
    p=sp.expit(-c);m,v=absolute.meanvar(c);M=n*m;V=n*v
    lo,hi=int(binom.ppf(1e-16,n,p)),int(binom.isf(1e-16,n,p))
    k=np.arange(lo,hi+1);wt=binom.pmf(k,n,p);wt/=sum(wt)
    def residual(offset):
        xx=M+offset-c*k;ell=-xx/n-n*np.log1p(-1/n);x=(1-1/n)*np.maximum(xx,0)
        if mode=='scipy':cdf=sp.gammainc(n,x)
        elif isinstance(mode,str) and mode.startswith('stable'):
            from stable_gamma import lower
            cdf=lower(n,x,int(mode[6:]))
        else:cdf=gamma_edgeworth(n,x,mode)
        return float(wt@(np.expm1(ell)-np.exp(ell)*cdf))
    q=brentq(residual,-min(M,4*np.sqrt(V*np.log(n))),0,xtol=2e-10)
    z=brentq(lambda z:sp.log_ndtr(-z)-np.log(np.sqrt(V)*z/n),.001,12,xtol=1e-14)
    R=np.exp(sp.log_ndtr(-z)+z*z/2+.5*np.log(2*np.pi));H=1/(1+R/z)
    an=(2+c**3*p*(1-p)*(1-2*p))/(6*v)
    T=H*(an*(z*z-1)+v*(1+(1-z*z)*R/(2*z)))
    return dict(offset=q,centered=q+n*(m-MU),residual=abs(residual(q)),
                refined_soft_error=q+np.sqrt(V)*z-T,terms=len(k),z_V=z,mean=M,variance=V)

def run(n):
    start=time.time();L=np.log(n);z,H,eta=absolute.normal(n)
    c=brentq(lambda c:absolute.derivative(c,eta),.05,4,xtol=1e-14)
    base=full_soft(n,c,'scipy');R=np.exp(sp.log_ndtr(-z)+z*z/2+.5*np.log(2*np.pi));N=z/R
    T=H*(A0*(z*z-1)+MU*(1+(1-z*z)*R/(2*z)))
    gain=(MU-1)*z*z*H*H/8
    slab=(MU*z*z/2-L/2-np.log(N)-.5*np.log(2*np.pi))/(1+N)-1/N
    reduced=T+gain+slab+D0/L
    expanded=absolute.CQ*z*z+absolute.BQ+(-4*np.log(z)+AQ)/(z*z)
    row=dict(n=n,c_eta=c,z=z,N=N,**base,refined_soft_scaled_error=L*base['refined_soft_error'],
             resummed_score_centered=-np.sqrt(MU*n)*z+reduced,
             expansion_score_centered=-np.sqrt(MU*n)*z+expanded,
             expansion_minus_resummed=(expanded-reduced),
             resummed_slab_gain=slab,asymptotic_AQ=AQ)
    if n>=16384:
        e3=full_soft(n,c,3);e4=full_soft(n,c,4)
        row.update(edgeworth3_centered=e3['centered'],edgeworth4_centered=e4['centered'],
            edgeworth_order_stability=e4['centered']-e3['centered'],
            scipy_minus_edgeworth4=base['centered']-e4['centered'],
            edgeworth_refined_soft_error=e4['refined_soft_error'])
    if n in [64,256,1024,4096,16384,65536,262144,1048576]:
        sys.path.insert(0,str(HERE.parent/'r700'));sys.path.insert(0,str(HERE.parent/'r692'));sys.path.insert(0,str(HERE.parent/'r704'))
        import joint,refine,profile
        for mod in [joint,refine,profile]:mod.MU=c;mod.P=sp.expit(-c)
        old={r['n']:r for r in json.loads((HERE.parent/'r700/refined.json').read_text())}
        ramp=refine.run(old[n]);ev,_,_,_=joint.evaluator(n,1e-16)
        rc=ev(ramp['a'],ramp['b'])[0]+n*(c-MU)
        profiles=[profile.run((n,order,18)) for order in [256,512]]
        sc=profiles[-1]['score']-MU*n
        row.update(ramp_centered=rc,ramp_success=ramp['success'],ramp_residual=ramp['residual'],
          nonlinear_centered=sc,quadrature_stability=abs(profiles[1]['score']-profiles[0]['score']),
          mean_residual=profiles[-1]['average_residual'],
          nonlinear_minus_resummed=sc-row['resummed_score_centered'],
          nonlinear_minus_expanded=sc-row['expansion_score_centered'],
          scaled_expansion_remainder=z*z*(sc-row['expansion_score_centered']),
          scaled_slab_remainder=L*((rc-base['centered'])-slab))
    row['seconds']=time.time()-start
    return row
if __name__=='__main__':
    jobs=[64,128,256,512,1024,4096,16384,65536,262144,1048576,4194304,16777216,33554432,67108864,134217728,268435456]
    out=[]
    with ProcessPoolExecutor(12) as pool:
        for future in as_completed([pool.submit(run,n) for n in jobs]):
            row=future.result();out.append(row);print(json.dumps(row),flush=True)
            (HERE/'refined.json').write_text(json.dumps(sorted(out,key=lambda r:r['n']),indent=2)+'\n')
