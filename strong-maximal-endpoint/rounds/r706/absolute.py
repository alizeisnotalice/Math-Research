"""Deterministic soft-root optimization; exact Gamma-binomial mixture, floating point."""
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS']:
    os.environ[key]='1'
import json,time,sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from scipy import special as sp
from scipy.stats import binom
from scipy.optimize import brentq,minimize_scalar
MU=float(1+sp.lambertw(1/np.e).real);SIGMA=np.sqrt(MU)
CQ=(9-MU)/24;B0=(5*MU-6)/6;BQ=(13*MU-15)/12
HERE=Path(__file__).resolve().parent

def meanvar(c):
    p=sp.expit(-c)
    return 1+c*p,1+c*c*p*(1-p)
def derivative(c,eta):
    p=sp.expit(-c);r=p*(1-p)
    return p-c*r-eta*(2*c*r-c*c*r*(1-2*p))
def normal(n):
    z=brentq(lambda z:sp.log_ndtr(-z)-np.log(SIGMA*z/np.sqrt(n)),.001,12,xtol=1e-14)
    H=1/(1+SIGMA/(np.sqrt(n)*np.exp(-z*z/2)/np.sqrt(2*np.pi)))
    eta=z*H/(2*SIGMA*np.sqrt(n))
    return z,H,eta

def soft(n,c,tail=1e-15):
    m,v=meanvar(c);p=sp.expit(-c);M=n*m
    lo,hi=int(binom.ppf(tail,n,p)),int(binom.isf(tail,n,p))
    k=np.arange(lo,hi+1);weights=binom.pmf(k,n,p);weights/=sum(weights)
    omitted=float(binom.cdf(lo-1,n,p)+binom.sf(hi,n,p))
    def residual(offset):
        x=M+offset-c*k;ell=-x/n-n*np.log1p(-1/n)
        return float(weights@(np.expm1(ell)-np.exp(ell)*sp.gammainc(n,(1-1/n)*np.maximum(x,0))))
    offset=brentq(residual,-min(M,4*np.sqrt(n*v*np.log(n))),0,xtol=2e-10)
    return dict(centered=float(offset+n*(m-MU)),offset=float(offset),root_residual=abs(residual(offset)),
                omitted_probability=omitted,terms=len(k),mean=M,variance=n*v)

def run(n):
    start=time.time();L=np.log(n);rho=np.sqrt(L/n);z,H,eta=normal(n)
    ce=brentq(lambda c:derivative(c,eta),.05,4,xtol=1e-14)
    bounds=(max(.05,MU-2*rho),MU+.5*rho)
    opt=minimize_scalar(lambda c:-soft(n,c)['centered'],bounds=bounds,method='bounded',
                        options={'xatol':2e-10,'maxiter':120})
    best=soft(n,float(opt.x));pred=soft(n,ce);fine=soft(n,float(opt.x),1e-17)
    dn=SIGMA*np.sqrt(n)*z
    row=dict(n=n,c_soft=float(opt.x),c_eta=ce,optimizer_success=bool(opt.success),evaluations=int(opt.nfev),
             optimizer_message=str(opt.message),search_bounds=bounds,z=z,H=H,
             H_constant_diagnostic=z*z*(H*H-1),soft_constant_diagnostic=best['centered']+dn-CQ*z*z,
             transformed_soft_diagnostic=best['centered']+dn-CQ*z*z+(MU-1)/2,
             asymptotic_soft_constant=B0-(MU-1)/4,asymptotic_full_constant=BQ,
             c_eta_score_loss=best['centered']-pred['centered'],
             tail_stability=abs(best['centered']-fine['centered']),
             c_eta_soft_centered=pred['centered'],**best)
    if n in [64,1024,16384,262144]:
        sys.path.insert(0,str(HERE.parent/'r700'));sys.path.insert(0,str(HERE.parent/'r692'))
        sys.path.insert(0,str(HERE.parent/'r704'))
        import joint,refine,profile
        c=float(opt.x);p=sp.expit(-c)
        for mod in [joint,refine,profile]:mod.MU=c;mod.P=p
        old={r['n']:r for r in json.loads((HERE.parent/'r700/refined.json').read_text())}
        ramp=refine.run(old[n]);ev,_,_,_=joint.evaluator(n,1e-16)
        centered_ramp=ev(ramp['a'],ramp['b'])[0]+n*(c-MU)
        values=[profile.run((n,order,18)) for order in [256,512]]
        row.update(ramp_a=ramp['a'],ramp_b=ramp['b'],ramp_success=ramp['success'],ramp_residual=ramp['residual'],
            ramp_constant_diagnostic=centered_ramp+dn-CQ*z*z,
            nonlinear_constant_diagnostic=values[-1]['score']-MU*n+dn-CQ*z*z,
            nonlinear_order_stability=abs(values[-1]['score']-values[0]['score']),
            nonlinear_root_residual=values[-1]['average_residual'])
    row['seconds']=time.time()-start
    return row
if __name__=='__main__':
    ns=[64,128,256,512,1024,4096,16384,65536,262144,1048576,4194304,16777216]
    out=[]
    with ProcessPoolExecutor(12) as pool:
        for fut in as_completed([pool.submit(run,n) for n in ns]):
            row=fut.result();out.append(row);print(json.dumps(row),flush=True)
            (HERE/'absolute.json').write_text(json.dumps(sorted(out,key=lambda r:r['n']),indent=2)+'\n')
