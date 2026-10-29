from pathlib import Path
import numpy as np,json,hashlib
from scipy.special import log_ndtr,ndtri_exp
from scipy.optimize import minimize_scalar
H=Path(__file__).resolve().parent
rows=[];N=16384;seed=202610295000
for v in [16.,64.]:
 best=minimize_scalar(lambda h:-(h+log_ndtr(-(h+v/2)/np.sqrt(v))),bounds=(-v,2*v),method='bounded');fixed=float(np.exp(-best.fun));hs=np.unique(np.r_[np.arange(-3,v+4,.5),v/2-2])
 for rho in [.5,1.]:
  for sg in [-1,1]:
   for u in [0.,.5,1.,2.]:
    rng=np.random.default_rng(seed+len(rows));r=1/3;t=u*np.sqrt(v);lp=float(log_ndtr(-t));mode=rng.random(N)>.3
    # Under P, Z has mean -sqrt(v), X|Z retains the Gaussian Q conditional law.
    z=rng.standard_normal(N)-np.sqrt(v);x=rho*z+np.sqrt(1-rho*rho)*rng.standard_normal(N)
    ids=np.flatnonzero(mode);uu=-sg*ndtri_exp(np.log(rng.random(len(ids)))+lp)
    xx=r*uu+np.sqrt(1-r*r)*rng.standard_normal(len(ids));zz=rho*xx+np.sqrt(1-rho*rho)*rng.standard_normal(len(ids));x[ids]=xx;z[ids]=zz
    lb=v/2+np.sqrt(v)*z+log_ndtr(-(t-sg*r*x)/np.sqrt(1-r*r))-lp
    lq=np.logaddexp(np.log(.3),np.log(.7)+lb);ww=np.exp(-lq);order=np.argsort(lb);ls=lb[order];weights=ww[order];tail=np.r_[np.cumsum(weights[::-1])[::-1],0.];k=np.searchsorted(ls,hs,side='right');W=np.exp(hs)*tail[k]/N;j=int(W.argmax());hi=int(np.flatnonzero(hs==v/2-2)[0])
    assert np.exp(lb-lq).max()<=1/.7*(1+1e-12)
    rows.append(dict(v=v,r=r,rho=rho,gate_sign=sg,t_over_sqrtv=u,log_source_gate_mass=lp,best_grid_h=float(hs[j]),best_grid_W=float(W[j]),gain_over_fixed_exact_profile=float(W[j]/fixed),high_h=float(hs[hi]),high_W=float(W[hi]),high_gain_over_fixed_same_h=float(W[hi]/np.exp(hs[hi]+log_ndtr(-(hs[hi]+v/2)/np.sqrt(v)))),RN_mass=float(ww.mean()),normalized_gate_output_mass=float(np.exp(lb-lq).mean())))
out=dict(rows=rows,N=N,base_seed=seed,script_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='Normalized Gaussian source-gate channel exploration, not a proven HL source/channel or interval certificate. Grid-selected gains descriptive, no independent holdout.')
(H/'RESULTS.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(dict(cases=len(rows),best_profile=max(rows,key=lambda r:r['gain_over_fixed_exact_profile']),best_high=max(rows,key=lambda r:r['high_gain_over_fixed_same_h'])),indent=2))
