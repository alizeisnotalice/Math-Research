"""Three-candidate exact bin search for Gaussian auxiliary entropy; no R-sized arrays."""
from pathlib import Path
import json,math,time,hashlib
import numpy as np
from scipy.special import ndtr,ndtri
HERE=Path(__file__).resolve().parent
def three_candidates(tau,R,y):
    s=math.sqrt(1-tau*tau)
    k=np.minimum(R-1,np.floor(R*ndtr(y/tau)).astype(np.int64))
    indices=np.stack([np.maximum(0,k-1),k,np.minimum(R-1,k+1)],axis=1)
    lo=(ndtri(indices/R)-tau*y[:,None])/s
    hi=(ndtri((indices+1)/R)-tau*y[:,None])/s
    masses=np.where(lo>=0,ndtr(-lo)-ndtr(-hi),ndtr(hi)-ndtr(lo))
    return R*masses.max(axis=1)
def main():
    target=HERE/'gaussian_order_results.json'
    if target.exists():raise RuntimeError('Refusing replacement')
    started=time.time();y=np.linspace(-12,12,9601)
    w=np.exp(-y*y/2)/math.sqrt(2*math.pi)
    reg={'tau':[.2,1/3,.6,.8],'R_powers_of_two':list(range(1,25)),
         'Y_grid':'[-12,12], step .0025, ordinary trapezoid',
         'scope':'Gaussian auxiliary only; deterministic integration, no MC, not interval certified',
         'prior_results_SHA256':hashlib.sha256((HERE/'results.json').read_bytes()).hexdigest()}
    (HERE/'gaussian_order_registration.json').write_text(json.dumps(reg,indent=2)+'\n')
    rows=[];old=json.loads((HERE/'results.json').read_text())['gaussian']
    for tau in reg['tau']:
        s=math.sqrt(1-tau*tau);coefficient=2*tau/(s*math.sqrt(math.pi))
        for k in reg['R_powers_of_two']:
            R=2**k;H=three_candidates(tau,R,y)
            EH=float(np.trapezoid(w*H,y));a=float(ndtri(1-1/R))
            match=[x for x in old if x['tau']==tau and x['R']==R]
            rows.append({'tau':tau,'R':R,'a_R':a,'E_H':EH,
                         'E_H_over_sqrt_logR':EH/math.sqrt(math.log(R)),
                         'proposed_limit_coefficient':coefficient,
                         'upper_3_plus_linear_a':3+2*tau*a/(s*math.sqrt(2*math.pi)),
                         'prior_full_bin_result_abs_difference':abs(EH-match[0]['moments']['1.0'])if match else None})
        print(tau,'last',rows[-1],flush=True)
    # Independent all-bin reference for a coarser finite test grid, pointwise.
    maxerr=0.
    yy=np.linspace(-8,8,3201)
    for tau in reg['tau']:
        for R in [2,4,16,256,4096]:
            edges=np.r_[-np.inf,ndtri(np.arange(1,R)/R),np.inf]
            ref=[];s=math.sqrt(1-tau*tau)
            for first in range(0,len(yy),128):
                yb=yy[first:first+128,None]
                lo=(edges[:-1]-tau*yb)/s;hi=(edges[1:]-tau*yb)/s
                mass=np.where(lo>=0,ndtr(-lo)-ndtr(-hi),ndtr(hi)-ndtr(lo))
                ref.extend((R*mass.max(axis=1)).tolist())
            err=float(abs(three_candidates(tau,R,yy)-ref).max())
            assert err<1e-10;maxerr=max(maxerr,err)
    refined=[]
    yf=np.linspace(-12,12,19201);wf=np.exp(-yf*yf/2)/math.sqrt(2*math.pi)
    for tau in reg['tau']:
        v=float(np.trapezoid(wf*three_candidates(tau,2**24,yf),yf))
        coarse=next(x['E_H']for x in rows if x['tau']==tau and x['R']==2**24)
        refined.append({'tau':tau,'R':2**24,'refined_E_H':v,'coarse_abs_difference':abs(v-coarse)})
    out={'rows':rows,'three_candidate_vs_full_bin_max_abs_error':maxerr,
         'refinement_half_step':refined,'seconds':time.time()-started,
         'qualification':reg['scope']}
    target.write_text(json.dumps(out,indent=2)+'\n')
    print('max_pointwise_error',maxerr,'seconds',out['seconds'],flush=True)
if __name__=='__main__':main()
