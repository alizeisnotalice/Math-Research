"""Numerical-first channel diagnostics; no true-source tail MC or Gaussian substitution."""
from pathlib import Path
import math,json,sys,hashlib,importlib.util,time
import numpy as np
from scipy.special import ndtr,ndtri

HERE=Path(__file__).resolve().parent
OTHER=HERE.parent/'round038_highdim/source_tables.py'
def gauss_channel(tau,R,y):
    edges=np.r_[-np.inf,ndtri(np.arange(1,R)/R),np.inf]
    s=math.sqrt(1-tau*tau);out=[]
    for first in range(0,len(y),256):
        yy=y[first:first+256,None]
        lo=(edges[:-1][None,:]-tau*yy)/s
        hi=(edges[1:][None,:]-tau*yy)/s
        # Use survival differences for two positive endpoints.
        masses=np.where(lo>=0,ndtr(-lo)-ndtr(-hi),ndtr(hi)-ndtr(lo))
        out.extend((R*masses.max(axis=1)).tolist())
    return np.array(out)
def binomial(n,p):
    return np.array([math.comb(n,k)*p**k*(1-p)**(n-k)for k in range(n+1)])
def cuts(n,R):
    cums=np.cumsum(np.array([math.comb(n,k)for k in range(n+1)],dtype=object))
    hi=[next(k for k,c in enumerate(cums)if R*c>=b*2**n)for b in range(1,R)]+[n]
    return list(zip([0]+[x+1 for x in hi[:-1]],hi))
def main():
    if(HERE/'results.json').exists():raise RuntimeError('Refusing replacement')
    started=time.time()
    plan={'gaussian_tau':[.2,1/3,.6,.8],
          'R':[2,4,8,16,32,64,128,256,1024,4096],
          'Y_quadrature':'trapezoid [-12,12], step .005; ordinary FP, not interval',
          'moment_q':[.5,1.,1.5,2.],
          'binary_n':[31,63,127,255],'binary_channel_p':[1/3,2/3],
          'true_source_L':'7/6','source_SHA256':hashlib.sha256(OTHER.read_bytes()).hexdigest(),
          'scope':'Numerical diagnostics first; no Gaussian replacement of true cube',
          'script_SHA256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE/'registration.json').write_text(json.dumps(plan,indent=2)+'\n')
    y=np.linspace(-12,12,4801);density=np.exp(-y*y/2)/math.sqrt(2*math.pi)
    grows=[]
    for tau in plan['gaussian_tau']:
        prev=np.ones_like(y)
        for R in plan['R']:
            H=gauss_channel(tau,R,y)
            assert np.min(H-prev)>-1e-9
            envelope=np.exp(y*y/2)/math.sqrt(1-tau*tau)
            grows.append({'tau':tau,'R':R,'moments':{str(q):float(np.trapezoid(density*H**q,y))for q in plan['moment_q']},
                          'H_at_Y_0_2_4_6':[float(H[np.argmin(abs(y-z))])for z in [0,2,4,6]],
                          'max_H_over_proposed_envelope':float((H/envelope).max())})
            prev=H
        print('gaussian',tau,'done',flush=True)
    binary=[]
    for n in plan['binary_n']:
        prior=binomial(n,.5)
        kernel=np.array([np.convolve(binomial(t,2/3),binomial(n-t,1/3))for t in range(n+1)])
        for R in [2,4,8,16,32,64,128,256,1024]:
            bins=cuts(n,R)
            likelihoods=[kernel[:,a:b+1].sum(axis=1)/prior[a:b+1].sum()for a,b in bins if a<=b]
            H=np.max(likelihoods,axis=0)
            binary.append({'n':n,'R':R,'active_bins':len(likelihoods),
                           'moments':{str(q):float(np.dot(prior,H**q))for q in plan['moment_q']}})
        H=(kernel/prior[None,:]).max(axis=1)
        binary.append({'n':n,'R':'all_singletons','active_bins':n+1,
                       'moments':{str(q):float(np.dot(prior,H**q))for q in plan['moment_q']},
                       'max_kernel_mass_residual':float(abs(kernel.sum(axis=1)-1).max()),
                       'E_H_via_reverse_channel':float(kernel.max(axis=0).sum())})
        print('binary',n,'done',flush=True)
    sys.dont_write_bytecode=True
    spec=importlib.util.spec_from_file_location('frozen38source',OTHER)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    src=mod.Source();real=[]
    for b in range(4):
        d=src.evaluate(7/6,b)
        a=d['mean_plus']+d['mean_minus'];norm=float(a.mean());p=d['mean_plus']/a
        eta=float(np.mean(a*p)/norm);var=float(np.mean(a*(p-eta)**2)/norm)
        real.append({'b':b,'L':'7/6','Q':src.Q,'normalizer':norm,
                     'eta_under_bin_tilt':eta,'source_eta':src.pchi,'E_p_minus_source_eta':eta-src.pchi,
                     'p_min':float(p.min()),'p_max':float(p.max()),
                     'variance_p_under_bin_tilt':var,'signal_fraction':var/(eta*(1-eta)),
                     'macro_range_signal_cap':eta*(1-eta)/((1+eta)*(2-eta)),
                     'analytic_class_mean_error':d['analytic_mean_error_abs'],
                     'empirical_FP_guard':d['floating_point_guard_abs'],
                     'qualification':'Bin-averaged true source kernels; ordinary FP+analytic midpoint error; not interval certified continuous signal.'})
        del d,a,p
    result={'python':sys.version,'numpy':np.__version__,'gaussian':grows,'binary':binary,'true_source':real,
            'seconds':time.time()-started,'scope':plan['scope']}
    (HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    print('real',[(r['b'],r['signal_fraction'])for r in real],flush=True)
    print('seconds',result['seconds'],flush=True)
if __name__=='__main__':main()
