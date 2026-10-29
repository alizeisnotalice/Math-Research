"""True periodic cube averages; empirical Haar-query profiles, not certificates."""
import json, math, platform, time
from pathlib import Path
import numpy as np
from scipy.special import gammaln, logsumexp
OUT=Path(__file__).resolve().parent

def weak(logz):
    s=np.sort(logz)[::-1]
    val=s+np.log(np.arange(1,len(s)+1)/len(s))
    j=int(np.argmax(val))
    return {'W':float(np.exp(val[j])), 'log_threshold':float(s[j]), 'hits':j+1}

def run():
    t0=time.perf_counter(); rows=[]; amp=.9; pi=.5; c=.05
    for n in (16,64,128):
        rng=np.random.default_rng(202610295100+n); x=rng.random((4096,n))
        ks=sorted(set((0,n//8,n//4,n//2)))
        logbase=[]; loggate={k:[] for k in ks}; logtilt={k:[] for k in ks}
        checks=[]
        for L in (10/9,7/6,11/9):
            y=1+amp*np.sinc(L)*np.cos(2*np.pi*x)
            p=pi+c*(y-1)/y
            dp=np.zeros((len(x),n+1)); dp[:,0]=1
            for j in range(n):
                old=dp[:,:j+1].copy()
                dp[:,1:j+2]=dp[:,1:j+2]*(1-p[:,j,None])+old*p[:,j,None]
                dp[:,0]*=1-p[:,j]
            checks.append({'L':L,'max_pmf_sum_error':float(np.max(np.abs(dp.sum(1)-1))), 'p_min':float(p.min()),'p_max':float(p.max())})
            ly=np.log(y).sum(1); logbase.append(ly)
            for k in ks:
                logmass=float(gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1)-n*math.log(2))
                loggate[k].append(ly+np.log(dp[:,k])-logmass)
                # Matched exponential source tilt; k=0 is its exact limiting product.
                z=k/n
                alpha=1+c*(z-pi)/(pi*(1-pi))
                ys=1+alpha*(y-1)
                logtilt[k].append(np.log(ys).sum(1))
        base=weak(np.max(logbase,axis=0))
        for k in ks:
            g=weak(np.max(loggate[k],axis=0)); til=weak(np.max(logtilt[k],axis=0))
            logmass=float(gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1)-n*math.log(2))
            rows.append({'n':n,'k':k,'N':len(x),'source_log_mass':logmass,'base':base,'gate':g,'matched_tilt':til,'gate_over_base':g['W']/base['W'],'gate_over_tilt':g['W']/til['W'], 'log_pointwise_tilt_cap':float(-(gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1)+ (k*math.log(k/n) if k else 0)+((n-k)*math.log(1-k/n) if k<n else 0))), 'checks':checks})
        # n-dimensional exact limiting identity: count=0 is product of rho_0/(1-pi).
        assert np.max(np.abs(np.array(loggate[0])-np.array(logtilt[0])))<1e-11
    return {'model':'rho=1+0.9cos(2pi y), rho1=.5rho+.05(rho-1); fixed count k', 'scales':[10/9,7/6,11/9], 'source_min_rho':.1, 'soft_c':c,'status':'descriptive Monte Carlo query profiles, exact conditional DP to floating precision; not asymptotic proof', 'python':platform.python_version(),'numpy':np.__version__,'elapsed_seconds':time.perf_counter()-t0,'rows':rows}
if __name__=='__main__':
    r=run(); (OUT/'RESULTS.json').write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'rows':len(r['rows']),'seconds':r['elapsed_seconds'],'max_gate_over_tilt':max(z['gate_over_tilt'] for z in r['rows']),'max_gate_over_base':max(z['gate_over_base'] for z in r['rows'])},indent=2))
