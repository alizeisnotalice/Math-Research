"""Independent heterogeneous PB DP vs normal-CDF/bin probabilities.

No Gaussian source substitution: exact finite DP algorithm on explicit p,
ordinary floating point. True-cube examples use the frozen round37 kernels.
"""
from pathlib import Path
import importlib.util,json,math,hashlib,sys
import numpy as np

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'round037_quantile_cell/probe.py'
spec=importlib.util.spec_from_file_location('frozen37',SOURCE)
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
SEED=202610290381
def normal_cdf(x):return .5*math.erfc(-x/math.sqrt(2))
def normal_mass(lo,hi):
    if lo==float('-inf'):return normal_cdf(hi)
    if hi==float('inf'):return normal_cdf(-lo)
    if lo>=0:return .5*(math.erfc(lo/math.sqrt(2))-math.erfc(hi/math.sqrt(2)))
    return normal_cdf(hi)-normal_cdf(lo)
def analyze(p,label,origin):
    p=np.array(p,dtype=float);n=len(p)
    dp=np.zeros(n+1);dp[0]=1
    for i,prob in enumerate(p):
        nxt=np.zeros_like(dp)
        nxt[:i+1]+=dp[:i+1]*(1-prob)
        nxt[1:i+2]+=dp[:i+1]*prob
        dp=nxt
    mu=float(p.sum());var=float((p*(1-p)).sum());sigma=math.sqrt(var)
    third=float((p*(1-p)*(p*p+(1-p)**2)).sum());lyap=third/sigma**3
    cdf=np.cumsum(dp);priorcdf=np.r_[0.,cdf[:-1]]
    gauss=np.array([normal_cdf((k-mu)/sigma)for k in range(n+1)])
    ks=float(np.maximum(abs(cdf-gauss),abs(priorcdf-gauss)).max())
    cc=np.array([normal_cdf((k+.5-mu)/sigma)for k in range(n+1)])
    low,high,nums,den=mod.quantiles(n,8)
    bins=[]
    for b,(a,z)in enumerate(zip(low,high)):
        truth=float(dp[a:z+1].sum())if a<=z else 0.
        if a>z:
            approx=0.;finite=0
        else:
            ll=float('-inf')if a<=0 else (a-.5-mu)/sigma
            uu=float('inf')if z>=n else(z+.5-mu)/sigma
            approx=normal_mass(ll,uu);finite=int(math.isfinite(ll))+int(math.isfinite(uu))
        bins.append({'b':b,'a':a,'z':z,'PB_probability':truth,
                     'normal_half_lattice_probability':approx,
                     'absolute_error':abs(truth-approx),'finite_CDF_endpoints':finite})
    return {'n':n,'label':label,'origin':origin,'p':p.tolist(),'p_min':float(p.min()),'p_max':float(p.max()),
            'mu':mu,'variance':var,'third_absolute_centered_sum':third,'Lyapunov_ratio':lyap,
            'DP_total_mass_residual':float(abs(dp.sum()-1)),
            'Kolmogorov_jump_endpoint_error':ks,'observed_KS_over_Lyapunov':ks/lyap,
            'half_lattice_CDF_error_at_integer_k':float(abs(cdf-cc).max()),
            'max_bin_error':max(x['absolute_error']for x in bins),'bins':bins}
def main():
    if(HERE/'results.json').exists():raise RuntimeError('Refusing repeat experiment.')
    plan={'seed':SEED,'n':[31,63,255],'standalone_p_range':[.25,.75],
          'cube_source_SHA256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
          'script_SHA256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          'scope':'Explicit PB parameters, finite true-box kernels, no source-tail MC.'}
    (HERE/'registration.json').write_text(json.dumps(plan,indent=2)+'\n')
    rng=np.random.default_rng(SEED);rows=[]
    for n in plan['n']:
        patterns={'constant_quarter':np.full(n,.25),'constant_half':np.full(n,.5),
                  'constant_three_quarters':np.full(n,.75),
                  'alternating_endpoints':np.where(np.arange(n)%2,.25,.75),
                  'linear':np.linspace(.25,.75,n),'uniform_heterogeneous':rng.uniform(.25,.75,n),
                  'random_endpoints':rng.choice([.25,.75],n)}
        for label,p in patterns.items():rows.append(analyze(p,label,'standalone explicit p'))
        src=mod.Source(n,12,4)
        for offset in [.137,.391]:
            x=((offset+np.arange(n)*(math.sqrt(5)-1)/2)%1)[None,:]
            for L in [10/9,7/6,11/9]:
                plus,minus=mod.class_weights(src,x,L)
                for b in [0,3]:
                    p=plus[0,:,b]/(plus[0,:,b]+minus[0,:,b])
                    rows.append(analyze(p,f'cube_offset{offset}_L{L}_table{b}',
                                        'round37 K12/R4 actual box formula; n255 is a new finite-query extension, not old MC data'))
    out={'python':sys.executable,'numpy':np.__version__,'rows':rows,
         'scope':'Exact finite DP algorithm in ordinary float; observed normal error only, no BE constant certified by numerics.'}
    (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
    for n in plan['n']:
        z=[r for r in rows if r['n']==n]
        print(n,'cases',len(z),'max_KS',max(r['Kolmogorov_jump_endpoint_error']for r in z),
              'max_KS_over_Lyap',max(r['observed_KS_over_Lyapunov']for r in z),
              'max_bin',max(r['max_bin_error']for r in z))
    print('rows',len(rows),'DP_mass_residual',max(r['DP_total_mass_residual']for r in rows))
if __name__=='__main__':main()
