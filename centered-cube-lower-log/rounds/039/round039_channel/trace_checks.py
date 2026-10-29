"""Independent exact tiny posterior-count trace checks and finite weak diagnostics."""
from pathlib import Path
from fractions import Fraction as F
import itertools,math,json
import numpy as np
HERE=Path(__file__).resolve().parent
def pb(ps):
    out=[F(1)]
    for p in ps:
        z=[F(0)]*(len(out)+1)
        for k,v in enumerate(out):z[k]+=v*(1-p);z[k+1]+=v*p
        out=z
    return out
exact=[]
for eta,probs in [(F(1,2),[F(1,3),F(2,3)]),(F(2,5),[F(1,5),F(3,5)])]:
    f=sum((p-eta)**2 for p in probs)/2/(eta*(1-eta))
    for n in [1,3,5]:
        prior=pb([eta]*n);trace=F(0)
        for obs in itertools.product([0,1],repeat=n):
            post=pb([probs[o]for o in obs])
            trace+=sum(v*v/pi for v,pi in zip(post,prior))/2**n
        target=sum((f**k for k in range(n+1)),F(0))
        assert trace==target
        exact.append({'eta':str(eta),'posterior_p':list(map(str,probs)),
                      'n':n,'signal_fraction':str(f),'trace':str(trace),
                      'geometric_trace':str(target),'equality_exact':True})
def binomial(n,p):return np.array([math.comb(n,k)*p**k*(1-p)**(n-k)for k in range(n+1)])
weak=[]
for n in [31,63,127,255]:
    prior=binomial(n,.5)
    kernel=np.array([np.convolve(binomial(t,2/3),binomial(n-t,1/3))for t in range(n+1)])
    L=kernel/prior[None,:];H=L.max(axis=1)
    order=np.argsort(-H);tail=np.cumsum(prior[order])
    weak_value=float((H[order]*tail).max()) # inclusive-tail equivalent supremum as lambda approaches atom from below
    trace=float(np.sum(prior[:,None]*L*kernel))
    expected=sum((1/9)**k for k in range(n+1))
    assert abs(trace-expected)<1e-12
    weak.append({'n':n,'singleton_H_weak_L1':weak_value,'E_H':float(prior@H),
                 'weighted_L2_trace':trace,'geometric_trace':expected})
out={'scope':'Exact tiny Fraction channel identities plus deterministic finite DP; no MC',
     'exact_cases':exact,'binary_float':weak}
target=HERE/'trace_results.json'
if target.exists():raise RuntimeError('Refusing replacement')
target.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
