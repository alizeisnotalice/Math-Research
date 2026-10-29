#!/usr/bin/env python3
"""Read-only independent replay of round042_selector's frozen experiment."""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal, localcontext
from itertools import product
import importlib.util
import hashlib
import json
import time
import numpy as np
import audit

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'round042_selector'

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def exact_tables(c):
    K,R=c['K'],c['R']
    eps,delta=F(str(c['epsilon'])),F(str(c['delta']))
    ans=[]
    for b in range(R):
        table=[]
        for k in range(K):
            phase=F(2*b+1,2*R) if k%2==0 else F(2*(R-b)-1,2*R)
            center=(k+phase)/K
            left,right=center-delta/2,center+delta/2
            table.append([(F(k,K),left,eps),
                          (left,right,eps+(1-eps)/(K*delta)),
                          (right,F(k+1,K),eps)])
        ans.append(table)
    audit.validate_tables(ans)
    return ans

def independent_values(x,c,selector,pi):
    """Dictionary polynomial expansion and independently formed overlap kernels."""
    R,K,n=c['R'],c['K'],c['n']
    N,J=len(x),len(c['lengths'])
    tables=exact_tables(c)
    selected=np.zeros((N,J));prior=np.zeros((N,J));fixed=[]
    for b in range(R):
        rows=[]
        for i in range(n):
            row=[]
            for k in range(K):
                values=np.zeros((N,J))
                for j,L in enumerate(c['lengths']):
                    for left,right,height in tables[b][k]:
                        for shift in range(-3,4):
                            values[:,j]+=float(height)*np.maximum(0,
                                np.minimum(x[:,i]+L/2,float(right)+shift)-
                                np.maximum(x[:,i]-L/2,float(left)+shift))/L
                row.append(values)
            rows.append(row)
        polynomial={(0,)*K:np.ones((N,J))}
        fixed_b=np.ones((N,J))
        for row in rows:
            fixed_b*=sum(row)
            nxt={}
            for alpha,coefficient in polynomial.items():
                for k,w in enumerate(row):
                    key=alpha[:k]+(alpha[k]+1,)+alpha[k+1:]
                    term=coefficient*w
                    if key in nxt:nxt[key]+=term
                    else:nxt[key]=term.copy()
            polynomial=nxt
        for a,v in polynomial.items():
            if selector[a]==b:selected+=v
        prior+=float(pi[b])*fixed_b
        fixed.append(fixed_b)
    return selected,prior,np.asarray(fixed)

def decimal_ci(k,N,alpha):
    with localcontext() as ctx:
        ctx.prec=80
        if isinstance(alpha,F):
            alpha=Decimal(alpha.numerator)/Decimal(alpha.denominator)
        p=Decimal(k)/N
        target=(1/alpha).ln()/N
        def kl(q):
            if p==0:return -(1-q).ln()
            if p==1:return -q.ln()
            return p*(p/q).ln()+(1-p)*((1-p)/(1-q)).ln()
        if k==0:lo=Decimal(0)
        else:
            a,b=Decimal(0),p
            for _ in range(260):
                m=(a+b)/2
                if kl(m)>target:a=m
                else:b=m
            lo=a
        if k==N:hi=Decimal(1)
        else:
            a,b=p,Decimal(1)
            for _ in range(260):
                m=(a+b)/2
                if kl(m)>target:b=m
                else:a=m
            hi=b
        return float(lo),float(hi)

def main():
    started=time.time()
    spec=importlib.util.spec_from_file_location('main_probe',SOURCE/'probe.py')
    p=importlib.util.module_from_spec(spec);spec.loader.exec_module(p)
    frozen=json.loads((SOURCE/'FROZEN_SELECTOR.json').read_text())
    result=json.loads((SOURCE/'RESULTS.json').read_text())
    plan=json.loads((SOURCE/'PLAN.json').read_text())
    c=result['config']
    assert plan['config']==c==p.CONFIG
    assert plan['script_SHA256']==sha(SOURCE/'probe.py')
    assert frozen['plan_SHA256']==result['plan_SHA256']==sha(SOURCE/'PLAN.json')
    assert result['selector_SHA256']==sha(SOURCE/'FROZEN_SELECTOR.json')
    selector=dict(zip(map(tuple,frozen['types']),frozen['selector']))
    assert set(selector)==set(audit.compositions(c['n'],c['K']))
    exact_pi=audit.prior_weights(c['n'],c['K'],c['R'],selector)
    assert list(map(float,exact_pi))==frozen['prior_mass']
    tables=exact_tables(c)
    full_weights=audit.cube_weights(tables,[F(0)]*c['n'],F(1,2))
    exact_mass=audit.selected_numerator(full_weights,selector)
    assert exact_mass==1

    # Replay selection using original trainer, without calling main or writing its files.
    x=np.random.default_rng(c['train_seed']).random((c['train_N'],c['n']))
    feat,pr,types=p.features(x,c)
    best,records=p.fit(feat,pr,types,c)
    assert best[1].tolist()==frozen['selector']
    assert best[2]==frozen['lambda_']
    assert best[3].tolist()==frozen['prior_mass']
    assert records==frozen['training']
    del x,feat,pr

    # True Fraction box integral vs original features, on 24 actual heldout queries.
    x=np.random.default_rng(c['test_seed']).random((24,c['n']))
    feat,pr,types=p.features(x,c)
    selected=p.selector_value(feat,np.array(frozen['selector']))
    max_kernel_error=0.;max_gated_error=0.;max_prior_error=0.;checks=0
    for row,xx in enumerate(x):
        xf=[F.from_float(float(v)) for v in xx]
        for j,length in enumerate(c['lengths']):
            radius=F(str(length))/2
            exact_w=audit.cube_weights(tables,xf,radius)
            gv=audit.selected_numerator(exact_w,selector)/(2*radius)**c['n']
            pv=audit.prior_numerator(exact_w,selector)/(2*radius)**c['n']
            max_gated_error=max(max_gated_error,abs(float(gv)-selected[row,j]))
            max_prior_error=max(max_prior_error,abs(float(pv)-float(np.dot(np.asarray(exact_pi,dtype=float),pr[:,row,j]))))
            main_w=p.kernels(xx[None,:],c)
            for b,i,k in product(range(c['R']),range(c['n']),range(c['K'])):
                max_kernel_error=max(max_kernel_error,abs(float(exact_w[b][i][k]/(2*radius))-main_w[b,i,0,j,k]))
            checks+=1
    assert max_kernel_error<1e-12 and max_gated_error<1e-10 and max_prior_error<1e-10
    del feat,pr,selected

    # All heldout statistics recomputed with no main features/kernel calls.
    rng=np.random.default_rng(c['test_seed']);counts=np.zeros(6,dtype=np.int64)
    sums=np.zeros(6);min_margins=np.full(6,np.inf)
    for start in range(0,c['test_N'],256):
        xx=rng.random((min(256,c['test_N']-start),c['n']))
        gv,pv,fixed=independent_values(xx,c,selector,exact_pi)
        maxima=np.concatenate([gv.max(axis=1)[None,:],pv.max(axis=1)[None,:],fixed.max(axis=2)])
        counts+=(maxima>frozen['lambda_']).sum(axis=1)
        sums+=maxima.sum(axis=1)
        min_margins=np.minimum(min_margins,np.min(abs(maxima-frozen['lambda_']),axis=1))
    print(json.dumps(dict(independent_counts=counts.tolist(),
        main_counts=[r['successes'] for r in result['rows']],
        min_threshold_distances=min_margins.tolist(),
        fraction_errors=[max_kernel_error,max_gated_error,max_prior_error]),indent=2),flush=True)
    means=sums/c['test_N']
    mean_errors=abs(means-np.array([r['mean_finite_max'] for r in result['rows']]))
    assert max(mean_errors)<1e-11
    ci_errors=[]
    for row in result['rows']:
        ci=decimal_ci(row['successes'],row['N'],Decimal(1)/Decimal(240))
        ci_errors.extend(abs(ci[i]-row['weak_CI'][i]) for i in (0,1))
    assert max(ci_errors)<2e-13
    out=dict(status='PASS' if counts.tolist()==[r['successes'] for r in result['rows']] else 'STRICT_BOUNDARY_FAILURE',source_script_SHA256=sha(SOURCE/'probe.py'),
        source_plan_SHA256=sha(SOURCE/'PLAN.json'),source_selector_SHA256=sha(SOURCE/'FROZEN_SELECTOR.json'),
        exact_mass=str(exact_mass),exact_prior_mass=[str(v) for v in exact_pi],
        training_replay='selector, chosen threshold, pi, all training records identical',
        fraction_query_count=24,fraction_query_scale_checks=checks,
        max_kernel_absolute_error=max_kernel_error,max_gated_absolute_error=max_gated_error,
        max_prior_absolute_error=max_prior_error,independent_holdout_successes=counts.tolist(),
        max_holdout_mean_error=float(max(mean_errors)),minimum_threshold_distances=min_margins.tolist(),
        max_CI_error_80digit_reference=max(ci_errors),
        main_CI_rounding_note='main lo=b and hi=a are inward by at most numerical root error; this is not an outward interval certificate',
        conclusion='Original strict-boundary counts are not reproducible under equivalent kernel evaluation; raw significance is not accepted as mathematical evidence. See BOUNDARY_DIAGNOSTIC.json.',
        numpy_version=np.__version__,seconds=time.time()-started)
    (HERE/'INDEPENDENT_SELECTOR_CHECK.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
