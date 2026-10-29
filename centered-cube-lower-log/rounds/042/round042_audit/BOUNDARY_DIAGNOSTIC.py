#!/usr/bin/env python3
"""Three frozen cases, independent true-rational source kernels, no retraining.

The 1e-10 guard is a diagnostic tolerance, NOT a certified FP error bound.
All original experiment files are read only.
"""
import json
import time
from pathlib import Path
from fractions import Fraction as F
import numpy as np
import audit
from INDEPENDENT_SELECTOR_CHECK import independent_values, exact_tables, sha, decimal_ci
from decimal import Decimal

HERE=Path(__file__).resolve().parent
SRC=HERE.parent/'round042_selector'

def main():
    start=time.time();out=[]
    for name,folder in [('n8',SRC),('n7',SRC/'balanced/n7'),('n15',SRC/'balanced/n15')]:
        frozen=json.loads((folder/'FROZEN_SELECTOR.json').read_text())
        saved=json.loads((folder/'RESULTS.json').read_text())
        plan=json.loads((folder/'PLAN.json').read_text())
        c=plan['config'];lam=frozen['lambda_'];tol=1e-10
        selector=dict(zip(map(tuple,frozen['types']),frozen['selector']))
        pi=audit.prior_weights(c['n'],c['K'],c['R'],selector)
        assert list(map(float,pi))==frozen['prior_mass']
        assert frozen['plan_SHA256']==saved['plan_SHA256']==sha(folder/'PLAN.json')
        assert saved['selector_SHA256']==sha(folder/'FROZEN_SELECTOR.json')
        orbits=None
        if c['n']%2:
            orbits=set()
            for a,b in selector.items():
                orbit=[a[r:]+a[:r] for r in range(4)]
                assert len(set(orbit))==4
                for r,rot in enumerate(orbit):assert selector[rot]==(b+r)%4
                orbits.add(min(orbit))
            assert pi==[F(1,4)]*4
        # Exact fixed/prior L=3/2 constant check on 20 real heldout queries.
        tables=exact_tables(c)
        x=np.random.default_rng(c['test_seed']).random((20,c['n']))
        exact_kernel_checks=0
        for xx in x:
            xf=[F.from_float(float(v)) for v in xx]
            w=audit.cube_weights(tables,xf,F(3,4))
            for b in range(4):
                for row in w[b]:
                    assert sum(row)==F(3,2)
                    exact_kernel_checks+=1
            assert audit.prior_numerator(w,selector)/F(3,2)**c['n']==1
        # Float diagnostics on every original heldout query, with independent code.
        rng=np.random.default_rng(c['test_seed'])
        raw=np.zeros(6,dtype=np.int64);corrected=np.zeros(6,dtype=np.int64)
        definite=np.zeros(6,dtype=np.int64);possible=np.zeros(6,dtype=np.int64)
        max_l15_residual=np.zeros(5)
        min_margin=np.full(6,np.inf)
        j=c['lengths'].index(1.5)
        for offset in range(0,c['test_N'],256):
            xx=rng.random((min(256,c['test_N']-offset),c['n']))
            gv,pv,fixed=independent_values(xx,c,selector,pi)
            mx=np.concatenate([gv.max(axis=1)[None,:],pv.max(axis=1)[None,:],fixed.max(axis=2)])
            raw+=(mx>lam).sum(axis=1)
            max_l15_residual=np.maximum(max_l15_residual,
                np.max(abs(np.concatenate([pv[:,j][None,:],fixed[:,:,j]])-1),axis=1))
            pv[:,j]=1.;fixed[:,:,j]=1.
            mx=np.concatenate([gv.max(axis=1)[None,:],pv.max(axis=1)[None,:],fixed.max(axis=2)])
            corrected+=(mx>lam).sum(axis=1)
            definite+=(mx>lam+tol).sum(axis=1)
            possible+=(mx>lam-tol).sum(axis=1)
            min_margin=np.minimum(min_margin,np.min(abs(mx-lam),axis=1))
        alpha=Decimal(1)/Decimal(240 if name=='n8' else 480)
        cis=[dict(lower=decimal_ci(int(lo),c['test_N'],alpha)[0],
                  upper=decimal_ci(int(hi),c['test_N'],alpha)[1])
             for lo,hi in zip(definite,possible)]
        row=dict(case=name,n=c['n'],N=c['test_N'],lambda_=lam,tolerance=tol,
            selector_SHA256=sha(folder/'FROZEN_SELECTOR.json'),prior_mass_exact=[str(v) for v in pi],
            orbit_count=None if orbits is None else len(orbits),
            exact_L15_coordinate_table_checks=exact_kernel_checks,
            original_saved_counts=[r['successes'] for r in saved['rows']],
            independent_raw_counts=raw.tolist(),
            independent_L15_fixed_prior_exact1_counts=corrected.tolist(),
            diagnostic_definite_counts=definite.tolist(),diagnostic_possible_counts=possible.tolist(),
            diagnostic_ambiguous_counts=(possible-definite).tolist(),
            diagnostic_CIs=cis,
            diagnostic_gated_lower_minus_prior_upper=cis[0]['lower']-cis[1]['upper'],
            L15_max_float_residuals_prior_fixed=max_l15_residual.tolist(),
            minimum_threshold_distances_after_exact1_correction=min_margin.tolist(),
            status='boundary diagnostic only; guard has no certified global FP guarantee')
        out.append(row)
        print(json.dumps(row,indent=2),flush=True)
    result=dict(status='STRICT_BOUNDARY_FAILURE_CONFIRMED',models=['gated','prior','fixed0','fixed1','fixed2','fixed3'],
        scope='Read-only three-case replay. No retraining. Guard tolerance is not a numerical certificate.',
        cases=out,seconds=time.time()-start)
    (HERE/'BOUNDARY_DIAGNOSTIC.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
