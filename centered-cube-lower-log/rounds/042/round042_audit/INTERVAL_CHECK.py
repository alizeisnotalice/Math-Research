#!/usr/bin/env python3
"""Independent Fraction checks of the frozen outward interval repair; no MC rerun."""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal
import importlib.util
import json
import sys
import time
import numpy as np
import audit
from INDEPENDENT_SELECTOR_CHECK import exact_tables, decimal_ci, sha

HERE=Path(__file__).resolve().parent
SRC=HERE.parent/'round042_selector'
sys.path.insert(0,str(SRC))
spec=importlib.util.spec_from_file_location('interval_source',SRC/'interval_recheck.py')
source=importlib.util.module_from_spec(spec);spec.loader.exec_module(source)

def inside(value,lo,hi):
    assert F.from_float(float(lo))<=value<=F.from_float(float(hi)),(value,lo,hi)

def main():
    started=time.time()
    results=json.loads((SRC/'INTERVAL_RESULTS.json').read_text())
    plan=json.loads((SRC/'INTERVAL_RESULTS_PLAN.json').read_text())
    assert results['plan_SHA256']==sha(SRC/'INTERVAL_RESULTS_PLAN.json')
    assert plan['script_SHA256']==sha(SRC/'interval_recheck.py')
    assert plan['core_SHA256']==sha(SRC/'probe.py')
    diagnostic=json.loads((HERE/'BOUNDARY_DIAGNOSTIC.json').read_text())
    records=[];allci=[]
    for case in results['cases']:
        n=case['n'];folder=SRC if n==8 else SRC/f'balanced/n{n}'
        cfg=json.loads((folder/'PLAN.json').read_text())['config']
        frozen=json.loads((folder/'FROZEN_SELECTOR.json').read_text())
        assert case['plan_SHA256']==sha(folder/'PLAN.json')
        assert case['selector_SHA256']==sha(folder/'FROZEN_SELECTOR.json')
        s=np.array(frozen['selector'])
        selector=dict(zip(map(tuple,frozen['types']),frozen['selector']))
        tables=exact_tables(cfg)
        pi=audit.prior_weights(n,4,4,selector)
        points=3 if n<10 else 2
        x=np.random.default_rng(cfg['test_seed']).random((points,n))
        kl,kh=source.kernels(x,cfg)
        low,high=source.bounds(x,cfg,s)
        kernel_checks=0;output_checks=0
        maxwidth=0.
        for qi,xx in enumerate(x):
            xf=[F.from_float(float(v)) for v in xx]
            for j,length in enumerate(cfg['lengths']):
                L=F(str(length));w=audit.cube_weights(tables,xf,L/2)
                for b in range(4):
                    for i in range(n):
                        for k in range(4):
                            inside(w[b][i][k]/L,kl[b,i,qi,j,k],kh[b,i,qi,j,k])
                            kernel_checks+=1
                gated=audit.selected_numerator(w,selector)/L**n
                prior=audit.prior_numerator(w,selector)/L**n
                fixed=[]
                for b in range(4):
                    f=F(1)
                    for row in w[b]:f*=sum(row)/L
                    fixed.append(f)
                exact=[gated,prior]+fixed
                for model,v in enumerate(exact):
                    inside(v,low[model,qi,j],high[model,qi,j])
                    output_checks+=1
                if L==F(3,2):
                    assert exact[1:]==[F(1)]*5
                    assert np.all(low[1:,qi,j]==1) and np.all(high[1:,qi,j]==1)
                maxwidth=max(maxwidth,float(np.max(high[:,qi,j]-low[:,qi,j])))
        # Decimal reference targets exactly alpha=1/720, not binary64 .05/36.
        for row in case['rows']:
            a=decimal_ci(row['lower_successes'],row['N'],F(1,720))[0]
            b=decimal_ci(row['upper_successes'],row['N'],F(1,720))[1]
            allci.extend([abs(a-row['weak_CI'][0]),abs(b-row['weak_CI'][1])])
        old=next(v for v in diagnostic['cases'] if v['n']==n)
        assert [r['lower_successes'] for r in case['rows']]==old['diagnostic_definite_counts']
        assert all(r['lower_successes']==r['upper_successes'] for r in case['rows'])
        records.append(dict(n=n,real_queries=points,kernel_fraction_enclosures=kernel_checks,
            exact_selected_prior_fixed_enclosures=output_checks,max_checked_width=maxwidth,
            exact_prior=[str(v) for v in pi],
            all_saved_counts_match_independent_boundary_replay=True))
    assert max(allci)<2e-13
    out=dict(status='PASS',source_interval_script_SHA256=sha(SRC/'interval_recheck.py'),
        source_core_SHA256=sha(SRC/'probe.py'),interval_results_SHA256=sha(SRC/'INTERVAL_RESULTS.json'),
        cases=records,CI_endpoints_checked=len(allci),
        max_CI_absolute_error_vs_80digit_Decimal=max(allci),
        CI_scope='Independent high precision numerical check, not outward libm or interval-log certification',
        scope='Outward local kernel/DP proof reviewed; eight real dyadic queries checked exactly; saved all-sample counts agree with earlier independent replay. No full MC rerun.',
        seconds=time.time()-started,numpy_version=np.__version__)
    (HERE/'INTERVAL_CHECK.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
