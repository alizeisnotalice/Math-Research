#!/usr/bin/env python3
"""Small exact query-cell tests and independent integer CI certificate replay."""
from pathlib import Path
from fractions import Fraction as F
import importlib.util
import json
import sys
import time
import numpy as np
import audit
from INDEPENDENT_SELECTOR_CHECK import exact_tables, sha

HERE=Path(__file__).resolve().parent
SRC=HERE.parent/'round042_selector'
sys.path.insert(0,str(SRC))
spec=importlib.util.spec_from_file_location('cell_source',SRC/'cell_recheck.py')
cell=importlib.util.module_from_spec(spec);spec.loader.exec_module(cell)

def contains(v,a,b):assert F.from_float(float(a))<=v<=F.from_float(float(b))

def main():
    started=time.time();ir=cell.ir
    results=json.loads((SRC/'CELL_RESULTS.json').read_text())
    plan=json.loads((SRC/'CELL_PLAN.json').read_text())
    exactci=json.loads((SRC/'EXACT_CI.json').read_text())
    assert plan['script_SHA256']==sha(SRC/'cell_recheck.py')
    assert plan['parent_SHA256']==sha(SRC/'interval_recheck.py')
    assert plan['core_SHA256']==sha(SRC/'probe.py')
    assert results['plan_SHA256']==sha(SRC/'CELL_PLAN.json')
    assert exactci['source_SHA256']==sha(SRC/'CELL_RESULTS.json')
    assert exactci['script_SHA256']==sha(SRC/'exact_ci.py')
    old_overlap=ir.overlap;ir.overlap=cell.cell_overlap
    step=F(1,2**53);checks=0;dpchecks=0;records=[]
    for case in results['cases']:
        n=case['n'];folder=SRC if n==8 else SRC/f'balanced/n{n}'
        cfg=json.loads((folder/'PLAN.json').read_text())['config']
        frozen=json.loads((folder/'FROZEN_SELECTOR.json').read_text())
        selector=dict(zip(map(tuple,frozen['types']),frozen['selector']))
        xx=np.random.default_rng(cfg['test_seed']).random((1,n))
        kl,kh=ir.kernels(xx,cfg)
        low,high=ir.bounds(xx,cfg,np.array(frozen['selector']))
        tables=exact_tables(cfg)
        # Opposite corners plus alternating corner, all rational and exactly in cell.
        for corner in (0,1,2):
            point=[F.from_float(float(v))+(step if corner==1 or (corner==2 and i%2) else 0)
                   for i,v in enumerate(xx[0])]
            for j,length in enumerate(cfg['lengths']):
                L=F(str(length));w=audit.cube_weights(tables,point,L/2)
                for b in range(4):
                    for i in range(n):
                        for k in range(4):
                            contains(w[b][i][k]/L,kl[b,i,0,j,k],kh[b,i,0,j,k]);checks+=1
                gated=audit.selected_numerator(w,selector)/L**n
                prior=audit.prior_numerator(w,selector)/L**n
                fixed=[]
                for b in range(4):
                    f=F(1)
                    for row in w[b]:f*=sum(row)/L
                    fixed.append(f)
                for model,v in enumerate([gated,prior]+fixed):
                    contains(v,low[model,0,j],high[model,0,j]);dpchecks+=1
        counts=[(r['lower_successes'],r['upper_successes']) for r in case['rows']]
        assert all(a==b for a,b in counts)
        records.append(dict(n=n,one_real_cell_three_corners=True,counts=[a for a,b in counts]))
    ir.overlap=old_overlap
    # Independently assemble the exact likelihood expression; no source certified().
    comparisons=0;maxbits=0
    for case in exactci['cases']:
        source_case=next(c for c in results['cases'] if c['n']==case['n'])
        for row,srow in zip(case['rows'],source_case['rows']):
            N=row['N'];assert row['lower_count']==srow['lower_successes']
            assert row['upper_count']==srow['upper_successes']
            for side,k in [('lower',row['lower_count']),('upper',row['upper_count'])]:
                endpoint=row[side];a=endpoint['numerator'];d=endpoint['denominator']
                assert 0<k<N and 0<a<d
                assert a*N<=k*d if side=='lower' else a*N>=k*d
                # Compare alpha-denominator * likelihood numerator to numerator
                # of empirical likelihood after cancelling (N*d)**N denominators.
                lhs=720*(a**k)*((d-a)**(N-k))*(N**N)
                rhs=(d**N)*(k**k)*((N-k)**(N-k))
                assert lhs<=rhs
                maxbits=max(maxbits,lhs.bit_length(),rhs.bit_length());comparisons+=1
    assert comparisons==36
    out=dict(status='PASS',cell_script_SHA256=sha(SRC/'cell_recheck.py'),
        exact_ci_script_SHA256=sha(SRC/'exact_ci.py'),cell_results_SHA256=sha(SRC/'CELL_RESULTS.json'),
        exact_ci_results_SHA256=sha(SRC/'EXACT_CI.json'),cases=records,
        exact_kernel_corner_checks=checks,exact_full_average_corner_checks=dpchecks,
        independent_integer_comparisons=comparisons,largest_integer_bits=maxbits,
        scope='Whole-cell algorithm reviewed; selected exact corners tested; every rational CI departure inequality independently certified. Assumes ideal iid uniform53bit cell indices, not a PRNG proof or deterministic continuum integral.',
        seconds=time.time()-started)
    (HERE/'CELL_CI_CHECK.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
