#!/usr/bin/env python3
"""Check output accounting and scalar reductions without recomputing any vertex quadratic."""
from pathlib import Path
import json,hashlib,math,numpy as np
p=Path(__file__).resolve().parent;sha=lambda x:hashlib.sha256(x.read_bytes()).hexdigest();m=json.loads((p/'manifest.json').read_text())
for f,h in m['files'].items():assert sha(p/f)==h,f
for line in (p/'prerun-sha256.txt').read_text().splitlines():
 h,f=line.split('  ',1);q=Path(f);assert sha(q if q.is_absolute() else p/q)==h,f
cache=Path('/Users/zhengzhihao/Documents/Codex/2026-10-27/n-n-5-1-1-verify/work/cache/gaussian-round001/original-8192-shapes.npz');assert sha(cache)=='6dfa48bf8954a7d02270af5bda9c1d33d4fcaf0ba4f5e486589186b7f7c09552'
old=json.loads((p/'r016-manifest.json').read_text());assert sha(p/'r016-result.json')==old['files']['result.json'];R=json.loads((p/'results.json').read_text());rc=json.loads((p/'receipt.json').read_text());z=np.load(p/'observables.npz');base=np.load(cache);v=z['L_vertex'];t=z['L_trace'];q=z['max_vertex_Q'];tr=z['trace_B'];d=base['det_A_minus_half'];k=R['kd'];assert rc['status']==R['status']=='completed' and len(v)==len(t)==len(q)==len(tr)==8192
assert R['vertex_quadratic_evaluations']==8192*128 and R['trace_evaluations']==8192 and len(R['same_input_cross_checks'])==8 and R['cross_check_vertex_evaluations']==1024
assert np.isfinite(v).all() and np.isfinite(t).all() and (v>0).all() and (t>=v*(1-1e-12)).all() and np.allclose(v,k*d/q**3,rtol=2e-15,atol=0) and np.allclose(t,k*d/tr**3,rtol=2e-15,atol=0)
assert max(t)<=R['uniform_bound']*(1+1e-10) and all(x['max_vertex_relative_discrepancy']<=1e-6 and x['trace_relative_discrepancy']<=1e-6 and x['detweight_relative_discrepancy']<=1e-6 for x in R['same_input_cross_checks'])
for name,x in [('vertex',v),('trace',t)]:
 s=R[name];assert math.isclose(float(x.mean()),s['mean'],rel_tol=1e-15) and np.allclose(x.reshape(8,1024).mean(axis=1),s['block_means'],rtol=1e-15,atol=0) and math.isclose(float(x.std(ddof=1)/np.sqrt(8192)),s['nominal_shape_SE'],rel_tol=1e-15)
 assert math.isclose(s['mean']/R['g_required_reference'],s['ratio_to_required_g'],rel_tol=1e-15)
assert [x['block'] for x in R['blocks']]==list(range(8)) and all(x['status']=='completed' and x['shapes']==1024 for x in R['blocks'])
print('PASS: frozen cache/source hashes, 8192 retained new observables, eight block summaries, same-input cross-check ledger and arithmetic hierarchy. No vertex recomputation; no expectation or global g certificate.')
