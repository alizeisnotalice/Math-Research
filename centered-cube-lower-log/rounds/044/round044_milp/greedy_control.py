"""Same-data round043 coordinate-search control. Frozen before evaluating heldout."""
from pathlib import Path
import sys
sys.dont_write_bytecode=True
import importlib.util,json,numpy as np,math
H=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('p44',H/'probe.py');p=importlib.util.module_from_spec(spec);spec.loader.exec_module(p);C=p.C.copy()
p.dump('GREEDY_PLAN.json',dict(config=C,script_SHA256=p.sha(__file__),scope='Same real source, formal train seed/queries/thresholds, round043 coordinate optimizer. Its frozen selector evaluated on same independent test seed. Added before reading MILP heldout; does not alter MILP selection.'))
feat,fixed,_=p.source.features(np.random.default_rng(C['train_seed']).random((C['train_N'],32)),C);best,records=p.source.fit(feat,fixed,C,32);score,s,lam,pi=best
p.dump('GREEDY_FROZEN.json',dict(lambda_=lam,selector=s.tolist(),prior_mass=pi.tolist(),train_weak=score,training=records,plan_SHA256=p.sha(H/'GREEDY_PLAN.json')))
rng=np.random.default_rng(C['test_seed']);counts=np.zeros(6,int)
for at in range(0,C['test_N'],C['batch']):
 ft,fx,_=p.source.features(rng.random((min(C['batch'],C['test_N']-at),32)),C);vals=np.concatenate([p.source.value(ft,s).max(-1)[None,:],np.einsum('b,bnj->nj',pi,fx).max(-1)[None,:],fx.max(-1)]);counts+=(vals>lam+C['guard']).sum(-1)
p.dump('GREEDY_SUMMARY.json',dict(lambda_=lam,train_weak=float(score),test_counts=counts.tolist(),test_weak=(lam*counts/C['test_N']).tolist(),model_order=['greedy_gate','actual_prior','fixed0','fixed1','fixed2','fixed3'],prior_mass=pi.tolist(),scope='Same train and test sample control; float exploration only. No extra inference from unadjusted model selection.'))
print(json.dumps(dict(lambda_=lam,train_weak=float(score),test_weak=(lam*counts/C['test_N']).tolist())))
