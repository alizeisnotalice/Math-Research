"""Independent frozen-selector rerun and exact cap/normalization checks."""
from pathlib import Path
import sys
sys.dont_write_bytecode=True
import importlib.util,json,math
from fractions import Fraction as F
from itertools import product
import numpy as np
H=Path(__file__).resolve().parent
sp=importlib.util.spec_from_file_location('milp_probe',H/'probe.py');p=importlib.util.module_from_spec(sp);sp.loader.exec_module(p)
C=p.C
assert p.sha(H/'probe.py')==json.loads((H/'PLAN.json').read_text())['script_SHA256']
assert p.sha(p.SOURCE)==json.loads((H/'PLAN.json').read_text())['source_SHA256']
checked=0
for ts in product([F(0),F(1,10),F(1,2),F(1),F(3,2),F(10)],repeat=3):
 for bits in product([0,1],repeat=3):
  assert (sum(t*b for t,b in zip(ts,bits))>=1)==(sum(min(F(1),t)*b for t,b in zip(ts,bits))>=1)
  checked+=1
rng=np.random.default_rng(C['train_seed']);feat,fixed,e=p.source.features(rng.random((C['train_N'],32)),C);training=json.loads((H/'TRAINING.json').read_text());rechecks=[]
for r in training:
 if not r['incumbent_available']:continue
 s=np.array(r['selector']);lam=r['lambda_'];assert len(s)==33 and all(0<=v<4 for v in s)
 vals=p.source.value(feat,s);v=vals.max(-1);k=int((v>lam+C['guard']).sum());assert k==r['recomputed_train_successes']
 nums=[sum(math.comb(32,j)for j in range(33)if s[j]==b)for b in range(4)];assert sum(nums)==2**32;pi=np.array(nums)/2**32
 rechecks.append(dict(lambda_=lam,recomputed_train_successes=k,gate_weak=lam*k/C['train_N'],actual_prior_train_weak=float(lam*np.mean(np.einsum('b,bnj->nj',pi,fixed).max(-1)>lam+C['guard'])),bestfixed_train_weak=max(float(lam*np.mean(q.max(-1)>lam+C['guard']))for q in fixed),margin_successes=int((v>=lam+C['milp_margin']).sum()),guard_success_witness_L=[None if v[i]<=lam+C['guard']else C['lengths'][int(vals[i].argmax())]for i in range(len(v))]))
f=json.loads((H/'FROZEN.json').read_text());assert sum(f['prior_exact_numerators'])==f['prior_exact_denominator']==2**32
# Independent heldout streaming replay using exactly the frozen seed and rounded selector.
s=np.array(f['selector']);lam=f['lambda_'];pi=np.array(f['prior_exact_numerators'])/2**32;rng=np.random.default_rng(C['test_seed']);counts=np.zeros(6,int)
for at in range(0,C['test_N'],C['batch']):
 x=rng.random((min(C['batch'],C['test_N']-at),32));ft,fx,_=p.source.features(x,C);v=np.concatenate([p.source.value(ft,s).max(-1)[None,:],np.einsum('b,bnj->nj',pi,fx).max(-1)[None,:],fx.max(-1)]);counts+=(v>lam+C['guard']).sum(-1)
result=json.loads((H/'SUMMARY.json').read_text());assert counts.tolist()==[r['k']for r in result['rows']]
p.dump('VERIFICATION.json',dict(passed=True,cap_exact_binary_tests=checked,all_incumbents_recomputed_train=rechecks,independent_heldout_replay_counts=counts.tolist(),integer_prior_mass_verified=True,scope='Numerical implementation recheck only; no interval error or continuous query cell proof.'))
print(json.dumps(dict(passed=True,cap_tests=checked,test_counts=counts.tolist())))
