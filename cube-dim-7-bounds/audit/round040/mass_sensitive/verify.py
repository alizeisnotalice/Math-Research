"""Minimal R040 scalar/payment checks; not a replay of the author's 1D toy.

No production history, geometry atomization, or LP solver is run.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,sys

HERE=Path(__file__).resolve().parent
eta=F(1,7)
w=[F(3,10),F(5,12),F(7,20)]
budget_ratios=[F(1),1+eta/3,1+eta]
weights=[v/h for v,h in zip(w,budget_ratios)]
assert all(v/(1+eta)<=p<=v for v,p in zip(w,weights))
assert sum(w)/(1+eta)<=sum(weights)<=sum(w)
assert sum(v/(1+eta) for v in w)==sum(w)/(1+eta)
# Exact disjoint planting identity, with homothety volume multiplier delta.
FB,FH,RB,RH,delta=F(9),F(5),F(7,8),F(2,3),F(1,1000)
mixture=(FB*RB+delta*FH*RH)/(FB+delta*FH)
assert mixture-RB==delta*FH*(RH-RB)/(FB+delta*FH)

# A mass-consistent weighted-cluster LP illustration. Nested open stars
# (-1,1) inside (-2,2), and a disjoint star (10,11); f=2 on their union.
# Thus actual costs are M=(4,8,2), not independently chosen prices.
M=[F(4),F(8),F(2)]
columns=[({0},0),({1},1),({0,1},1),({2},2)]
theta=[F(0),F(0),F(1),F(1)]
s=[F(0),F(8),F(2)]
assert all(sum(t for t,(J,j) in zip(theta,columns) if i in J)>=1 for i in range(3))
assert all(sum(s[i] for i in J)<=M[j] for J,j in columns)
primal=sum(t*M[j] for t,(J,j) in zip(theta,columns))
assert primal==sum(s)==10
# The two separated components have costs 8 and 2; no cross cluster exists.
assert primal==F(8)+F(2)

receipt={'status':'passed minimal new scalar and mass-consistent cost checks',
         'source_proof_sha256':hashlib.sha256(Path(sys.argv[1]).read_bytes()).hexdigest(),
         'budget_bracket':{'sum_w':str(sum(w)),'endpoint':str(sum(w)/(1+eta)),'hybrid':str(sum(weights))},
         'planting_ratio':str(mixture),'weighted_cluster_cost':str(primal),
         'weighted_cluster_primal_and_dual_equal':True,
         'interpretation':'Primal feasibility gives an upper certificate; dual feasibility alone gives only a cost lower bound. No uniform estimate follows.',
         'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'scope':'No R037 history computation, author toy replay, production LP/search, MATLAB, n7, or Git mutation.'}
(HERE/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
