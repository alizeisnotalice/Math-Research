"""Exact finite proposal/RN/cap probes before the continuous-source proof."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json,hashlib
HERE=Path(__file__).resolve().parent
def main():
 target=HERE/'FINITE.json';assert not target.exists();rows=[];alpha=F(1,2)
 us=[F(3,5),F(1,5),F(4,5)];R=2;J=len(us);pi=[F(1,2)]*2
 for n in [1,2,3,4,6]:
  points=list(product([0,1],repeat=n));records=[]
  for x in points:
   A=[[F(1)for l in range(J)]for b in range(R)]
   for b,l in product(range(R),range(J)):
    u=us[l]*(-1 if b else 1)
    for xi in x:A[b][l]*=1+u*(1 if xi==0 else -1)
   q=alpha+(1-alpha)*sum(sum(v)for v in A)/(R*J)
   # Product of per-coordinate mixtures is a DIFFERENT proposal density.
   qwrong=F(1)
   for xi in x:qwrong*=alpha+(1-alpha)*sum(1+us[l]*(-1 if b else 1)*(1 if xi==0 else -1)for b,l in product(range(R),range(J)))/(R*J)
   G=[sum(pi[b]*A[b][l]for b in range(R))for l in range(J)]
   records.append((q,qwrong,max(G)))
  assert sum(q for q,qw,g in records)/2**n==1
  for lam in [F(11,10),F(3,2),F(2),F(5)]:
   true=lam*sum(g>lam for q,qw,g in records)/2**n
   actual=sum(q*(lam/q if g>lam else 0)for q,qw,g in records)/2**n
   wrong=sum(q*(lam/qw if g>lam else 0)for q,qw,g in records)/2**n
   assert actual==true
   # Uniform anchor weights and uniform input pi: sharper J/(1-alpha).
   cap=min(lam/alpha,F(J)/(1-alpha))
   contributions=[lam/q if g>lam else F(0)for q,qw,g in records]
   assert max(contributions)<=cap
   rows.append(dict(n=n,threshold=str(lam),true_weak=str(true),correct_RN_estimator_mean=str(actual),wrong_per_coordinate_RN_mean=str(wrong),cap_exact=str(cap),max_observed_contribution=str(max(contributions))))
 assert any(r['wrong_per_coordinate_RN_mean']!=r['true_weak']for r in rows)
 out=dict(script_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),rows=rows,
 scope='Exact abstract finite-probability RN/cap audit, not an HL source or C_n example. Product proposals can be valid if actual sampling and density agree; mismatch is the counterexample.')
 target.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(dict(cases=len(rows),mismatch_examples=sum(r['wrong_per_coordinate_RN_mean']!=r['true_weak']for r in rows))))
if __name__=='__main__':main()
