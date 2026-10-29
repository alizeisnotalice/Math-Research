"""Exact mass and marginal checks can pass a WRONG RN for a whole-vector mixture."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json,hashlib
HERE=Path(__file__).resolve().parent

def main():
 p=HERE/'MASS_TRAP.json';assert not p.exists();rows=[];us=[F(3,5),F(1,5),F(4,5)];lam=F(11,10)
 mass=F(0);actual=F(0);wrong=F(0);wrong_rn_mass=F(0);marg=[F(0),F(0)];correct_anchor=F(0);wrong_anchor=F(0)
 for x in product([0,1],repeat=2):
  A=[]
  for b in [0,1]:
   A.append([((1+u*(-1 if b else 1)*(1 if x[0]==0 else -1))*(1+u*(-1 if b else 1)*(1 if x[1]==0 else -1)))for u in us])
  q=F(1,2)+sum(sum(v)for v in A)/12
  target=max((A[0][j]+A[1][j])/2 for j in range(3));event=target>lam
  mass+=q/4;wrong_rn_mass+=q/4
  actual+=F(event)*lam/4;wrong+=F(event)*q*lam/4
  for i in range(2):marg[i]+=q*F(x[i]==0)/4
  correct_anchor+=A[0][0]/4;wrong_anchor+=q*A[0][0]/4
  rows.append(dict(x=x,q_exact=str(q),wrong_per_coordinate_q_exact='1',model_max_exact=str(target),hit=event))
 assert mass==wrong_rn_mass==1 and marg==[F(1,2)]*2 and correct_anchor==1 and wrong!=actual
 out=dict(rows=rows,threshold_exact=str(lam),proposal_mass_exact=str(mass),wrong_RN_mass_exact=str(wrong_rn_mass),one_coordinate_marginals=list(map(str,marg)),true_weak_exact=str(actual),wrong_RN_weak_exact=str(wrong),correct_RN_anchor_integral_exact=str(correct_anchor),wrong_RN_anchor_integral_exact=str(wrong_anchor),
 scope='Abstract finite-probability diagnostic counterexample, not an HL lower bound.',script_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
 p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
