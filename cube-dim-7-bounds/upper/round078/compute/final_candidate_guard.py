from fractions import Fraction as F
from math import factorial
from pathlib import Path
import hashlib,json,os
p=Path(__file__).resolve().parent
r=json.loads((p/'final_candidate_registration.json').read_text())
L=F(8192); sigma=F(25,256); N=L**64/factorial(64)
def S(y):
 term=total=F(1)
 for j in range(1,65):
  term*=y/j; total+=term
 return total
remainder=1/N+1/S(sigma*L)+182/S(sigma*L/2)
eB=eN=eT=F(1,2**19); eJ=eR=eS=eD=F(1,2**20); mass=F(3043,10000); ret=F(999999,1000000)
P0=1-eB-F(6,768); lo=P0*(1-eN)/(1+eJ); hi=(1+eN)/(P0*(1-eJ)); M=2*(1+eR)*hi
budget=(1-eD)*mass*(2*ret*(1-eR)*lo-eS-eT)/M
gamma=F(r['gamma']); price=F(r['price_upper']); c=F(r['c']); allowance=F(r['scalar_allowance'])
checks={'gamma_below_complete_budget':gamma<budget,'slack_exact':gamma-price-c==allowance,'scalar_remainder_below_new_slack':remainder<allowance}
assert all(checks.values())
out={'status':'PASS exact Fraction checks, conditional candidate only','pid':os.getpid(),'checks':checks,'registration_sha256':hashlib.sha256((p/'final_candidate_registration.json').read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'complete_budget_exact':str(budget),'gamma_margin_exact':str(budget-gamma),'remainder_decimal':float(remainder),'remainder_upper_simple':'1/10000','remainder_exact':str(remainder),'slack_exact':str(allowance),'tail':'Use R073 analytically established monotonicity of N and scalar prefix bounds; original geometry unchanged.'}
(p/'final_candidate_guard.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:out[k] for k in ['status','pid','checks','remainder_decimal','slack_exact']},indent=2))
