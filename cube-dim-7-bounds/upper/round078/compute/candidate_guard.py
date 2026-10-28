from fractions import Fraction as F
from pathlib import Path
import hashlib,json,os
p=Path(__file__).resolve().parent
r=json.loads((p/'candidate_registration.json').read_text())
eB=eN=eT=F(1,2**19); eJ=eR=eS=eD=F(1,2**20); mass=F(3043,10000); ret=F(999999,1000000)
P0=1-eB-F(6,768); lo=P0*(1-eN)/(1+eJ); hi=(1+eN)/(P0*(1-eJ)); M=2*(1+eR)*hi
Gamma=mass*(2*ret*(1-eR)*lo-eS-eT)/M; budget=(1-eD)*Gamma
gamma=F(r['gamma']); c=F(r['c']); price=F(r['price_upper'])
assert gamma<budget and gamma-price-c==F(1,5000)
out={'pid':os.getpid(),'status':'PASS exact Fraction guards, conditional on registered Bernstein assumptions','registration_sha256':hashlib.sha256((p/'candidate_registration.json').read_bytes()).hexdigest(),'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'budget_exact':str(budget),'gamma_margin_exact':str(budget-gamma),'gamma_minus_price_minus_c':str(gamma-price-c),'budget_decimal':float(budget),'gamma_margin_decimal':float(budget-gamma),'scope':'The previously proved scalar remainder < 1/5000 is not rerun; no new Cn theorem asserted.'}
(p/'candidate_guard.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
