"""The sole predeclared coefficient scaling candidate; exact rationals only."""
from pathlib import Path
from fractions import Fraction as Q
import json,hashlib
root=Path(__file__).resolve().parents[1]
src=root/'upstream-v4/cert/cert_h32_a4740_R3320_3p59123.json'
d=json.loads(src.read_text());fac=Q(100000,100009)
report=json.loads((root/'second_public.json').read_text())
assert report['c2']['min']=='[1.000099345 +/- 3.68e-10]'
m=Q('1.000099345')-Q('3.68e-10')
assert m*fac>1
out={key:d[key] for key in ['alpha','R','h','N']}
out['c']=str(Q(d['c'])*fac);out['a']={key:str(Q(val)*fac) for key,val in d['a'].items()}
for key,val in d['a'].items():assert Q(out['a'][key])==fac*Q(val)
base=json.loads((root/'exact_cost.json').read_text());L=Q(base['lower']);U=Q(base['upper'])
assert fac*U<L and fac*U<Q(359091,100000)
p=Path(__file__).with_name('scaled_certificate.json')
if p.exists():raise SystemExit('Candidate already exists, inspect instead of regenerate')
p.write_text(json.dumps(out,indent=2)+'\n')
proof={'scale':str(fac),'public_certificate_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'new_certificate_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'public_dominance_lower':str(m),'scaled_dominance_lower':str(m*fac),'new_cost_lower':str(fac*L),'new_cost_upper':str(fac*U),'strict_gain_lower':str(L-fac*U),'target_rational':'359091/100000','status':'Exact algebraic scaling; full new-certificate checks pending.'}
p.with_name('scaling_arithmetic.json').write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof,indent=2))
