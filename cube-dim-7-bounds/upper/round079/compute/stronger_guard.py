from pathlib import Path
from fractions import Fraction as F
import json,os,hashlib
p=Path(__file__).resolve().parent
r=json.loads((p/'stronger_registration.json').read_text()); previous=json.loads((p/'certificate_result.json').read_text()); k=previous['key_rational_bounds']
gamma=F(r['gamma']); price=F(r['price']); c=F(r['c']); slack=F(r['slack'])
budget=F(k['gamma_budget']); price_enclosure=F(k['price_upper']); remainder=F(k['full_root_remainder'])
gates={'same_tuple_certificate_pass':previous['status']=='PASS','gamma_below_budget':gamma<budget,'exact_price_upper':price_enclosure<price,'scalar_slack':gamma-price-c==slack,'full_scalar_remainder':remainder<slack}
assert all(gates.values())
out={'status':'PASS exact Fraction arithmetic on locked previous certificate; conditional candidate only','pid':os.getpid(),'gates':gates,'budget_margin':str(budget-gamma),'budget_margin_decimal':float(budget-gamma),'price_margin_decimal':float(price-price_enclosure),'root_remainder_decimal':float(remainder),'slack':str(slack),'registration_sha256':hashlib.sha256((p/'stronger_registration.json').read_bytes()).hexdigest(),'certificate_sha256':hashlib.sha256((p/'certificate_result.json').read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(p/'stronger_guard.json').write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps(out,indent=2))
