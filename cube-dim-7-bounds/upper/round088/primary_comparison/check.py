from fractions import Fraction as Q
from pathlib import Path
import json,hashlib,os

def logb(z):
 assert z>1
 u=(z-1)/(z+1)
 a=2*sum((u**(2*j+1)/Q(2*j+1) for j in range(24)),Q(0))
 return a,a+2*u**49/(49*(1-u*u))
old_hi=Q(619,400)-logb(Q(100000,37073))[0]
cap=Q(200000,37073*65536)
mu_hi=4-6*logb(Q(5,3))[0]
new_lo=logb(5/(1+mu_hi))[0]
checks={'old_complete_coefficient':old_hi<Q(556,1000),'old_cap_error':cap<Q(1,10000),'full_cost_gap':Q(9,10)-old_hi-cap>Q(3439,10000),'new_family_s5':new_lo>Q(9493,10000)}
assert all(checks.values())
out={'status':'PASS','pid':os.getpid(),'checks':checks,'old_coefficient_upper':str(old_hi),'cap_error_upper':str(cap),'s5_free_coefficient_lower':str(new_lo),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'scalar comparisons conditional on R088 full proof; no threshold'}
Path(__file__).with_name('receipt.json').write_text(json.dumps(out,indent=2))
print('PASS: 4 fresh exact complete-cost comparisons')
