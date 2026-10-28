from fractions import Fraction as Q
from pathlib import Path
import json,os,hashlib

def lb(z):
 assert z>1
 u=(z-1)/(z+1);v=2*sum((u**(2*j+1)/Q(2*j+1) for j in range(24)),Q(0))
 return v,v+2*u**49/(49*(1-u*u))
l2=lb(Q(2));linf=lb(1/(1-l2[0]))[0]
l53=lb(Q(5,3));c5hi=lb(5/(5-6*l53[1]))[1]
c88hi=Q(3,2)-lb(Q(20,11))[0]
checks={'c_infty_gt_1_18':linf>Q(59,50),'old_family_c5_lt_point95':c5hi<Q(19,20),'fixed_old_c88_lt_point903':c88hi<Q(903,1000),'new_full_cost_gap':Q(59,50)-c88hi-Q(40,11*4096)>Q(275,1000)}
assert all(checks.values())
Path(__file__).with_name('receipt.json').write_text(json.dumps({'status':'PASS','pid':os.getpid(),'checks':checks,'c_infty_lower':str(linf),'c5_upper':str(c5hi),'c88_upper':str(c88hi),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'complete scalar comparison conditional on accepted R090 theorem, no dimension threshold'},indent=2))
print('PASS: 4 new exact full-cost family comparisons')
