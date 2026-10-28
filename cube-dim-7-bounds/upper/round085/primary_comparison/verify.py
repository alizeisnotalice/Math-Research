from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import hashlib,json,os

def log_bounds(z):
    assert z>1
    u=(z-1)/(z+1)
    v=2*sum((u**(2*j+1)/Q(2*j+1) for j in range(24)),Q(0))
    return v,v+2*u**49/(49*(1-u*u))
a_hi=sum((Q(-1,2)**k/factorial(k) for k in range(25)),Q(0))
a_lo=a_hi+Q(-1,2)**25/factorial(25)
cs_lo=log_bounds(1/(3-4*a_lo))[0]
old_hi=Q(1501,1000)-log_bounds(Q(1000,387))[0]
checks={
 'candidate_coefficient_below_family_limit':cs_lo>Q(5553,10000),
 'old_full_cost_coefficient_below_5517_over_10000':old_hi<Q(5517,10000),
 'old_cap_loss_error_below_1_over_10000':Q(2000,387*65536)<Q(1,10000),
 'eventual_full_cost_gap_gt_35_over_10000':Q(5553-5517,10000)-Q(2000,387*65536)>Q(35,10000)
}
assert all(checks.values())
out={'status':'PASS','pid':os.getpid(),'checks':checks,'family_limit_lower':str(cs_lo),'old_coefficient_upper':str(old_hi),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'exact scalar comparison only; parameter-family theorem audited separately'}
Path(__file__).with_name('receipt.json').write_text(json.dumps(out,indent=2))
print('PASS: 4 new exact full-cost comparison checks')
