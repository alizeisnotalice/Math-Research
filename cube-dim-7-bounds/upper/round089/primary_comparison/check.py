from fractions import Fraction as Q
from pathlib import Path
import json,hashlib,os

def lb(z):
 u=(z-1)/(z+1);a=2*sum((u**(2*j+1)/Q(2*j+1) for j in range(24)),Q(0))
 return a,a+2*u**49/(49*(1-u*u))
c81_hi=Q(13,10)-lb(Q(5,2))[0]
c88_lo=Q(3,2)-lb(Q(20,11))[1]
c89_hi=Q(3193,2000)-lb(Q(10000,4987))[0]
checks={'c81_upper':c81_hi<Q(384,1000),'full_comparison_gap':Q(9,10)-c81_hi-Q(5,2**27)>Q(5159,10000),'do_not_claim_eventual_dominance_over_H88':c88_lo>c89_hi}
assert all(checks.values())
Path(__file__).with_name('receipt.json').write_text(json.dumps({'status':'PASS','pid':os.getpid(),'checks':checks,'c81_upper':str(c81_hi),'c88_lower':str(c88_lo),'c89_upper':str(c89_hi),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'scalar branch comparison; R089 operator certification remains separate'},indent=2))
print('PASS: three exact full-cost branch comparisons')
