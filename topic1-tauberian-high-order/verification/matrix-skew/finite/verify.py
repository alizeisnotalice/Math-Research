from fractions import Fraction as F
from decimal import Decimal
from pathlib import Path
import json,time
from compute import case
base=Path(__file__).resolve().parent
if __name__=='__main__':
 t=time.perf_counter();data=json.loads((base/'results'/'certificates.json').read_text());items=data['cases'];assert len(items)==15
 for item in items:
  assert len(item['attempts'])>=1
  for r in item['attempts']:
   u=F(r['u']);assert u.numerator==1 and u.denominator&(u.denominator-1)==0;j=u.denominator.bit_length()-1
   assert r.get('roots') is not None or r['status']=='guard_failed'
   fresh,exact=case(r['n'],j,roots=r.get('roots'))
   for k in ('n','u','kappa','lambda','status','guards'):assert fresh[k]==r[k],(r['n'],k)
   if exact is None:continue
   for k in ('roots','exact_parameters','N','center','finite_halo_minus_Ln','definitions'):assert fresh[k]==r[k],(r['n'],k)
   lo,hi=map(lambda x:F(Decimal(x)),r['D_minus_B']);assert lo<=exact[0]<=exact[1]<=hi
   for name,i in (('A_bounds',4),('B_bounds',6)):
    lo,hi=map(lambda x:F(Decimal(x)),r[name]);assert lo<=exact[i]<=exact[i+1]<=hi
 print(json.dumps({'verify':'PASS','cases':len(items),'certified':sum(x['attempts'][-1]['status']=='certified' for x in items),'seconds':time.perf_counter()-t,'scope':'reads roots; reuses generator formulas; no bisection'}))
