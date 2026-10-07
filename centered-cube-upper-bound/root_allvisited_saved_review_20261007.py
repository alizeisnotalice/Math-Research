from pathlib import Path
from decimal import Decimal,localcontext
from fractions import Fraction as F
import hashlib,json,math
base=Path(__file__).parent
src=base/'frozen_jump_allvisited_hormander_guard_results_20261007.json'
d=json.loads(src.read_text());count=0;max_elastic=0.0;max_inverse=0.0

def ck(v):
 global count
 assert v;count+=1

def bd(v):
 return v/(1+v).ln()-1

with localcontext() as ctx:
 ctx.prec=80
 for rnd in d['rounds']:
  n=rnd['n'];ck(rnd['cutoff']==1024*n*n)
  for rec in rnd['symbol']:
   v=Decimal(str(rec['lambda']));l=(1+v).ln()
   e=v*(l-v/(1+v))/(l*l*bd(v))
   gap=abs(float(e)-rec['elasticity']);max_elastic=max(max_elastic,gap)
   ck(gap<1e-10);ck(Decimal(2)/3<=e<=1)
  for rec in rnd['density']:
   c=Decimal(str(rec['c']));a=Decimal(str(rec['t_over_c']))
   for key,target in [('beta_t',1/(c*a)),('beta_c',1/c)]:
    rel=abs(bd(Decimal(str(rec[key])))/target-1)
    max_inverse=max(max_inverse,float(rel));ck(rel<Decimal('1e-10'))
   ck(0<=rec['normalized_peak_fine']+rec['analytic_positive_tail_upper']<=5)
  for rec in rnd['contours']:
   rr=complex(*rec['real']);ss=complex(*rec['shifted'])
   ck(abs(abs(rr-ss)-rec['contour_discrepancy'])<1e-25)
  ck(F(rnd['exact_H_exponent_upper'])==12+3*n-24*n*n)
  ck(F(rnd['exact_L2_exponent_upper'])==F(2)+F(7*n,2)-F(3072*n,5))
for kind,name in [('registration','frozen_jump_allvisited_hormander_registration_20261007.json'),('script','frozen_jump_allvisited_hormander_guard_20261007.py')]:
 ck(hashlib.sha256((base/name).read_bytes()).hexdigest()==d[kind+'_sha256'])
out={'status':'PASS','checks':count,'max_80digit_elasticity_difference':max_elastic,'max_80digit_inverse_relative_difference':max_inverse,'method':'independent 80-digit Decimal original-symbol and inverse-clock verification, exact cutoff reconstruction, saved contour discrepancy arithmetic and file hashes','limits':'Decimal cross-checks are not interval certificates; no independent GL density or contour quadrature rerun; no actual cube transfer','source_sha256':hashlib.sha256(src.read_bytes()).hexdigest()}
(base/'root_allvisited_saved_review_20261007.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out))
