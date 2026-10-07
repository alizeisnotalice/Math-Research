from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
import json, numpy as np
p=Path(__file__).resolve().parent
s=p/'joint_low_future_contract_probe_20261007_v2.py'
j=json.loads((p/'joint_low_future_contract_probe_20261007_v2.json').read_text())
assert sha256(s.read_bytes()).hexdigest()==j['script_sha256']
assert 'copy_negate()' in s.read_text()
checks=[]
for row in j['rounds']:
 n=row['n'];R=F(row['radius_over_a']);d=F(row['spacing_over_a'])
 vals=[]
 for index,offset in enumerate(row['offsets_over_radius']):
  t=float(F(offset)); lo,hi=map(F,row['original_phi_intervals'][index])
  qs=[]
  for N in (128,256):
   nodes,weights=np.polynomial.legendre.leggauss(N);v=(nodes+1)/2
   density=v*(-np.expm1(-(.5+t)/v)-np.expm1(-(.5-t)/v))
   q=float(np.dot(weights,density)/2);assert float(lo)<q<float(hi);qs.append(q)
  vals.append({'phi_gauss_128':qs[0],'phi_gauss_256':qs[1],'difference':abs(qs[0]-qs[1])})
 sl,sh=map(F,row['sigma_interval']);assert F(1,n)<sl<sh<F(1,16)
 avg=sum(v['phi_gauss_256'] for v in vals)/2
 sig=-np.expm1(-np.log(3)/n)/(1-avg)
 assert float(sl)<sig<float(sh)
 # Direct exact hard-arrival values and nonconcentration certificate.
 assert R-2*d<1<R<2
 assert F(1,2**n)<1/(3*R**n)
 assert F(1,2**n)<F(49,65536*n)
 checks.append({'n':n,'original_phi_independent_floating_checks':vals,'sigma_independent_float':float(sig),'exact_hard_arrival_and_microbox_guards':True})
out={'status':'PASS','source':'Only v2 with exact Decimal copy_negate endpoints reviewed. v1 invalid for interval certification.','scope':'Independent floating Gauss quadrature inside reviewed directed Riemann intervals, exact finite hard arrival checks. No birth/CPGP/history certification.','rounds':checks}
(p/'root_joint_contract_review_20261007.json').write_text(json.dumps(out,indent=2))
print(json.dumps({'status':'PASS','rounds':len(checks),'max_GL_difference':max(v['difference'] for r in checks for v in r['original_phi_independent_floating_checks'])}))
