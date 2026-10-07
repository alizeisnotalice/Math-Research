from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
import json
p=Path(__file__).resolve().parent
j=json.loads((p/'actual_gate_contract_exact_guard_20261007_results.json').read_text())
def verify(x,r):
 if isinstance(r,str):assert x==F(r);return
 num=x.numerator.to_bytes((x.numerator.bit_length()+7)//8,'big');den=x.denominator.to_bytes((x.denominator.bit_length()+7)//8,'big')
 raw=len(num).to_bytes(8,'big')+num+len(den).to_bytes(8,'big')+den
 assert sha256(raw).hexdigest()==r['sha256_framed_binary']
 assert x.numerator.bit_length()==r['numerator_bits'] and x.denominator.bit_length()==r['denominator_bits']
count=0
for row in j['rounds']:
 n=row['n']; k=row['k'];assert 2**(k-1)<n+2<=2**k
 assert F(row['Hmark_over_a_certified_lower'])>k
 for ix,r in enumerate(row['receiver_box_components']):
  # Reconstruction from box delta multiplier (0,1,2,4)/(64n).
  mult=(0,1,2,4)[ix]
  ratio=1+F(mult,64*(n+16))
  qRn=ratio**n/3
  verify(ratio,r['R_over_R0']);verify(qRn,r['qR_to_n']);verify(n*qRn*qRn,r['squared_CP_threshold_lower'])
  assert n*qRn*qRn>1
  for mass,m in zip((F(0),F(1,2**n),F(1,3),F(1)),r['mass_components']):
   verify(mass,m['h']);assert m['strict_CP_possible'] is False and m['GP_strict_for_this_h']==(mass>qRn)
   count+=1
out={'status':'PASS','reconstructed_parent_components':count,'source':'Independent reconstruction of all saved parameter fractions and binary receipts; original gate definitions reviewed in tex. No new original-kernel tests.'}
(p/'root_actual_gate_saved_review_20261007.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out))
