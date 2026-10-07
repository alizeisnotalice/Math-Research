from pathlib import Path
from fractions import Fraction as F
from math import comb
from hashlib import sha256
import json
root=Path('/Users/zhengzhihao/Desktop/T/output/cube_general_20261003/geom_20261006')
data=json.loads((root/'low_future_source_payment_guard_20261007_results.json').read_text())
count=0
for r in data['rounds']:
 for c in r['cases']:
  assert F(c['traffic'])<=F(c['posterior_bound'])<=F(c['kernel_bound']);count+=1
s=json.loads((root/'low_future_replacement_exact_guard_20261007_results.json').read_text())
def check(q,rec):
 for label,k in [('numerator',q.numerator),('denominator',q.denominator)]:
  assert k.bit_length()==rec[label+'_bits']
  assert sha256(k.to_bytes(max(1,(k.bit_length()+7)//8),'big')).hexdigest()==rec[label+'_sha256_binary']
checks=0
for row in s['rounds']:
 n,a=row['n'],row['alpha']; den=comb(n+a,a)
 for t in row['tail_checks']:
  M=t['M'];q=F(1)
  for j in range(a):q*=F(n-M+j,n+1+j)
  check(q,t['tail']);checks+=1
 # Truncated forward/reversed integer polynomial recurrence, independent
 # of grouped binomial convolution. Exact same saved coefficient targets.
 def polynomial(constants,linear,degree):
  co=[1]+[0]*degree
  for index,c in enumerate(constants):
   for j in range(min(degree,index+1),0,-1):co[j]=c*co[j]+linear*co[j-1]
   co[0]*=c
  return co
 vals=[3]*row['mixed_inside']+[8]*row['mixed_outside']
 reverse=polynomial(vals,4,a)
 for c in row['mixed_coefficients']:
  k=c['k']
  if k==0:q=F(1)
  elif k==1:q=F(sum(vals),4*n)
  else:q=F(reverse[n-k],4**n*comb(n,k))
  check(q,c['B_k']);checks+=1
out={'saved_source_inequality_chains':count,'independent_spectrum_receipts':checks,'status':'PASS','execution_note':'Full Fraction coefficient recurrence stopped for cost; exact truncated forward/reverse integer recurrence used, thresholds unchanged','scope':'Saved rational results and independent coefficient recurrence only; no FIRST certification.'}
(root/'root_low_future_saved_review_20261007.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out))
