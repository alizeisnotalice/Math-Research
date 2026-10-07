from pathlib import Path
from fractions import Fraction as F
from math import comb
import hashlib,json
D=Path(__file__).parent
src=D/'low_count_future_capacity_exact_guard_results_20261007.json'
data=json.loads(src.read_text()); checks=0; sample_count=0

def check(x):
 global checks
 checks+=1; assert x

def digest(x):
 chunks=[]
 for v in (x.numerator,x.denominator):
  b=abs(v).to_bytes(max(1,(abs(v).bit_length()+7)//8),'big')
  chunks.append((b'-' if v<0 else b'+')+len(b).to_bytes(8,'big')+b)
 return hashlib.sha256(b''.join(chunks)).hexdigest()

def cert(x,c):
 check(digest(x)==c['framed_pair_sha256'])
 check(abs(x.numerator).bit_length()==c['numerator_bits'])
 check(x.denominator.bit_length()==c['denominator_bits'])
 if 'exact' in c: check(x==F(c['exact']))

for rnd in data['rounds']:
 n,alpha=rnd['n'],rnd['alpha']; eps=F(rnd['checker_epsilon'])
 for model in rnd['models']:
  sig=F(model['sigma']); rho=(1-sig)/sig
  limit=max(x['k'] for x in model['sample_weights'])
  ws=[F(1,comb(n+alpha,alpha))]; delta=F(0)
  # Independent all-positive integration-by-parts recurrence, no Beta sum.
  for k in range(limit):
   delta=(alpha/sig*ws[-1]+k*rho*delta)/(n-k)
   check(delta>0); ws.append(ws[-1]+delta)
  for x in model['sample_weights']:
   cert(ws[x['k']],x['w']); check((ws[x['k']]>=eps)==x['w_ge_checker_epsilon']); sample_count+=1
  cc=model['cutoff']; hits=[k for k,w in enumerate(ws) if w>=eps]
  if 'exact_within_registered_range' in cc:
   check(hits[0]==cc['exact_within_registered_range'])
  else:
   check(not hits); check(limit==cc['strictly_greater_than']); check(cc['analytic_upper_bound']==n)
  rn=model['independent_RN_identity']; d=rn['supported_soft_coordinates']; p=F(rn['current_count_p'])
  Z=sum((comb(d,k)*p**k*(1-p)**(d-k)*ws[k] for k in range(d+1)),F(0))
  cert(Z,rn['future_RN_ratio'])
  for c in model['joint_fee_models']:
   old,new,paid=map(F,(c['old_paid_mass'],c['new_label_paid_mass'],c['combined_paid_mass']))
   check(old+new==paid); check(0<=paid<=1)
   check(old==int(c['source_high']))
  ep=rnd['endpoints']; cert(ws[0],ep['sigma0_w0'])
  check(F(ep['sigma0_new_future_allsoft_prior_mass'])==F(alpha,n+alpha))
 sp=F(1,alpha*((1<<n)-1)); cert(sp,rnd['source_entropy']['source_p'])
 check(sp*(1-F(1,1<<n))==F(1,alpha*(1<<n)))
 for x in rnd['source_entropy']['binary_depths']: check(F(x['full_mass_layer_reuse'])==x['depth'])
for key,name in [('script_sha256','low_count_future_capacity_exact_guard_20261007.py'),('registration_sha256','low_count_future_capacity_guard_registration_20261007.json')]:
 check(hashlib.sha256((D/name).read_bytes()).hexdigest()==data[key])
result={'status':'PASS_INDEPENDENT_POSITIVE_RECURRENCE','checks':checks,'saved_sample_weights':sample_count,'cutoff_models':9,'RN_ratios':9,'method':'positive integration-by-parts recurrence; author Beta-sum script neither imported nor run','source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'scope':'coefficient identities and registered rational epsilon only; not actual cube/FIRST or theoretical log threshold'}
(D/'root_low_count_saved_review_20261007.json').write_text(json.dumps(result,indent=2)+'\n'); print(json.dumps(result))
