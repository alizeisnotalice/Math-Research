from pathlib import Path
from fractions import Fraction as F
from math import comb
import hashlib,json
base=Path(__file__).parent
src=base/'frozen_jump_cross_layer_face_guard_results_20261007.json'
data=json.loads(src.read_text())
checks=0
for rnd in data['rounds']:
 n=rnd['dimension']
 for row in rnd['coupon']:
  k=row['k']
  # Independent inclusion-exclusion; the author uses positive occupancy recurrence.
  full=sum((-1)**j*comb(n,j)*(n-j)**k for j in range(n+1))
  missing=F(n**k-full,n**k)
  assert missing==F(row['missing_mass']);checks+=1
  assert hashlib.sha256(str(full).encode()).hexdigest()==row['all_visited_count_sha256'];checks+=1
  assert 0<=missing<=min(F(1),F(n*(n-1)**k,n**k));checks+=1
 assert F(rnd['proper_face_tail_over_W_analytic_upper'])<F(1,8);checks+=1
out={'status':'PASS','checks':checks,'method':'independent inclusion-exclusion reconstruction of all 21 saved coupon records and hashes','not_checked':'floating Gamma diagnostic values are read-audited only; no full high-tail or cube claim','source_sha256':hashlib.sha256(src.read_bytes()).hexdigest()}
(base/'root_face_coupon_saved_review_20261007.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out))
