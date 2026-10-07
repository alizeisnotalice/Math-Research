from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import json,hashlib
D=Path(__file__).parent
src=D/'double_cutoff_actual_geometry_exact_guard_20261007_results.json'
data=json.loads(src.read_text()); count=0

def check(v):
 global count
 count+=1; assert v

def prefix(node): return bin(node)[3:]
def leaves(node): return [i for i in range(8) if format(i,'03b').startswith(prefix(node))]
def ancestor(node,leaf): return format(leaf,'03b').startswith(prefix(node))
patterns=[{1,2,4,8,9,5},{2,5,10,11,12,14},set(range(8,16))]
m=[F(16+i*i,17) if i%2==0 else F(1,64*(i+1)) for i in range(8)]
h=[m[i]*F(2+i%3,4) for i in range(8)]; W=sum(m)
for rnd in data['rounds']:
 n=rnd['n']; zeta=F(1,isqrt(n)); cores={}
 for node in range(1,16):
  ids=leaves(node); keep=[i for i in ids if i%2==0]
  if node!=1: keep += [i for i in ids if i%2][:1]
  if sum(m[i] for i in ids if i not in keep)>zeta*sum(m[i] for i in ids): keep=ids
  cores[node]=set(keep)
 S=[]; H=[]
 for x in range(12):
  a=[F(1+(x+2*i)%7,3+i) for i in range(8)]
  b=[F(1+(2*x+i)%5,2+i) for i in range(8)]
  S.append([a[i]/(sum(a)*m[i]) for i in range(8)])
  H.append([F(5,2)*b[i]/(sum(b)*h[i]) for i in range(8)])
 den=max(sum(max(S[x][i],H[x][i]) for x in range(12)) for i in range(8)); rx=1/den
 for j,marked in enumerate(patterns):
  top=sorted(p for p in marked if not any(p!=q and prefix(p).startswith(prefix(q)) for q in marked))
  bad={i for p in top for i in leaves(p) if i not in cores[p]}; mb=sum(m[i] for i in bad)
  rec=rnd['antichain_receipts'][j]
  check(top==rec['top']); check(sorted(bad)==rec['bad_ids']); check(mb==F(rec['bad_mass'])); check(mb<=zeta*W)
  cp=gp=F(0); reclass=0
  for y in range(8):
   for z in range(8):
    if y==z: continue
    pa,pb=format(y,'03b'),format(z,'03b'); common=''
    for a,b in zip(pa,pb):
     if a!=b: break
     common+=a
    parent=int('1'+common,2)
    selected=[p for p in top if ancestor(p,y) and ancestor(p,z)]
    check(len(selected)<=1)
    if not selected: continue
    A=selected[0]
    for v in (y,z): reclass+=int(v in cores[parent] and v not in cores[A])
    for x in range(12):
     weight=m[y]*S[x][y]*h[z]*H[x][z]/(2+(x+y+z)%9)
     gp+=weight*int(z not in cores[A])
     if x%3!=1: cp+=weight*int(z not in cores[A])
     else: cp+=weight*int(y not in cores[A])/3
  cp*=rx; gp*=rx
  check(reclass==rec['child_good_ancestor_bad_occurrences'])
  saved=rnd['marginal_receipts'][j]
  check(cp==F(saved['CP_bad_traffic'])); check(gp==F(saved['GP_bad_traffic']))
  check(gp<=mb); check(cp<=F(5,2)*mb)
# Hash of saved predicate list checks data integrity; it is not a rerun.
check(hashlib.sha256(json.dumps(data['checks'],sort_keys=True).encode()).hexdigest()==data['checks_sha256'])
check(all(c['pass'] for c in data['checks'])); check(len(data['checks'])==1530)
result={'status':'PASS_INDEPENDENT_PREFIX_AND_MARGINAL_RECONSTRUCTION','checks':count,'rounds':3,'antichain_patterns':9,'marginal_receipts':9,'method':'binary-prefix ancestry and independent pair sums; author script neither imported nor run','not_reconstructed':'box-union volumes and radius predicates read-audited only; no actual FIRST data','source_sha256':hashlib.sha256(src.read_bytes()).hexdigest()}
(D/'root_core_capacity_saved_review_20261007.json').write_text(json.dumps(result,indent=2)+'\n'); print(json.dumps(result))
