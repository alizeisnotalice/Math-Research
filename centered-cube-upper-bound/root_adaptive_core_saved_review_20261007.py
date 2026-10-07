from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import hashlib,json
base=Path(__file__).parent
src=base/'adaptive_core_source_budget_exact_guard_20261007_results.json'
d=json.loads(src.read_text()); count=0
def ck(v):
 global count
 assert v;count+=1
leaves=[2,12,12,13,13,13,14,14,30,30,31]
paths=[]
for leaf in leaves:
 path=[]
 while leaf:
  path.append(leaf);leaf//=2
 paths.append(path)
nodes=sorted(set().union(*map(set,paths)))
m=[F(1)]+[F(19+(i-1)%3,100) for i in range(1,10)]+[F(1,1000)]
for row in d['rounds']:
 n=row['n'];zeta=F(1,isqrt(n));T=isqrt(n)*(n+1).bit_length()**4
 q=F(99,100)*T/F(2**n);D=max(map(len,paths))
 weights={p:F(1) if p==2 else F(0) if p in paths[0] else F(1,D) for p in nodes}
 mp={p:sum(m[i] for i in range(11) if p in paths[i]) for p in nodes}
 charge=sum(weights[p]*mp[p] for p in nodes)
 ck(charge==F(row['charged_mass']));ck(D==row['D'])
 for path in paths:ck(sum(weights[p] for p in path)<=1)
 grid=[F(1)]
 while grid[-1]<2:grid.append(min(F(2),grid[-1]*F(n+1,n)))
 predicted={2:F(2)}
 for p in nodes:
  if p==2 or not weights[p]:continue
  # Minimum cardinality determined independently by descending mass.
  masses=sorted([m[i] for i in range(11) if p in paths[i]],reverse=True)
  required=mp[p]*(1-zeta*weights[p]);cum=F(0)
  for num,x in enumerate(masses,1):
   cum+=x
   if cum>=required:break
  largest=None
  for r in grid:
   if q*num*r**n<=T*weights[p]*mp[p]:largest=r
   else:break
  if largest is not None:predicted[p]=largest
 saved={int(p):F(r) for p,r in row['radii'].items()}
 ck(predicted==saved)
 c=row['hard_child_only'];y,z=c['y'],c['z']
 lca=max(set(paths[y])&set(paths[z]),key=int.bit_length)
 ck(lca==c['original_LCA']==1);ck(2 in paths[z] and 2 not in paths[y])
 ck(F(row['positive_marginals']['strict_new_hard_child_traffic'])>0)
 big=row['compressed_large_dimension'];N=int(big['n']);K=big['k_n'];t=isqrt(N)*K**4
 ck(F(big['T_over_n'])==F(t,N));ck(big['T_less_n_over_2']==(2*t<N))
ck(len(d['checks'])==d['total_checks']==1381)
ck(hashlib.sha256(json.dumps(d['checks'],sort_keys=True).encode()).hexdigest()==d['checks_sha256'])
out={'status':'PASS','checks':count,'method':'heap-parent paths and descending-mass cardinality certificates; independent radius selection reconstruction','scope':'all three saved source charge/radius/eligibility records and large-n budget ratios; arrays read-audited, saved positive traffic sign checked but not independently recalculated; no actual FIRST verification','source_sha256':hashlib.sha256(src.read_bytes()).hexdigest()}
(base/'root_adaptive_core_saved_review_20261007.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out))
