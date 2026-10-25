"""Exact finite extension audit; does not claim full lattice optimality."""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from pathlib import Path
import sys,json
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'baseline'))
from verify_bounds import floor_root
@lru_cache(None)
def root(k,n):
 D=10**40;q=floor_root(k*D**n,n)
 assert q**n<=k*D**n<(q+1)**n
 return F(q,D),F(q+1,D)
def add(a,b):return a[0]+b[0],a[1]+b[1]
def scale(a,q):return a[0]*q,a[1]*q
def mul(a,b):return a[0]*b[0],a[1]*b[1]
def power(a,m):
 r=(F(1),F(1))
 for _ in range(m):r=mul(r,a)
 return r
out=[]
for n in [3,4]:
 u=root(2,n);v=root(3,n)
 a=[add((F(1),)*2,scale(power(u,n-2),F(1,2))),scale(add((F(2),)*2,mul(power(u,n-1),v)),F(1,4))]+[scale(add((F(1),)*2,u),F(1,2))]*(n-2)
 amin=min(aa[0] for aa in a)
 tail=F(1)
 for aa in a:tail*=1/aa[0]+1/(8*amin)
 assert tail<1
 empty=0;contained=[];gaps=[]
 for l in product(range(1,9),repeat=n):
  if max(l)<=4:continue
  N=1
  for x in l:N*=x
  r=scale(root(N,n),F(1,2))
  # Any nonempty center interval requires r>(l_j-1)*a_j/2.
  ds=[r[1]-F(x-1,2)*aa[0] for x,aa in zip(l,a)]
  if min(ds)<0:empty+=1;gaps.append(-min(ds));continue
  # Remaining all-odd box has lower corner 0, contained in singleton box [0,1/2]^n.
  assert all(x%2==1 for x in l)
  assert all(z<F(1,2) for z in ds)
  assert all(r[0]-F(x-1,2)*aa[1]>0 for x,aa in zip(l,a))
  contained.append(list(l))
 out.append({'dimension':n,'cap':8,'new_count_tuples':8**n-4**n,'strictly_empty':empty,'contained_in_singleton':contained,'smallest_empty_margin':str(min(gaps)),'tail_product_upper':str(tail),'all_counts_above_8_excluded':True})
print(json.dumps(out,indent=2))
if '--save' in sys.argv:(Path(__file__).resolve().parents[1]/'results/round002_extension.json').write_text(json.dumps(out,indent=2)+'\n')
