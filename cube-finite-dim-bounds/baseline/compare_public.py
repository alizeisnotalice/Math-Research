"""Strict rational-root interval comparisons for the public and internal baselines."""
from fractions import Fraction as F
from verify_bounds import floor_root
from pathlib import Path
import json,sys
p=Path(__file__).resolve().parents[1]
def add(a,b):return a[0]+b[0],a[1]+b[1]
def neg(a):return -a[1],-a[0]
def sub(a,b):return add(a,neg(b))
def mul(a,b):
 v=[x*y for x in a for y in b];return min(v),max(v)
def scalar(a,c):return mul(a,(F(c),F(c)))
def power(a,n):
 r=(F(1),F(1))
 for _ in range(n):r=mul(r,a)
 return r
def root(numer,denom,n):
 D=10**40;k=floor_root(numer*D**n//denom,n)
 assert denom*k**n<=numer*D**n<denom*(k+1)**n
 return F(k,D),F(k+1,D)
old=json.loads((p/'baseline/original_certificates.json').read_text())
pub=json.loads((p/'results/public_lower_run.json').read_text())['result']
t=F(int(pub['ratio_num']),int(pub['ratio_den']))
out=[]
for c in old:
 n=c['dimension'];u=root(2,1,n);w=root(3,2,n);a=scalar(add((F(1),F(1)),u),F(1,2))
 x=sub(scalar(a,3),scalar(w,2));y=scalar(sub((F(3),F(3)),u),F(1,2))
 A=sub(scalar(power(a,n-1),F(3,2)),scalar(sub(power(x,n-1),power(y,n-1)),F(1,2)))
 Rlo=F(c['lower_bound_interval'][0]);assert Rlo>max(A[1],t)
 if n==3:assert t>A[1]
 else:assert A[0]>t
 out.append({'dimension':n,'Aldaz_interval':list(map(str,A)),'public_winner':'Tomar_monotonicity' if n==3 else 'Aldaz','internal_strictly_above_public':True})
print(json.dumps(out,indent=2))
if '--save' in sys.argv:(p/'results/public_comparison.json').write_text(json.dumps(out,indent=2)+'\n')
