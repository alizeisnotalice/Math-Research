from fractions import Fraction as Q
from math import comb
from pathlib import Path
import json
zero=[Q(0)]*4
def add(*p):return [sum(q[i] for q in p) for i in range(4)]
def scale(p,c):return [c*t for t in p]
def power(z,k):
 out=[Q(1),Q(0),Q(0),Q(0)]
 for _ in range(k):out=[sum(out[j]*z[i-j] for j in range(i+1)) for i in range(4)]
 return out
s=Q(1,9);a=Q(10,27);al=Q(729,29);g=Q(27,29)
def c(z,b,mid):
 v=z[0]+z[1]*mid
 if v<0:z=scale(z,-1);v=-v
 if v<b:return add([2*b**3/3,0,0,0],scale(power(z,2),-b),scale(power(z,3),Q(1,2)))
 if v<2*b:return scale(power(add([2*b,0,0,0],scale(z,-1)),3),Q(1,6))
 return zero
breaks=sorted({Q(0),s,2*s,a,2*a,1-s,1-2*s,1-a,1-2*a,Q(1)})
rows=[]
for l,r in zip(breaks[:-1],breaks[1:]):
 m=(l+r)/2;x=[Q(0),Q(1),Q(0),Q(0)];xm=[Q(-1),Q(1),Q(0),Q(0)]
 f=add([1+al*(2*a**3/3+2*s**3),0,0,0],scale(c(x,a,m),-al),scale(c(xm,a,m),-al),scale(c(x,s,m),-3*al),scale(c(xm,s,m),-3*al))
 d=add(f,scale(power([Q(1),Q(1),Q(0),Q(0)],2),-g))
 t=zero[:]
 for j in range(4):t=add(t,scale(power([l,r-l,Q(0),Q(0)],j),d[j]))
 bern=[sum(t[j]*Q(comb(k,j),comb(3,j)) for j in range(k+1)) for k in range(4)]
 if r<=s:assert all(z>=0 for z in bern)
 elif l>=2*s:assert all(z<=0 for z in bern)
 else:assert all(z==0 for z in bern)
 rows.append(dict(left=str(l),right=str(r),polynomial=[str(z) for z in d],bernstein=[str(z) for z in bern]))
Path(__file__).with_name('exact.json').write_text(json.dumps(rows,indent=2));print(json.dumps(rows,indent=2))
