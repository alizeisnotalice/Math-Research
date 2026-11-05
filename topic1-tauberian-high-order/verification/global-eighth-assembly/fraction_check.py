from fractions import Fraction as F
from math import comb,factorial
import json
N=10
zero=lambda:[F(0)]*(N+1)
def const(x):
 p=zero();p[0]=F(x);return p
def add(*ps):return [sum(p[i] for p in ps) for i in range(N+1)]
def sc(p,k):return [x*k for x in p]
def mul(p,q):
 r=zero()
 for i,x in enumerate(p):
  for j,y in enumerate(q[:N+1-i]):r[i+j]+=x*y
 return r
def inv(p):
 r=zero();r[0]=1/p[0]
 for n in range(1,N+1):r[n]=-sum(p[k]*r[n-k] for k in range(1,n+1))/p[0]
 return r
def pw(p,k):
 if k<0:return pw(inv(p),-k)
 r=const(1)
 for _ in range(k):r=mul(r,p)
 return r
v=zero();v[1]=F(1); A=pw(add(const(1),v),2);a=pw(v,2);K=F(4,27);c8=-F(13,18)+F(400,243)
Q=sc(add(mul(pw(v,6),pw(add(const(1),v),-4)),sc(mul(pw(v,7),pw(add(const(1),v),-5)),F(10,9)),sc(mul(pw(v,8),pw(add(const(1),v),-2)),c8)),K)
def coef(q):
 d=mul(a,mul(A,mul(q,mul(inv(add(A,sc(a,-1))),inv(add(A,sc(a,-1),sc(mul(a,q),-1)))))))
 return d[10]
normal=coef(Q);inside=coef(add(Q,sc(pw(v,8),F(8,27))))
base=sum((-1)**(8-(j-1)*2)*comb(9,8-(j-1)*2) for j in range(1,6))
b8=2**10*(base+inside)
values=[]
for q in (F(0),F(1,2),F(1)):
 p=F(1)
 for j in range(1,5):p*=j-q
 values.append({'q':str(q),'sigma_e5_coefficient':str(p/120)})
print(json.dumps({'normal_J2':str(normal),'inside_J2':str(inside),'baseline_j8':base,'ordinary_b28':str(b8),'sigma_examples':values,'uniform_remainder_bound':str(F(120*64,factorial(6))),'all_exact':normal==F(15694,6561) and inside==F(17638,6561) and b8==F(1737988096,6561)},indent=2))
