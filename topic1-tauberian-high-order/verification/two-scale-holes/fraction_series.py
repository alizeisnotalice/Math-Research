from fractions import Fraction as F
import json
N=10
zero=lambda:[F(0)]*(N+1)
def const(x):
 p=zero();p[0]=F(x);return p
def add(*ps):return [sum(p[i] for p in ps) for i in range(N+1)]
def scale(p,k):return [x*k for x in p]
def mul(p,q):
 r=zero()
 for i,x in enumerate(p):
  for j,y in enumerate(q[:N+1-i]):r[i+j]+=x*y
 return r
def inv(p):
 assert p[0]!=0
 r=zero();r[0]=1/p[0]
 for n in range(1,N+1):r[n]=-sum(p[k]*r[n-k] for k in range(1,n+1))/p[0]
 return r
def power(p,k):
 if k<0:return power(inv(p),-k)
 r=const(1)
 for _ in range(k):r=mul(r,p)
 return r
u=zero();u[1]=F(1);g=inv(add(const(1),scale(u,2)));a=scale(power(u,2),4)
k=add(scale(power(g,2),F(16,3)),scale(mul(power(g,3),u),F(640,81)))
d=add(const(1),scale(mul(power(k,2),power(u,6)),-1))
l0=scale(mul(power(k,3),inv(g)),F(1,32))
l1=scale(mul(power(k,4),inv(power(g,2))),F(-5,768))
l2=scale(mul(power(k,5),inv(power(g,3))),F(13,8192))
old=add(l0,mul(l1,u),mul(l2,power(u,2)))
def ratio(lam):
 w=add(mul(u,g),mul(lam,power(u,6)));rho=add(const(1),scale(w,-2));return mul(d,inv(add(d,scale(mul(a,power(rho,2)),-1))))
oldratio=ratio(old);Ln=inv(add(const(1),scale(mul(a,power(g,2)),-1)))
oldcoef=add(oldratio,scale(Ln,-1))[10]/1024
out=[]
for c in (0,1,7):
 new=add(old,scale(mul(a,l0),-1),scale(power(u,3),c));diff=add(ratio(new),scale(oldratio,-1));gain=diff[10]/1024
 out.append({'c':c,'coefficients_0_to_9_zero':all(x==0 for x in diff[:10]),'normalized_gain':str(gain),'old_av8':str(oldcoef),'new_av8':str(oldcoef+gain),'all_exact':gain==F(8,27) and oldcoef+gain==F(17638,6561)})
print(json.dumps(out,indent=2))
