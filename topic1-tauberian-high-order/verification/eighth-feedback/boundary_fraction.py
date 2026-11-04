from fractions import Fraction as F
from math import comb
import json

def add(*polys):
 out=[F(0)]*max(map(len,polys))
 for p in polys:
  for i,c in enumerate(p):out[i]+=c
 return out
def mul(p,q):
 r=[F(0)]*(len(p)+len(q)-1)
 for i,x in enumerate(p):
  for j,y in enumerate(q):r[i+j]+=x*y
 return r
def sc(a,p):return [a*x for x in p]
def integ(p,l,r):return sum(c*(r**(i+1)-l**(i+1))/F(i+1) for i,c in enumerate(p))
out=[]
for d in (F(1,4),F(1,10),F(1,100)):
 theta=d/(1-d); mlo=[F(0),F(1),-1/(2*theta)]
 # integral of low branch to d + theta integral_d^q (1-y)
 mhi=[d-d*d/(2*theta)-theta*(d-d*d/2),theta,-theta/2]
 def model(m):return add(sc(2*d,mul([1,1],m)),sc(-2*d*d,add(sc(2,mul([0,1],m)),sc(-1,mul(m,m)))))
 actual=integ(model(mlo),F(0),d)+integ(model(mhi),d,F(1))
 closed=d*d*(155-155*d+86*d*d-22*d**3-4*d**4)/60
 out.append({'delta':str(d),'M0_continuous':sum(c*d**i for i,c in enumerate(mlo))==sum(c*d**i for i,c in enumerate(mhi)),'integral_exact':str(actual),'closed_equal':actual==closed})
for n in (3,4):
 for u in (F(1,16),F(1,64)):
  k=F(8*n,3);w=u/(1+2*u);e=k*u**3;a=(2*u)**n;f=4*n*u*u/(1+2*u);B=comb(n,2)
  inside=n*comb(n-1,1)*2*(2*w)**(n-2)*e**3/(2*a*f)
  outside=n*comb(n-1,2)*2*2*(2*w)**(n-3)*e**4/(4*a*f)
  P=2*B*(2*w)**(n-2)*e**3/(a*f); T=B*(n-2)*(2*w)**(n-3)*e**4/(a*f)
  out.append({'n':n,'u':str(u),'inside_P_equal':inside==P,'outside_T3_equal':outside==T,'positive':inside>0 and outside>0,'scope':'outer box integral counts; no actual deleted set'})
print(json.dumps(out,indent=2))
