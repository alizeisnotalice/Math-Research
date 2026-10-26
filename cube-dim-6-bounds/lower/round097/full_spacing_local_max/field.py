"""Pure accepted E058 field formulas only; old review main/run is never executed."""
from fractions import Fraction as Q
from decimal import Decimal, localcontext
from math import prod
DEN=10**100
Z=(Q(0),)*36
def constant(x):return (Q(x),)+Z[1:]

def mono(i,j,c=Q(1)):
 a=list(Z);a[6*(i%6)+j%6]=Q(c)*2**(i//6)*3**(j//6);return tuple(a)

def parse(d):
 a=list(Z)
 for ij,c in d.items():
  i,j=map(int,ij.split(','));assert 0<=i<6 and 0<=j<6 and type(c) in (str,int)
  a[6*i+j]=Q(c)
 return tuple(a)

def encode(a):return {f'{i//6},{i%6}':str(c) for i,c in enumerate(a) if c}

def add(a,b):return tuple(x+y for x,y in zip(a,b))

def scale(a,c):return tuple(x*c for x in a)

def sub(a,b):return add(a,scale(b,-1))

def mul(a,b):
 r=list(Z)
 for i,x in enumerate(a):
  if x:
   for j,y in enumerate(b):
    if y:
     p,q=i//6+j//6,i%6+j%6
     r[6*(p%6)+q%6]+=x*y*2**(p//6)*3**(q//6)
 return tuple(r)

def rootbox(m):
 lo,hi=DEN,2*DEN;target=m*DEN**6
 assert lo**6<=target<hi**6
 while hi-lo>1:
  mid=(lo+hi)//2
  if mid**6<=target:lo=mid
  else:hi=mid
 assert lo**6<=target<(lo+1)**6
 return Q(lo,DEN),Q(lo+1,DEN)

def bounds(a):
 lo=hi=Q(0)
 for c,(l,h) in zip(a,POWER_BOUNDS):
  lo+=c*(l if c>=0 else h);hi+=c*(h if c>=0 else l)
 return lo,hi

def sign(a):
 if a==Z:return 0
 lo,hi=bounds(a)
 assert lo>0 or hi<0,'Unresolved exact sign'
 return 1 if lo>0 else -1

def cmp(a,b):return sign(sub(a,b))

def minimum(a,b):return a if cmp(a,b)<=0 else b

def maximum(a,b):return a if cmp(a,b)>=0 else b

def dec(q):
 with localcontext() as c:
  c.prec=65;return str(Decimal(q.numerator)/Decimal(q.denominator))

def view(a):
 lo,hi=bounds(a)
 return {'coefficients':encode(a),'exact_interval':[str(lo),str(hi)],'display_lower':dec(lo)}

def rad(ell):
 n=prod(ell);p=q=0
 while n%2==0:p+=1;n//=2
 while n%3==0:q+=1;n//=3
 assert n==1
 return mono(p,q,Q(1,2))
UL,UH=rootbox(2)
VL,VH=rootbox(3)
POWER_BOUNDS=[(UL**(k//6)*VL**(k%6),UH**(k//6)*VH**(k%6)) for k in range(36)]
