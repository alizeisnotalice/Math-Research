from fractions import Fraction as Q
from pathlib import Path
import json
N=128
def at(x):
 s=2*sum((x**(2*j+1)/Q(2*j+1)for j in range(N)),Q())
 return s,s+2*x**(2*N+1)/(Q(2*N+1)*(1-x*x))
l2,u2=at(Q(1,3));rows=[]
for a in [32,36,40,64]:
 k=a.bit_length()-1;l,u=at(Q(a-2**k,a+2**k));ll=k*l2+l;uu=k*u2+u
 A=(1+Q(2,a))**7;B=(1+Q(1,a))**7
 for t in [6,7]:
  lo=(1+A/t)*(t+B*(1+7*ll));hi=(1+A/t)*(t+B*(1+7*uu))
  rows.append(dict(a=a,tau=t,lower=str(lo),upper=str(hi)))
w=min(rows,key=lambda z:Q(z['upper']));upper=Q(w['upper']);assert all(upper<Q(z['lower'])for z in rows if z is not w)
U=Q((upper.numerator*10**12+upper.denominator-1)//upper.denominator,10**12)
assert upper<U<Q(77433340322253,10**12)
print(json.dumps(dict(status='PASS independent exact rational eight-pair enclosures',N=N,winner=w,rational_C7_strict_upper=str(U),rows=rows,finite_set_optimum_only=True),indent=2))
