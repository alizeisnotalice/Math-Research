from pathlib import Path
from fractions import Fraction as F
import math,json,hashlib
P=Path(__file__).parent;B=1024;D=1<<B
class I:
 def __init__(self,a,b=None,raw=False):
  if raw:self.a,self.b=a,b
  else:
   a=F(a);b=a if b is None else F(b)
   self.a=(a*D).numerator//(a*D).denominator
   q=b*D;self.b=-((-q.numerator)//q.denominator)
  assert self.a<=self.b
 def __add__(self,o):
  o=o if isinstance(o,I) else I(o);return I(self.a+o.a,self.b+o.b,True)
 __radd__=__add__
 def __neg__(self):return I(-self.b,-self.a,True)
 def __sub__(self,o):return self+(-o if isinstance(o,I) else -I(o))
 def __mul__(self,o):
  o=o if isinstance(o,I) else I(o);v=[self.a*o.a,self.a*o.b,self.b*o.a,self.b*o.b]
  return I(min(v)//D,-((-max(v))//D),True)
 __rmul__=__mul__
 def inv(self):
  assert self.a>0;return I(D*D//self.b,-((-D*D)//self.a),True)
 def __pow__(self,n):
  v=I(1);a=self
  while n:
   if n&1:v=v*a
   a=a*a;n//=2
  return v
 def data(self):return {'lo_numerator':str(self.a),'hi_numerator':str(self.b),'denominator_power2':B}
def expminus(q):
 scale=1
 while q/scale>F(1,2):scale*=2
 x=q/scale;term=lo=F(1);M=256
 for j in range(1,M+1):term*=x/j;lo+=term
 hi=lo+(term*x/(M+1))/(1-x/(M+2))
 return I(lo,hi).inv()**scale
checks=0;rows=[]
for n in [32,512,8192]:
 K=0
 while (K+1)**3<=n:K+=1
 cost=F(K*(K+1)*(K+2),3*(n+1))
 assert cost==sum((F(k*(k+1),n+1) for k in range(1,K+1)),F(0));checks+=1
 assert cost<=1;checks+=1
 for k in range(1,K+1):
  star=F(k+1,n+1);c=math.comb(n,k)
  peak=I(star)**k*I(1-star)**(n-k)*c
  assert peak.b<=D;checks+=1
  major=peak*(k*star)
  assert major.b<=I(F(k*(k+1),n+1)).b;checks+=1
  for r in sorted(set([F(k,2*n),F(k,n),star,min(F(2*k,n),F(1,2))])):
   base=I(r)**k;bern=I(1-r)**(n-k)
   diff=base*(bern-expminus(n*r))*c
   upper=base*bern*(c*k*r)
   # All tested nodes have strictly positive comparison margin.
   assert max(0,diff.b)<=upper.a;checks+=1
   if r==star:
    assert max(upper.a,major.a)<=min(upper.b,major.b)
   else:
    assert upper.b<=major.a
   checks+=1
  rows.append({'n':n,'k':k,'critical_r':str(star),'aggregated_positive_envelope_upper':major.data(),'simplified_fee':str(F(k*(k+1),n+1))})
 rows.append({'n':n,'K':K,'total_envelope_fee':str(cost)})
reg=P/'direct_difference_low_support_registration_20261007.json'
out={'status':'PASS','checks':checks,'rows':rows,'scope':'original difference universal coefficient envelope only, not actual source response or complete weak estimate','registration_sha256':hashlib.sha256(reg.read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(P/'direct_difference_low_support_results_20261007.json').write_text(json.dumps(out,indent=2))
print(json.dumps({'status':'PASS','checks':checks,'round_costs':[x for x in rows if 'K'in x]}))
