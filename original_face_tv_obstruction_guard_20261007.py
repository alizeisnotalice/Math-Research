from fractions import Fraction as F
from pathlib import Path
import json,hashlib,math
P=Path(__file__).parent; B=160; D=1<<B

def floorf(q):return q.numerator//q.denominator
def ceilf(q):return -((-q.numerator)//q.denominator)
class I:
 def __init__(self,a,b=None,raw=False):
  if raw:self.a,self.b=a,b
  else:
   a=F(a);b=a if b is None else F(b)
   self.a,self.b=floorf(a*D),ceilf(b*D)
  assert self.a<=self.b
 def __add__(self,o):
  o=o if isinstance(o,I) else I(o)
  return I(self.a+o.a,self.b+o.b,True)
 __radd__=__add__
 def __neg__(self):return I(-self.b,-self.a,True)
 def __sub__(self,o):return self+(-o if isinstance(o,I) else -I(o))
 def __rsub__(self,o):return I(o)+(-self)
 def __mul__(self,o):
  o=o if isinstance(o,I) else I(o)
  q=[self.a*o.a,self.a*o.b,self.b*o.a,self.b*o.b]
  return I(min(q)//D,-((-max(q))//D),True)
 __rmul__=__mul__
 def inv(self):
  assert self.a>0
  return I(D*D//self.b,-((-D*D)//self.a),True)
 def __truediv__(self,o):return self*(o.inv() if isinstance(o,I) else I(o).inv())
 def __pow__(self,n):
  v=I(1);a=self
  while n:
   if n&1:v=v*a
   a=a*a;n//=2
  return v
 def save(self):return {'lower_numerator':str(self.a),'upper_numerator':str(self.b),'denominator_power2':B}
class C:
 def __init__(self,a,b=0):self.x=a if isinstance(a,I) else I(a);self.y=b if isinstance(b,I) else I(b)
 def __add__(self,o):
  o=o if isinstance(o,C) else C(o);return C(self.x+o.x,self.y+o.y)
 __radd__=__add__
 def __neg__(self):return C(-self.x,-self.y)
 def __sub__(self,o):return self+(-o if isinstance(o,C) else -C(o))
 def __rsub__(self,o):return C(o)+(-self)
 def __mul__(self,o):
  o=o if isinstance(o,C) else C(o)
  return C(self.x*o.x-self.y*o.y,self.x*o.y+self.y*o.x)
 __rmul__=__mul__
 def __pow__(self,n):
  v=C(1);a=self
  while n:
   if n&1:v=v*a
   a=a*a;n//=2
  return v
 def norm2lower(self):
  def dist(a):return min(abs(a.a),abs(a.b)) if a.a*a.b>0 else 0
  return F(dist(self.x)**2+dist(self.y)**2,D*D)
 def save(self):return {'real':self.x.save(),'imag':self.y.save()}
def exp_pos(x,M=48):
 term=F(1);lo=term
 for k in range(1,M+1):term*=x/k;lo+=term
 first=term*x/(M+1);ratio=x/(M+2)
 return I(lo,lo+first/(1-ratio))
def trig(x,sine,M=48):
 start=1 if sine else 0
 lo=sum(((-1)**k*x**(2*k+start)/math.factorial(2*k+start) for k in range(M)),F())
 nxt=(-1)**M*x**(2*M+start)/math.factorial(2*M+start)
 return I(min(lo,lo+nxt),max(lo,lo+nxt))
def char_N_part(p,z,n):
 h=C(1)+ (z-1)*C(p)
 return (h**(n-1))*(h+(1-z)*C(n*(I(1)-p)))
rows=[];checks=0
for m in [4,16,64]:
 n=m**3;r=F(1,m);z=C(trig(r,False),trig(r,True));er=exp_pos(r).inv();a=I(1)-er
 assert 0<a.a<=a.b<D;checks+=1
 kr=char_N_part(I(r),z,n);ha=char_N_part(a,z,n)
 diff=(kr-ha)*C(F(1,m*m))
 atom=(n+1)*(I(1-r)**n-er**n)
 nonlocaldiff=diff-C(atom/(m*m))
 assert diff.norm2lower()>F(1,100);checks+=1
 assert nonlocaldiff.norm2lower()>F(1,100);checks+=1
 # Symbolic binomial-coordinate score identity after cancellation.
 for k in [0,1,n//2,n]:
  assert F(n-k)-F(k)*(1-r)/r==n-F(k)/r;checks+=1
 rows.append({'m':m,'n':n,'r':str(r),'N_characteristic_div_m2':diff.save(),'Nnl_characteristic_div_m2':nonlocaldiff.save(),'certified_normalized_TV_lower':'1/10','atom_interval':atom.save(),'display_N_abs_lower':math.sqrt(float(diff.norm2lower())),'display_Nnl_abs_lower':math.sqrt(float(nonlocaldiff.norm2lower()))})
reg=P/'original_face_tv_obstruction_registration_20261007.json'
out={'status':'PASS','checks':checks,'rows':rows,'binary_fraction_bits':B,'method':'Every operation outward dyadic rounded; trig alternating rational brackets; positive exp rational tail; binary integer powers','registration_sha256':hashlib.sha256(reg.read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'scope':'Original face-marginal TV lower only, not weak endpoint or actual gated counterexample'}
(P/'original_face_tv_obstruction_results_20261007.json').write_text(json.dumps(out,indent=2))
print(json.dumps({'status':out['status'],'checks':checks,'lower_bounds':[(x['n'],x['display_N_abs_lower'],x['display_Nnl_abs_lower']) for x in rows]}))
