"""Rational interval certificate for the Phi_2 Hessian sign reversal.
Only integer/Fraction arithmetic is used; no floating-point transcendental calls.
"""
from fractions import Fraction as F
from math import isqrt
from functools import lru_cache
import json,sys
D=10**32

def require(condition,message):
 if not condition:raise ValueError(message)

def floor_grid(x):return F((x.numerator*D)//x.denominator,D)
def ceil_grid(x):return F(-((-x.numerator*D)//x.denominator),D)
class I:
 def __init__(self,lo,hi=None):
  self.lo=F(lo);self.hi=F(lo if hi is None else hi);require(self.lo<=self.hi,"reversed interval")
 @staticmethod
 def rounded(lo,hi):return I(floor_grid(lo),ceil_grid(hi))
 def __add__(self,o):
  o=iv(o);return I.rounded(self.lo+o.lo,self.hi+o.hi)
 __radd__=__add__
 def __neg__(self):return I(-self.hi,-self.lo)
 def __sub__(self,o):return self+-iv(o)
 def __rsub__(self,o):return iv(o)+-self
 def __mul__(self,o):
  o=iv(o);p=[x*y for x in (self.lo,self.hi) for y in (o.lo,o.hi)];return I.rounded(min(p),max(p))
 __rmul__=__mul__
 def __truediv__(self,o):
  o=iv(o);require(o.lo>0 or o.hi<0,"zero-straddling divisor")
  return self*I.rounded(1/o.hi,1/o.lo)
 def __rtruediv__(self,o):return iv(o)/self
 def sq(self):return I.rounded(0 if self.lo<=0<=self.hi else min(self.lo**2,self.hi**2),max(self.lo**2,self.hi**2))
 def export(self):return {'lower':str(self.lo),'upper':str(self.hi)}
def iv(x):return x if isinstance(x,I) else I(x)
def expi(x,terms=40):
 x=iv(x);require(0<=x.lo<=x.hi<=1,"exp domain must be [0,1]")
 term=I(1);total=I(1)
 for j in range(1,terms+1):term=term*x/j;total+=term
 nextterm=term*x/(terms+1)
 tail=nextterm.hi/(1-x.hi/F(terms+2))
 return I.rounded(total.lo,total.hi+tail)
def logi(x,terms=40):
 x=iv(x);require(1<=x.lo<=x.hi<=3,"log domain must be [1,3]")
 z=(x-1)/(x+1);z2=z.sq();term=z;total=I(0)
 for j in range(terms):
  total+=term/F(2*j+1);term*=z2
 tail=2*term.hi/(F(2*terms+1)*(1-z2.hi))
 return I.rounded(2*total.lo,2*total.hi+tail)
def sqrti(x):
 x=iv(x);require(x.lo>=0,"negative sqrt input")
 def lower(y):return F(isqrt((y.numerator*D*D)//y.denominator),D)
 lo=lower(x.lo);h=lower(x.hi);hi=h if h*h==x.hi else h+F(1,D)
 return I(lo,hi)
ELL=logi(I(2))
@lru_cache(None)
def endpoint(x):
 y=iv(x)*ELL;a=2+y;g=expi(y)/sqrti(a);lg=logi(a)
 return [ELL*g,ELL*g*lg,ELL*g*lg.sq()]
# For j=0,1,2 the transformed integrand has positive logarithmic
# derivative log(2)*(1-1/(2a)+j/(a*log(a))) on a in[2,2+log(2)].
# This mathematical monotonicity is what licenses the endpoint sums.
def enclose_integrals(K):
 fs=[endpoint(F(j,K)) for j in range(K+1)]
 return [I.rounded(sum((fs[j][k].lo for j in range(K)),F(0))/K,sum((fs[j][k].hi for j in range(1,K+1)),F(0))/K) for k in range(3)]
def evaluate(K):
 j0,j1,j2=enclose_integrals(K);C=1+j0;aa=logi(C);r1=j1/C;r2=j2/C;U=2*aa+r1
 Ftrue=F(81,8)*(2*U.sq()-F(8,3)*U-r2+r1.sq())
 P=1/C.sq();lam=3*P;m=3-3*aa;v=F(3,2)-F(9,4)*U
 B11=F(1,2)*(m.sq()-3*m*r1+F(9,4)*r2)
 omit_core=-B11+2*v.sq();T=F(9,4)*(4*aa.sq()+4*aa*r1+r2);W=F(1,2)/C;M1=-F(3,2)*U
 omit_soft_core=-T-W*m.sq()+(M1+W*m).sq()/W
 actual=P*Ftrue;omitT=lam*omit_core;omitSoft=lam*omit_soft_core
 return {'K':K,'I0':j0.export(),'I1':j1.export(),'I2':j2.export(),'C':C.export(),'log2':ELL.export(),'actual_P_second':actual.export(),'omit_kernel_prediction':omitT.export(),'omit_soft_prediction':omitSoft.export(),'strict_signs_certified':actual.hi<0<omitT.lo and omitSoft.lo>0}
if __name__=='__main__':
 attempts=[]
 for K in [64,128,256,512,1024,2048]:
  r=evaluate(K);attempts.append({'K':K,'strict_signs_certified':r['strict_signs_certified']})
  if r['strict_signs_certified']:break
 require(r['strict_signs_certified'],"no strict sign certificate found")
 result={'status':'certified using exact rational outward interval arithmetic','model':'Phi_2, cap2, mean target1, Exp(rate theta), theta=1/3','fixed_scale':str(D),'Taylor_terms':40,'quadrature':'monotone left/right Riemann sums after a=2+x log2, x in[0,1]','attempts':attempts,'certificate':r,'scope':'signs only in the stated continuous convex model; no strong maximal C2 bound'}
 s=json.dumps(result,indent=2)
 if len(sys.argv)>1:open(sys.argv[1],'w').write(s+'\n')
 print(s)
