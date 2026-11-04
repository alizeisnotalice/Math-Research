"""Exact atom/positive-certificate checks and high-precision rectangle regressions."""
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import json,platform,time,argparse
import mpmath as mp
import sympy as sp
from witness_capacity import atoms,solve
ROOT=Path(__file__).resolve().parents[1]
def val(c,ell,L):return max(mp.mpf(0),min(L/2,c+ell/2)-max(-L/2,c-ell/2))
def oracle(w,x,y):
 L=1/w;contacts={abs(L-2*x),L+2*x,abs(w-2*y),w+2*y}
 return max(val(x,l,L)*val(y,l,w)/l for l in contacts if l>0)
def threshold(a):
 c=2*a*a-1;d=(1+a)/(a+mp.sqrt(c));u=min(mp.mpf(1),a*(1+a)/(a+mp.sqrt(c))**2);w=a+(1-a)*u
 sa=u*(1+a-a*u)/(w*(1+a));sd=(d*d-u*u)/(2*a*d)
 return max(sa,sd),w,u

def positive_certificates():
 src=json.loads((ROOT/'results/eight_full_search.json').read_text())['rows'];rows=[]
 for ts in ['1/100','1/200','1/1000']:
  t=F(ts);a=1-t;sig=F(39,100);h=1-sig*t;c=2*a*a-1;D=10**24;lo=F(isqrt(c.numerator*D*D//c.denominator),D);hi=lo+F(1,D)
  assert lo*lo<=c<hi*hi
  A=(1/a-a)*h/2;B=(a-(lo+hi)/2)*h/2;centers=[(A,0),(-A,0),(0,A),(0,-A),(B,B),(-B,B),(B,-B),(-B,-B)]
  row=min(src,key=lambda r:abs(r['t']-float(t))+abs(r['sigma']-float(sig)));ell=[F(round(v*10**12),10**12) for v in row['scales']]
  boxes=[(x-l/2,x+l/2,y-l/2,y+l/2) for (x,y),l in zip(centers,ell)];caps=atoms(boxes)
  result=solve(caps,[a*l for l in ell],denominator_limit=10**12,solver_options={'primal_feasibility_tolerance':1e-10,'dual_feasibility_tolerance':1e-10})
  masses={int(s):F(v) for s,v in result['primal_masses'].items()};slack=min(sum(m for mask,m in masses.items() if mask>>i&1)-a*l for i,l in enumerate(ell));eps=hi-lo
  assert slack>4*max(ell)*eps
  rows.append({'t':str(t),'sigma':str(sig),'root_lower':str(lo),'root_upper':str(hi),'boxes':[[str(x) for x in b] for b in boxes],'scales':[str(v) for v in ell],'masses':{str(k):str(v) for k,v in masses.items()},'original_center_error_bound':str(eps),'minimum_capture_slack':str(slack),'scope':'exact positive masses; square-root center uncertainty absorbed by translation bound'})
 (ROOT/'results/eight_positive_certificates.json').write_text(json.dumps(rows,indent=2)+chr(10))
 check_certificates()

def check_certificates():
 rows=json.loads((ROOT/'results/eight_positive_certificates.json').read_text())
 for r in rows:
  t=F(r['t']);a=1-t;c=2*a*a-1;lo=F(r['root_lower']);hi=F(r['root_upper']);assert lo*lo<=c<hi*hi
  boxes=[tuple(map(F,b)) for b in r['boxes']];caps=atoms(boxes);m={int(k):F(v) for k,v in r['masses'].items()};ell=list(map(F,r['scales']));eps=F(r['original_center_error_bound'])
  assert all(0<=v<=caps[k] for k,v in m.items()) and sum(m.values())<=1
  for i,l in enumerate(ell):assert a<l<1/a and sum(v for k,v in m.items() if k>>i&1)-a*l>4*l*eps
 return len(rows)

class Interval:
 def __init__(self,lo,hi=None):self.lo=F(lo);self.hi=F(lo if hi is None else hi)
 @staticmethod
 def cast(v):return v if isinstance(v,Interval) else Interval(v)
 def __add__(self,v):
  v=self.cast(v);return Interval(self.lo+v.lo,self.hi+v.hi)
 __radd__=__add__
 def __neg__(self):return Interval(-self.hi,-self.lo)
 def __sub__(self,v):return self+-self.cast(v)
 def __rsub__(self,v):return self.cast(v)+-self
 def __mul__(self,v):
  v=self.cast(v);values=[self.lo*v.lo,self.lo*v.hi,self.hi*v.lo,self.hi*v.hi];return Interval(min(values),max(values))
 __rmul__=__mul__
 def __truediv__(self,v):
  v=self.cast(v);assert v.lo>0;return self*Interval(1/v.hi,1/v.lo)
 def __pow__(self,n):
  assert n in [0,1,2]
  if n==0:return Interval(1)
  if n==1:return self
  if self.lo>=0:return Interval(self.lo*self.lo,self.hi*self.hi)
  if self.hi<=0:return Interval(self.hi*self.hi,self.lo*self.lo)
  return Interval(0,max(self.lo*self.lo,self.hi*self.hi))
 def sqrt(self):
  assert self.lo>=0;D=10**80
  ll=F(isqrt(self.lo.numerator*D*D//self.lo.denominator),D);hh=F(isqrt(self.hi.numerator*D*D//self.hi.denominator)+1,D)
  assert ll*ll<=self.lo and hh*hh>=self.hi
  return Interval(ll,hh)
 def encode(self):return [str(self.lo),str(self.hi)]

def l_corner_intervals(write=True):
 rows=[]
 for ts in ['1/1000','1/10000','1/100000']:
  tt=F(ts);a=Interval(1-tt);t=Interval(tt);sig=F(7,20);q=(Interval(4)-Interval(7).sqrt())/6;p=2*q-q*q
  discr=a*a-F(3,4)*(1-a*a);r0=(1-a*a)/(a+discr.sqrt());r=r0*(1-sig*t);K=F(1,4)+4*p;ell=(1+K*r*r).sqrt()-r;ell2=ell+2*(1-q)*r
  low=ell*(ell-a)-p*r*r;radial=r-(1-a*a)/(2*a)
  if tt==F(1,1000):assert low.hi<0
  else:
   assert low.lo>0 and radial.lo>0 and (ell-a).lo>0 and (ell-3*r/2).lo>0
   assert (1-a*ell).lo>0 and (1-a*ell2).lo>0
  rows.append({'t':ts,'sigma':'7/20','sqrt_decimal_grid':80,'q':q.encode(),'r':r.encode(),'ell':ell.encode(),'ell2':ell2.encode(),'low_slack':low.encode(),'radial_excess_over_rectangle_bound':radial.encode(),'strict_counterexample':low.lo>0 and radial.lo>0})
 if write:(ROOT/'results/l_corner_exact_certificates.json').write_text(json.dumps(rows,indent=2)+chr(10))
 else:
  saved=json.loads((ROOT/'results/l_corner_exact_certificates.json').read_text());assert saved==rows
 return sum(row['strict_counterexample'] for row in rows)

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--check-only',action='store_true');args=parser.parse_args()
 if args.check_only:print({'exact_positive_certificates_rechecked':check_certificates(),'exact_L_counterexamples_rechecked':l_corner_intervals(False)});raise SystemExit
 started=time.monotonic();count=0
 for t in [F(1,10),F(1,20),F(1,100)]:
  for lx in [F(j,4) for j in range(-4,5)]:
   for ly in [F(j,4) for j in range(-4,5)]:
    x=1+lx*t;y=1+ly*t;boxes=[(t-x/2,t+x/2,-x/2,x/2),(-t-x/2,-t+x/2,-x/2,x/2),(-y/2,y/2,t-y/2,t+y/2),(-y/2,y/2,-t-y/2,-t+y/2)];caps=atoms(boxes)
    layers=[sum(v for mask,v in caps.items() if mask.bit_count()>=j) for j in [2,3,4]]
    assert layers==[x*y,x*x+y*y-x*y-4*t*t,(x-2*t)*(y-2*t)];count+=1
 mp.mp.dps=80;ac=mp.findroot(lambda a:4*a**3-a-1,(mp.mpf('.75'),mp.mpf('.8')));samples=[]
 for a in [mp.mpf('.71'),mp.mpf('.72'),mp.mpf('.75'),ac-mp.mpf('1e-8'),ac,ac+mp.mpf('1e-8'),mp.mpf('.8'),mp.mpf('.9'),mp.mpf('.95'),mp.mpf('.99'),mp.mpf('.999'),mp.mpf('.999999')]:
  s,w,u=threshold(a);t=1-a;A0=(1/a-a)/2;B0=(a-mp.sqrt(2*a*a-1))/2;h=1-s*t
  mm=[oracle(w,A0*h,0),oracle(w,0,A0*h),oracle(w,B0*h,B0*h)];assert abs(min(mm)-a)<mp.mpf('1e-65')
  ep=mp.mpf('1e-4');hm=1-(s+ep)*t;hp=1-(s-ep)*t
  assert min(oracle(w,A0*hm,0),oracle(w,0,A0*hm),oracle(w,B0*hm,B0*hm))>a
  assert min(oracle(w,A0*hp,0),oracle(w,0,A0*hp),oracle(w,B0*hp,B0*hp))<a
  samples.append({'alpha':mp.nstr(a,70),'optimal_short_side':mp.nstr(w,70),'u':mp.nstr(u,70),'sigma':mp.nstr(s,70),'oracle_min_minus_alpha':mp.nstr(min(mm)-a,8)})
 # Independent rational formal recursion, cross-checked with actual Wolfram output.
 z=sp.Symbol('z');w=sp.Integer(1)
 for n in range(1,15):
  expr=(2*z-z*z)*w*w-4*((1-z)**3-1)*w+4*(1-z)**4
  cn=sp.expand(expr).coeff(z,n)/4;w+=cn*z**n
 a=1-z;expr=(w-a)*(1-a*w)/(w*(1-a*a)*z);series=sp.series(expr,z,0,13).removeO();coefs=[str(sp.expand(series).coeff(z,n)) for n in range(13)]
 cas=json.loads((ROOT/'results/eight_witness_wolfram.json').read_text());assert coefs==cas['sigma_coefficients']
 # Uniform explicit Cauchy remainder at t<1/8, checked numerically at selected orders.
 tails=[]
 for ts in ['.001','.01','.05','.1']:
  tt=mp.mpf(ts);ss=threshold(1-tt)[0]
  for N in [0,1,2,4,8,12]:
   approx=sum(mp.mpf(str(sp.Rational(v).p))/int(sp.Rational(v).q)*tt**n for n,v in enumerate(coefs[:N+1]));bd=4*(8*tt)**(N+1)/(1-8*tt);assert abs(ss-approx)<=bd;tails.append({'t':ts,'N':N,'absolute_error':mp.nstr(abs(ss-approx),25),'cauchy_bound':mp.nstr(bd,25)})
 positive_certificates();lcount=l_corner_intervals()
 out={'python':platform.python_version(),'sympy':sp.__version__,'mpmath':mp.__version__,'precision':80,'exact_layer_cases':count,'transition_alpha':mp.nstr(ac,70),'sigma_coefficients':coefs,'rectangle_samples':samples,'taylor_tail_samples':tails,'positive_certificates':3,'L_counterexamples':lcount,'seconds':time.monotonic()-started,'scope':'exact finite layer and positive certificates; regression of exact rectangle model, not a global eight-scale upper bound'}
 (ROOT/'results/eight_witness_python.json').write_text(json.dumps(out,indent=2)+chr(10));print({k:v for k,v in out.items() if k not in ['rectangle_samples','taylor_tail_samples']})
