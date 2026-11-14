"""Standard-library rational certificate: Gamma information-spectrum rate violates actual Lp bound at t=1/10."""
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal,localcontext,ROUND_FLOOR,ROUND_CEILING
import json,sys,hashlib,math
HERE=Path(__file__).resolve().parent
source=Path(sys.argv[1]) if len(sys.argv)>1 else HERE.parent/'r708'/'constant_certificate.json'
d=json.loads(source.read_text());ml,mh=map(F,d['mu_interval'])
def expi(x):
 assert 0<x<2
 term=F(1);s=term
 for k in range(1,101):term*=x/k;s+=term
 nxt=term*x/101;return s,s+nxt/(1-F(2,102))
assert (ml-1)*expi(ml)[1]<1<(mh-1)*expi(mh)[0]
def logi(x):
 e=0
 while x>=2:x/=2;e+=1
 while x<1:x*=2;e-=1
 def unit(y):
  z=(y-1)/(y+1);s=sum((2*z**(2*j+1)/F(2*j+1) for j in range(80)),F(0));return s,s+2*z**161/(161*(1-z*z))
 a,b=unit(x);c,d=unit(F(2));return(a+e*c,b+e*d) if e>=0 else(a+e*d,b+e*c)
P=lambda x:10*x**11-11*x**10-1
lo=F(1);hi=F(2)
for j in range(220):
 m=(lo+hi)/2
 if P(m)<0:lo=m
 else:hi=m
assert P(lo)<0<P(hi)
ll,lh=logi(F(11,10));cl=logi(lo)[0];ch=logi(hi)[1];gaplo=F(11,10)*(ml*ll-ch);gaphi=F(11,10)*(mh*lh-cl);assert gaplo>0

def dec(x,up=False):
 with localcontext() as c:
  c.prec=55;c.rounding=ROUND_CEILING if up else ROUND_FLOOR;return str(Decimal(x.numerator)/Decimal(x.denominator))
out=dict(status='complete',t='1/10',q=11,polynomial='10*c^11-11*c^10-1',c_bracket=[str(lo),str(hi)],mu_bracket=d['mu_interval'],mu_source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),exp_terms=100,log_terms=80,rate_gap_lower=dec(gaplo),rate_gap_upper=dec(gaphi,True),scope='Strict asymptotic information-moment violation by the specified Gamma rank construction; no general endpoint closure')
(HERE/'certify_lp.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k] for k in ['status','rate_gap_lower','rate_gap_upper']}))

# Exact sufficient conditions for the three-point all-positive-moment comparison.
l2,h2=logi(F(2));checks=[]
for k in range(12,35,2):
 n=2**k;dl=F(k,2)*l2+4;dh=F(k,2)*h2+4
 for v in [F(1,10),F(1),F(4),None]:
  if v is None:al=ah=F(2**(k//2),2)
  else:al=dl+v*k*l2;ah=dh+v*k*h2
  flags=[ah*ah<=ml*n/2,ah<=ml*n/4,ml*n-dh-ah>h2,ml*n>h2,ml*n>2*h2,dl+al>3*h2]
  assert all(flags)
  checks.append(dict(n=n,variant=str(v) if v is not None else 'sqrt',quadratic_condition=flags[0],linear_condition=flags[1],positive_low_threshold=flags[2],rare_probability_below_half=flags[3],rare_probability_below_quarter=flags[4],mode_gap_log8=flags[5]))
(HERE/'three_point_conditions.json').write_text(json.dumps(dict(status='complete',arithmetic='exact rational logarithm and mu enclosures',scope='sufficient all-positive-t comparison conditions for abstract three-step spectra only',cases=len(checks),results=checks),indent=2)+'\n')
print(json.dumps(dict(three_point_exact_conditions='pass',cases=len(checks))))
