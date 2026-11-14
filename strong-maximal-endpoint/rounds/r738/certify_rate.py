from fractions import Fraction as F
from pathlib import Path
from decimal import Decimal,localcontext,ROUND_FLOOR,ROUND_CEILING
import json
HERE=Path(__file__).resolve().parent
P=lambda t:-5*t**19-18*t**15+45*t**9+30*t**4+36
lo=F(1);hi=F(2)
for _ in range(220):
 mid=(lo+hi)/2
 if P(mid)>0:lo=mid
 else:hi=mid
assert P(lo)>0>P(hi)
num=lambda t:t**15+t**10+2*t**6
den=lambda t:t**15+3
al=num(lo)/den(hi);ah=num(hi)/den(lo)
def logi(x):
 e=0
 while x>=2:x/=2;e+=1
 while x<1:x*=2;e-=1
 def unit(y):
  z=(y-1)/(y+1);s=sum((2*z**(2*j+1)/F(2*j+1) for j in range(80)),F(0));return s,s+2*z**161/(161*(1-z*z))
 a,b=unit(x);c,d=unit(F(2));return (a+e*c,b+e*d) if e>=0 else(a+e*d,b+e*c)
cl=logi(al)[0];ch=logi(ah)[1];assert ch<F(414,1000)
def dec(x,up=False):
 with localcontext() as c:
  c.prec=60;c.rounding=ROUND_CEILING if up else ROUND_FLOOR;return str(Decimal(x.numerator)/Decimal(x.denominator))
out=dict(status='complete',method='exact rational polynomial bisection plus 80-term atanh log enclosure',scope='three-cell specified grid witness exponential rate only; not full source halo or general endpoint coefficient',root_bracket=[str(lo),str(hi)],root_signs_exact=True,log_rate_lower=dec(cl),log_rate_upper=dec(ch,True),strict_log_rate_less_than='0.414',general_mu_greater_than=1,log_terms=80)
(HERE/'certify_rate.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k] for k in ['status','log_rate_lower','log_rate_upper']}))
