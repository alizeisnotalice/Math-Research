"""Actual parameter arithmetic for a conditional flat-variance mixture budget.
The finite-n flatness assumption D is NOT numerically certified here.
"""
from fractions import Fraction as F
from decimal import Decimal as D,localcontext
from pathlib import Path
import json
out=Path(__file__).resolve().parent
rows=[]
for n,d in ((10**32,47),(10**42,63),(10**60,93)):
 K=2**d-1;p=F(99,100);v=9*p*p*n/(29*K*K);h=p*F(58,9)/(K*F(10,9));eps=1/v
 J=(K*v).numerator//(K*v).denominator+2;qbar=n*h*h
 for A in (F(51,50),F(2)):
  mean_bias=A*h*qbar;t=eps-mean_bias;assert t>0 and A*h<1
  varcap=(1+A*h)*h*h*qbar;E=t*t/(2*(varcap+h*h*t/3))
  lgJ=(J-1).bit_length();binary_tail_exponent=1+2*lgJ-E.numerator//E.denominator
  cubic=n*(A*h)**3/(3*(1-A*h))
  for flatD in (F(3,100),F(1)):
   delta=A*A*(flatD+eps)/2+cubic
   with localcontext()as ctx:
    ctx.prec=50
    dec=lambda x:D(x.numerator)/D(x.denominator)
    rows.append({'n':str(n),'d':d,'v_decimal':str(dec(v)),'A':str(A),'h':str(h),'epsilon':str(eps),'qbar_unconditional':str(qbar),'anchor_mean_bias':str(mean_bias),'Bernstein_exponent':str(E),'log2_bad_output_mass_upper_integer':binary_tail_exponent,'cubic_remainder':str(cubic),'flatness_D_ASSUMED':str(flatD),'good_delta':str(delta),'positive_amplitude_multiplier_decimal':str((2*dec(delta)).exp()),'signed_amplitude_multiplier_decimal':str(2*(2*dec(delta)).exp())})
r={'status':'Exact Fraction budget; Decimal displays. Concentration unconditional from h. Flat-variance D is an explicit conditional assumption at finite n; only eventual bounded D is taken from accepted source theorem.','rows':rows}
(out/'RESULTS.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps({'rows':len(rows),'first':{k:rows[0][k] for k in ['n','A','flatness_D_ASSUMED','log2_bad_output_mass_upper_integer','positive_amplitude_multiplier_decimal','signed_amplitude_multiplier_decimal']}},indent=2))
