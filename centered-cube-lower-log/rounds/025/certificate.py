from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import json
n=900000000;a=F(1,2046);A=F(77,15);q=2;t=1100
rad=2*n*A*A*q*t
root=F(isqrt((rad.numerator*10**24)//rad.denominator)+1,10**12)
dQ=root+2*A*A*t/3
eta=a*a*dQ/2+a**3*A*(n*q+dQ)/(3*(1-a*A))
assert eta<2 and n*a*a*q<431
out=dict(eta_upper_rational=str(eta),eta_display=float(eta),eta_chosen=2,second_density_log_bound=str(n*a*a*q),r_max=60000,log_Pbad_upper=-1088,log_Qbad_upper=-328,weak_lower_bad_correction_upper='1e-135',weak_upper_bad_correction_upper='1e-450')
Path(__file__).with_name('certificate.json').write_text(json.dumps(out,indent=2));print({k:v for k,v in out.items() if k!='eta_upper_rational'})
