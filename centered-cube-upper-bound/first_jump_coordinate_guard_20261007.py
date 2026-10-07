"""Exact coefficients in the raw-kernel coordinate concentration proof.
Not an implementation of the spectral kernel; no full FIRST qualification.
"""
from pathlib import Path
from fractions import Fraction as F
from math import exp,log1p
import json,hashlib
rows=[]
for n in (512,4096,32768):
 sigma=F(2,n)
 assert sigma>F(1,n) and sigma<=F(1,128)
 for t in (F(0),sigma/2,sigma):
  phi_lower=F(1,64)
  # Lower resolvent bound divided by phi; worst phi at the known lower endpoint.
  ratio=(1-t/phi_lower)/(1-t)
  assert ratio>=F(1,2)
  rows.append(dict(n=n,sigma=str(sigma),t=str(t),resolvent_ratio_lower=str(ratio)))
 b=1-sigma
 assert 16*pow(b.numerator,n)>=pow(b.denominator,n)
 assert F(1,32)>F(1,n)
 rows.append(dict(n=n,coordinate_ratio_certified_lower='1/32',uniform_coordinate_bound=str(F(1,n)),illustrative_sharper_ratio=.5*exp(n*log1p(-float(sigma)))))
p=Path(__file__);out=p.with_name(p.stem+'_results.json');assert not out.exists()
out.write_text(json.dumps(dict(status='PASS_EXACT_EFFECTIVE_SOFTNESS_KERNEL_GUARD',rows=rows,scope='Exact coefficient checks for analytical raw-kernel lower bound at sigma=2/n. Kernel identities and clipping proved separately. No full FIRST/GOOD/CP/GP input.',script_sha256=hashlib.sha256(p.read_bytes()).hexdigest()),indent=2)+'\n')
print('Three rounds: 9 resolvent coefficient and 3 coordinate concentration guards passed.')
