"""Exact rational envelopes for U02 tails, not a quadrature error certificate.

Domain a>=47/100, delta=1/1000, L=32, 3/25<=split<=1/5,
image M=8, spectral N=40. Uses only pi>3, e>2, erfc(z)<=exp(-z*z).
"""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json

rows=[]
for n in (3,4):
    b=Q(n)*Q(47,100)**2/(4*Q(1,1000))
    exponent=b.numerator//b.denominator
    # delta*(4*pi*delta)^(-n/2) <= delta^(1-n/2).
    prefactor=32 if n==3 else 1000
    low_green=Q(prefactor,2**exponent)
    low_mass=Q(2*n,1000*2**250)
    high_green=Q(2**n,2**(72*n))
    high_mass=Q(4**n,2**(72*n))
    # Power perturbation: |x^n-y^n|<=n*2^(n-1)|x-y| for |x|,|y|<=2.
    image_green=Q(n*2**(n-1)*8,2**1445)
    image_mass=Q(n*2**(n-1)*4,2**451)
    spectral_green=Q(n*2**n,2**1771)
    spectral_mass=spectral_green
    total_green=low_green+high_green+image_green+spectral_green
    total_mass=low_mass+high_mass+image_mass+spectral_mass
    assert total_green<Q(1,10**48)
    assert total_mass<Q(1,10**62)
    rows.append({'n':n,'low_time_green_exponent_floor':exponent,
                 'low_time_green_upper':str(low_green),
                 'low_time_green_upper_decimal':float(low_green),
                 'all_analytic_green_omissions_lt':'1e-48',
                 'all_analytic_mass_omissions_lt':'1e-62',
                 'scope':'time/series omissions only; finite-interval quadrature error NOT included'})
print(json.dumps({'rows':rows,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2))
