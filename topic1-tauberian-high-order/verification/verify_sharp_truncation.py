"""Exact finite periodic inputs, no enumeration of their N^n cells.
E=[0,1]^n minus (N-1)^n explicitly specified cubes as in prop:lower.
Checks samples only; the uniform error estimate is proved in the paper.
"""
from fractions import Fraction as F
import json
rows=[]
for n in range(2,9):
    for k in (4,6,8):
        u=F(1,2**k); t=u**n; rho=u*(1-u*u); N=2**(3*k)
        h=(2*u/(1+2*u))**n
        hf=(1-F(1,N))**n*(2*rho/(1+2*rho))**n
        finite=1/(1-hf); envelope=1/(1-h)
        margin=t-rho**n
        # A direct rational bound: h-hf <= n*h*(u^2+1/N).
        assert margin>0 and 2*rho<1 and hf<h< F(1,2)
        assert h-hf <= n*h*(u*u+F(1,N))
        assert envelope-finite <= 4*n*(2**n*t)*(u*u+F(1,N))
        rows.append({'n':n,'u':str(u),'N':N,'witness_defect_margin':str(margin),
                     'finite_gain':float(finite-1),
                     'normalized_second':float((finite-1-2**n*t)/(t*u)),
                     'envelope_gap_over_t_u2':float((envelope-finite)/(t*u*u))})
print(json.dumps({'arithmetic':'fractions.Fraction exact; decimals displayed only',
                  'all_checks_pass':True,'samples':rows},ensure_ascii=False,indent=2))
