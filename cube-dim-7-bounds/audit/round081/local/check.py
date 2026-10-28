#!/usr/bin/env python3
"""Bounded independent s=4 one-coordinate check, not a search or proof.

Exact Fraction checks support the elementary constant inequalities.
SciPy quadrature checks the three integral formulas in full, edge and
central-deletion regimes. Floating evidence is separate from proof.md.
"""
from fractions import Fraction as F
from pathlib import Path
import json,math,os,platform
from scipy.integrate import quad

HERE=Path(__file__).resolve().parent
exact={
 'Q_less_half_via_positive_exp_prefix': F(1)+F(1,2)+F(1,8)>F(8,5),
 'F_partial_constants': F(3,5)*4+F(3,5)*4+2*4*F(6,5)+1==F(77,5),
 'log_derivative_bound': (F(77,5)+32)/2/F(2,3)==F(711,20)<100,
 'first_profile_derivative_bound': F(28,5)+F(1,1000)==F(5601,1000)<40,
 'second_profile_derivative_bound': F(14,5)+F(1,3000)<3,
}
assert all(exact.values())
h=.001;epsilon=h/(4*.99);t=64.

def formula(q,w):
    R=4*q-h;a0=(1-epsilon)*q;g=R-a0;rho=t*q
    m=max(abs(w)-a0,0.);v=max(g-m,0.)
    Q=lambda z:-math.log1p(-z/(2*R))
    FF=lambda z:2*R*Q(z)-z
    return (1+(FF(g)+FF(v)-m)/(2*rho),
      R*(Q(g)**2+Q(v)**2)/(2*rho),
      R*(Q(g)**3+Q(v)**3)/(3*rho))

samples=[];max_error=0.;max_scaled_derivative=[0.,0.,0.]
for q in ((1+h)/4,.5,.99):
    R=4*q-h;a0=(1-epsilon)*q;g=R-a0;rho=t*q
    for w in (0.,a0/2,a0,a0+g/2,R,R+h,R+2*h,-R-h):
        expected=formula(q,w)
        lo=max(-rho,w-rho-a0);hi=min(rho,w+rho+a0)
        points=[z for z in (-rho+g,0.,rho-g) if lo<z<hi]
        direct=[]
        for j in range(3):
            def integrand(z):
                ell=min(1.,(R+rho+a0-abs(z))/(2*R))
                phi=-math.log(ell)
                return phi**j/ell/(2*rho)
            value,error=quad(integrand,lo,hi,points=points,epsabs=2e-13,epsrel=2e-13)
            direct.append(value)
        discrepancy=max(abs(a-b) for a,b in zip(expected,direct))
        max_error=max(max_error,discrepancy)
        samples.append({'q':q,'w':w,'regime':'central' if abs(w)>R else ('edge' if abs(w)>a0 else 'full'),
          'formula':list(expected),'direct_quadrature':direct,'max_absolute_error':discrepancy})
        dq=1e-7
        plus,minus=formula(q+dq,w),formula(q-dq,w)
        derivatives=[(math.log(plus[0])-math.log(minus[0]))/(2*dq),
          (plus[1]-minus[1])/(2*dq),(plus[2]-minus[2])/(2*dq)]
        for j,val in enumerate(derivatives):
            max_scaled_derivative[j]=max(max_scaled_derivative[j],rho*abs(val))
assert max_error<1e-10
out={'status':'PASS bounded arithmetic and numerical formula diagnostics; rigorous proof separate',
 'pid':os.getpid(),'python':platform.python_version(),'exact_constant_checks':exact,
 'numerical_scope':'24 samples, direct coordinate quadrature; finite differences are diagnostics only',
 'h':h,'epsilon':epsilon,'t':t,'max_formula_quadrature_error':max_error,
 'max_rho_scaled_finite_difference':[float(v) for v in max_scaled_derivative],
 'samples':samples}
(HERE/'result.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='samples'},indent=2))
