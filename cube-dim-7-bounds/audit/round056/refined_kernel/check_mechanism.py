#!/usr/bin/env python3
"""New local moment/count identities and exact rates, not the full finite guard suite."""
import hashlib
import json
import platform
from fractions import Fraction as F
from pathlib import Path

import sympy as sy

HERE=Path(__file__).resolve().parent
SOURCE=Path('/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round056/refined_inverse/PROOF.md')
EXPECTED='6986ed9ec48b46f3136a9a56722d0092db4f2942d52e35c9f4de37eb279f399a'
assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()==EXPECTED

R,g,q,D,eta,r,h=sy.symbols('R g q D eta r h',positive=True)
lp=-sy.log(1-g/(2*R))
symbolic={
    'free_first_moment_antiderivative':sy.simplify(sy.diff(R*lp**2,g)-lp/(1-g/(2*R)))==0,
    'free_second_moment_antiderivative':sy.simplify(sy.diff(2*R*lp**3/3,g)-lp**2/(1-g/(2*R)))==0,
    'ellpin_monotone_derivative':sy.simplify(sy.diff((1+(1-eta)*q/(q+D))/2,q)-(1-eta)*D/(2*(q+D)**2))==0,
    'cross_endpoint_count_identity':sy.simplify((1-eta)*(r-D)/r-(1-eta)*(r-D-2*h)/(r-2*h)-2*h*(1-eta)*D/(r*(r-2*h)))==0,
}
assert all(symbolic.values())

e=F(1,1024); p0=F(65275,65536); qmin=F(31,40);qmax=F(37,40)
ulo=p0*(1-e)/(1+e);uhi=(1+e)/(p0*(1-e))
u=F(500,3083); Q=u-u*u/2+u**3/3
fac=F(43,31)/(1-e);mrate=fac*Q**2;vrate=F(2,3)*fac*Q**3
x0=F(93,100)*2**20;chi=F(7,200);sigma=F(7,100)
rational={
    'Qpin_Taylor_upper_below_one_sixth':Q<F(1,6),
    'sourceface_mean_rate':mrate<F(4,125),
    'sourceface_variance_rate':vrate<F(1,300),
    'pin_over_x_budget':Q/x0<F(1,2000),
    'Chebyshev_gap':chi-F(4,125)-F(1,2000)==F(1,400),
    'Chebyshev_coefficient':1/(F(300)*F(1,400)**2)==F(1600,3),
    'perface_tail_below_retention_error':F(1600,3)/x0<e,
    'row_upper_below_21_over_8':2*uhi/qmin<F(21,8),
    'sourcecap_full_n_minus_one_rate':F(74,75)*F(3,74)==F(1,25),
    'sourcecap_theta_gap':F(1,25)-chi==F(1,200),
    'tau_theta_exponent':sigma/2-chi*F(93,100)==F(49,20000),
}
assert all(rational.values())
files=[SOURCE,HERE/'mechanism_review.md',Path(__file__).resolve()]
out={'status':'PASS_LOCAL_MECHANISM_CHECKS','python':platform.python_version(),'sympy':sy.__version__,
    'symbolic_identities':symbolic,'rational_checks':rational,
    'values':{'Qpin_upper':str(Q),'mean_rate_upper':str(mrate),'variance_rate_upper':str(vrate),
        'Ulo':str(ulo),'Uhi':str(uhi),'row_upper':str(2*uhi/qmin),
        'perface_tail_L0':str(F(1600,3)/x0)},
    'sha256':{str(path):hashlib.sha256(path.read_bytes()).hexdigest() for path in files},
    'scope':'Only new local symbolic/count moment identities and rational mechanism rates. Does not rerun the local-J proof, full 29 finite guards, scalar remainder, old grids, production batches or MATLAB. Scoped analytic review and separate full-chain audit are required before any Cn conclusion.'}
(HERE/'local_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':out['status'],'symbolic_identities':symbolic,'rational_checks':rational},indent=2))
