#!/usr/bin/env python3
"""Independent compressed one-tuple audit; no author helper or grid replay."""
import hashlib
import json
import platform
from fractions import Fraction as F
from pathlib import Path

import sympy as s

HERE=Path(__file__).resolve().parent
SOURCE=Path('/Users/zhengzhihao/Desktop/T/cube-r040-cover-selector/cube-dim-7-bounds/upper/round059/wide_source_band/PROOF.md')
AUTHOR_CERT=SOURCE.with_name('certificate.json')
LOCAL=Path('/Users/zhengzhihao/Desktop/T/cube-r030-roots-audit/cube-dim-7-bounds/audit/round056/one_coordinate/LEMMA.md')

def alternating_log_ratio(lower,upper):
    """log(upper/lower), 40-even lower / 41-odd upper; not author's atanh recipe."""
    u=(upper-lower)/lower
    assert 0<u<1
    even=sum(((-1)**(k+1)*u**k/k for k in range(1,41)),F())
    odd=even+u**41/41
    return even,odd

d=F(391,1905); aa=F(283,2032); bb=F(63,64)
e=F(1,2**20); hc=F(1,2**48); eret=F(1,4096); ea=F(1,65536)
chi=F(147,10000); sigma=F(3,100)
vmin=aa*d;vmax=bb*d;qmin=1-vmax;qmax=1-vmin
base=qmin+d;rlo=base+e;rhi=qmax+d-e
p0=1-e-1280*hc;ulo=p0*(1-e)/(1+e);uhi=(1+e)/(p0*(1-e))
degree=2*uhi/qmin;bracket=(1-eret)*2*ulo/qmax-2*e
mass=F(159477,10**6); Gamma=mass*bracket/degree
gamma_budget=(1-e)*Gamma;gamma=F(1637,12500)
Qbar=F(26979,250000)
target_low,target_high=alternating_log_ratio(rlo,rhi)
pin_low,pin_high=alternating_log_ratio(2*base-d,2*base)
mean_rate=base*Qbar**2/(qmin*(1-e))
var_rate=F(2,3)*base*Qbar**3/(qmin*(1-e))
eta=vmin/(2*qmax);csrc=(qmax-qmin)/qmin
penalty=sigma/(1-sigma);craw=gamma-penalty

checks={
    'd_local_domain':0<d<F(1,4)<F(1,3),
    'source_band_valid':0<vmin<vmax<1 and F(1,2)<qmin<qmax<1,
    'target_and_active_base':1<base<rlo<rhi<F(5,4),
    'inset_leaves_nonempty_band':2*e<qmax-qmin,
    'positive_product_and_row_budgets':0<p0<1 and degree>2 and bracket>0,
    'target_mass_new_lower':(1-ea)*target_low>mass,
    'active_pinlog_upper':pin_high<Qbar<F(1,9),
    'active_face_mean_rate':mean_rate<F(293,20000),
    'active_face_variance_rate':var_rate<F(1,900),
    'Chebyshev_mean_gap':chi-F(293,20000)==F(1,20000),
    'Chebyshev_pin_gap_arithmetic':F(1,9)/(F(40000,9))==F(1,40000),
    'Chebyshev_tail_coefficient':1/(F(900)*F(1,40000)**2)==F(40000**2,900),
    'wide_band_Csrc_is_not_d':csrc>d,
    'cap_exponent_exact':eta-chi==F(499871,37603070000)>F(1,100000),
    'tau_exponent_exact':sigma/2-chi*(1-sigma)==F(741,10**6)>0,
    'chosen_gamma_paid':0<gamma<gamma_budget<Gamma<1,
    'complete_scalar_raw':craw==F(121289,1212500),
    'strict_one_tenth_margin':craw-F(1,10)==F(39,1212500)>0,
}
assert all(checks.values())

# Symbolic identities check new bands/counts without choosing n or replaying geometry.
r,h,D,eps,q=s.symbols('r h D eps q',positive=True)
ap=(1-eps)*(r-D);am=(1-eps)*(r-D-2*h)
symbolic={
    'count_cross_endpoint_identity':s.simplify(ap/r-am/(r-2*h)-2*h*(1-eps)*D/(r*(r-2*h)))==0,
    'pin_length_monotone':s.simplify(s.diff((1+(1-eps)*q/(q+D))/2,q)-(1-eps)*D/(2*(q+D)**2))==0,
    'active_radius_lower_difference':s.expand((qmin+e-h)+d-h-base)==e-2*h,
}
assert all(symbolic.values())
values={'vmin':vmin,'vmax':vmax,'qmin':qmin,'qmax':qmax,'Rbase':base,
    'rlo':rlo,'rhi':rhi,'P0':p0,'Ulo':ulo,'Uhi':uhi,'M':degree,'b':bracket,
    'mL':mass,'Gamma_L':Gamma,'gamma_budget':gamma_budget,'gamma':gamma,
    'Qbar':Qbar,'mean_rate_upper':mean_rate,'variance_rate_upper':var_rate,
    'eta':eta,'Csrc':csrc,'cap_gap':eta-chi,'tau_gap':sigma/2-chi*(1-sigma),
    'c_raw':craw,'one_tenth_slack':craw-F(1,10),
    'target_log_lower_40':target_low,'target_log_upper_41':target_high,
    'mass_lower_slack':(1-ea)*target_low-mass,'pinlog_upper_41':pin_high,
    'gamma_budget_slack':gamma_budget-gamma}
paths=[SOURCE,AUTHOR_CERT,LOCAL,Path(__file__).resolve()]
out={'status':'PASS_NEW_FIXED_TUPLE_AND_BAND_CHECKS','python':platform.python_version(),
    'sympy':s.__version__,'arithmetic':'independent Fraction and SymPy',
    'log_recipe':'40/41-term alternating log1p enclosure, distinct from author 64-term positive atanh helper',
    'rational_checks':checks,'symbolic_checks':symbolic,
    'values':{k:str(v) for k,v in values.items()},
    'sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
    'scope':'One new fixed tuple/band and key margins only. No author helper imported, no grid or accepted baseline rerun, no MATLAB, no explicit n threshold certified. Full operator activation requires the separate analytic audit.'}
(HERE/'new_margins_receipt.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':out['status'],'rational_checks':checks,'symbolic_checks':symbolic,
    'one_tenth_slack':str(craw-F(1,10))},indent=2))
