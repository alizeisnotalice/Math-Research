"""Minimal exact gates for the endpoint true-operator threshold; no grid."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math

root=Path(__file__).resolve().parent
L0=2**52
d=F(1,2**24); eps=F(1,2**20); sigma=32*d*d
h0=eps/(L0*d)
nmin=F(L0**2,2); tmin=F(L0,2)
rmin=1+d/4; rmax=1+3*d/4
qbridge_min=1-3*d/4-h0; q0max=1-d/4+h0
Rbridge_min=rmin-2*h0
rho_min=tmin*qbridge_min
zeta_t_upper=F(2,L0*nmin*nmin)
M=math.factorial(8)*(24*d**3/F(L0**6)+4*d/F(L0**7))
V=48*math.factorial(8)*d**3/F(L0**6)
cap=F(24*2**44)/(L0**3*d**2)
tau=F(2**74*math.factorial(64),2**384)
scalar=F(1,2**103)+F(math.factorial(16),2**144)+F(2**27*math.factorial(16),2**128)
checks={
 'parameter_sigma':sigma==F(1,2**43)<F(1,2),
 'd_coordinate_domain':d<=F(1,16),
 'h_over_d':h0/d==F(2**28,L0)<=F(1,8),
 'zeta_t_over_h':zeta_t_upper/h0<=F(1,2),
 'matching_shift_below_half_d':h0+zeta_t_upper<d/2,
 't_at_least_16':tmin>=16,
 'all_bridge_q_above_half':qbridge_min>F(1,2),
 'baseline_q0_below_one':q0max<1,
 'Rbridge_above_one':Rbridge_min>1,
 'target_r_below_two':rmax<2,
 'rho_above_4r':rho_min>4*rmax,
 'g_lower_positive':d-h0>0,
 'g_upper_at_most_d':d-h0+zeta_t_upper<=d,
 'baseline_J_delta':d/rho_min<=F(1,8),
 'ell_pinned_at_least_threequarters':1-d/(2*Rbridge_min)>=F(3,4),
 'baseline_abs_mean':M<=F(1,32),
 'baseline_variance':V<=F(9,8192),
 'baseline_Cheb_failure':F(9,8192)/F(3,32)**2<=F(1,8),
 'nd_at_least_eight':nmin*d>=8,
 'Kw_Cheb_failure':F(8)/(9*nmin*d)<=F(1,8),
 'n_at_least_seventeen':nmin>=17,
 'target_cap_face_fraction':F(3,4)**16<=F(1,8),
 'target_cap_coordinate_gap':1/rmax-F(2,L0)>F(1,2),
 'target_mass_fraction':(1-F(3,8))*F(1,4)>=F(1,16),
 'bridge_log_bound':1280*eps==F(5,4096)<F(1,8),
 'full_inverse_product_log':F(1,8)+F(5,4096)<F(1,4),
 'Ld_at_least_one':L0*d>=1,
 'source_cap_rate':12*d<=F(1,64),
 'source_cap_row_loss':cap<F(1,36)<F(3,64),
 'tau_Taylor64_base_exponent':4*d*d*L0==64,
 'tau_Taylor64_base_prefactor':2**20*d*(2**78)==2**74,
 'tau_condition':tau<F(3,128),
 'energy_deficit':F(3,128)*F(1,16)/32==F(3,65536)>=F(1,65536),
 'budget_positivity_gate':L0-52-F(4,L0)>8,
 'x_larger_than_e_coarse':F(L0,2)>3,
 'delta_and_scalar_margin':d/65536==F(1,2**40) and (d/65536)/2-2*sigma==F(1,2**42),
 'scalar_Taylor16':scalar<F(1,2**42),
}
assert all(checks.values()),[k for k,v in checks.items() if not v]
source=Path('/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round045/soft_inverse_column/PROOF.md')
sourcehash=hashlib.sha256(source.read_bytes()).hexdigest()
assert sourcehash=='add435fdf2f55993eba9064bd50c493e4c42769c1c48463329c4015be2e66ab9'
receipt={
 'status':'PASS: all exact endpoint-operator gates',
 'threshold':'log(n) >= 2^52',
 'parameters':{'d':str(d),'epsilon':str(eps),'sigma':str(sigma)},
 'checks':checks,
 'tail_proof':'baseline L^-6/L^-7; cap L^-3; tau L^-125/2; scalar L^-2/L^-16/L^-31/2, plus uniform coordinate guards',
 'frozen_mechanism_source':{'path':str(source),'sha256':sourcehash},
 'source_sha256':{name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in ('THRESHOLD_PROOF.md','verify_threshold.py')},
 'no_production_search':True,'no_MATLAB':True,'no_author_source_changes':True,
}
(root/'verification.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'status':receipt['status'],'checks':len(checks),
                  'proof_sha256':receipt['source_sha256']['THRESHOLD_PROOF.md'],
                  'mechanism_source_unchanged':sourcehash,
                  'operator_threshold':'L>=2^52, subject to independent threshold audit'},ensure_ascii=False,indent=2))
