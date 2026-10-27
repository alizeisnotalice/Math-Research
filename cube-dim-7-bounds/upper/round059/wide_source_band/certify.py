#!/usr/bin/env python3
"""One exact fixed-tuple certificate. No grid, MATLAB, or numerical solve."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import platform

HERE = Path(__file__).resolve().parent
KTERMS = 64


def log_ratio_bounds(lo, hi):
    assert F(0) < lo < hi
    z = (hi-lo)/(hi+lo)
    low = 2*sum((z**(2*j+1)/F(2*j+1) for j in range(KTERMS)), F(0))
    high = low + 2*z**(2*KTERMS+1)/((2*KTERMS+1)*(1-z*z))
    return low, high


d, aband, bband = F(391,1905), F(283,2032), F(63,64)
e = F(1,2**20)
eret, ea = F(1,4096), F(1,65536)
chi, sigma = F(147,10000), F(3,100)
delta, hc = e, F(1,2**48)
vmin, vmax = aband*d, bband*d
qmin, qmax = 1-vmax, 1-vmin
Rbase = qmin+d
rlo, rhi = Rbase+delta, qmax+d-delta
P0 = 1-e-1280*hc
Ulo = P0*(1-e)/(1+e)
Uhi = (1+e)/(P0*(1-e))
M = 2*Uhi/qmin
b = (1-eret)*2*Ulo/qmax-2*e
mcert = F(159477,10**6)
Gamma = mcert*b/M
gamma_budget = (1-e)*Gamma
gamma = F(1637,12500)
craw = gamma-sigma/(1-sigma)
massloglow, massloghigh = log_ratio_bounds(rlo,rhi)
pinloglow, pinloghigh = log_ratio_bounds(2*Rbase-d,2*Rbase)
Qcert = F(26979,250000)
mubar, varbar = F(293,20000), F(1,900)
muup = Rbase*Qcert**2/(qmin*(1-e))
varup = F(2,3)*Rbase*Qcert**3/(qmin*(1-e))
eta = vmin/(2*qmax)
Csource = (qmax-qmin)/qmin

assert 0 < d < F(1,4) < F(1,3)
assert 0 < vmin < vmax < 1 and F(1,2) < qmin < qmax < 1
assert 1 < Rbase < rlo < rhi < F(5,4)
assert 2*delta < qmax-qmin
assert 0 < P0 < 1 and M > 2 and b > 0
assert (1-ea)*massloglow > mcert
assert F(1079155,10**7) < pinloglow < pinloghigh < Qcert < F(1,9)
assert muup < mubar < chi < eta
assert varup < varbar
assert chi-mubar == F(1,20000)
assert eta-chi > F(1,100000)
assert sigma/2-chi*(1-sigma) == F(741,10**6) > 0
assert 0 < gamma < gamma_budget < Gamma < 1
assert craw == F(121289,1212500)
assert craw-F(1,10) == F(39,1212500) > 0

source_paths = {
    'local_coordinate_lemma': Path('/Users/zhengzhihao/Desktop/T/cube-r030-roots-audit/cube-dim-7-bounds/audit/round056/one_coordinate/LEMMA.md'),
    'R058_model_certificate': Path('/Users/zhengzhihao/Desktop/T/cube-r030-compute/cube-dim-7-bounds/upper/round058/compute/model_certificate.json'),
    'R051_bridge_reference': Path('/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round051/flexible_budget/PROOF.md'),
}
model = json.loads(source_paths['R058_model_certificate'].read_text())
assert model['input']['d'] == str(d)
assert model['input']['a'] == str(aband)
assert model['input']['b'] == str(bband)

fixed = dict(d=d, a_band=aband, b_band=bband, chi=chi, sigma=sigma,
             inset=delta, h_coefficient=hc, epsilon_J=e, epsilon_N=e,
             epsilon_S=e, epsilon_T=e, epsilon_D=e, epsilon_Ret=eret,
             epsilon_A=ea, baseline_log_tolerance=e)
derived = dict(vmin=vmin,vmax=vmax,qmin=qmin,qmax=qmax,Rbase=Rbase,
               rlo=rlo,rhi=rhi,P0=P0,Ulo=Ulo,Uhi=Uhi,M=M,b=b,
               target_mass_lower=mcert,Gamma_lower=Gamma,
               gamma_budget_lower=gamma_budget,gamma_chosen=gamma,
               pinlog_upper=Qcert,meanrate_upper=muup,meanrate_simple=mubar,
               varrate_upper=varup,varrate_simple=varbar,caprate=eta,
               source_cap_exponent_gap=eta-chi,
               tau_exponent_gap=sigma/2-chi*(1-sigma),
               source_band_logmass_upper=Csource,c_raw_chosen=craw,
               margin_above_one_tenth=craw-F(1,10))
result = {
    'status': 'Fixed-tuple exact scalar/local-moment certificate; operator proof remains proposed pending independent review.',
    'python': platform.python_version(), 'arithmetic': 'fractions.Fraction',
    'log_method': '64 positive atanh-series terms with exact geometric tail upper bound',
    'parameters': {k:str(v) for k,v in fixed.items()},
    'derived_rationals': {k:str(v) for k,v in derived.items()},
    'checks': ['full target mass lower includes epsilon_A', 'active-source pin and inverse moments',
               'source cap exponent positive', 'complete tau exponent positive',
               'full degree budget M>2', 'chosen deficit below complete energy budget',
               'complete scalar margin above 1/10 positive'],
    'source_sha256': {k:hashlib.sha256(p.read_bytes()).hexdigest() for k,p in source_paths.items()},
    'all_assertions_passed': True,
}
(HERE/'certificate.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({
    'all_assertions_passed': True, 'gamma_chosen':str(gamma),
    'c_raw_chosen':str(craw), 'margin_above_one_tenth':str(craw-F(1,10)),
    'tau_exponent_gap':str(sigma/2-chi*(1-sigma)),
    'source_cap_exponent_gap':str(eta-chi),
    'gamma_budget_display':float(gamma_budget),
    'source_identity':result['source_sha256'],
},indent=2))
