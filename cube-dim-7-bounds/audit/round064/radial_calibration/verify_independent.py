"""Independent R064 fixed-tuple checks; no author helper/imports.

Target log enclosure uses 40/41 alternating log1p terms, rather than
the author's 64-term atanh recipe. All assertions use Fraction/SymPy.
Decimal displays are diagnostics only. This script does not certify N.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import platform
import os
import sys
from datetime import datetime, timezone
import sympy as sp

HERE = Path(__file__).resolve().parent
PROOF = Path('/Users/zhengzhihao/Desktop/T/cube-r040-cover-selector/cube-dim-7-bounds/upper/round064/radial_calibration/PROOF.md')
LEMMA = Path('/Users/zhengzhihao/Desktop/T/cube-r030-roots-audit/cube-dim-7-bounds/audit/round064/relative_coordinate/LEMMA.md')
expected = {PROOF: '17ce38ab5af46c0358046bdcc2aa57a05844c0799e9a97280d959887f59615d6', LEMMA: 'c60d48d116863c9c9cac177e19a9fa4158f12fd8dcef0d5239a5ad0a5487d9f1'}
def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()
checks = {}
def check(name, predicate):
    assert bool(predicate), name
    checks[name] = True
def rat(x):
    return {'exact': str(x), 'display': float(x)}
for p, digest in expected.items():
    check('source_identity_' + p.name, sha(p) == digest)

e, s = F(1,1024), F(3,2)
qlo, qhi = F(2,3), F(909,1000)
trim = F(1,2**30)
rlo, rhi = 1+trim, s*qhi-trim
sig, chi, hc = F(91,1000), F(1,20), F(1,2**64)
P = 1-e-F(1,2**50)
ulo, uhi = P*(1-e)/(1+e), (1+e)/(P*(1-e))
degree = 2*(1+e)*uhi
rowminus = 2*(1-e)**2*ulo-2*e
massfloor = F(1239,4000)
Gamma = massfloor*rowminus/degree
budget, gamma = (1-e)*Gamma, F(61,200)
check('radii_and_band', 0 < qlo < qhi < 1 < rlo < rhi < F(3,2))
check('s_qlo_lower_trim', s*qlo == 1 and rlo-1 == trim)
check('s_qhi_upper_trim', s*qhi-rhi == trim)
check('P_bounds', 0 < P < 1)
check('positive_budget', rowminus > 0 and degree > 2)
check('Gamma_order', 0 < gamma < budget < Gamma < 1)

# log(rhi/rlo)=log(1+u). Alternating series: S40 < log < S41.
u = (rhi-rlo)/rlo
check('alternating_series_domain', 0 < u < 1)
loglo = sum(((-1)**(j+1)*u**j/F(j) for j in range(1,41)), F(0))
loghi = loglo + u**41/F(41)
check('log_enclosure_nonempty', loglo < loghi)
masslo, masshi = (1-e)*loglo, (1-e)*loghi
check('author_massfloor_strict', masslo > massfloor)
check('distinct_massfloor_identity', massfloor-F(3097,10000) == F(1,20000))

z = F(1,5)
Q = sum(((-1)**(j+1)*z**j/F(j) for j in range(1,6)),F(0))
meanrate = s*Q**2/(1-e)
varrate = (2*s/F(3))*Q**3/(1-e)
check('pin_log_polynomial', Q == F(34187,187500) and Q < F(1,5))
check('source_meanrate', meanrate < F(999,20000))
check('source_varrate', varrate < F(1,150))
check('retention_exponent_gap', chi-F(999,20000) == F(1,20000))
check('pin_paid_when_x_ge8000', F(1,5) <= F(8000,40000))
check('pin_relative_bound', 1-F(1,3)/2 == F(5,6))
check('bridge_value', 512*hc == F(1,2**55))
check('bridge_slack', 512*hc < F(1,2**50))
check('nrho_bridge', (1-sig)/qlo < F(3,2))

Csrc = (qhi-qlo)/qlo
eta = (1-qhi)/(2*qhi)
check('wide_source_mass_constant', Csrc == F(727,2000))
check('cap_rate', eta == F(91,1818))
check('cap_positive_gap', eta-chi == F(1,18180))
check('eventual_safe_cap_gap', (eta-chi)/2 == F(1,36360))
check('tau_positive_gap', sig/2-chi*(1-sig) == F(1,20000))
raw = gamma-sig/(1-sig)
check('full_scalar_coefficient', raw == F(37249,181800))
check('strict_one_fifth_slack', raw-F(1,5) == F(889,181800) > 0)

# New symbolic geometry/count/measure identities, not sample evaluation.
q,h,r,eps,t,n = sp.symbols('q h r eps t n', positive=True)
R = sp.Rational(3,2)*q-h
a0 = (1-eps)*q
g = R-a0
ell = (1+a0/R)/2
check('new_R_derivative', sp.diff(R,q) == sp.Rational(3,2))
check('new_g_derivative', sp.diff(g,q) == sp.Rational(1,2)+eps)
check('relative_gap_identity', sp.simplify(g-R/3-(eps*q-2*h/3)) == 0)
check('negative_pin_derivative', sp.simplify(sp.diff(ell,q)+(1-eps)*h/(2*R**2)) == 0)
qm, qp = (r-h)/sp.Rational(3,2), (r+h)/sp.Rational(3,2)
check('right_endpoint_mean_one_radius', sp.simplify(R.subs(q,qp)-r) == 0)
check('left_endpoint_reference', sp.simplify(R.subs(q,qm)-(r-2*h)) == 0)
check('ordinary_source_width', sp.simplify(qp-qm-4*h/3) == 0)
check('logwidth_s_cancels', sp.simplify(qp/qm-(r+h)/(r-h)) == 0)
amin=(1-eps)*qm
deficit=1-(1+(1-1/n)*amin/r)/(1+amin/(r-2*h))
displayed=(2*h*amin/(r*(r-2*h))+amin/(n*r))/(1+amin/(r-2*h))
check('count_lower_mean_exact', sp.simplify(deficit-displayed)==0)
ap=(1-eps)*qp
upperdef=1-(1+(1-1/n)*ap/r)/(1+ap/r)
check('count_upper_mean_exact', sp.simplify(upperdef-ap/(n*(r+ap)))==0)
rh=sp.symbols('rhi',positive=True)
check('target_radial_measure_cancel', sp.simplify((r/rh)/r-1/rh)==0)
beta, dilation, m = sp.symbols('beta dilation m',positive=True)
alpha=dilation/(dilation+1)
cc=(dilation*beta-1)/(dilation+1)
check('finite_a_envelope_identity', sp.simplify(alpha-cc*m-(1-beta*m)-((beta+1)*m-1)/(dilation+1))==0)
zeta=sp.symbols('zeta',positive=True)
kappa=beta-zeta
check('finite_a_capture_guard', sp.simplify(alpha*beta-kappa-2/(dilation+1)-(zeta-(beta+2)/(dilation+1)))==0)

result = {
 'status':'PASS', 'scope':'new independent fixed-tuple/symbolic checks; no explicit threshold; operator arguments audited separately',
 'utc':datetime.now(timezone.utc).isoformat(), 'pid':os.getpid(),
 'python':sys.executable,'python_version':platform.python_version(),'sympy':sp.__version__,
 'source_hashes':{str(p):sha(p) for p in expected}, 'script_sha256':sha(Path(__file__)),
 'recipe':{'target_log':'40-even/41-odd alternating log1p sums in exact Fraction', 'pin_log':'odd-five log1p upper polynomial', 'symbolic':'new matching/count/radial and finite-a identities'},
 'check_count':len(checks), 'checks':checks,
 'margins':{k:rat(v) for k,v in {
   'mass_lower_minus_author_floor':masslo-massfloor,
   'target_log_enclosure_width':loghi-loglo,
   'gamma_budget_minus_61_200':budget-gamma,
   'Gamma_L':Gamma,'gamma_budget':budget,'M':degree,
   'meanrate_margin':F(999,20000)-meanrate,
   'varrate_margin':F(1,150)-varrate,
   'raw_minus_one_fifth':raw-F(1,5)}.items()}}
(HERE/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'status':result['status'],'check_count':len(checks),'script_sha256':result['script_sha256'],'margins':{k:v['display'] for k,v in result['margins'].items()}},indent=2))
