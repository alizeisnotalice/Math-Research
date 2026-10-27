"""New R069 independent explicit-threshold check, exact rational arithmetic.

This uses the weaker n>=L^8/2^16, cap/tau Taylor12, and scalar Taylor5.
It imports no author helper and recomputes no accepted logarithm series.
All statements about infinite tails are justified in REVIEW.md, not sampled.
"""
from fractions import Fraction as Q
from pathlib import Path
from math import factorial
from datetime import datetime, timezone
import hashlib, json, os, platform, sys

HERE=Path(__file__).resolve().parent
BASE=Path('/Users/zhengzhihao/Desktop/T')
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
inputs={
 'R069_threshold':('cube-r030-parameters/cube-dim-7-bounds/upper/round069/profile_threshold/THRESHOLD_PROOF.md','d1e0420ba6ce509d2a1e9d5e0c6c01fd8aedd62c52adc3bc187b584fd646da10'),
 'R066_operator':('cube-r030-geometry/cube-dim-7-bounds/upper/round066/target_profile/PROOF.md','4b7b3210f063312527c08adab5b7f81dc4b2c0ecf879e466dc715eb56ea6ebd1'),
 'R066_full_review':('cube-r030-uniform-audit/cube-dim-7-bounds/audit/round066/target_profile/review.md','5a13e5a4364a2fbcc4b9d89fb2536ebeded2569033041017a7355237fd110ea8'),
 'R064_relative':('cube-r030-roots-audit/cube-dim-7-bounds/audit/round064/relative_coordinate/LEMMA.md','c60d48d116863c9c9cac177e19a9fa4158f12fd8dcef0d5239a5ad0a5487d9f1'),
 'R066_profile':('cube-r030-roots-audit/cube-dim-7-bounds/audit/round066/first_profile/LEMMA.md','67bdf01d9aa65c70bf058df08b527b904e954d19f3327ecbb43e5b5c316b5030'),
 'R066_price':('cube-r030-parameters/cube-dim-7-bounds/audit/round066/exact_scalar/SCALAR_PROOF.md','d058c87d90db9748ac05a6d42fd711f68546daa217eee5ef8106b9c6e1d689f0')}
checks={}
def ok(name,p):
 assert bool(p),name
 checks[name]=True
sources={}
for name,(rel,expected) in inputs.items():
 p=BASE/rel
 ok('identity_'+name,sha(p)==expected)
 sources[name]={'path':str(p),'sha256':sha(p)}

L=Q(2**30); sqrtL=Q(2**15)
N=L**8/2**16   # deliberately weaker than the author's L^8/8!
ql,qu,s=Q(2,3),Q(917,1000),Q(3,2)
sig,chi=Q(83,1000),Q(9,200)
x=(1-sig)*L
delta=Q(1,2**30); rlo=1+delta; rhi=s*qu-delta
hc=Q(1,2**80); h=hc/L
e,epsA,epsRet=Q(1,2**20),Q(1,2**16),Q(1,2**12)
epsmax=1/(x*N**2)
ok('weaker_exponential_lower_recipe',factorial(8)<2**16 and N>L**2/2)
ok('base_domain',L>=64 and sqrtL**2==L and 0<sig<Q(1,2) and L/2<x<L and x>1000)
ok('trim_and_positive_source_reference',0<h<min(delta/2,Q(1,1000)) and s*ql-h>0 and ql/2-h>0)
ok('epsilon_relative_active_guards',0<epsmax<min(Q(1,2),h) and epsmax*qu<2*h/3 and rlo-2*h>1)
ok('matching_trim_both_ends',s*ql==1 and s*qu-rhi==delta and 2*h<delta)
ok('separated_supports_and_rho_domain',N/L>16 and ql*N/L>4*rhi>4 and 1<rlo<rhi<Q(3,2))
ok('enough_n_for_cap_and_probabilities',N>=400 and N>=65 and N*ql>=200000)
ok('relative_gap_domain',rhi/3<Q(1,2) and qu>ql and 0<qu<1)

mean_log=(Q(3,4)*L**2/ql**2+Q(1,2)*L/ql)/N
var_log=Q(3,2)*L**2/(ql**2*N)
fail_log=var_log/(e-mean_log)**2
fail_count=4/(e**2*N)
fail_profile=Q(200000**2)/(ql**2*N)
fail_cap=Q(3,4)**64
ok('baseline_mean_available',mean_log<e/2)
ok('baseline_log_Chebyshev_direct',fail_log<epsA/8)
ok('both_count_failure_bounds',fail_count<epsA/8)
ok('count_endpoint_drift',2*h+1/N<3*e/4)
ok('target_cap_coordinate_margin',1/rhi-L/N>Q(1,2))
ok('target_cap_each_signed_face',fail_cap<epsA/8)
ok('profile_pin_budget',1/(N*ql)<=Q(1,200000))
ok('profile_failure_coefficient',Q(200000**2)/ql**2==90000000000)
ok('profile_event_Chebyshev',fail_profile<epsA/8)
ok('five_event_union',5*epsA/8<epsA and fail_log+2*fail_count+fail_cap+fail_profile<epsA)

ok('product_log_bridge',512*L*h==Q(1,2**71)<Q(1,2**60))
ok('individual_J_separate_precision',L/(3*ql*N)<e)
ok('profile_transport_budget',h<Q(1,1600000) and 16*h<Q(1,100000))
ok('whole_q_profile_coefficient',Q(221,5000)+2*Q(1,100000)==Q(2211,50000))
ok('radial_relative_width',0<h/rlo<1 and 1/(1-(h/rlo)**2)<1+e)

ret=Q(5000**2)/(150*x)
ok('source_pin_and_retention_gap',Q(1,5)/x<Q(1,5000) and chi-Q(223,5000)-Q(1,5000)==Q(1,5000))
ok('retention_exact_base_value',ret==Q(1953125,11538530304))
ok('every_source_face_retention',ret<epsRet)

Csrc=(qu-ql)/ql; eta=(1-qu)/(2*qu)
safe=Q(1,10000); cS=(1-sig)*safe
nu=sig/2-chi*(1-sig)
cap12=Q(factorial(12))*Csrc*rhi/(hc*cS**12*L**11)
tau12=Q(factorial(12))*rhi*sqrtL/(hc*nu**12*L**11)
ok('wide_source_prefactor',Csrc==Q(751,2000) and 0<Csrc*rhi<1)
ok('cap_rate_identity',eta==Q(83,1834) and eta-chi==Q(47,183400))
ok('finite_n_source_cap_gap',eta-chi-eta/400>safe and eta-chi-eta/N>safe)
ok('new_cap_Taylor12_payment',cap12<e)
ok('tau_rate_identity',nu==Q(47,200000)>0)
ok('new_tau_Taylor12_payment',tau12<e)

gamma=Q(317,1000)
ok('B_lower_and_sqrt_positive_guard',L>=16 and L/4-4/L-1>0 and 0<gamma<1 and x>4)
ok('finite_a_all_dilations',8*N**4>2*N and 1/N**3-Q(4)/(8*N**4)>0)
ok('target_item_and_source_item_order',rhi+1<3 and ql*N/L>4)

# Existing accepted price/budget receipts are read and identified, not recomputed.
pricepath=BASE/'cube-r030-parameters/cube-dim-7-bounds/audit/round066/exact_scalar/verification.json'
price=json.loads(pricepath.read_text())
ok('accepted_price_interface',Q(price['parameters']['gamma'])==gamma and Q(price['parameters']['sigma'])==sig and Q(price['parameters']['price_upper'])==Q(4333,50000) and Q(price['net_slack'])==Q(17,50000))
sources['historical_price_receipt']={'path':str(pricepath),'sha256':sha(pricepath),'read_only':'accepted price interface; no old log sum rerun; historical conditional identity unchanged'}
fixedpath=BASE/'cube-r030-geometry/cube-dim-7-bounds/upper/round066/target_profile/verification.json'
fixed=json.loads(fixedpath.read_text())
ok('accepted_fixed_budget_interface',Q(fixed['gamma'])==gamma and Q(fixed['mass_lower'])==Q(637,2000) and Q(fixed['budget'])>gamma and Q(fixed['profile_mu_upper'])==Q(221,5000) and Q(fixed['source_mean_upper'])==Q(223,5000) and Q(fixed['source_variance_upper'])==Q(1,150))
sources['historical_fixed_receipt']={'path':str(fixedpath),'sha256':sha(fixedpath),'read_only':'accepted fixed constants only; 16-gate computation not rerun'}
slack=gamma-Q(4333,50000)-Q(23,100)
scalar5=2/L**2+Q(factorial(5))/(sig*L)**5+2*sqrtL*Q(factorial(5))/(sig*L/2)**5
ok('same_full_Phi_exact_price_slack',slack==Q(17,50000)>0)
ok('new_full_scalar_Taylor5_payment',scalar5<slack)
ok('old_full_cost_buffer_comparison',Q(61,200)-Q(91,1000)+Q(4,2**20)<Q(9,40)<Q(23,100))

values={k:{'rational':str(v),'display_only':float(v)} for k,v in {
 'N_lower':N,'baseline_mean':mean_log,'baseline_variance':var_log,
 'baseline_failure_direct':fail_log,'each_count_failure':fail_count,
 'profile_failure':fail_profile,'cap_perface_failure':fail_cap,
 'inverse_perface_failure':ret,'retention_margin':epsRet-ret,
 'cap_Taylor12':cap12,'tau_Taylor12':tau12,
 'scalar_Taylor5':scalar5,'scalar_slack':slack}.items()}
out={'status':'PASS','scope':'new independent full-operator threshold guards for integer n with natural log n >= 2^30; no old receipt revision',
 'utc':datetime.now(timezone.utc).isoformat(),'pid':os.getpid(),'interpreter':sys.executable,'python':platform.python_version(),
 'recipe':'n>=L^8/2^16; direct baseline Chebyshev; cap/tau Taylor12; full scalar Taylor5; exact Fraction arithmetic',
 'code_sha256':sha(Path(__file__)),'sources':sources,'checks':checks,'check_count':len(checks),'values':values,
 'tails':{'N':'increasing L^8','N/L':'increasing L^7','epsilon/h':'decreasing L^-16',
 'meanlog':'sum of decreasing L^-6,L^-7','variance':'decreasing L^-6',
 'baseline_failure':'decreasing numerator/increasing positive denominator',
 'counts/profile':'decreasing L^-8','count_drift':'decreasing sum L^-1,L^-8',
 'J_individual':'decreasing L^-7','retention':'decreasing L^-1',
 'cap_Taylor12':'decreasing L^-11','tau_Taylor12':'decreasing L^-21/2',
 'scalar_Taylor5':'sum decreasing L^-2,L^-5,L^-9/2'},
 'no_author_helper':True,'no_old_log_or_baseline_replay':True,'no_production_or_MATLAB_or_Git':True}
(HERE/'verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':'PASS','checks':len(checks),'code_sha256':out['code_sha256'],'values':{k:v['display_only'] for k,v in values.items()}},indent=2))
