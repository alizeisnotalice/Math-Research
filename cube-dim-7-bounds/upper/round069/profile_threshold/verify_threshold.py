"""A new fixed L0 threshold certificate; no old verifier/log/grid replay."""
from fractions import Fraction as F
from pathlib import Path
import datetime,hashlib,json,math,os,sys
root=Path(__file__).resolve().parent
L=2**30; sqrtL=2**15; fac8=math.factorial(8); fac10=math.factorial(10); fac4=math.factorial(4)
N=F(L**8,fac8); ql=F(2,3); qu=F(917,1000); s=F(3,2)
delta=F(1,2**30); rl=1+delta; ru=s*qu-delta
sigma=F(83,1000); chi=F(9,200); x=(1-sigma)*L; gamma=F(317,1000)
hc=F(1,2**80); h=hc/L; e=F(1,2**20); eA=F(1,65536); eRet=F(1,4096)
epsilon=1/(x*N*N); gmax=F(1,2)
Mlog=fac8*(6*gmax**3/(ql*ql*L**6)+gmax/(ql*L**7))
Vlog=12*fac8*gmax**3/(ql*ql*L**6)
profile_bad=F(200000**2)/(ql*ql*N)
radial_ratio=1/(1-(h/rl)**2)
eta=(1-qu)/(2*qu); capgap=eta-chi; cS=(1-sigma)/10000
Csrc=(qu-ql)/ql
cap10=F(fac10)*Csrc*ru/(hc*cS**10*L**9)
nu=sigma/2-chi*(1-sigma)
tau10=F(fac10)*ru*sqrtL/(hc*nu**10*L**9)
R4=2/F(L*L)+F(fac4)/(sigma*L)**4+2*sqrtL*F(fac4)/(sigma*L/2)**4
Rpower=F(1,2**59)+F(fac4,2**104)+F(2**16*fac4,2**100)
slack=F(17,50000); priceupper=F(4333,50000)
checks={
 'fixed_L0':L==2**30 and sqrtL**2==L,
 'band_domain':F(1,2)<ql<qu<1,
 'target_domain':1<rl<ru<F(3,2),
 'trimmed_matching_nonempty':2*delta<s*(qu-ql),
 'sigma_domain':0<sigma<F(1,2) and F(L,2)<x<L,
 'n_Taylor8_dominates_n_Taylor2':N>F(L*L,2),
 'n_at_least_400':N>=400,
 'n_at_least_65':N>=65,
 'h_below_trim':h<delta/2,
 'h_coordinate_domain':h<F(1,1000),
 'epsilon_below_half':epsilon<F(1,2),
 'relative_gap_guard':epsilon*qu<2*h/3,
 'epsilon_below_h':epsilon<h,
 'whole_band_R_positive':s*ql-h>0,
 'whole_band_gap_positive':F(1,2)*ql-h>0,
 'whole_band_relative_gap_endpoint':epsilon*qu-2*h/3<0,
 'active_R_above_one':rl-2*h>1,
 't_above_16':N/L>16,
 'source_rho_ge_4target':ql*N/L>4*ru,
 'source_target_disjoint':ql*N/L>4 and ru<2,
 'target_inside_t':ru<N/L,
 'kappa_positive_gap':epsilon<1,
 'g_below_half_active':ru/3<gmax,
 'pinned_baseline_log_domain':gmax/(ql*N/L)<F(1,8),
 'baseline_mean':Mlog<e/2,
 'baseline_variance_failure':Vlog<(e/2)**2*eA/8,
 'minus_count_failure':4/(e*e*N)<eA/8,
 'plus_count_failure':4/(e*e*N)<eA/8,
 'count_whole_interval_drift':2*h+1/N<3*e/4,
 'target_cap_coordinate':1/ru-L/N>F(1,2),
 'target_cap_perface_probability':F(3,4)**64<eA/8,
 'profile_pin_mean':N*ql>=200000,
 'profile_bad_identity':F(200000**2)/(ql*ql)==90000000000,
 'profile_failure':profile_bad<eA/8,
 'five_event_total':5*eA/8<eA,
 'product_log_bridge':512*L*h==F(1,2**71)<F(1,2**60),
 'P0_positive':0<1-e-F(1,2**60)<1,
 'individual_J_precision':F(L)/(3*ql*N)<e,
 'profile_transport':16*h<F(1,100000),
 'profile_transport_h_guard':h<F(1,1600000),
 'profile_transport_sum':F(221,5000)+F(2,100000)==F(2211,50000),
 'radial_logwidth_domain':0<h/rl<1,
 'radial_logwidth_upper':radial_ratio<1+e,
 'inverse_pin_budget':x>=1000 and F(1,5)/x<=F(1,5000),
 'inverse_mean_gap':chi-F(223,5000)-F(1,5000)==F(1,5000),
 'inverse_face_retention':F(5000**2)/(150*x)<eRet,
 'source_Csrc_radial_prefactor':0<Csrc*ru<1,
 'source_eta_identity':eta==F(83,1834),
 'source_cap_rate_gap':capgap==F(47,183400)>0,
 'source_cap_n400_payment':capgap-eta/400>F(1,10000),
 'source_cap_actual_n_payment':capgap-eta/N>F(1,10000),
 'source_prefactor_power':Csrc*ru/h<2**110,
 'source_argument_power':cS*L>2**16,
 'source_Taylor10_power':F(fac10,2**50)<e,
 'source_Taylor10_exact':cap10<e,
 'scalar_positive_root_domain':L>=64 and F(L,4)-F(4,L)-1>0 and 0<gamma<1,
 'scalar_x_domain':x>4 and N/L>1,
 'tau_rate_identity':nu==F(47,200000)>0,
 'tau_radial_prefactor_power':ru*sqrtL/h<2**126,
 'tau_argument_power':nu*L>2**17,
 'tau_Taylor10_power':F(fac10,2**44)<e,
 'tau_Taylor10_exact':tau10<e,
 'finite_a_radius_guard':8*N**4>2*N,
 'finite_a_center_guard':1/N**3-F(4)/(8*N**4)>0,
 'exact_scalar_slack':gamma-priceupper-F(23,100)==slack>0,
 'scalar_argument_powers':sigma*L>2**26 and sigma*L/2>2**25,
 'scalar_Taylor4_enclosure':R4<Rpower,
 'scalar_Taylor4_full_margin':Rpower<slack,
 'old_cost_positive_Bbuffer':F(61,200)-F(91,1000)+F(4,2**20)<F(9,40)<F(23,100),
 'old_comparison_algebra_domain':L>=2**20,
}
assert all(checks.values()),[k for k,v in checks.items() if not v]
base=Path('/Users/zhengzhihao/Desktop/T')
hashof=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
sources={
 'R066_full_proof':(base/'cube-r030-geometry/cube-dim-7-bounds/upper/round066/target_profile/PROOF.md','4b7b3210f063312527c08adab5b7f81dc4b2c0ecf879e466dc715eb56ea6ebd1'),
 'R066_full_review':(base/'cube-r030-uniform-audit/cube-dim-7-bounds/audit/round066/target_profile/review.md','5a13e5a4364a2fbcc4b9d89fb2536ebeded2569033041017a7355237fd110ea8'),
 'relative_coordinate_lemma':(base/'cube-r030-roots-audit/cube-dim-7-bounds/audit/round064/relative_coordinate/LEMMA.md','c60d48d116863c9c9cac177e19a9fa4158f12fd8dcef0d5239a5ad0a5487d9f1'),
 'R066_exact_scalar_proof':(base/'cube-r030-parameters/cube-dim-7-bounds/audit/round066/exact_scalar/SCALAR_PROOF.md','d058c87d90db9748ac05a6d42fd711f68546daa217eee5ef8106b9c6e1d689f0'),
 'R064_full_proof':(base/'cube-r040-cover-selector/cube-dim-7-bounds/upper/round064/radial_calibration/PROOF.md','17ce38ab5af46c0358046bdcc2aa57a05844c0799e9a97280d959887f59615d6'),
}
source_receipts={}
for name,(path,expected) in sources.items():
 actual=hashof(path); assert actual==expected,(name,actual)
 source_receipts[name]={'path':str(path),'sha256':actual}
pricepath=base/'cube-r030-parameters/cube-dim-7-bounds/audit/round066/exact_scalar/verification.json'
price=json.loads(pricepath.read_text())
assert F(price['parameters']['gamma'])==gamma and F(price['parameters']['sigma'])==sigma
assert F(price['parameters']['price_upper'])==priceupper and F(price['net_slack'])==slack
assert F(price['log_upper_S7'])<priceupper
source_receipts['R066_exact_scalar_historical_receipt']={'path':str(pricepath),'sha256':hashof(pricepath),'scope':'read existing rational log certificate; seven-term recipe not rerun; original conditional-run identity preserved'}
oldpath=base/'cube-r030-geometry/cube-dim-7-bounds/upper/round066/target_profile/verification.json'
old=json.loads(oldpath.read_text());assert F(old['gamma'])==gamma and F(old['profile_mu_upper'])==F(221,5000) and F(old['source_mean_upper'])==F(223,5000) and F(old['source_variance_upper'])==F(1,150)
source_receipts['R066_author_historical_receipt']={'path':str(oldpath),'sha256':hashof(oldpath),'scope':'reads accepted fixed budget interfaces only; author 16 gates not rerun'}
review64=base/'cube-r030-soft/cube-dim-7-bounds/audit/round064/radial_calibration/REVIEW.md'
source_receipts['R064_full_review']={'path':str(review64),'sha256':hashof(review64),'scope':'eventual operator validity retained; no explicit N64 imported'}
receipt={
 'status':'PASS: new L0=2^30 threshold gates; soft independent threshold audit still required',
 'actual_pid':os.getpid(),'actual_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'python':sys.version,
 'domain':'integer n with log n >= 2^30; original centered common-radius cube operator',
 'L0':str(L),'arithmetic':'standard-library Fraction/integer, no floating inference',
 'checks':checks,'check_count':len(checks),
 'values':{'N':str(N),'x0':str(x),'h0':str(h),'baseline_mean_upper':str(Mlog),'baseline_variance_upper':str(Vlog),'profile_failure_upper':str(profile_bad),'radial_width_upper':str(radial_ratio),'perface_inverse_failure_upper':str(F(5000**2)/(150*x)),'source_cap_Taylor10_upper':str(cap10),'tau_Taylor10_upper':str(tau10),'scalar_Taylor4_upper':str(R4),'scalar_power_upper':str(Rpower),'exact_scalar_slack':str(slack)},
 'tails':{'N':'L^8 increasing','N_over_L':'L^7 increasing','h':'L^-1','epsilon_over_h':'L^-16 via N','baseline_mean':'positive sum L^-6,L^-7','baseline_variance':'L^-6','counts_and_profile_failure':'L^-8','count_drift':'positive sum L^-1,L^-8','individual_J':'L^-7','radial_ratio':'1/(1-constant/L^2) decreases','inverse_retention':'L^-1','source_cap_Taylor10':'L^-9','tau_Taylor10':'L^-17/2','scalar_Taylor4':'positive sum L^-2,L^-4,L^-7/2','old_cost_buffer':'constant+4/L decreases'},
 'sources':source_receipts,
 'own_sources':{name:hashof(root/name) for name in ('REGISTRATION.md','THRESHOLD_PROOF.md','verify_threshold.py')},
 'no_search':True,'no_old_program_replay':True,'no_log_recalculation':True,'no_MATLAB':True,'no_Git':True,
 'scope':'Same Phi66 and gamma; explicit threshold is new and not retroactively assigned to R066 receipt. R064 comparison only in each mechanism own valid domain; not a new Phi optimization.',
}
(root/'verification.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'status':receipt['status'],'checks':len(checks),'actual_pid':os.getpid(),'proof_sha256':receipt['own_sources']['THRESHOLD_PROOF.md'],'inverse_failure':str(F(5000**2)/(150*x)),'profile_failure':str(profile_bad),'slack':str(slack)},indent=2))
