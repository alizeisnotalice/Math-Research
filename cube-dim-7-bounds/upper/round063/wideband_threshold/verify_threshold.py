"""One new L0 certificate; no author helper, grid, log-series, or old-gate replay."""
from fractions import Fraction as F
from pathlib import Path
import datetime, hashlib, json, math, os, sys

root=Path(__file__).resolve().parent
L=2**33; L56=2**20; sqrt_upper=2**17
fac=math.factorial(8); N=F(L**8,fac); n2=F(L**2,2)
d=F(391,1905); ab=F(283,2032); bb=F(63,64)
sigma=F(3,100); chi=F(147,10000); e=F(1,2**20); eA=F(1,65536); eRet=F(1,4096)
hc=F(1,2**48); delta=e; h=hc/(L*d); x=(1-sigma)*L
vmin=ab*d; vmax=bb*d; ql=1-vmax; qu=1-vmin; Rbase=ql+d
rl=Rbase+delta; ru=qu+d-delta
eta=vmin/(2*qu); gap=eta-chi; nu=sigma/2-chi*(1-sigma)
Csrc=(qu-ql)/ql
mean=fac*(6*d**3/(ql*ql*L**6)+d/(ql*L**7))
var=12*fac*d**3/(ql*ql*L**6)
Jerror=5*d*L/(8*ql*N)
zt=1/(x*N*N)
cS=(1-sigma)/200000
cap8=fac*Csrc*d/(hc*cS**8*L**7)
# Exact sqrt(L) is irrational here; upper sqrt_upper is deliberately rational.
tau8=F(fac)*d*sqrt_upper/(hc*nu**8*L**7)
scalar8=2/F(L*L)+F(fac)/(sigma*L)**8+2*sqrt_upper*F(fac)/(sigma*L/2)**8
scalar_power=F(1,2**65)+F(fac,2**216)+F(2**18*fac,2**208)
gamma=F(1637,12500); slack=F(39,1212500)
gamma56=F(172531397335235056462877787662549380239,1643616601687301779182642153860479582208)
checks={
 'fixed_L0':L==2**33,
 'local_d_domain':0<d<F(1,4)<F(1,3),
 'source_band_domain':F(1,2)<ql<qu<1,
 'target_band_domain':1<Rbase<rl<ru<F(5,4),
 'target_interval_nonempty':2*delta<qu-ql,
 'scalar_sigma_domain':0<sigma<F(1,2) and F(L,2)<x<L,
 'sqrt_rational_upper':L<sqrt_upper**2,
 'n_Taylor8_dominates_n_Taylor2':N>n2,
 'n_at_least_200000':N>=200000,
 'n_at_least_65':N>=65,
 'nd_at_least_one':N*d>=1,
 'h_inside_inset':h<delta/2,
 'h_inside_coordinate_domain':h<d/8,
 'whole_band_R_above_one':h<(Rbase-1)/2,
 'zeta_t_small':zt<h/2,
 'zeta_t_coordinate_domain':zt<d,
 'kappa_positive_gap':zt<1,
 't_above_16':N/L>16,
 'source_rho_ge_4target':ql*N/L>4*ru,
 'source_target_disjoint':ql*N/L>4 and ru<2,
 'target_inside_t':ru<N/L,
 'local_gap_range':0<d-h and d-h+zt*qu<d<Rbase/2,
 'pinned_J_domain':d/(ql*N/L)<F(1,8),
 'baseline_mean_half_tolerance':mean<e/2,
 'baseline_variance_failure_budget':var<(e/2)**2*eA/8,
 'Kw_failure':1/(4*N*d*d)<eA/8,
 'minus_count_failure':4/(e*e*N)<eA/8,
 'plus_count_failure':4/(e*e*N)<eA/8,
 'count_lower_drift':2*h+1/N<3*e/4,
 'count_upper_drift':2*h*d<3*e/4,
 'target_cap_coordinate':1/ru-L/N>F(1,2),
 'target_cap_probability':F(3,4)**64<eA/8,
 'five_bad_events_total':5*eA/8<eA,
 'product_lower_positive':0<1-e-1280*hc<1,
 'J_uniform_precision':Jerror<e,
 'inverse_mean_gap_argument':x>=F(40000,9),
 'inverse_gap_identity':chi-F(293,20000)==F(1,20000),
 'inverse_retention_failure':F(40000**2)/(900*x)<eRet,
 'source_Csrc_domain':0<Csrc<1,
 'source_eta_domain':0<eta<1,
 'source_cap_gap_identity':gap==F(499871,37603070000)>F(1,100000),
 'source_cap_finite_n_gap':gap-eta/N>F(1,200000),
 'source_prefactor_power':Csrc/h<2**81,
 'source_argument_power':cS*L>2**15,
 'source_Taylor8_power_budget':F(2**81*fac,2**120)<e,
 'source_Taylor8_exact':cap8<e,
 'scalar_positive_domain':L>=64 and F(L,4)-F(4,L)-1>0,
 'scalar_t_positive_domain':x>4 and N/L>1,
 'tau_exponent_identity':nu==F(741,1000000)>0,
 'tau_prefactor_power':F(sqrt_upper)/h<2**98,
 'tau_argument_power':nu*L>2**22,
 'tau_Taylor8_power_budget':F(2**98*fac,2**176)<e,
 'tau_Taylor8_exact_upper':tau8<e,
 'finite_a_radius_guard':8*N**4>2*N,
 'finite_a_center_guard':1/N**3-F(4)/(8*N**4)>0,
 'scalar_slack_identity':gamma-sigma/(1-sigma)-F(1,10)==slack>0,
 'scalar_full_argument':sigma*L>2**27 and sigma*L/2>2**26,
 'scalar_power_enclosure':scalar8<scalar_power,
 'scalar_full_margin':scalar_power<slack,
 'old_full_cost_buffer':gamma56-F(7,100)+F(4,L56)<F(1,10),
 'valid_domain_intersection':L>=L56,
}
assert all(checks.values()),[k for k,v in checks.items() if not v]

base=Path('/Users/zhengzhihao/Desktop/T')
sources={
 'R059_full_proof':(base/'cube-r040-cover-selector/cube-dim-7-bounds/upper/round059/wide_source_band/PROOF.md','fe95f369634bbc753f5e44a32d6a2b74b288748ebd4bfe1f6c3cac43b086bd8f'),
 'R059_full_review':(base/'cube-r030-soft/cube-dim-7-bounds/audit/round059/wide_source_band/review.md','c8aff39bf68749d25f9978cef11be82ccb7907ec7e72c5678ac66b06a65a33d4'),
 'local_coordinate_lemma':(base/'cube-r030-roots-audit/cube-dim-7-bounds/audit/round056/one_coordinate/LEMMA.md','97430af9980c2de5d105bc4ea86667573c5c28e0ff40a92f1c613cd95a243eec'),
 'R056_full_proof':(base/'cube-r030-geometry/cube-dim-7-bounds/upper/round056/refined_inverse/PROOF.md','6986ed9ec48b46f3136a9a56722d0092db4f2942d52e35c9f4de37eb279f399a'),
}
hashof=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
source_receipts={}
for name,(path,expected) in sources.items():
 actual=hashof(path);assert actual==expected,(name,actual)
 source_receipts[name]={'path':str(path),'sha256':actual}
certpath=base/'cube-r040-cover-selector/cube-dim-7-bounds/upper/round059/wide_source_band/certificate.json'
cert=json.loads(certpath.read_text())
# Structural read of accepted exact inputs only; no 64-term log or old checks replay.
for key,value in {'d':d,'a_band':ab,'b_band':bb,'chi':chi,'sigma':sigma,'inset':delta,'h_coefficient':hc,'epsilon_J':e,'epsilon_N':e,'epsilon_S':e,'epsilon_T':e,'epsilon_D':e,'epsilon_Ret':eRet,'epsilon_A':eA,'baseline_log_tolerance':e}.items():
 assert F(cert['parameters'][key])==value,key
for key,value in {'qmin':ql,'qmax':qu,'Rbase':Rbase,'rlo':rl,'rhi':ru,'gamma_chosen':gamma,'caprate':eta,'source_cap_exponent_gap':gap,'tau_exponent_gap':nu}.items():
 assert F(cert['derived_rationals'][key])==value,key
source_receipts['R059_historical_certificate']={'path':str(certpath),'sha256':hashof(certpath),'scope':'exact input/derived-field compatibility only; original log certificate was not rerun or relabeled'}
oldpath=base/'cube-r030-geometry/cube-dim-7-bounds/upper/round056/refined_inverse/verification.json'
old=json.loads(oldpath.read_text());assert F(old['gamma'])==gamma56
source_receipts['R056_historical_receipt']={'path':str(oldpath),'sha256':hashof(oldpath),'scope':'reads accepted exact gamma only; no old verifier replay'}
receipt={
 'status':'PASS: new R063 explicit threshold certificate; independent threshold audit still required',
 'L0':str(L),'domain':'log n >= 2^33; n integer',
 'arithmetic':'Python standard-library Fraction/integer; no floating inference',
 'actual_pid':os.getpid(),'actual_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'python':sys.version,
 'checks':checks,'check_count':len(checks),
 'values':{'N_lower':str(N),'x0':str(x),'h0':str(h),'baseline_mean_upper':str(mean),'baseline_variance_upper':str(var),'J_error_upper':str(Jerror),'retention_probability_upper':str(F(40000**2)/(900*x)),'source_Taylor8_upper':str(cap8),'tau_Taylor8_upper_rational_sqrt':str(tau8),'scalar_Taylor8_upper':str(scalar8),'scalar_power_upper':str(scalar_power),'scalar_slack':str(slack),'R056_cost_buffer_at_L20':str(gamma56-F(7,100)+F(4,L56))},
 'tails':{'N':'L^8 increasing','N_over_L':'L^7 increasing','h':'L^-1','zeta_t_over_h':'L^-16 via N','baseline_mean':'positive sum L^-6 and L^-7','baseline_variance':'L^-6','Kw_and_counts':'L^-8','drift':'positive sum L^-1 and L^-8','Jerror':'L^-7','retention':'L^-1','source_cap_Taylor8':'L^-7','tau_Taylor8':'L^-13/2','scalar':'positive sum L^-2,L^-8,L^-15/2','R056_buffer':'constant + 4/L decreasing'},
 'sources':source_receipts,
 'own_sources':{name:hashof(root/name) for name in ('REGISTRATION.md','THRESHOLD_PROOF.md','verify_threshold.py')},
 'no_search':True,'no_old_verifier_replay':True,'no_MATLAB':True,'no_Git':True,
 'activation':'Do not enter accepted Cn ledger before independent new threshold audit. R059 full ACCEPT alone was eventual.'}
(root/'verification.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'status':receipt['status'],'checks':len(checks),'actual_pid':os.getpid(),'proof_sha256':receipt['own_sources']['THRESHOLD_PROOF.md'],'retention_upper':str(F(40000**2)/(900*x)),'sourcecap_upper':str(cap8),'tau_upper':str(tau8),'scalar_slack':str(slack)},indent=2))
