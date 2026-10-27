#!/usr/bin/env python3
"""New independent R063 finite-threshold recipe; no author code or old verifier imports."""
import hashlib
import json
import math
import platform
from fractions import Fraction as Q
from pathlib import Path

HERE=Path(__file__).resolve().parent
BASE=Path('/Users/zhengzhihao/Desktop/T')
SRC=BASE/'cube-r030-parameters/cube-dim-7-bounds/upper/round063/wideband_threshold/THRESHOLD_PROOF.md'
EXPECTED='70fc76d864f233900baadd7f177b7532b654355c3ce3b515a822fe6b65a19133'
assert hashlib.sha256(SRC.read_bytes()).hexdigest()==EXPECTED

L=Q(2**33); d=Q(391,1905); ab=Q(283,2032); bb=Q(63,64)
sig=Q(3,100);chi=Q(147,10000);x=(1-sig)*L
e=Q(1,2**20);ea=Q(1,65536);eret=Q(1,4096);hc=Q(1,2**48)
vmin=ab*d;vmax=bb*d;qlo=1-vmax;qhi=1-vmin
rb=qlo+d;rlo=rb+e;rhi=qhi+d-e;h=hc/(L*d)
nb=L**8/2**16  # e^L >= L^8/8! > this weaker exact power bound.
sqrt_upper=Q(2**17)
eta=vmin/(2*qhi);gap=eta-chi;nu=sig/2-chi*(1-sig)
csrc=(qhi-qlo)/qlo;gamma=Q(1637,12500);slack=Q(39,1212500)

base_mean=6*d**3*L**2/(qlo**2*nb)+d*L/(qlo*nb)
base_var=12*d**3*L**2/(qlo**2*nb)
zt=1/(x*nb**2)
tail_per_face=Q(40000**2,900)/x
cS=(1-sig)/200000

# Different orders from the author: source cap Taylor-10, tau Taylor-6,
# full scalar Taylor-4. All true tail functions are decreasing powers.
cap10=Q(math.factorial(10))*csrc/(h*(cS*L)**10)
tau6=Q(math.factorial(6))*sqrt_upper/(h*(nu*L)**6)
scalar4=1/nb+Q(math.factorial(4))/(sig*L)**4+2*sqrt_upper*math.factorial(4)/(sig*L/2)**4

checks={
 'weaker_exponential_lower_bound':math.factorial(8)<2**16 and nb==2**248,
 'tuple_and_band_domain':0<d<Q(1,4)<Q(1,3) and Q(1,2)<qlo<qhi<1 and 1<rb<rlo<rhi<Q(5,4),
 'x_and_t_base':Q(L,2)<x<L and x>Q(40000,9) and nb/L>16,
 'large_n_and_Kw_threshold':nb>=200000 and nb*d>1,
 'matching_trim_and_global_R':h<min(e/2,d/8,(rb-1)/2),
 'zeta_local_and_gap':zt<h/2 and zt<d and zt<1 and d-h+zt*qhi<d<rb/2,
 'source_target_scale_separation':qlo*nb/L>4*rhi and rhi<nb/L,
 'baseline_failure_budget':base_mean<e/2 and base_var<(e/2)**2*ea/8,
 'Kw_fifth_event':1/(4*nb*d*d)<ea/8,
 'both_endpoint_count_events':4/(e*e*nb)<ea/8,
 'whole_interval_count_drift':2*h+1/nb<3*e/4 and 2*h*d<3*e/4,
 'target_cap_face_budget':1/rhi-L/nb>Q(1,2) and Q(3,4)**64<ea/8,
 'five_bad_event_union':5*ea/8<ea,
 'individual_J_and_separate_product':5*d*L/(8*qlo*nb)<e and 0<1-e-1280*hc<1,
 'inverse_face_retention':tail_per_face<eret and chi-Q(293,20000)==Q(1,20000),
 'wide_Csrc_and_eta':0<csrc<1 and csrc>d and 0<eta<1,
 'finite_n_source_cap_exponent':gap>Q(1,100000) and gap-eta/nb>Q(1,200000),
 'cap10_power_enclosure':csrc/h<2**79 and cS*L>2**15 and cap10<Q(math.factorial(10),2**71)<e,
 'sqrt_base_upper':L<sqrt_upper**2,
 'tau6_power_enclosure':sqrt_upper/h<2**96 and nu*L>2**22 and tau6<Q(math.factorial(6),2**36)<e,
 'tau_exponent_identity':nu==Q(741,10**6)>0,
 'root_positive_domain':L>=64 and x>4 and L/4-4/L-1>0 and 0<gamma<1,
 'finite_algorithmic_dilation':nb>1 and 8*nb**4>2*nb and 2*nb-1>0,
 'complete_scalar_slack':gamma-sig/(1-sig)-Q(1,10)==slack>0,
 'scalar4_power_enclosure':sig*L>2**27 and sig*L/2>2**26 and scalar4<(Q(1,2**248)+Q(24,2**108)+Q(24,2**86))<slack,
}
assert all(checks.values()),[k for k,v in checks.items() if not v]

# Read only already accepted fixed data; do not re-evaluate log series or old budgets.
paths={
 'R063_source':SRC,
 'R059_final_source':BASE/'cube-r040-cover-selector/cube-dim-7-bounds/upper/round059/wide_source_band/PROOF.md',
 'R059_full_audit':BASE/'cube-r030-soft/cube-dim-7-bounds/audit/round059/wide_source_band/review.md',
 'R056_local_lemma':BASE/'cube-r030-roots-audit/cube-dim-7-bounds/audit/round056/one_coordinate/LEMMA.md',
 'R059_historical_certificate':BASE/'cube-r040-cover-selector/cube-dim-7-bounds/upper/round059/wide_source_band/certificate.json',
 'R056_final_source':BASE/'cube-r030-geometry/cube-dim-7-bounds/upper/round056/refined_inverse/PROOF.md',
 'R056_historical_values':BASE/'cube-r030-geometry/cube-dim-7-bounds/upper/round056/refined_inverse/verification.json',
}
expected={
 'R059_final_source':'fe95f369634bbc753f5e44a32d6a2b74b288748ebd4bfe1f6c3cac43b086bd8f',
 'R059_full_audit':'c8aff39bf68749d25f9978cef11be82ccb7907ec7e72c5678ac66b06a65a33d4',
 'R056_local_lemma':'97430af9980c2de5d105bc4ea86667573c5c28e0ff40a92f1c613cd95a243eec',
 'R056_final_source':'6986ed9ec48b46f3136a9a56722d0092db4f2942d52e35c9f4de37eb279f399a',
}
hashes={k:hashlib.sha256(p.read_bytes()).hexdigest() for k,p in paths.items()}
assert all(hashes[k]==v for k,v in expected.items())
old59=json.loads(paths['R059_historical_certificate'].read_text())
for key,value in {'d':d,'a_band':ab,'b_band':bb,'chi':chi,'sigma':sig,
 'inset':e,'h_coefficient':hc,'epsilon_J':e,'epsilon_N':e,'epsilon_S':e,
 'epsilon_T':e,'epsilon_D':e,'epsilon_Ret':eret,'epsilon_A':ea,'baseline_log_tolerance':e}.items():
    assert Q(old59['parameters'][key])==value
assert Q(old59['derived_rationals']['gamma_chosen'])==gamma
g56=Q(json.loads(paths['R056_historical_values'].read_text())['gamma'])
comparison={
 'accepted_gamma56_below_coarse_upper':g56<Q(1,8),
 'positive_B56_buffer_separates_full_costs':Q(1,8)-Q(7,100)+Q(4,2**20)<Q(1,10),
 'actual_domain_intersection':L>=2**20,
}
assert all(comparison.values())

receipt={'status':'PASS_INDEPENDENT_COMPRESSED_R063_THRESHOLD',
 'python':platform.python_version(),'arithmetic':'standard-library Fraction/integer',
 'threshold':'log(n)>=2^33, integer dimension n; this is new R063 domain, not a retrospective R059 receipt edit',
 'recipe':'nbar=L^8/2^16; source Taylor10; tau Taylor6; scalar Taylor4; no author imports or log/grid/baseline replay',
 'checks':checks,'comparison_checks':comparison,
 'base_values':{k:str(v) for k,v in {'L0':L,'nbar0':nb,'x0':x,'h0':h,
 'baseline_mean_bound':base_mean,'baseline_variance_bound':base_var,
 'inverse_face_tail':tail_per_face,'cap10_bound':cap10,'tau6_bound':tau6,
 'scalar4_bound':scalar4,'scalar_slack':slack,'cap_effective_gap_lower':gap-eta/nb,
 'gamma56_accepted_readonly':g56}.items()},
 'tail_powers':{'nbar':'L^8','nbar_over_L':'L^7','h':'L^-1','zeta_t_over_h':'L^-16',
 'baseline_mean':'positive sum L^-6,L^-7','baseline_variance':'L^-6','Kw_counts':'L^-8',
 'count_drift':'positive sum L^-1,L^-8','individual_J_error':'L^-7',
 'inverse_face_tail':'L^-1','source_cap_T10':'L^-9','tau_T6':'L^-9/2',
 'scalar_T4':'positive sum L^-8,L^-4,L^-7/2'},
 'sources':{k:{'path':str(paths[k]),'sha256':v} for k,v in hashes.items()},
 'own_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'scope':'Exact base inequalities and recorded all-tail analytic powers. The accompanying review supplies their mathematical implication and actual-history/general-L1 bridge. No author verifier, old log certificate, old gates, production scan, MATLAB or Git operation was executed.'}
(HERE/'compressed_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':receipt['status'],'check_groups':len(checks),
 'comparison_groups':len(comparison),'inverse_face_tail':str(tail_per_face),
 'R063_source_sha256':hashes['R063_source']},indent=2))
