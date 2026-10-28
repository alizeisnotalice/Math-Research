"""One fixed geometry, 12 registered payment rows. All BE/H2 use is conditional."""
from fractions import Fraction as F
from pathlib import Path
import datetime,hashlib,json,math,os,sys,time
started=time.perf_counter(); root=Path(__file__).resolve().parent
L=4096; sqrtL=64; fac8=math.factorial(8); N=F(L**8,fac8)
s=F(3,2); ql=F(2,3); qu=F(91,100); sigma=F(19,200); chi=F(9,200)
x=(1-sigma)*L; delta=F(1,2**19); rl=1+delta; ru=s*qu-delta
hc=F(1,512); h=hc/L; eB=F(1,2**16); eJ=F(1,2**20); eN=F(1,2**19); eR=eJ; eA=F(1,65536); eD=eJ; eRet=F(1,3)
P0=1-eB-6*hc; ulo=P0*(1-eN)/(1+eJ); uhi=(1+eN)/(P0*(1-eJ)); M=2*(1+eR)*uhi; mL=F(311,1000)
eta=(1-qu)/(2*qu); Csrc=(qu-ql)/ql; cS=(1-sigma)/225; nu=sigma/2-chi*(1-sigma)

def prefix64(z):
 terms=[F(1)]; term=F(1)
 for k in range(1,65):
  term*=z/k; terms.append(term)
 total=sum(terms,F()); mu=sum((k*t for k,t in enumerate(terms)),F())/total
 return total,mu
PS,muS=prefix64(cS*L); PT,muT=prefix64(nu*L)
cap=(Csrc*ru/h)/PS; tau=(ru*sqrtL/h)/PT
# New current tuple target radial log enclosure and new price, not old certificates.
ulog=(ru-rl)/rl
mass20=sum(((-1)**(k+1)*ulog**k/k for k in range(1,21)),F())*(1-eA)
uprice=sigma/(1-sigma)
price7=sum(((-1)**(k+1)*uprice**k/k for k in range(1,8)),F())
price_upper=F(999,10000)
mean=(F(183,1000)+F(2211,50000)/(1-eJ)*x)
gap=chi*x-mean; varlo=F(13,2500)*x; varhi=F(531,100000)*x
normal=F(1,2)+F(199,500)*(F(3,5)-F(3,5)**3/6)
be_error=F(21,500); kept=normal-be_error
Mlog=fac8*(6*F(1,2)**3/(ql*ql*L**6)+F(1,2)/(ql*L**7))
Vlog=12*fac8*F(1,2)**3/(ql*ql*L**6)
Rpower=2/F(L*L)+F(fac8)/(sigma*L)**8+2*sqrtL*F(fac8)/(sigma*L/2)**8
common={
 'fixed_domain':0<sigma<F(1,2) and F(1,2)<ql<qu<1 and 1<rl<ru<F(3,2),
 'root_positive_scalar_domain':L>=64 and F(L,4)-F(4,L)-1>0,
 'h_trim':h<delta/2,
 'active_R':rl-2*h>1,
 'relative_g_guard':qu/(x*N*N)<2*h/3,
 'rho_guard':ql*N/L>4*ru,
 'beta_kappa':N/L>16 and 1/(x*N*N)<1,
 'baseline_mean':Mlog<eB/2,
 'baseline_failure':Vlog<(eB/2)**2*eA/8,
 'minus_count':4/(eN*eN*N)<eA/8,
 'plus_count':4/(eN*eN*N)<eA/8,
 'count_drift':2*h+1/N<3*eN/4,
 'target_cap_coordinate':1/ru-L/N>F(1,2),
 'target_cap_probability':N>=65 and F(3,4)**64<eA/8,
 'profile1_pin':N*ql>=200000,
 'profile1_failure':F(200000**2)/(ql*ql*N)<eA/8,
 'profile1_transport':16*h<=F(1,100000),
 'H2_interface_h':h<=F(1,2**20),
 'H2_interface_n':N>=2000 and x/N<=F(1,1000),
 'H2_interface_bad_probability':F(5760000)/N<eA/8,
 'H2_interface_transport':h/5<F(1,1000000),
 'six_bad_event_budget':6*eA/8<eA,
 'product_bridge':6*L*h==6*hc and 0<P0<1,
 'individual_J':F(L)/(3*ql*N)<eJ,
 'relative_radial_width':1/(1-(h/rl)**2)<1+eR,
 'target_mass_Taylor20_lower':0<ulog<1 and mass20>mL,
 'mean_gap_positive':gap>0,
 'standardized_z_at_least_three_fifths':gap*gap>=F(9,25)*varhi,
 'source_variance_at_least_19_assuming_H2':varlo>=19,
 'BE_error_assuming_constant_at_most_one':F(183,1000)**2<19*be_error**2,
 'pi_density_lower':F(199,500)**2<F(7,44),
 'normal_integral_identity':normal==F(90559,125000),
 'normal_after_BE_retention':kept>F(2,3)==1-eRet,
 'price_Taylor7':0<uprice<1 and price7<price_upper,
 'cap_exponent':eta==F(9,182) and eta-chi-eta/N>F(1,225),
 'cap_n10000_suffices':eta-chi-eta/10000>F(1,225),
 'cap_rate_identity':cS==F(181,45000),
 'tau_rate_identity':nu==F(271,40000)>0,
 'cap_prefix_monotonicity':muS>1,
 'tau_prefix_monotonicity':muT>F(3,2),
 'full_scalar_small':Rpower<F(1,10000),
 'finite_a_guard':8*N**4>2*N and 1/N**3-F(4)/(8*N**4)>0,
}
# Analytic H2 and BE theorem hypotheses remain open regardless of these booleans.

def row(eS,eT,gamma=None,c=None):
 b=2*(1-eRet)*(1-eR)*ulo-eS-eT
 Gamma=mL*b/M; budget=(1-eD)*Gamma
 if gamma is None:
  scaled=budget*100000; gamma=F(scaled.numerator//scaled.denominator,100000)
  if gamma==budget: gamma-=F(1,100000)
 if c is None: c=gamma-F(1,10)
 slack=gamma-price_upper-c
 gates={'cap_payment':cap<eS,'tau_payment':tau<eT,'energy_budget':0<gamma<budget<1,'scalar_complete':slack>Rpower,'geometry_common_arithmetic':all(common.values())}
 normalized={'cap':(eS-cap)/eS,'tau':(eT-tau)/eT,'energy':(budget-gamma)/budget,'scalar':(slack-Rpower)/slack}
 return {'epsilon_S':str(eS),'epsilon_T':str(eT),'Gamma':str(Gamma),'gamma_budget':str(budget),'gamma':str(gamma),'display_c_conditional':str(c),'exact_price_upper':str(price_upper),'scalar_slack':str(slack),'gates':gates,'passed_arithmetic':all(gates.values()),'normalized_margins':{k:str(v) for k,v in normalized.items()},'minimum_normalized_margin':str(min(normalized.values())),'status':'conditional parameter gates only: R070 H2 final identity and qualified BE theorem/full operator chain not supplied'}
base_row=row(F(3,40),F(1,4096),F(19,100),F(9,100))
rows=[row(F(sv,1000),tv) for sv in [75,74,73,72] for tv in [F(1,4096),F(1,5000),F(1,6000)]]
feasible=[r for r in rows if r['passed_arithmetic']]
best=max(feasible,key=lambda r:F(r['display_c_conditional'])) if feasible else None
hashof=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
srcbase=Path('/Users/zhengzhihao/Desktop/T')
sources={}
for name,path,expected in [
 ('R066_full_proof',srcbase/'cube-r030-geometry/cube-dim-7-bounds/upper/round066/target_profile/PROOF.md','4b7b3210f063312527c08adab5b7f81dc4b2c0ecf879e466dc715eb56ea6ebd1'),
 ('relative_coordinate_lemma',srcbase/'cube-r030-roots-audit/cube-dim-7-bounds/audit/round064/relative_coordinate/LEMMA.md','c60d48d116863c9c9cac177e19a9fa4158f12fd8dcef0d5239a5ad0a5487d9f1')]:
 actual=hashof(path); assert actual==expected
 sources[name]={'path':str(path),'sha256':actual,'scope':'formula/domain source only; original parameter theorem does not certify new tuple or BE step'}
receipt={
 'status':'BASE '+('PASS' if base_row['passed_arithmetic'] else 'FAIL')+' for conditional arithmetic only; no new Cn theorem',
 'actual_pid':os.getpid(),'actual_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'python':sys.version,'wall_seconds':time.perf_counter()-started,
 'registered_rows':12,'passed_rows':len(feasible),'common_checks':common,'common_check_count':len(common),
 'base':base_row,'rows':rows,'best_registered_row':best,
 'constants':{'L0':str(L),'x0':str(x),'h0':str(h),'P0':str(P0),'Ulo':str(ulo),'Uhi':str(uhi),'M':str(M),'target_mass_lower_Taylor20':str(mass20),'price_upper_Taylor7':str(price7),'cap_bound_S64':str(cap),'tau_bound_S64':str(tau),'cap_prefix_log_derivative_at_base':str(muS),'tau_prefix_log_derivative_at_base':str(muT),'mean_upper_assumed_profile':str(mean),'variance_lower_assumed_H2':str(varlo),'variance_upper_assumed_H2':str(varhi),'z_squared_margin':str(gap*gap-F(9,25)*varhi),'normal_lower':str(normal),'normal_minus_BE':str(kept),'scalar_remainder_upper':str(Rpower)},
 'unproved_dependencies':['R070 final H2 lemma/source identity and full new target profile/domain implementation, pending at run time','Qualified independent-array Berry-Esseen theorem with applicable constant C<=1, not yet supplied','Original history/general L1/all-positive-scale operator proof for this exact new tuple; arithmetic cannot replace it'],
 'sources':sources,'R070_H2_source':'pending; communicated draft interface only; no hash invented','BE_source':'missing; conditional inequality, not an accepted theorem',
 'own_sources':{name:hashof(root/name) for name in ['REGISTRATION.md','verify_candidate.py']},
 'no_old_certificate_replay':True,'single_process':True,'no_MATLAB':True,'no_Git':True,'no_operator_activation':True,
}
(root/'verification.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'status':receipt['status'],'actual_pid':os.getpid(),'common_checks':len(common),'failed_common':[k for k,v in common.items() if not v],'passed_rows':len(feasible),'base_gamma':base_row['gamma'],'base_Gamma':base_row['Gamma'],'base_pass':base_row['passed_arithmetic'],'best':{k:best[k] for k in ['epsilon_S','epsilon_T','gamma','display_c_conditional','minimum_normalized_margin']} if best else None,'wall_seconds':receipt['wall_seconds']},indent=2))
