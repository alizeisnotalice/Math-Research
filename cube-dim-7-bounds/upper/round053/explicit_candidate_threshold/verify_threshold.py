"""One retained exact tuple and finite threshold gates, no enumeration replay."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,math

root=Path(__file__).resolve().parent
L0=2**32; Lcmp=2**40
A=F(73,8); S=F(293,16); d=F(1023,2**20)
e=F(1,1024); eA=F(1,64); e0=F(1,2**20)
P=F(3579,4096); qmin=1-3*d/4
ul=P/(2*(1+e)); uh=1/(P*(1-d/2)*(1-e))
Mstar=max(F(4),2*uh/qmin)
m=(1-eA)*d/(2*(1+3*d/4)); bracket=(1-3/A)*ul-2*e
Gamma=m*bracket/Mstar; gamma=(1-e)*Gamma
sigma=S*d*d; price=sigma/(1-sigma); craw=gamma-price; c=F(1,2**16)
assert craw==F(6468305130094850733787980147,370780099533215320770594052505600)
compute=Path('/Users/zhengzhihao/Desktop/T/cube-r030-compute/cube-dim-7-bounds/upper/round052/compute')
results=compute/'results.json'; cert=json.loads(results.read_text())['top64'][0]
fields={'a':A,'s':S,'d':d,'u_lo':ul,'u_hi':uh,'mass':m,'M':Mstar,
        'bracket':bracket,'delta_lim':Gamma,'delta':gamma,'sigma':sigma,'penalty':price,'c_raw':craw}
for name,val in fields.items():
 assert F(cert[name])==val,name
h0=e0/(L0*d); nmin=F(L0**2,2); tmin=F(L0,2)
rmin=1+d/4; rmax=1+3*d/4; qb=1-3*d/4-h0
Rmin=rmin-2*h0; rhomin=tmin*qb
zt=F(2,L0*nmin*nmin)
mean=math.factorial(8)*(24*d**3/F(L0**6)+4*d/F(L0**7))
var=48*math.factorial(8)*d**3/F(L0**6)
Jerr=4*d*math.factorial(8)/F(L0**7)
eta=d/32-A*d*d; nu=(S/2-A)*d*d
cap=F(24*2**20)*d*d/(eta**4*L0**3)
tau=F(2**58*math.factorial(64),2**384)
scalar=F(1,2**63)+F(math.factorial(16),2**256)+F(2**17*math.factorial(16),2**240)
checks={
 'domain_a_s':A>3 and S>2*A,
 'domain_d_sigma':0<d<=F(1,16) and d<F(1,32)/A and 0<sigma<F(1,2),
 'all_error_domains':0<e<1 and 0<eA<1,
 'positive_bracket':bracket>0,
 'exact_M_branch':Mstar==4,
 'craw_above_chosen_c':craw>c,
 'h_over_d':h0/d<F(1,8),
 'zeta_t_over_h':zt/h0<F(1,2),
 'shift_below_half_d':h0+zt<d/2,
 'q_bridge_lower':qb>F(1,2),
 'q0_upper':1-d/4+h0<1,
 'R_bridge':Rmin>1,
 'rho_geometry':rhomin>4*rmax and rmax<2,
 'g_range':0<d-h0 and d-h0+zt<=d,
 'baseline_abs_mean':mean<=F(1,32),
 'tight_baseline_variance':var<=F(3,32)**2*eA/3,
 'log_Cheb_probability':(F(3,32)**2*eA/3)/F(3,32)**2==eA/3,
 'Kw_Cheb_probability':nmin*d>=2 and F(8)/(9*nmin*d)<=eA/3,
 'target_cap_probability':nmin>=33 and F(3,4)**32<eA/3,
 'target_cap_geometry':1/rmax-F(2,L0)>F(1,2),
 'individual_J_precision':Jerr<=e,
 'product_exact_P0':1-F(1,8)-1280*e0==P,
 'Ld':L0*d>=1,
 'source_decay_positive':eta>0,
 'source_cap_epsS':cap<e,
 'tau_decay_exact':nu==d*d/32,
 'tau_base_argument':nu*L0>64,
 'tau_base_prefactor':2**20*d*2**48<2**58,
 'tau_epsT':tau<e,
 'fixed_epsD_margin':Gamma-gamma==e*Gamma>0,
 'budget_positivity':0<gamma<1 and L0-32-F(4,L0)>8,
 'scalar_base_argument':sigma*L0>2**16,
 'full_scalar_slack':scalar<craw-c,
 'comparison_new_domain':Lcmp>=L0,
 'comparison_old_delta42':Lcmp>=2**27,
 'comparison_old_R049_R042':F(5,2**40)<c,
 'comparison_old_sqrt_log_base':40<7**2 and 40*F(2,3)>1,
 'comparison_old_q_full_cost':F(1,2**38)+F(14,2**20)<c,
}
assert all(checks.values()),[k for k,v in checks.items() if not v]
source=Path('/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round051/flexible_budget/PROOF.md')
numaudit=Path('/Users/zhengzhihao/Desktop/T/cube-r030-num-audit/cube-dim-7-bounds/audit/round052/README.md')
paths={'R051_generic_pending':source,'R052_results':results,'R052_retained_check':compute/'verify_top64.json','R052_independent_audit':numaudit}
sources={k:{'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for k,p in paths.items()}
assert sources['R051_generic_pending']['sha256']=='aeddd67b7a2109d402ea76d9f6fb3d93accf44bea80240340bd163759a7c9df6'
receipt={
 'status':'PASS: exact fixed-tuple and threshold gates; operator conditional on R051 full ACCEPT',
 'operator_status':'conditional draft; generic R051 full chain acceptance required',
 'candidate_threshold':'log(n)>=2^32',
 'comparison_formula_domain':'log(n)>=2^40; each old operator keeps its own valid domain',
 'exact_fields':{k:str(v) for k,v in fields.items()},
 'epsilon_A':str(eA),'other_epsilons':str(e),'chosen_c':str(c),
 'checks':checks,
 'sources':sources,
 'source_sha256':{name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in ('THRESHOLD_PROOF.md','verify_threshold.py')},
 'no_grid_replay':True,'no_production_computation':True,
}
(root/'verification.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'status':receipt['status'],'checks':len(checks),
                  'proof_sha256':receipt['source_sha256']['THRESHOLD_PROOF.md'],
                  'c_raw':str(craw),'threshold':'L>=2^32 conditional'},ensure_ascii=False,indent=2))
