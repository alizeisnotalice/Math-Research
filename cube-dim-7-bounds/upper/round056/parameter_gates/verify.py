"""One registered proposed tuple: exact constants, finite gates, and local identities."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,math
import sympy as sp

root=Path(__file__).resolve().parent
L0=2**20;d=F(3,10);sigma=F(7,100);chi=F(7,200);e0=F(1,2**24)
e=F(1,1024);eA=F(1,64);P0=F(65275,65536);M=F(21,8)
ql=F(31,40);qu=F(37,40);trim=F(1,2**20)
rl=F(43,40)+trim;ru=F(49,40)-trim
u=(ru-rl)/rl
log_lower=u-u*u/2+u**3/3-u**4/4
mL=(1-eA)*log_lower
ulo=P0*(1-e)/(1+e);uhi=(1+e)/(P0*(1-e))
b=(1-e)*2*ulo/qu-2*e
Gamma=mL*b/M;gamma=(1-e)*Gamma;craw=gamma-F(7,93);c=F(1,40)
z=F(500,3083);Q=z-z*z/2+z**3/3
fac=F(43,31)/(1-e);mr=fac*Q*Q;vr=F(2,3)*fac*Q**3
x0=(1-sigma)*L0;h0=e0/(L0*d);nmin=F(L0**2,2)
Rcoarse=F(10749,10000);tmin=F(L0,2);rhomin=tmin*ql
zt=F(2,L0*nmin*nmin)
mean=math.factorial(8)*(6*d**3/(ql*ql*L0**6)+d/(ql*L0**7))
var=12*math.factorial(8)*d**3/(ql*ql*L0**6)
Jerror=2*d*L0/(ql*nmin)
cap=F(2**44*math.factorial(8),2**96)
tau=F(2**54*math.factorial(16),2**176)
rem=F(1,2**39)+F(math.factorial(8),2**128)+F(2**11*math.factorial(8),2**120)
checks={
 'target_log_domain':0<u<1,
 'target_m_lower_positive':mL>0,
 'M_row_budget':2*uhi/ql<M,
 'M_column_budget':2/Rcoarse<M,
 'row_bracket_positive':b>0,
 'gamma_positive_under_one':0<gamma<1,
 'craw_above_1over40':craw>c,
 'h_inside_trim':h0<trim,
 'h_R_coarse':h0<F(1,10000),
 'h_small_gap':h0<d/8,
 'zeta_t_small':zt<h0/2,
 'reference_R_coarse':ql+d-h0>Rcoarse>1,
 'reference_g_range':0<d-h0 and d-h0+zt<d<Rcoarse/2,
 'kappa_positive_gap':zt<1,
 'rho_ge_4target':rhomin>4*ru,
 'target_below_5over4':ru<F(5,4),
 'matched_source_bounds':rl-d-h0>=ql and ru-d+h0<=qu,
 'baseline_target_J_domain':d/rhomin<F(1,8),
 'baseline_mean':mean<=F(1,512),
 'baseline_variance':var<=F(1,2**26),
 'baseline_failure':F(1,2**26)/F(1,512)**2==eA/4,
 'two_count_failures':nmin>=2**30 and 4/(e**2*nmin)<=eA/4,
 'count_drift_gate':2*h0+1/nmin<=3*e/4,
 'count_upper_drift':2*h0*d+e/4<=e,
 'target_cap_coordinate':1/ru-F(2,L0)>F(1,2),
 'target_cap_probability':nmin>=33 and F(3,4)**32<eA/4,
 'bridge_tail_deterministic':1<=4*d,
 'bridge_product_P0':1-F(1,256)-1280*e0==P0,
 'individual_J_error':Jerror<=e,
 'pin_ratio_coarse':2*Rcoarse/(2*Rcoarse-d)==F(3583,3083),
 'pin_Q_upper':Q==F(13226208500,87910718361)<F(1,6),
 'mean_rate':mr<F(4,125),
 'variance_rate':vr<F(1,300),
 'R_over_q':1+d/ql==F(43,31),
 'pin_budget':x0>=1000 and Q/x0<F(1,2000),
 'correct_Cheb_gap':chi-F(4,125)-F(1,2000)==F(1,400),
 'retention_probability':F(1600)/(3*x0)<e,
 'source_cap_n75':nmin>=75 and F(74,75)*F(3,74)==F(1,25),
 'source_exponent_difference':F(1,25)-chi==F(1,200),
 'source_prefactor':d/h0<2**44,
 'source_argument':x0/200>2**12,
 'source_cap_epsS':cap<e,
 'tau_exponent':sigma/2-chi*(1-sigma)==F(49,20000)>0,
 'tau_prefactor':F(2**10)/h0<2**54,
 'tau_argument':F(49,20000)*L0>2**11,
 'tau_epsT':tau<e,
 'scalar_argument':sigma*L0>2**16,
 'scalar_slack':rem<craw-c,
 'epsilonD_margin':Gamma-gamma==e*Gamma>0,
}
assert all(checks.values()),[k for k,v in checks.items() if not v]
# Local symbolic equalities, not a numerical integration batch.
w,R,rho,Qsym=sp.symbols('w R rho Q',positive=True)
assert sp.integrate(2*R*w/rho,(w,0,Qsym))==R*Qsym**2/rho
assert sp.integrate(2*R*w*w/rho,(w,0,Qsym))==2*R*Qsym**3/(3*rho)
r,hs,ds,zz=sp.symbols('r h d z')
qm=r-ds-hs;qp=r-ds+hs
assert sp.factor((1-zz)*qp/r-(1-zz)*qm/(r-2*hs)-2*hs*(1-zz)*(ds-hs)/(r*(r-2*hs)))==0
sources={}
for name,path in {
 'inverse_coordinate_formula_only':Path('/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round045/soft_inverse_column/PROOF.md'),
 'old_generic_formula_only':Path('/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round051/flexible_budget/PROOF.md')}.items():
 sources[name]={'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'scope':'formula source only; old operator parameter domain not extended'}
contract={'d':str(d),'sigma':str(sigma),'chi':str(chi),'L0':str(L0),'P0':str(P0),'epsilon':str(e),'epsilon_A':str(eA),'rlo':str(rl),'rhi':str(ru),'qmin':str(ql),'qmax':str(qu),'M':str(M)}
receipt={
 'status':'PASS: proposed fixed-tuple constants and finite gates only; no new Cn claimed',
 'input_contract':contract,'input_contract_sha256':hashlib.sha256(json.dumps(contract,sort_keys=True).encode()).hexdigest(),
 'm_lower':str(mL),'Gamma_lower':str(Gamma),'gamma_lower':str(gamma),'c_raw_lower':str(craw),'chosen_safe_c':str(c),
 'c_margin_lower':str(craw-c),'pin_Q_upper':str(Q),'inverse_mean_rate_upper':str(mr),'inverse_variance_rate_upper':str(vr),
 'checks':checks,'symbolic_checks':['conditional inverse first moment','conditional inverse second moment','count upper endpoint drift identity'],
 'sympy_version':sp.__version__,
 'registration_sha256':hashlib.sha256((root/'REGISTRATION.md').read_bytes()).hexdigest(),
 'sources':sources,'new_author_operator_source':'pending; no hash invented',
 'source_sha256':{name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in ('PROOF.md','verify.py')},
 'no_search':True,'no_old_grid_replay':True,'no_operator_activation':True,
}
(root/'verification.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n')
print(json.dumps({'status':receipt['status'],'checks':len(checks),'proof_sha256':receipt['source_sha256']['PROOF.md'],
                  'c_raw_lower':str(craw),'safe_c':'1/40','Cheb_pin_allocation':'pin/x<=1/2000 => gap1/400'},ensure_ascii=False,indent=2))
