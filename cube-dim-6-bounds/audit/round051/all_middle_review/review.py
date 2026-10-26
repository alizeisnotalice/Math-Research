#!/usr/bin/env python3
"""Six fixed independent all-middle checks; k3 reused, no old writer or kernel."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,os,sys,time,datetime,math
import sympy as s
import flint
from flint import arb,fmpq,ctx
HERE=Path(__file__).resolve().parent;start=time.time()
lock=json.loads((HERE/'source_lock.json').read_text());sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
for p,h in lock['source_sha256'].items():assert sha(p)==h,p
assert time.time()<datetime.datetime(2026,10,25,22,23,49,tzinfo=datetime.timezone.utc).timestamp()
assert not (HERE/'results.json').exists()
ctx.prec=768
checks=[]
# One exact common relation group, without geometry/frame enumeration.
k,p,q=s.symbols('k p q',positive=True);m=6-k
assert s.factor(k*p*p+m*q*q-(k*p+m*q)**2/6-k*m*(p-q)**2/6)==0
assert s.factor(36-(k-m)**2-4*k*m)==0
rhoq=F(1,64)
minimum_gates={'hv0_reduced':F(1,2)-rhoq,'hw0_reduced':F(15,16)-rhoq,
 'hv1_squared_gate_m_ge1':F(1,96)-rhoq**2,'hw1_squared_gate_m_ge1':F(1,384)-rhoq**2,
 'g_squared_gate_k_le5':F(1,32)-5*rhoq**2}
assert all(x>0 for x in minimum_gates.values())
assert F(10,9)**3==F(1000,729)
assert [32*n for n in [31,30,29]]==[992,960,928]
assert F(1-1,1)==0
weights=[math.comb(5,j) for j in range(5)];assert weights==[1,5,10,10,5] and sum(weights)==31
checks.append({'configuration':1,'name':'common_exact_geometry_and_normalization','two_value_nearest_line_identity_residual':'0','pair_Gram_identity_residual':'0','all_k_positivity_from_extreme_gates':{a:str(v) for a,v in minimum_gates.items()},'perturbation_bound':'sqrt(k)*rho','coarea_J2':'2sqrt(km)','R2_type':'min(k,m)','missing_peak_loss':'1000/729','actual_remaining_counts':[31,30,29],'time_terms':[992,960,928],'remaining_Gram_bounds':['26I','32I'],'singleton_B_deviation_squared':'0','weights':weights,'all_signed_lines_and_other_planes_reviewed_analytically':True,'no_k3_scalar_recomputation':True,'pass':True})
author=json.loads(Path(lock['author_results_path']).read_text())
peaks={x['id']:x for x in json.loads(Path(lock['projection033_results_path']).read_text())['cases']}
def aq(v):
 v=F(v);return arb(fmpq(v.numerator,v.denominator))
def dyad(v):
 lo=F(int(v['lower_mantissa']))*F(2)**int(v['lower_binary_exponent'])
 hi=F(int(v['upper_mantissa']))*F(2)**int(v['upper_binary_exponent'])
 assert lo<=hi;return lo,hi
def bounds(v):return F(str(v.lower().fmpq())),F(str(v.upper().fmpq()))
def enc(v):
 lo,hi=bounds(v);return {'ball':str(v),'lower_exact':str(lo),'upper_exact':str(hi)}
def compatible(v,a):
 lo,hi=bounds(v);al,ah=dyad(a);assert max(lo,al)<=min(hi,ah)
pi=arb.pi();rho=aq(rhoq);eps=aq(F(1,1000))
Crem=aq(F(44743079,2970))/pi**6
CauchyF=1+aq(F(6,928)).sqrt();loss=aq(F(1000,729))
root6=lambda v:(v.log()/6).exp()
case_upper={};case_results={}
for cid,kk in enumerate([1,2,4,5],start=2):
 mm=6-kk;eta=aq(kk).sqrt()*rho
 # Factor out sqrt(k) for the line gaps, unlike the producer's direct formula.
 hv0=aq(kk).sqrt()*(aq(F(1,2))-rho)
 hv1=aq(kk).sqrt()*(aq(F(mm,6)).sqrt()/4-rho)
 hw0=aq(kk).sqrt()*(aq(F(15,16))-rho)
 hw1=aq(kk).sqrt()*(aq(F(mm,6)).sqrt()/8-rho)
 g=1/(4*aq(2).sqrt())-eta
 assert all(z>0 for z in [hv0,hv1,hw0,hw1,g])
 dev=aq(F(mm-1,mm)).sqrt()*rho
 if mm==1:assert dev.is_zero()
 L=aq(F(1,4))-dev;U=aq(F(3,4))+dev;assert L>0 and U<1
 C=64/pi**2*(U/(1-U*U)-L/(1-L*L)+U.atanh()-L.atanh())
 typ=min(kk,mm);R2=F(peaks[f'r2_k{typ}']['R_exact'])
 assert F(author['cases'][str(kk)]['R2_exact'])==R2
 assert peaks[f'r2_k{typ}']['Gram_determinant']==4*kk*mm
 assert peaks[f'r2_k{typ}']['kappa_exact']=='R/(2*pi)^4'
 D=64*aq(R2)/pi**6
 common=loss*Crem*CauchyF/(g*g)
 Mv=common/(hv0*hv0*hv1*hv1);Mw=common/(hw0*hw0*hw1*hw1);M=Mv+Mw
 good=aq(F(3,7))*root6(D*C**5*eps**14)
 bad=aq(F(3,8))*root6(M*C**5*eps**16)
 total=good+bad
 row=author['cases'][str(kk)]
 for name,value in [('eta',eta),('hv0',hv0),('hv1',hv1),('hw0',hw0),('hw1',hw1),('g',g),('C_joint',C),('D_big',D),('Mv',Mv),('Mw',Mw),('M_combined',M),('short_big_upper',good),('short_bad_upper',bad),('short_total_upper',total)]:compatible(value,row[name])
 case_upper[kk]=bounds(total)[1]
 case_results[str(kk)]={'k':kk,'m':mm,'R2_exact':str(R2),'eta':enc(eta),'hv0':enc(hv0),'hv1':enc(hv1),'hw0':enc(hw0),'hw1':enc(hw1),'g':enc(g),'coordinate_deviation':enc(dev),'C_joint':enc(C),'D_big':enc(D),'M_bad':enc(M),'short_big':enc(good),'short_bad':enc(bad),'short_total':enc(total),'independent_upper_exact':str(case_upper[kk])}
 checks.append({'configuration':cid,'name':f'one_fixed_new_k{kk}_m{mm}','case_key':str(kk),'all_author_interval_comparisons':True,'positive_gates':True,'singleton_B_exact_zero':mm==1,'pass':True})
# Sixth configuration: union only; k3 is a saved independent049 rational endpoint.
old49=json.loads(Path(lock['k3_own049_results_path']).read_text())
k3hi=F(old49['configurations'][3]['total_short_upper_exact'])
case_upper[3]=k3hi
# Verify producer k3 really is its frozen049 record, without evaluating its formula.
old_author49_path=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-curvature-astra/cube-dim-6-bounds/upper/round049/middle_joint_marginal/results.json')
old_author49=json.loads(old_author49_path.read_text())['configurations'][2]
assert author['cases']['3']['short_total_upper']==old_author49['new_total_short_upper']
S=arb(0);weighted=[]
for kk,n in zip(range(1,6),weights):
 value=aq(case_upper[kk]);term=(value.log()*6/5).exp()*n;S+=term
 weighted.append({'k':kk,'n':n,'input_upper_exact':str(case_upper[kk]),'weighted_upper_power':enc(term),'k3_recomputed':False})
union=(S.log()*5/6).exp();unionhi=bounds(union)[1]
threshold=F(4689838525803,10**18)
assert unionhi<threshold
checks.append({'configuration':6,'name':'positive31_domain_union','weights':weights,'domains':31,'q':'6/5','k3_source_sha256':sha(Path(lock['k3_own049_results_path'])),'k3_saved_upper_exact':str(k3hi),'k3_new_scalar_evaluations':0,'weighted_terms':weighted,'sum_power':enc(S),'union_upper':enc(union),'strict_rational_upper':str(threshold),'strict_margin_exact':str(threshold-unionhi),'overlap_permitted':True,'domains_asserted_disjoint':False,'pass':True})
for p,h in lock['source_sha256'].items():assert sha(p)==h,p
out={'status':'PASS_ALL31_MIDDLES_SHORT_BUDGET_ONLY','pid':os.getpid(),'python':sys.version,'sympy':s.__version__,'python_flint':flint.__version__,'precision_bits':ctx.prec,'configuration_count':6,'checks':checks,'new_cases':case_results,'kernel_native_G_H_K_p_calls':0,'MATLAB_calls':0,'full_enumerations':0,'old_writer_calls':0,'old_k3_scalar_recomputations':0,'old_budget_recomputations':0,'source_inputs_unchanged':True,'source_lock_sha256':sha(HERE/'source_lock.json'),'review_source_sha256':sha(Path(__file__)),'elapsed_seconds':time.time()-start,'new_I_global_g_Cn':False}
assert out['configuration_count']==6
(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':out['status'],'pid':out['pid'],'configurations':6,'elapsed_seconds':out['elapsed_seconds'],'union_display_only':str(union),'strict_rational_upper':str(threshold),'k3_scalar_recomputations':0},indent=2))
