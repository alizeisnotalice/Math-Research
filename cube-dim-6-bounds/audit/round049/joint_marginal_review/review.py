#!/usr/bin/env python3
"""A049 four fixed independent checks after final source lock; no author writer."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib,time,os,sys
import sympy as s
import flint
from flint import arb,fmpq,ctx
HERE=Path(__file__).resolve().parent;start=time.time()
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
lock=json.loads((HERE/'source_lock.json').read_text())
assert start>=lock['audit_started_unix'] and start-lock['audit_started_unix']<720
for p,h in lock['source_sha256'].items():assert sha(p)==h,p
assert not (HERE/'results.json').exists()
ctx.prec=768
checks=[]
# Configuration1: a retained-coordinate pivot plus an independent Schur check.
v=s.ones(6,1)
flip=lambda j:s.Matrix([-1 if i==j else 1 for i in range(6)])
B=s.Matrix.hstack(v,flip(3),flip(1),flip(2),flip(4),flip(5))
rows=[0,3,1,2,4,5];C=B[rows,:]
P=C[:2,:2];Q=C[:2,2:];R=C[2:,:2];S=C[2:,2:]
Schur=S-R*P.inv()*Q
assert abs(B.det())==32 and abs(P.det())==2 and abs(Schur.det())==16
j=s.Matrix([[1,1],[1,-1]]);assert abs(j.det())==2
z,t=s.symbols('z t',real=True,positive=True)
f=lambda a:16*t/(s.pi*((16*t)**2+a*a))
actual=f((1+z)/2)*f((1-z)/2)/2
formula=2048*t*t/(s.pi**2*((1+z)**2+1024*t*t)*((1-z)**2+1024*t*t))
assert s.factor(actual-formula)==0
checks.append({'configuration':1,'retained_j':4,'block_counts':[2**4,2**4],'body_frame_enumeration':False,'basis_abs_det':32,'pivot_abs_det':2,'Schur_abs_det':16,'inverse_joint_Jacobian':'1/2','remaining_four_positive_orthant_factor':'1/16','density_identity_residual':'0','reflection_argument':'independent sign-bit flips outside coordinates1,4 permute the body column set; applied only to full p','pass':True})
# Configuration2: independent rational/log primitive and exact exponents.
P3=s.eye(3)-s.ones(3)/3;e=s.Matrix([1,0,0]);w=P3*e
assert (w.T*w)[0]==s.Rational(2,3)
primitive=(z/(1-z*z)+s.atanh(z))/2
assert s.factor(s.diff(primitive,z)-(1-z*z)**(-2))==0
q=s.Rational(6,5);theta=1-1/q
Ppower=-2*theta+2/q;Bpower=2/q
assert theta==s.Rational(1,6) and Ppower==s.Rational(4,3) and Bpower==s.Rational(5,3)
assert 1/(Ppower+1)==s.Rational(3,7) and 1/(Bpower+1)==s.Rational(3,8)
checks.append({'configuration':2,'B_deviation_functional_norm_squared':'2/3','primitive_derivative_residual':'0','interpolation_infinity_exponent':'1/6','interpolation_L1_exponent':'5/6','P_time_power':'4/3','B_time_power':'5/3','integrated_time_powers':['7/3','8/3'],'integrated_prefactors':['3/7','3/8'],'branch_orthant_symmetry_assumed':False,'actual_complement_L1_used_once':True,'pass':True})
# Configuration3: a single fixed joint-marginal constant, using atanh rather than log quotient.
aq=lambda q:arb(fmpq(q.numerator,q.denominator))
def save(a):return {'ball':str(a),'lower_exact':str(a.lower().fmpq()),'upper_exact':str(a.upper().fmpq())}
rho=arb(1)/64;dev=(arb(2)/3).sqrt()*rho
L=arb(1)/4-dev;U=arb(3)/4+dev
assert L>0 and U<1
# F(U)-F(L), evaluated as one difference with a symmetric atanh primitive.
J=(U/(1-U*U)-L/(1-L*L)+U.atanh()-L.atanh())/2
Cjoint=128*J/arb.pi()**2
assert Cjoint>0
checks.append({'configuration':3,'rho':'1/64','fixed_mean_interval':['1/4','3/4'],'fixed_j':4,'coordinate_deviation':save(dev),'L':save(L),'U':save(U),'integral':save(J),'C_joint':save(Cjoint),'fixed_Cjoint_scalar_only':True,'pass':True})
# Configuration4: monotone sixth-root evaluation using saved043 complete endpoints.
def dyadic(a):
 lo=F(int(a['lower_mantissa']))*(F(2)**int(a['lower_binary_exponent']))
 hi=F(int(a['upper_mantissa']))*(F(2)**int(a['upper_binary_exponent']))
 assert 0<lo<=hi
 return lo,hi
oldpath=Path(lock['old043_results_path']);old=json.loads(oldpath.read_text())
Dinterval=dyadic(old['configurations'][3]['C2'])
Minterval=dyadic(old['configurations'][4]['M_combined'])
oldinterval=dyadic(old['configurations'][5]['total_short_norm_upper'])
eps=arb(1)/1000
root6=lambda a:(a.log()/6).exp()
def item_bounds(interval,epower,pref):
 lo=root6(aq(interval[0])*Cjoint**5*eps**epower)*aq(pref)
 hi=root6(aq(interval[1])*Cjoint**5*eps**epower)*aq(pref)
 return F(str(lo.lower().fmpq())),F(str(hi.upper().fmpq()))
pl,pu=item_bounds(Dinterval,14,F(3,7));bl,bu=item_bounds(Minterval,16,F(3,8))
nlo,nup=pl+bl,pu+bu
threshold=F(66002272277,250000000000000000)
assert 0<nlo<=nup<threshold<oldinterval[0]
new=json.loads(Path(lock['author_results_path']).read_text())
newscalar=new['configurations'][2]
for key,lo,hi in [('C_joint',F(str(Cjoint.lower().fmpq())),F(str(Cjoint.upper().fmpq()))),('new_big_short_upper',pl,pu),('new_bad_short_upper',bl,bu),('new_total_short_upper',nlo,nup)]:
 al,au=dyadic(newscalar[key]);assert max(al,lo)<=min(au,hi),key
checks.append({'configuration':4,'source_D_big_key':'configurations[3].C2','source_M_bad_key':'configurations[4].M_combined','source_old_budget_key':'configurations[5].total_short_norm_upper','source_D_big_endpoints':list(map(str,Dinterval)),'source_M_bad_endpoints':list(map(str,Minterval)),'source_old_budget_endpoints':list(map(str,oldinterval)),'independent_evaluation':'sixthroot(D*Cjoint^5*epsilon^14)*3/7 and sixthroot(M*Cjoint^5*epsilon^16)*3/8, monotone source endpoints','big_short_lower_exact':str(pl),'big_short_upper_exact':str(pu),'bad_short_lower_exact':str(bl),'bad_short_upper_exact':str(bu),'total_short_lower_exact':str(nlo),'total_short_upper_exact':str(nup),'strict_rational_upper':str(threshold),'strict_margin_exact':str(threshold-nup),'new_upper_below_saved_old_lower':True,'author_interval_comparisons_pass':True,'old043_budget_recomputed':False,'actual_norm_reduction_factor_claimed':False,'pass':True})
for p,h in lock['source_sha256'].items():assert sha(p)==h,p
out={'status':'PASS_JOINT_MARGINAL_ACTUAL_BRANCH_SHORT_BUDGET','pid':os.getpid(),'python':sys.version,'sympy':s.__version__,'python_flint':flint.__version__,'precision_bits':ctx.prec,'configurations':checks,'configuration_count':4,'source_inputs_unchanged':True,'source_sha256':sha(Path(__file__)),'source_lock_sha256':sha(HERE/'source_lock.json'),'elapsed_seconds':time.time()-start,'elapsed_since_source_lock':time.time()-lock['audit_started_unix'],'old_writer_calls':0,'old043_budget_recomputed':False,'frame_fan_triple_enumerations':0,'kernel_native_G_H_K_p_calls':0,'MATLAB_calls':0,'new_I_certificate':False,'new_global_g':False,'new_C6_bound':False}
assert out['elapsed_since_source_lock']<720
(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':out['status'],'pid':out['pid'],'configurations':4,'elapsed_seconds':out['elapsed_seconds'],'elapsed_since_lock':out['elapsed_since_source_lock'],'Cjoint':str(Cjoint),'new_total_display_only':float(nup),'strict_rational_upper':str(threshold)},indent=2))
