#!/usr/bin/env python3
"""Finite exact symbolic prototype for ONE C130 cone/pair; no k_h evaluation."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,os,sys,time,traceback
sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent
S=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper/round130/true_second_count_boxes')
PAT=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-brascamp108/cube-dim-6-bounds/coordination/terra_20261027/pair_pushforward_patterns/patterns.json')
read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(n,d):
 with(P/n).open('x')as f:json.dump(d,f,indent=2);f.write('\n')
sourcepins={str(S/n):sha(S/n)for n in ['box_core.py','cone_table.json','proof.md','inputs.json']}
assert sourcepins[str(S/'cone_table.json')]=='86f81fd6c69d8c9a6b7194d34201732c528ac4cde320aefbb19570d854097040'
assert sha(PAT)=='01a05b2633dfdc7994227b24618847e4b1e7da98b0cde0704acf72ec9753876a'
cones=read(S/'cone_table.json')['cones'];c=cones[0];assert c['id']==0
b=list(map(F,c['nodes_by_body'][16]));d=list(map(F,c['nodes_by_body'][17]));assert b==[F(0)]*4+[F(1,60),F(1,30)] and d==[F(0)]*4+[F(1,60),F(0)]
write('finite_registration.json',{'status':'REGISTERED_BEFORE_SYMBOLIC_PROTOTYPE','source_hashes':sourcepins,'terra_pattern_path':str(PAT),'terra_pattern_sha256':sha(PAT),'real_cones':1,'cone_id':0,'real_body_pairs':1,'body_pair_zero_based':[16,17],'source_b':list(map(str,b)),'source_c':list(map(str,d)),'symbolic_jobs':{'two_by_two_linear_change':1,'joint_Laplace_integral':1,'normalization_integral':1,'first_moment_integrals':2,'mixed_moment_integrals':1,'generic_antiderivative_monomial_degrees':list(range(5))},'new_k_h_or_special_function_evaluations':0,'old_6D_grid_evaluations':0,'optimizer_native_calls':0,'workers':1,'nice_minimum':10,'time_limit':None,'script_sha256':sha(__file__)})
os.nice(max(0,10-os.getpriority(os.PRIO_PROCESS,0)));start=time.monotonic();cpu=time.process_time()
try:
 import sympy as sp
 a,z,s,t=sp.symbols('a z s t',positive=True);r=sp.symbols('r',real=True);eta=sp.symbols('eta',nonzero=True)
 M=sp.Matrix([[sp.Rational(b[j].numerator,b[j].denominator)for j in (4,5)],[sp.Rational(d[j].numerator,d[j].denominator)for j in (4,5)]])
 assert M.det()==-sp.Rational(1,1800);inverse=M.inv()*sp.Matrix([a,z]);assert inverse==sp.Matrix([60*z,30*a-30*z])
 density=sp.simplify(sp.exp(-sum(inverse))/abs(M.det()));assert density==1800*sp.exp(-30*a-30*z)
 def integral(poly):return sp.integrate(sp.integrate(poly*density,(a,z,sp.oo)),(z,0,sp.oo))
 norm=integral(1);EA=integral(a);EB=integral(z);EAB=integral(a*z);cov=sp.simplify(EAB-EA*EB)
 assert norm==1 and EA==sp.Rational(1,20) and EB==sp.Rational(1,60) and EAB==sp.Rational(1,900) and cov==sp.Rational(1,3600)>0
 LT=integral(sp.exp(-s*a-t*z));product=sp.prod(1/(1+sp.Rational(bb.numerator,bb.denominator)*s+sp.Rational(cc.numerator,cc.denominator)*t)for bb,cc in zip(b,d))
 assert sp.simplify(LT-product)==0 and sp.simplify(LT-1800/((s+30)*(s+t+60)))==0
 # Verify finite integration-by-parts formula used in the general proof, incl resonance.
 checks=[]
 for degree in range(5):
  poly=r**degree;Q=sum((-1)**k*sp.diff(poly,r,k)/eta**(k+1)for k in range(degree+1))
  assert sp.simplify(sp.diff(sp.exp(eta*r)*Q,r)-sp.exp(eta*r)*poly)==0
  zero=sp.Rational(1,degree+1)*r**(degree+1);assert sp.diff(zero,r)==poly;checks.append(degree)
 for path,h in sourcepins.items():assert sha(path)==h
 out={'status':'PASS_ONE_REAL_PAIR_EXACT_POSITIVE_2D_PUSHFORWARD','sympy_version':sp.__version__,'source_hashes':sourcepins,'registration_sha256':sha(P/'finite_registration.json'),'cone_id':0,'body_pair_zero_based':[16,17],'alpha':c['alpha'],'four_zero_columns_integrated_out_exactly':True,'active_column_matrix':[[str(v)for v in row]for row in M.tolist()],'determinant':str(M.det()),'inverse_active_exponentials':list(map(str,inverse)),'density_in_open_support':str(density),'support':'0 < B < A; boundary values immaterial','joint_Laplace':str(sp.factor(LT)),'product_Laplace_identity_verified':True,'normalization':str(norm),'mean_A':str(EA),'mean_B':str(EB),'mixed_moment_AB':str(EAB),'covariance_AB':str(cov),'A_B_independent':False,'nonnegative_density':True,'kernel_pair_identity':'E[k_h(A) k_h(B)] = 900*(integral_0^infinity exp(-30*u)*k_h(u) du)^2; follows symmetric integrand, not independence','generic_nonresonant_and_resonant_antiderivative_degrees_checked':checks,'new_k_h_or_special_function_evaluations':0,'old_6D_grid_evaluations':0,'native_optimizer_calls':0,'new_C6_bound':False,'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu}
 write('prototype_result.json',out);write('exit_receipt.json',{'status':out['status'],'exit_code':0,'result_sha256':sha(P/'prototype_result.json'),'wall_seconds':out['wall_seconds'],'CPU_seconds':out['CPU_seconds'],'new_k_h_or_special_function_evaluations':0,'old_6D_grid_evaluations':0,'native_optimizer_calls':0,'time_limit':None});print(json.dumps(out,indent=2))
except BaseException:
 write('exit_receipt.json',{'status':'FAILED_STOP_NO_RETRY','failure':traceback.format_exc(),'new_k_h_evaluations':0});raise
