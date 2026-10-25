#!/usr/bin/env python3
"""Fresh fixed-delta C009 review; existing stable peak balls only."""
from pathlib import Path
from time import perf_counter
from fractions import Fraction as Q
import hashlib,json,math,os
import sympy as s
from flint import arb,ctx

HERE=Path(__file__).resolve().parent
C9=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-audit/cube-dim-6-bounds/upper/round009/alpha15_off_strata')
C8=C9.parents[1]/'round008/alpha15_coercivity'
F8=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-literature/cube-dim-6-bounds/audit/round008/alpha15_product')
ctx.prec=384
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def A(q):return arb(q.numerator)/q.denominator
def main():
 start=perf_counter();c9=json.loads((C9/'off_strata_results.json').read_text())
 c8=json.loads((C8/'coercivity_results.json').read_text());f8=json.loads((F8/'check_results.json').read_text())
 assert c9['source_code_sha256']==sha(C9/'off_strata_check.py')
 for path,h in c9['input_file_hashes'].items():assert sha(Path(path))==h
 tailA=arb(f8['final_constants']['A']);M=arb(f8['final_constants']['M'])
 P4=arb(c8['eight_peak_cases']['d4_dimension6']['peak_upper'])
 P5=arb(c8['eight_peak_cases']['d5_dimension6']['peak_upper'])
 b32=tailA*32**(arb(3)/2);b31=tailA*32**(arb(5)/4)*31**(arb(1)/4)
 geometries=[];J2=arb(0);v=s.ones(6,1);I=s.eye(6)
 P1=I-v*v.T/6
 assert P1*(32*I-v*v.T)*P1==32*P1
 for k in range(1,6):
  w=s.Matrix([1]+[-1]*k+[1]*(5-k));B=s.Matrix.hstack(v,w);G=B.T*B;D=G.det()
  assert D==4*k*(6-k)
  projected=w.T*P1*w;assert projected[0]==s.Rational(D,6)
  P2=I-B*G.inv()*B.T
  assert P2*(32*I-B*B.T)*P2==32*P2
  count=32*math.comb(5,k);J2+=count/arb(int(D)).sqrt()
  geometries.append({'Hamming_k':k,'ordered_pair_count':count,'Gram_det':int(D),
   'second_projected_length_squared':str(projected[0]),'pair_coarea':'1/sqrt('+str(D)+')',
   'projected_tight_frame_constant':32,'matrix_residual_zero':True})
 assert sum(row['ordered_pair_count']for row in geometries)==992
 B0=32/arb(6).sqrt()*P5*b32;B2=J2*P4*b32*b31
 prodD=arb(15)/4*B0**(arb(1)/5)*B2**(arb(4)/5)
 # Independent symbolic minimum-time integration, retaining the r factor.
 u,a,b,r,rho=s.symbols('u a b r rho',positive=True)
 tstar=(a/b)**s.Rational(3,5)*rho**s.Rational(3,2)
 integ=r**s.Rational(-5,2)*(3*b*rho**s.Rational(-5,2)*tstar**s.Rational(1,3)+s.Rational(3,4)*a*tstar**s.Rational(-4,3))
 target=s.Rational(15,4)*a**s.Rational(1,5)*b**s.Rational(4,5)*r**s.Rational(-5,2)*rho**(-2)
 assert s.simplify(integ-target)==0
 # Exact face metric and cone geometry. Order probability is analytic, not sampled.
 Dface=s.eye(5)-s.ones(5,5)/6
 assert Dface.det()==s.Rational(1,6)and Dface.inv()==s.eye(5)+s.ones(5,5)
 assert all(Dface.inv()[j,j]==2 for j in range(5))
 zz=s.symbols('z0:5');z=s.Matrix(zz)
 assert s.expand((z.T*Dface*z)[0]-(sum(t*t for t in zz)-sum(zz)**2/6))==0
 mean,var=s.symbols('mean var',real=True)
 assert s.expand(1+5*mean**2+var-(1+5*mean)**2/6-s.Rational(5,6)*(1-mean)**2-var)==0
 assert s.integrate(u**5,(u,0,1))==s.Rational(1,6)
 assert max(s.Rational(k*(6-k),6)for k in range(1,7))==s.Rational(3,2)
 delta=Q(1,1000);assert 2*delta*delta<1
 K=arb(4)*arb(6).sqrt()*arb.pi()**2/45
 R=(1+5*(1-(arb(6)/5).sqrt()*A(delta))**2).sqrt()
 capI=arb(15)/7*K*prodD**(arb(4)/3)*R**(-arb(10)/3)*A(delta)**(arb(7)/3)
 capNorm=capI**(arb(3)/4)
 assert capI<arb(749)/10**10 and capNorm<arb(4525)/10**9
 for key,val in [('b32',b32),('b31',b31),('J2',J2),('B0',B0),('B2',B2),('D_product',prodD)]:
  assert(val-arb(c9['rank2_constants'][key])).contains(0)
 assert(capI-arb(c9['fixed_true_S1_cap']['G_I_upper'])).contains(0)
 assert(capNorm-arb(c9['fixed_true_S1_cap']['G_Lq_upper'])).contains(0)
 inputs=[C9/'report.md',C9/'off_strata_check.py',C9/'off_strata_results.json',C8/'coercivity_results.json',F8/'check_results.json',HERE/'registration.md']
 out={'status':'PASS_FRESH_INDEPENDENT_C009_STABLE_GEOMETRY_REVIEW','pid':os.getpid(),'arb_bits':384,
  'fixed_delta':'1/1000','fixed_cap_cases':1,'canonical_pair_Gram_checks':geometries,
  'new_stable_peak_evaluations':0,'kernel_evaluations':0,'MC_samples':0,'matlab_calls':0,'twelve_layer_budgets_rerun':False,
  'stable_inverse_density':{'frame_constant_after_deletion':32,'N_second_step':31,
   'b31_exact_form':'A*32^(5/4)*31^(1/4)','overlaps_treated_as_positive_sum':True},
  'time_min_integration':{'short_density_exponent':'-2/3','long_density_exponent':'-7/3',
   'G_r_exponent':'-5/2','G_rho_exponent':'-2','exact_symbolic_residual':'0'},
  'face_geometry':{'nearest_S1_line':'all-positive body line','D':'I_5-11T/6','det_D':'1/6',
   'inverse_D':'I_5+11T','common_Gaussian_order_probability':'1/6','not_orthant_reflection_fraction':'1/32',
   'rho_ball_volume_K':'4sqrt(6)*pi^2/45','delta_formula_valid_through':'1/sqrt(2)',
   'coordinate_support_bound':'sqrt(2)*delta','rho_max_squared':'3/2'},
  'constants':{key:val.str(100,more=True)for key,val in [('b32',b32),('b31',b31),('J2',J2),('B0',B0),('B2',B2),('D_product',prodD),('K',K),('R_delta',R)]},
  'full_G_cap_I_upper_arb384':capI.str(100,more=True),'full_G_cap_norm_upper_arb384':capNorm.str(100,more=True),
  'strict_rational_cap_I_upper':'749/10000000000','strict_rational_cap_norm_upper':'4525/1000000000',
  'proof_gaps_found':[],'new_C6_bound':False,'whole_face_I_certified_by_this_review':False,
  'scope':'Full-time G on the true rho<1/1000 S1 cap only; not the entire face or the short-potential complement.',
  'input_sha256':{str(p):sha(p)for p in inputs},'source_sha256':sha(Path(__file__)),
  'elapsed_seconds':perf_counter()-start}
 (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
 (HERE/'input_manifest.json').write_text(json.dumps({'input_sha256':out['input_sha256'],'source_sha256':out['source_sha256']},indent=2)+'\n')
 print('cap I',out['full_G_cap_I_upper_arb384']);print('cap norm',out['full_G_cap_norm_upper_arb384']);print(out['status'],'seconds',out['elapsed_seconds'])
if __name__=='__main__':main()
