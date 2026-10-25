#!/usr/bin/env python3
"""Four finite analytic body32 alpha=3/2 global lower schemes, Arb384."""
from fractions import Fraction as Q
from pathlib import Path
from time import perf_counter
import hashlib,json,math,os
import sympy as s
from flint import arb,ctx

HERE=Path(__file__).resolve().parent
CROOT=HERE.parents[2]/'upper/round008/alpha15_coercivity'
FROOT=HERE.parents[2]/'audit/round008/alpha15_product'
INTERFACE=HERE.parents[2]/'audit/round006/alpha15_review/report.md'
ctx.prec=384
def A(q):return arb(q.numerator)/q.denominator
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def power(x,q):return x**A(q)
def coeffs(K,rate0,c):
 out=[]
 for k in range(K+1):
  q=Q(4)+Q(4*k,3);rate=rate0 if k==0 else(arb(32)if k%2==0 else c)
  rising=math.prod(range(3,3+k))
  out.append(A(q).gamma()/(96*arb.pi()**3*4**k*math.factorial(k)*rising*power(rate,q)))
 return out
def polynomial(z,cs):
 out=arb(0)
 for k in range(len(cs)-1,-1,-1):out=out*z+((-1)**k)*cs[k]
 return out
def F(z,cs):
 total=sum((((-1)**k)*a*z**k/A(Q(3)+Q(4*k,3))for k,a in enumerate(cs)),arb(0))
 return power(z,Q(9,4))*total
def optimize(cs):
 lo=Q(0);hi=Q(1)
 while not polynomial(A(hi),cs)<0:
  hi*=2;assert hi<2**40
 for _ in range(220):
  mid=(lo+hi)/2;value=polynomial(A(mid),cs)
  if value>0:lo=mid
  elif value<0:hi=mid
  else:raise RuntimeError('Need more precision for strict root bisection')
 assert polynomial(A(lo),cs)>0 and polynomial(A(hi),cs)<0
 root=arb(A((lo+hi)/2),A((hi-lo)/2))
 optimum=power(arb(6),Q(-9,4))*F(root,cs)
 constructed=power(arb(6),Q(-9,4))*F(A(lo),cs)
 assert constructed>0
 # A decimal rational certificate with scale 1e18 is independently compared.
 qlo,qhi=0,10**18
 while qhi-qlo>1:
  mid=(qlo+qhi)//2
  if constructed>arb(mid)/10**18:qlo=mid
  else:qhi=mid
 lower=Q(qlo,10**18);assert constructed>A(lower)>0
 return {'root_interval_exact':[str(lo),str(hi)],'constructive_z':str(lo),
  'unique_optimal_z_arb384':root.str(100,more=True),
  'global_g_constructive_arb384':constructed.str(100,more=True),
  'optimal_method_g_arb384':optimum.str(100,more=True),
  'strict_rational_global_g_lower':str(lower)},constructed,optimum,lower
def symbolic_checks():
 z=s.symbols('z',positive=True);out={}
 for K in [1,3,5]:
  cs=s.symbols('a0:'+str(K+1));R=sum((-1)**k*cs[k]*z**k/(s.Rational(3)+s.Rational(4*k,3))for k in range(K+1))
  residual=s.simplify(s.diff(z**s.Rational(9,4)*R,z)-s.Rational(3,4)*z**s.Rational(5,4)*sum((-1)**k*cs[k]*z**k for k in range(K+1)))
  assert residual==0;out['tail_derivative_K'+str(K)]=str(residual)
 x=s.symbols('x');polys=[sum((-1)**k*x**(2*k)/s.factorial(2*k)for k in range(K+1))for K in range(6)]
 for k in range(1,6):assert s.expand(s.diff(polys[k],x,2)+polys[k-1])==0
 out['global_cosine_induction_second_derivative_identities']='PASS for all degrees through 10'
 return out
def main():
 started=perf_counter();c=32*power(arb(8)/3,Q(-1,4))
 mu=64*power(arb(6),Q(3,4))*(arb(5)/4).gamma()/(arb.pi().sqrt()*(arb(15)/4).gamma())
 assert c<mu<32
 fdata=json.loads((FROOT/'check_results.json').read_text());U=arb(fdata['final_constants']['full_G_norm_upper'])
 threshold=U/(4*arb(2).sqrt());schemes=[];closed_checks=[]
 for label,K,rate0 in [('quadratic_radial',1,arb(32)),('quadratic_Jensen',1,mu),('degree6_Jensen',3,mu),('degree10_Jensen',5,mu)]:
  cs=coeffs(K,rate0,c);mono={}
  if K==3:
   assert cs[2]**2<3*cs[1]*cs[3]
   mono['A2_squared_over_3_A1_A3']=(cs[2]**2/(3*cs[1]*cs[3])).str(95,more=True)
  if K==5:
   assert cs[2]**2<cs[1]*cs[3] and 2*cs[4]**2<5*cs[3]*cs[5]
   mono['A2_squared_over_A1_A3']=(cs[2]**2/(cs[1]*cs[3])).str(95,more=True)
   mono['2_A4_squared_over_5_A3_A5']=(2*cs[4]**2/(5*cs[3]*cs[5])).str(95,more=True)
  row,constructed,optimum,lower=optimize(cs)
  if K==1:
   closed=4*cs[0]/39*power(cs[0]/(6*cs[1]),Q(9,4))
   assert(closed-optimum).contains(0);closed_checks.append({'scheme':label,'closed_form_g':closed.str(100,more=True)})
   row['optimal_cube_T_arb384']=power(6*cs[1]/cs[0],Q(3,4)).str(100,more=True)
  rootlo=Q(row['root_interval_exact'][0]);row['constructive_cube_tail_T']=power(arb(6),Q(3,4))*power(A(rootlo),Q(-3,4))
  row['constructive_cube_tail_T']=row['constructive_cube_tail_T'].str(100,more=True)
  allocated_cost=power(4*U/optimum,Q(4,3))
  valid_large_cost_upper=power(4*U.upper()/A(lower),Q(4,3))
  assert optimum<threshold and allocated_cost>64
  row.update({'scheme':label,'cosine_degree':2*K,'positive_zero_rate':'32'if label=='quadratic_radial'else'mu',
   'coefficients_Ak_arb384':[a.str(100,more=True)for a in cs],
   'unique_root_monotonicity_checks':mono,'global_optimum_of_this_method_below_fixed_norm_threshold':True,
   'allocated_cost_arb384':allocated_cost.str(95,more=True),
   'valid_noncompetitive_cost_upper_using_rational_g':valid_large_cost_upper.str(95,more=True),
   'allocated_cost_strictly_exceeds64':True})
  schemes.append(row)
 # Determine the largest optimized certificate by strict interval comparison.
 best=schemes[0]
 for row in schemes[1:]:
  if arb(row['optimal_method_g_arb384'])>arb(best['optimal_method_g_arb384']):best=row
  else:assert arb(row['optimal_method_g_arb384'])<arb(best['optimal_method_g_arb384'])
 inputs=[CROOT/'report.md',CROOT/'coercivity_results.json',FROOT/'report.md',FROOT/'check_results.json',INTERFACE,HERE/'registration.md',HERE/'theorem.md']
 out={'status':'PASS_ALPHA15_GLOBAL_CUBE_LOWER_NONCOMPETITIVE','pid':os.getpid(),'precision_bits':384,
  'analytic_schemes':4,'scalar_budget_scenarios':8,'kernel_evaluations':0,'MC_samples':0,
  'alpha':'3/2','dimension':6,'cube_Euclidean_radius_squared':6,
  'c6_arb384':c.str(100,more=True),'angular_mean_mu_arb384':mu.str(100,more=True),
  'F_norm_budget_input':U.str(95,more=True),'fixed_norm_required_g_for_cost64':threshold.str(95,more=True),
  'schemes':schemes,'best_scheme':best['scheme'],'best_strict_rational_global_g_lower':best['strict_rational_global_g_lower'],
  'best_optimal_method_g_arb384':best['optimal_method_g_arb384'],
  'global_g_proof_status':'whole cube including boundaries and singular strata as extended-valued inequality',
  'new_C6_bound':False,'generic_C6_upper_retained':'64',
  'budget_obstruction_scope':'Only these four optimized lower-certificate schemes paired with the fixed F norm-envelope value; not an upper bound on actual g or lower bound on actual cost.',
  'closed_quadratic_formula_checks':closed_checks,'symbolic_checks':symbolic_checks(),
  'input_sha256':{str(p):sha(p)for p in inputs},'source_sha256':sha(Path(__file__)),
  'elapsed_seconds':perf_counter()-started}
 (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
 (HERE/'input_manifest.json').write_text(json.dumps({k:out[k]for k in ['alpha','dimension','analytic_schemes','scalar_budget_scenarios','input_sha256','source_sha256']},indent=2)+'\n')
 print('mu',mu.str(30,more=True),'c6',c.str(30,more=True),flush=True)
 for row in schemes:print(row['scheme'],'g',row['optimal_method_g_arb384'][:62],'cost',row['allocated_cost_arb384'][:45],flush=True)
 print('best',out['best_scheme'],'strict global g >',out['best_strict_rational_global_g_lower'])
 print(out['status'],'seconds',out['elapsed_seconds'])
if __name__=='__main__':main()
