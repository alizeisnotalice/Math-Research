#!/usr/bin/env python3
"""Independent coefficient normalization and stored-certificate reconstruction.

No new analytic scenario: use Gamma angular moments instead of the original
Pochhammer coefficient implementation, and reconstruct the four stored bounds.
"""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json
from flint import arb,ctx
ctx.prec=384;HERE=Path(__file__).resolve().parent
def A(x):return arb(x.numerator)/x.denominator
def pw(x,q):return x**A(q)
def poly(z,cs):return sum((((-1)**k)*v*z**k for k,v in enumerate(cs)),arb(0))
def tail(z,cs):return pw(z,Q(9,4))*sum((((-1)**k)*v*z**k/A(Q(3)+Q(4*k,3))for k,v in enumerate(cs)),arb(0))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 d=json.loads((HERE/'results.json').read_text());c=32*pw(arb(8)/3,Q(-1,4));mu=arb(d['angular_mean_mu_arb384'])
 U=arb(d['F_norm_budget_input']);rows=[]
 for row in d['schemes']:
  K=row['cosine_degree']//2;cs=[]
  for k in range(K+1):
   angular=(arb(k)+arb(1)/2).gamma()*arb(3).gamma()/(arb.pi().sqrt()*arb(k+3).gamma()*arb(2*k+1).gamma())
   rate=(arb(32)if row['positive_zero_rate']=='32'else mu)if k==0 else(arb(32)if k%2==0 else c)
   q=Q(4)+Q(4*k,3)
   value=A(q).gamma()*angular/(96*arb.pi()**3*pw(rate,q));cs.append(value)
   assert(value-arb(row['coefficients_Ak_arb384'][k])).contains(0)
  lo,hi=map(Q,row['root_interval_exact']);assert poly(A(lo),cs)>0 and poly(A(hi),cs)<0
  g=pw(arb(6),Q(-9,4))*tail(A(lo),cs);lower=Q(row['strict_rational_global_g_lower'])
  assert g>A(lower)>0
  if row['scheme']=='quadratic_radial':
   assert(cs[0]-1/(16*32**4*arb.pi()**3)).contains(0)
   explicit_b=(arb(16)/3).gamma()/(1152*arb.pi()**3*pw(c,Q(16,3)))
   assert(cs[1]-explicit_b).contains(0)
  optimum=arb(row['optimal_method_g_arb384']);assert optimum<arb(109)/10**9
  upper=lower+Q(1,10**18);assert A(lower)<optimum<A(upper)
  cost=pw(4*U/optimum,Q(4,3));assert cost>1200000
  rows.append({'scheme':row['scheme'],'Gamma_angular_coefficients_match':True,
   'strict_root_bracket_reconstructed':True,'strict_global_g_lower':str(lower),
   'optimal_method_strict_rational_interval':[str(lower),str(upper)],
   'method_optimum_strict_upper':'109/1000000000','allocated_cost_strict_lower':'1200000'})
 threshold=U/(4*arb(2).sqrt());assert threshold>arb(177)/10**6
 out={'status':'PASS_INDEPENDENT_NORMALIZATION_AND_ROOT_RECONSTRUCTION','new_scenarios':0,
  'whole_cube_radius_squared':6,'checked_schemes':rows,
  'fixed_norm_required_g_strict_lower':'177/1000000',
  'best_global_g_strict_lower':d['best_strict_rational_global_g_lower'],
  'actual_g_upper_bound_claimed':False,'actual_cost_lower_bound_claimed':False,'new_C6_bound':False,
  'results_sha256':sha(HERE/'results.json'),'checker_sha256':sha(Path(__file__))}
 (HERE/'checked_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
