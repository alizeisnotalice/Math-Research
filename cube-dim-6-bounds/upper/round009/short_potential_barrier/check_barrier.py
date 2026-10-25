#!/usr/bin/env python3
"""Five exact pair geometries and eight scalar lower-barrier budgets."""
from fractions import Fraction as Q
from pathlib import Path
from time import perf_counter
import hashlib,json,math,os
import sympy as s
from flint import arb,ctx

HERE=Path(__file__).resolve().parent
MAIN=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds')
FROOT=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-literature/cube-dim-6-bounds/audit/round009/short_potential')
ctx.prec=384
LS=[Q(1,4),Q(1,2),Q(5,8),Q(3,4),Q(7,8),Q(1),Q(3,2),Q(2)]
def A(q):return arb(q.numerator)/q.denominator
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rational_floor(value,den=10**12):
 lo,hi=0,den
 while hi-lo>1:
  mid=(lo+hi)//2
  if value>arb(mid)/den:lo=mid
  else:hi=mid
 return Q(lo,den)
def main():
 started=perf_counter();pi=arb.pi();rt6=arb(6).sqrt();c=8*rt6;d=7*rt6
 a=3/(4*pi*pi*32**4);b=15/(8*pi*pi*c**6)
 rest=3/(4*pi*pi*d**4);tail=1800*rt6/pi*rest;kappa=(arb(5)/27).sqrt()
 T=Q(4,125);theta=T/4096;assert theta==Q(1,128000)
 assert(b/a-arb(5)/108).contains(0)
 geometries=[];v=s.ones(6,1)
 for k in range(1,6):
  w=s.Matrix([1]+[-1]*k+[1]*(5-k));B=s.Matrix.hstack(v,w);G=B.T*B
  D=G.det();assert D==4*k*(6-k)and G.inv()[0,0]==s.Rational(6,D)and G.inv()[1,1]==s.Rational(6,D)
  eig=G.eigenvals();assert min(eig)==2*min(k,6-k)
  geometries.append({'Hamming_k':k,'pair_count':16*math.comb(5,k),'Gram':[[int(x)for x in G.row(i)]for i in range(2)],
   'det_Gram':int(D),'min_Gram_eigenvalue':int(min(eig)),'pair_coarea_factor':'1/sqrt('+str(D)+')',
   'pseudoinverse_row_norm_squared':str(s.Rational(6,D))})
 assert sum(g['pair_count']for g in geometries)==496
 # Independent symbolic checks of constants and determinant monotonicity.
 P=s.symbols('P',positive=True);aa=3/(4*P**2*32**4);bb=15/(8*P**2*(8*s.sqrt(6))**6)
 assert s.simplify(bb/aa-s.Rational(5,108))==0
 D0,tau0,AA=s.symbols('D tau A',positive=True)
 expr=D0**s.Rational(3,2)/(tau0**2*D0+6*AA**2)**2
 logder=s.simplify(s.diff(expr,D0)/expr)
 assert s.simplify(logder-(18*AA**2-tau0**2*D0)/(2*D0*(tau0**2*D0+6*AA**2)))==0
 assert Q(36)*T*T<Q(108)
 rows=[]
 for L in LS:
  rawtau=a*A(L)/(4*tail)
  if rawtau<A(T):tau=rawtau
  else:assert rawtau>A(T);tau=A(T)
  Bsq=arb(3)/10*(rt6+A(L))**2 # worst determinant D=20, proved analytically
  lower_density=a/(2*pi*pi*arb(20).sqrt()*(tau*tau+Bsq)**2)
  rhomax=1/(kappa*(1/tau+A(theta)/lower_density))
  coeff=lower_density/kappa;offset=kappa/tau
  rows.append({'L':str(L),'tau_arb384':tau.str(100,more=True),'tau_is_raw_cutoff':bool(rawtau<A(T)),
   'B_worst_squared_arb384':Bsq.str(100,more=True),'density_coefficient_d_arb384':lower_density.str(100,more=True),
   'barrier_c_arb384':coeff.str(100,more=True),'barrier_C_arb384':offset.str(100,more=True),
   'rho_star_arb384':rhomax.str(100,more=True),'positive_time_tube_radius_tau_over_kappa':(tau/kappa).str(100,more=True)})
 best=rows[0]
 for row in rows[1:]:
  if arb(row['rho_star_arb384'])>arb(best['rho_star_arb384']):best=row
  else:assert arb(row['rho_star_arb384'])<arb(best['rho_star_arb384'])
 rho=arb(best['rho_star_arb384']);delta=Q(1,80000)
 if not A(delta)<rho:delta=rational_floor(rho,10**9)
 assert 0<A(delta)<rho
 tau=arb(best['tau_arb384']);coef=arb(best['density_coefficient_d_arb384'])
 lowerK=coef*(1/(kappa*A(delta))-1/tau);assert lowerK>A(theta)
 strictK=rational_floor(lowerK);assert A(theta)<A(strictK)<lowerK
 inputs=[MAIN/'upper/theory/001_n6_body_frame_and_singular_strata.md',FROOT/'registration.md',FROOT/'check_results.json',HERE/'registration.md',HERE/'theorem.md']
 out={'status':'PASS_UNIFORM_ALL_PAIR_SHORT_POTENTIAL_LOWER_BARRIER','pid':os.getpid(),
  'model':'alpha=1 equal-weight body32 Cauchy','precision_bits':384,'canonical_pair_geometries':geometries,
  'scalar_budget_cases':8,'kernel_evaluations':0,'MC_samples':0,'matlab_calls':0,
  'constants':{'spectral_lower_c':'8sqrt(6)','deleted_quotient_lower':'7sqrt(6)',
   'a_arb384':a.str(100,more=True),'b_arb384':b.str(100,more=True),'b_over_a_exact':'5/108',
   'Crest_arb384':rest.str(100,more=True),'Ctail_arb384':tail.str(100,more=True),
   'kappa_exact':'sqrt(5/27)','kappa_arb384':kappa.str(100,more=True)},
  'T':str(T),'theta_exact':str(theta),'rows':rows,'best_registered_L':best['L'],
  'best_barrier':best,'certified_distance_threshold_exact':str(delta),
  'K_T_uniform_strict_rational_lower_on_closed_tube':str(strictK),
  'K_T_lower_at_distance_threshold_arb384':lowerK.str(100,more=True),
  'includes_S1_and_all_plane_intersections':True,'K_T_on_S2':'positive infinity',
  'pair_Jacobian_retained':True,'zero_probability_conditioning_used':False,
  'symbolic_b_over_a_residual':'0','symbolic_det_monotonicity_residual':'0',
  'unproved_compact_complement':'{x in [-1,1]^6: dist(x,S2) >= '+str(delta)+'}',
  'whole_cube_m_T_certified':False,'new_C6_bound':False,
  'input_sha256':{str(p):sha(p)for p in inputs},'source_sha256':sha(Path(__file__)),
  'elapsed_seconds':perf_counter()-started}
 (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
 (HERE/'input_manifest.json').write_text(json.dumps({k:out[k]for k in ['model','scalar_budget_cases','input_sha256','source_sha256']},indent=2)+'\n')
 for row in rows:print('L',row['L'],'rho_star',row['rho_star_arb384'][:65],flush=True)
 print('best L',best['L'],'delta',delta,'K_T >',strictK,'> theta',theta)
 print(out['status'],'seconds',out['elapsed_seconds'],flush=True)
if __name__=='__main__':main()
