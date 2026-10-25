#!/usr/bin/env python3
"""Independent Gamma normalization and exact rational corollary reconstruction."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json
from flint import arb,ctx
ctx.prec=384;HERE=Path(__file__).resolve().parent
def A(q):return arb(q.numerator)/q.denominator
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def peak4(rate):
 sphere=2*arb.pi()**2/arb(2).gamma()
 return sphere*arb(4).gamma()/((2*arb.pi())**4*rate**4)
def main():
 data=json.loads((HERE/'results.json').read_text());pi=arb.pi();rt6=arb(6).sqrt()
 a=peak4(arb(32));b=(2*pi*pi)*arb(6).gamma()/((2*pi)**4*8*(8*rt6)**6)
 rest=peak4(7*rt6);tail=30*rest*(2/pi)*(30*rt6);kappa=2*(b/a).sqrt()
 for val,key in [(a,'a_arb384'),(b,'b_arb384'),(rest,'Crest_arb384'),(tail,'Ctail_arb384'),(kappa,'kappa_arb384')]:
  assert(val-arb(data['constants'][key])).contains(0)
 rows=[]
 for row in data['rows']:
  L=Q(row['L']);tau=a*A(L)/(4*tail);assert tau<arb(4)/125
  # Direct row norm of the k=1 Gram inverse: 6/(36-16)=3/10.
  Bsq=arb(6)/20*(rt6+A(L))**2
  d=a/(2*pi*pi*arb(20).sqrt()*(tau*tau+Bsq)**2)
  rho=1/(kappa*(1/tau+arb(1)/128000/d))
  assert(rho-arb(row['rho_star_arb384'])).contains(0)
  rows.append({'L':str(L),'constant_and_rho_reconstruction':True})
 best=data['best_barrier'];L=Q(best['L']);tau=a*A(L)/(4*tail)
 Bsq=arb(3)/10*(rt6+A(L))**2;d=a/(2*pi*pi*arb(20).sqrt()*(tau*tau+Bsq)**2)
 c=d/kappa;C=kappa/tau
 assert c>arb(1)/5000000000 and C<40000
 delta=Q(1,80000);simpleK=Q(1,5000000000)*(1/delta-40000)
 assert simpleK==Q(1,125000)>Q(1,128000)
 sharpK=d*(1/(kappa*A(delta))-1/tau)
 assert sharpK>A(Q(data['K_T_uniform_strict_rational_lower_on_closed_tube']))>arb(1)/128000
 out={'status':'PASS_INDEPENDENT_ALL_PAIR_BARRIER_AND_RATIONAL_COROLLARY','new_parameter_cases':0,
  'normalized_Gamma_peak_checks':True,'reconstructed_eight_L_budgets':rows,
  'best_registered_L':str(L),'barrier_c_strict_lower':'1/5000000000',
  'barrier_C_strict_upper':'40000','barrier_formula':'K_(4/125)(x) >= (1/5000000000)*(1/dist(x,S2)-40000)_+',
  'on_S2':'positive infinity','distance_threshold_exact':'1/80000',
  'uniform_K_lower_on_closed_tube_exact':'1/125000','required_threshold_exact':'1/128000',
  'strict_margin_exact':str(Q(1,125000)-Q(1,128000)),
  'sharp_K_lower_arb384':sharpK.str(100,more=True),
  'whole_cube_m_T_certified':False,'new_C6_bound':False,
  'remaining_compact_set':'{x in [-1,1]^6: dist(x,S2)>=1/80000}',
  'results_sha256':sha(HERE/'results.json'),'checker_source_sha256':sha(Path(__file__))}
 (HERE/'checked_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
