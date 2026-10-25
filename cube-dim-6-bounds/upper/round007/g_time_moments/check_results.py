#!/usr/bin/env python3
"""No kernel replay: exact geometry and independent acceptance postcheck."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json
import sympy as s
from flint import arb,ctx
import run_moments as m
import postprocess_sparse as sp

HERE=Path(__file__).resolve().parent;ctx.prec=384
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 raw=json.loads((HERE/'results.json').read_text());new=json.loads((HERE/'sparse_shift_results.json').read_text())
 assert raw['status']=='COMPLETE_LOCAL_ONLY'and len(raw['time_nodes'])==8 and len(raw['gradients'])==5
 assert raw['elapsed_seconds']<1200
 assert sum(g['complete_group_terms']==967680 for g in raw['time_nodes']+raw['gradients'])==13
 vals,masses=m.bands(raw['time_nodes']);assert all(lo>0 and hi>0 for lo,hi in masses)
 for row in raw['gradients']:
  assert(arb(row['H_arb384'])-vals[row['time_index']]).contains(0)
 m.U2=sp.sparse_u2;m.band_lower=sp.sparse_band_lower
 accepted=[]
 for box,gr in zip(new['taylor_boxes'],new['gradients']):
  r=Q(box['radius']);l=gr['time_index'];g=[arb(v)for v in gr['gradient_x2_to_x6']]
  value,linear,loss=m.taylor_lower(r,l,new['center_numerators'],new['denominator'],vals,masses,g)
  assert value>arb(1)/60000 and(value-arb(box['G_lower_arb384'])).contains(0)
  assert value>m.A(Q(box['G_lower_exact_rational']))>arb(1)/60000
  # The next radius increment fails this sufficient criterion, not the true G bound.
  assert not m.taylor_lower(r+Q(1,1000000),l,new['center_numerators'],new['denominator'],vals,masses,g)[0]>arb(1)/60000
  accepted.append({'epsilon':str(m.TAUS[l]),'radius':str(r),'G_lower_endpoint':value.lower().str(90,more=True),
   'G_lower_exact_rational':box['G_lower_exact_rational']})
 # Fair diagnostic only: retain the same strict gradient, geometry, Hessian
 # sign, and sparse spatial factor, replacing J2 by epsilon^-2 H.
 def uniform_u2(r,l,vals,masses,details=False):
  return m.A(m.TAUS[l])**-2*sp.phi(m.A(r)/m.A(m.TAUS[l]))*vals[l].upper()
 uniform=[];m.U2=uniform_u2
 for gr in new['gradients']:
  l=gr['time_index'];g=[arb(v)for v in gr['gradient_x2_to_x6']]
  rr=m.largest(lambda rr:m.taylor_lower(rr,l,new['center_numerators'],new['denominator'],vals,masses,g)[0]>arb(1)/60000)
  uniform.append({'epsilon':str(m.TAUS[l]),'radius':str(rr)})
 m.U2=sp.sparse_u2
 best_uniform=max(Q(v['radius'])for v in uniform)
 tenfold=10*Q(new['best_original_radius']);tenfold_checks=[]
 for gr in new['gradients']:
  l=gr['time_index'];g=[arb(v)for v in gr['gradient_x2_to_x6']]
  value=m.taylor_lower(tenfold,l,new['center_numerators'],new['denominator'],vals,masses,g)[0]
  assert value.upper()<arb(1)/60000
  tenfold_checks.append({'epsilon':gr['tau'],'Taylor_criterion_upper':value.upper().str(90,more=True)})
 best=new['best_taylor_box'];r=Q(best['radius']);bs=[(Q(a),Q(b))for a,b in best['box_bounds_x2_to_x6']]
 assert all(bs[i][1]<bs[j][0]or bs[j][1]<bs[i][0]for i in range(5)for j in range(i+1,5))
 witness=list(map(Q,['19/20','9/10','17/20','4/5','3/4']))
 assert any(all(v<a or v>b for v in witness)for a,b in bs)
 vol=Q(120)
 for a,b in bs:vol*=b-a
 assert vol<1
 # Independent symbolic check of the one-factor lower logarithmic curvature.
 z,t=s.symbols('z t',real=True,nonzero=True)
 logsecond=-2*(t*t-z*z)/(t*t+z*z)**2
 remainder=2*z*z*(z*z+3*t*t)/(t*t*(t*t+z*z)**2)
 residual=s.simplify(logsecond+2/(t*t)-remainder);assert residual==0
 l=next(g['time_index']for g in new['gradients']if g['tau']==best['epsilon'])
 moment=sp.sparse_u2(r,l,vals,masses)
 old_global_budget=m.A(m.TAUS[l])**-2*vals[l]
 old_box_budget=uniform_u2(r,l,vals,masses)
 # This scalar is a diagnostic comparison only; it is absent from the certificate construction.
 out={'status':'PASS_LOCAL_TAYLOR_ACCEPTANCE_AND_EXACT_COVERAGE','kernel_replayed':False,
  'new_H_nodes':8,'same_center_gradient_evaluations':5,'new_spatial_centers':1,'accepted':accepted,
  'best_epsilon':best['epsilon'],'best_radius':best['radius'],
  'best_original_multiplicative_radius':new['best_original_radius'],
  'best_radius_over_best_original':best['radius_over_best_original'],
  'best_radius_over_best_exponential_taylor':new['best_sparse_over_best_exponential'],
  'matched_gradient_uniform_epsilon_inverse_square_Taylor_radii':uniform,
  'best_matched_uniform_Taylor_radius':str(best_uniform),
  'best_time_moment_radius_over_best_matched_uniform_radius':str(r/best_uniform),
  'tenfold_reference_radius':str(tenfold),'tenfold_rejected_by_current_sufficient_budget':tenfold_checks,
  'best_G_lower_endpoint':next(a['G_lower_endpoint']for a in accepted if a['epsilon']==best['epsilon']),
  'best_G_lower_exact_rational':best['G_lower_exact_rational'],
  'best_J2_box_upper':moment.str(95,more=True),
  'comparison_epsilon_inverse_square_H_center':old_global_budget.str(95,more=True),
  'J2_box_over_that_comparison':(moment/old_global_budget).str(95,more=True),
  'comparison_epsilon_inverse_square_H_whole_box_upper':old_box_budget.str(95,more=True),
  'J2_box_over_matched_uniform_box_upper':(moment/old_box_budget).str(95,more=True),
  'exact_positive_face_volume_covered':str(vol),'exact_uncovered_volume':str(1-vol),
  'uncovered_rational_point':['1','19/20','9/10','17/20','4/5','3/4'],
  'lower_log_curvature_symbolic_residual':str(residual),'global_g_target_certified':False,'new_C6_bound':False,
  'raw_results_sha256':sha(HERE/'results.json'),'sparse_results_sha256':sha(HERE/'sparse_shift_results.json'),
  'checker_sha256':sha(Path(__file__))}
 (HERE/'checked_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
