#!/usr/bin/env python3
"""Independent postcheck of acceptance arithmetic and exact coverage bookkeeping.

This does not re-evaluate the full kernel. A primary/independent reviewer must
still replay the kernel formula before promoting these model-local results.
"""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,itertools,json
from flint import arb,ctx
import sympy as s

HERE=Path(__file__).resolve().parent;ctx.prec=384
def A(q):return arb(q.numerator)/q.denominator
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 d=json.loads((HERE/'local_results.json').read_text());target=arb(1)/60000
 assert d['status']=='COMPLETE_LOCAL_ONLY' and len(d['outputs'])==len(d['boxes'])==8
 assert sum(len(row['values'])for row in d['outputs'])==24
 assert d['elapsed_seconds']<1200
 volumes=[];permuted=[];witness=tuple(map(Q,['19/20','9/10','17/20','4/5','3/4']))
 for row,box in zip(d['outputs'],d['boxes']):
  assert row['name']==box['name'] and row['complete_cone_terms']==967680
  e=Q(box['epsilon']);r=Q(box['radius'])
  val=arb(next(v['H_arb384']for v in row['values']if v['epsilon']==str(e)))
  lower=val*(-A(r)/A(e)).exp();assert lower>target
  assert (lower-arb(box['G_lower_arb384'])).contains(0)
  bs=tuple((Q(a),Q(b))for a,b in box['bounds_x2_to_x6'])
  for c,(a,b)in zip(row['point_numerators'][1:],bs):
   assert (a,b)==(max(Q(0),Q(c,1000)-r),min(Q(1),Q(c,1000)+r))
  images=set(itertools.permutations(bs));vol=Q(1)
  for a,b in bs:vol*=b-a
  assert not any(all(a<=v<=b for v,(a,b)in zip(witness,img))for img in images)
  volumes.append({'name':row['name'],'unique_positive_face_permutations':len(images),'one_box_volume':str(vol),
   'sum_image_volume_upper':str(vol*len(images))})
  permuted.extend(images)
 total=sum((Q(v['sum_image_volume_upper'])for v in volumes),Q(0))
 assert total<Q(1,1000000)
 # Exact symbolic completion-of-squares verification for the valid new ODE.
 n=s.symbols('N',positive=True);b=s.symbols('b',real=True)
 residual=s.expand((n+1)*(b+1/(2*(n+1)))**2-1/(4*(n+1))-((n+1)*b*b+b))
 assert s.simplify(residual)==0
 c=Q(4225,4224);assert c==1+Q(1,4*32*33) and c<Q(129,128)
 # One already-certified S2 point is a calibration overlap, not new evidence.
 oldpath=HERE.parents[2]/'upper/round004/truncated_kernel/arb_reference_results.json'
 old=json.loads(oldpath.read_text());overlap=next(row for row in d['outputs']if row['name']=='S2_one_four')
 val=arb(next(v['H_arb384']for v in overlap['values']if v['epsilon']=='1/100'))
 assert (val-arb(old['H_epsilon_arb_256'])).contains(0)
 out={'status':'PASS_LOCAL_ACCEPTANCE_AND_EXACT_COVERAGE_BOOKKEEPING','new_point_epsilon_pairs':23,
  'calibration_overlap_pairs':1,'distinct_geometric_points':8,'local_boxes':8,
  'box_coverage':volumes,'total_permuted_volume_upper':str(total),
  'total_permuted_volume_upper_arb':A(total).str(95,more=True),
  'uncovered_volume_lower':str(1-total),'uncovered_rational_witness':['1']+list(map(str,witness)),
  'global_g_target_certified':False,'new_C6_bound':False,
  'new_ODE_coefficient_exact':str(c),'new_ODE_coefficient_arb384':A(c).str(95,more=True),
  'old_over_new_side_length_factor':((arb(129)/128)/A(c)).sqrt().str(95,more=True),
  'symbolic_square_residual':str(s.simplify(residual)),
  'kernel_replayed':False,'kernel_complete_group_summation_in_results':True,
  'calibration_source_sha256':sha(oldpath),'result_sha256':sha(HERE/'local_results.json'),
  'checker_source_sha256':sha(Path(__file__))}
 (HERE/'checked_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
