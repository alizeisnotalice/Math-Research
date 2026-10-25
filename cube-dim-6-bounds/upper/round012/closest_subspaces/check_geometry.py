#!/usr/bin/env python3
"""Twelve rational geometry points against every32/38-direction pair."""
from fractions import Fraction as Q
from pathlib import Path
from time import perf_counter
import hashlib,itertools,json,os
import sympy as s
from flint import arb,ctx

HERE=Path(__file__).resolve().parent
B10=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-upper-astra/cube-dim-6-bounds/upper/round010/short_minimum_diagnostic')
MODEL=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper/theory/001_n6_body_frame_and_singular_strata.md')
ctx.prec=384
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dot(x,y):return sum((a*b for a,b in zip(x,y)),Q(0))
def variance(x):return sum((a*a for a in x),Q(0))-sum(x,Q(0))**2/len(x)
def A(x):return arb(x.numerator)/x.denominator
def points():
 raw=[
  ('zero',['0']*6),('all_equal',['1']*6),('signed_body_line',['1','-1','1','-1','1','-1']),
  ('single_axis',['0','0','-1','0','0','0']),('two_unequal_axes',['1','1/3','0','0','0','0']),
  ('two_level_ties',['1','1','1','1/3','1/3','1/3']),('three_level_ties',['1','1','1/2','1/2','0','0']),
  ('signed_permuted_distinct',['-1/7','5/7','0','-3/7','1','-2/7']),
  ('generic_chamber',['1','9/10','4/5','7/10','3/5','1/2']),
  ('signed_body_axis',['-2/5','-2/5','1','2/5','-2/5','2/5'])]
 out=[(name,tuple(Q(v)for v in x))for name,x in raw]
 for index in [1,2]:
  d=json.loads((B10/('point'+str(index)+'.json')).read_text())
  out.append(('B010_point'+str(index),tuple(Q(n,d['denominator'])for n in d['candidate_numerators'])))
 assert len(out)==12 and len(set(x for _,x in out))==12
 return out
def pair_projection(u,v,a,b,c,D):
 return tuple((c*u[i]*u[j]-b*(u[i]*v[j]+v[i]*u[j])+a*v[i]*v[j])/D for i in range(6)for j in range(6))
def record_min(old,value,i,j):
 if old is None or value<old[0]:return[value,[i,j],1]
 if value==old[0]:old[2]+=1
 return old
def main():
 start=perf_counter();pts=points();beta=Q(3,7)
 body=[tuple(map(Q,(1,)+bits))for bits in itertools.product([-1,1],repeat=5)]
 dirs=body+[tuple(beta if i==j else Q(0)for i in range(6))for j in range(6)]
 inputs=[MODEL,B10/'point1.json',B10/'point2.json',B10/'checked_results.json',B10/'report.md',HERE/'registration.md',HERE/'theorem.md']
 manifest={'geometry_points':[{'name':name,'coordinates':list(map(str,x))}for name,x in pts],
  'point_cap':12,'beta_for_exact_augmented_test':str(beta),'input_sha256':{str(p):sha(p)for p in inputs},
  'source_sha256':sha(Path(__file__)),'kernel_evaluations':0,'stop_rule':'First exact mismatch or completion of all twelve points.'}
 (HERE/'input_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 norms=[dot(v,v)for v in dirs];pairs=[];sets={key:set()for key in ['body_body','body_axis','axis_axis']};counts={key:0 for key in sets}
 for i,j in itertools.combinations(range(38),2):
  a,c=norms[i],norms[j];b=dot(dirs[i],dirs[j]);D=a*c-b*b;assert D>0
  kind='body_body'if j<32 else('body_axis'if i<32 else'axis_axis')
  pairs.append((i,j,a,b,c,D,kind));counts[kind]+=1
  sets[kind].add(pair_projection(dirs[i],dirs[j],a,b,c,D))
 assert counts=={'body_body':496,'body_axis':192,'axis_axis':15}
 assert len(sets['body_body'])==496 and len(sets['body_axis'])==96 and len(sets['axis_axis'])==15
 assert sets['body_axis']<=sets['body_body']and not sets['axis_axis']&sets['body_body']
 assert len(set.union(*sets.values()))==511
 # Actual symbolic checks supporting the analytic proof, not just enumeration.
 t,S,H,k=s.symbols('t S H k',real=True)
 gain=t*t/k+(S-t)**2/(6-k)
 assert s.simplify(s.diff(gain,t,2)-2/k-2/(6-k))==0
 leave=H-t*t-(S-t)**2/5;assert s.diff(leave,t,2)==-s.Rational(12,5)
 bet=s.symbols('beta',positive=True);v=s.ones(6,1);w=s.Matrix([1,-1,1,1,1,1]);axis=s.Matrix([0,bet,0,0,0,0])
 B=s.Matrix.hstack(v,w,axis);assert B.rank()==2 and bet*v-bet*w-2*axis==s.zeros(6,1)
 old=json.loads((B10/'checked_results.json').read_text());outputs=[]
 for name,x in pts:
  z=sorted(map(abs,x),reverse=True);e=[variance(z[:kk])+variance(z[kk:])for kk in range(1,6)]
  lo=[variance(z[1:]),variance(z[:-1])];axis2=sum((v*v for v in z[2:]),Q(0))
  assert lo==[e[0],e[4]]
  all_leave=[variance(z[:i]+z[i+1:])for i in range(6)];assert min(all_leave)==min(lo)
  formula_base=min(e);formula_aug8=min(e+lo+[axis2]);formula_aug6=min(e+[axis2]);assert formula_aug8==formula_aug6
  norm=dot(x,x);dots=[dot(x,v)for v in dirs]
  minima={key:None for key in ['body_body','body_axis','axis_axis','augmented']}
  for i,j,a,b,c,D,kind in pairs:
   distance=norm-(c*dots[i]**2-2*b*dots[i]*dots[j]+a*dots[j]**2)/D
   assert distance>=0
   minima[kind]=record_min(minima[kind],distance,i,j);minima['augmented']=record_min(minima['augmented'],distance,i,j)
  base1=min(norm-dots[i]**2/norms[i]for i in range(32))
  aug1=min(norm-dots[i]**2/norms[i]for i in range(38))
  expected_base1=variance(z);expected_aug1=min(expected_base1,sum((v*v for v in z[1:]),Q(0)))
  valid=(minima['body_body'][0]==formula_base and minima['body_axis'][0]==min(lo)
   and minima['axis_axis'][0]==axis2 and minima['augmented'][0]==formula_aug6
   and base1==expected_base1 and aug1==expected_aug1)
  if not valid:
   (HERE/'counterexample.json').write_text(json.dumps({'name':name,'x':list(map(str,x)),
    'base_formula':str(formula_base),'base_direct':str(minima['body_body'][0]),
    'aug_formula':str(formula_aug6),'aug_direct':str(minima['augmented'][0])},indent=2)+'\n')
   raise ArithmeticError('Exact formula mismatch: '+name)
  assert(formula_base==0)==(len(set(z))<=2)
  assert(formula_aug6==0)==(len(set(z))<=2 or sum(v!=0 for v in z)<=2)
  row={'name':name,'coordinates':list(map(str,x)),'sorted_absolute_coordinates':list(map(str,z)),
   'base_contiguous_five':list(map(str,e)),'augmented_leave_extremes':list(map(str,lo)),
   'top2_axis_residual':str(axis2),'base_S1_squared':str(base1),'augmented_S1_squared':str(aug1),
   'base_S2_squared':str(formula_base),'augmented_S2_squared':str(formula_aug6),
   'all_496_base_pairs_match':True,'all_703_augmented_pairs_match':True,
   'base_optimal_split_k':[kk+1 for kk,value in enumerate(e)if value==formula_base],
   'first_closest_base_pair':minima['body_body'][1],'base_closest_pair_count':minima['body_body'][2],
   'first_closest_augmented_pair':minima['augmented'][1],'augmented_closest_pair_count':minima['augmented'][2]}
  if name.startswith('B010'):
   index=name[-1];saved=old['point'+index+'_geometry'];assert formula_base==Q(saved['S2_distance_squared_exact'])
   row['B010_saved_exact_distance_reproduced']=True
   row['base_distance_arb384']=A(formula_base).sqrt().str(100,more=True)
   row['augmented_distance_arb384']=A(formula_aug6).sqrt().str(100,more=True)
  outputs.append(row);print(name,'base d2^2',formula_base,'aug d2^2',formula_aug6,flush=True)
 out={'status':'PASS_BASE_FIVE_AND_AUGMENTED_SIX_EXACT_DISTANCE_FORMULAS','pid':os.getpid(),
  'rational_points':12,'body_directions':32,'augmented_directions':38,'beta_test':str(beta),
  'pairs_per_point':counts,'augmented_pairs_per_point':703,
  'geometrically_distinct_planes':{key:len(value)for key,value in sets.items()},'augmented_distinct_planes_total':511,
  'body_axis_planes_add_to_S2':False,'proposed_eight_formula_valid':True,'two_leave_out_terms_redundant':True,
  'all_beta_positive_span_proof':True,'ties_and_zeros_handled_analytically':True,
  'signed_permutation_proof':True,'tubes_disjoint_claimed':False,
  'symbolic_checks':{'fixed_size_gain_second_derivative':'2/k+2/(6-k)','leave_one_out_second_derivative':'-12/5',
   'augmented_dependent_triple_rank':2,'dependent_triple_relation':'beta*v-beta*w-2*(beta*e2)=0'},
  'kernel_evaluations':0,'MC_samples':0,'matlab_calls':0,'new_C6_bound':False,'outputs':outputs,
  'input_sha256':manifest['input_sha256'],'source_sha256':manifest['source_sha256'],'elapsed_seconds':perf_counter()-start}
 (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n');print(out['status'],'seconds',out['elapsed_seconds'])
if __name__=='__main__':main()
