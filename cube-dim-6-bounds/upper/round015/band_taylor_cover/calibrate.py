#!/usr/bin/env python3
"""Eight existing-data band Taylor tests; zero new kernel or gradient calls."""
from fractions import Fraction as Q
from pathlib import Path
from time import perf_counter
import hashlib,itertools,json,os
import sympy as s
from flint import arb,ctx

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parents[2]
A7=PROJECT/'upper/round007/g_time_moments/results.json'
BROOT=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-upper-astra/cube-dim-6-bounds/upper')
B10=BROOT/'round010/short_minimum_diagnostic'
B9=BROOT/'round009/short_potential'
ctx.prec=384
STARTS=[Q(1,1000),Q(1,500),Q(1,250),Q(1,125)]
RADII=[Q(1,250),Q(1,125)]
def A(x):return arb(x.numerator)/x.denominator
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def endpoint(x,upper=False):
 v=x.upper()if upper else x.lower();assert v.is_exact();f=v.fmpq();return Q(int(f.numerator),int(f.denominator))
def phi(u):return(12*(u/12).asinh()).exp()
def normsq(v):return sum((x*x for x in v),Q(0))
def within(v):return sum((x*x for x in v),Q(0))-sum(v,Q(0))**2/len(v)
def residual(v,k):return within(v[:k])+within(v[k:])
def cell(c,r):
 bounds=[(max(Q(0),x-r),min(Q(1),x+r))for x in c[1:]]
 vertices=[(Q(1),)+v for v in itertools.product(*bounds)]
 assert all(all(v[i]>=v[i+1]for i in range(5))and v[-1]>=0 for v in vertices)
 weights=[]
 for bits in itertools.product([0,1],repeat=5):
  w=Q(1)
  for x,(lo,hi),bit in zip(c[1:],bounds,bits):
   p=(x-lo)/(hi-lo);w*=p if bit else 1-p
  weights.append(w)
 assert sum(weights)==1 and all(w>=0 for w in weights)
 assert tuple(sum((w*v[i]for w,v in zip(weights,vertices)),Q(0))for i in range(6))==c
 return bounds,vertices
def toy_checks():
 root=[tuple([Q(1)]+[Q(1)]*j+[Q(0)]*(5-j))for j in range(6)]
 bary=tuple(sum((v[i]for v in root),Q(0))/6 for i in range(6))
 assert all(min(residual(v,k)for k in range(1,6))==0 for v in root)
 assert min(residual(bary,k)for k in range(1,6))==Q(1,9)
 delta=Q(1,100000)
 assert not any(all(residual(v,k)<=delta*delta for v in root)for k in range(1,6))
 lengths={(i,j):normsq(tuple(root[i][h]-root[j][h]for h in range(6)))for i,j in itertools.combinations(range(6),2)}
 edge=min(pair for pair,value in lengths.items()if value==max(lengths.values()));assert edge==(0,5)
 mid=tuple((root[0][i]+root[5][i])/2 for i in range(6))
 left=list(root);right=list(root);left[5]=mid;right[0]=mid
 def volume(vertices):
  M=s.Matrix([[s.Rational(vertices[j][i+1]-vertices[0][i+1])for j in range(1,6)]for i in range(5)])
  return abs(M.det())/s.factorial(5)
 assert volume(root)==s.Rational(1,120)and volume(left)==volume(right)==s.Rational(1,240)
 l=s.symbols('l0:6');vv=[s.Matrix(list(v))for v in root];mm=s.Matrix(list(mid))
 parent=sum((l[i]*vv[i]for i in range(6)),s.zeros(6,1))
 child1=(l[0]-l[5])*vv[0]+2*l[5]*mm+sum((l[i]*vv[i]for i in range(1,5)),s.zeros(6,1))
 child2=(l[5]-l[0])*vv[5]+2*l[0]*mm+sum((l[i]*vv[i]for i in range(1,5)),s.zeros(6,1))
 assert parent==child1==child2
 return {'toy_count':2,'root_vertices':[list(map(str,v))for v in root],
  'root_all_vertices_in_S2':True,'root_barycenter_S2_squared':'1/9','per_vertex_closest_plane_rule_rejected':True,
  'same_plane_tube_root_rejected_at_delta':'1/100000','longest_edge':[0,5],
  'root_volume':'1/120','child_volumes':['1/240','1/240'],'split_union_barycentric_identity':'PASS'}
def main():
 started=perf_counter();old=json.loads(A7.read_text());bg=json.loads((B10/'gradient_results.json').read_text())
 base=json.loads((B9/'postprocess_results.json').read_text());refuted=json.loads((B10/'checked_results.json').read_text())
 assert old['center_numerators']==bg['center_numerators']and old['denominator']==bg['denominator']
 c=tuple(Q(n,old['denominator'])for n in old['center_numerators']);T=Q(4,125);m0=T/4000;bare=T/4096
 assert str(T)in refuted['strict_refutation_times']
 times={Q(row['tau']):arb(row['H_arb384'])for row in old['time_nodes']}
 grads={Q(row['tau']):row for row in old['gradients']}
 end=next(row for row in bg['gradient_passes']if Q(row['tau'])==T)
 HT=arb(end['H_arb384']);assert(HT-times[T]).contains(0)
 gradT=[arb(v)for v in end['gradient_x2_to_x6']]
 oldrow=next(row for row in base['rows']if Q(row['T'])==T)
 oldr=Q(oldrow['local_band_radius_m0']);assert oldr==Q(991,500000)
 assert Q(oldrow['m0_for_cost_125_over_2'])==m0==Q(1,125000)
 inputs=[A7,B10/'gradient_results.json',B10/'checked_results.json',B9/'postprocess_results.json',B9/'report.md',
  PROJECT/'upper/round009/short_potential_barrier/checked_results.json',
  PROJECT/'upper/round012/closest_subspaces/theorem.md',HERE/'registration.md',HERE/'theorem.md',HERE/'checker_schema.md']
 meta={'center':list(map(str,c)),'T':str(T),'starts':list(map(str,STARTS)),'radii':list(map(str,RADII)),
  'matched_target':str(m0),'bare64_target':str(bare),'old_Phi_radius':str(oldr),
  'input_sha256':{str(p):sha(p)for p in inputs},'source_sha256':sha(Path(__file__))}
 (HERE/'input_manifest.json').write_text(json.dumps(meta,indent=2)+'\n')
 rows=[]
 for a in STARTS:
  Ha=times[a];ga=grads[a];assert(arb(ga['H_arb384'])-Ha).contains(0)
  grad=[arb(v)-w for v,w in zip(ga['gradient_x2_to_x6'],gradT)]
  LB=endpoint(Ha-HT);grid=sorted(t for t in times if a<=t<=T);assert grid[0]==a and grid[-1]==T
  bands=[]
  for t,u in zip(grid,grid[1:]):
   delta=times[t]-times[u];low=max(Q(0),endpoint(delta));high=endpoint(delta,True);assert high>=low>=0
   bands.append((t,u,low,high))
  for radius in RADII:
   bounds,vertices=cell(c,radius);hs=[tuple(v[i]-c[i]for i in range(6))for v in vertices]
   r=max(max(map(abs,h))for h in hs);D=max(map(normsq,hs));assert r==radius and D==5*r*r
   linear=[endpoint(sum((g*A(h[i+1])for i,g in enumerate(grad)),arb(0)))for h in hs]
   lin=min(linear);idx=linear.index(lin)
   contributions=[A(high)/(A(t)*A(t))*phi(A(r)/A(t))for t,u,low,high in bands]
   U=sum(contributions,arb(0));Uupper=endpoint(U,True)
   lower=LB+lin-D*Uupper/32
   row={'a':str(a),'T':str(T),'radius':str(r),'cell_free_bounds':[[str(lo),str(hi)]for lo,hi in bounds],
    'cell_vertices':len(vertices),'center_in_cell_exact':True,'cell_in_ordered_face_exact':True,
    'max_displacement_norm_squared':str(D),'B_center_lower_exact':str(LB),'min_vertex_linear_lower_exact':str(lin),
    'min_linear_vertex':list(map(str,vertices[idx])),'J2_cell_upper_arb384':U.str(100,more=True),
    'J2_cell_upper_exact':str(Uupper),'Taylor_curvature_loss_upper_exact':str(D*Uupper/32),
    'K_cell_lower_exact':str(lower),'K_cell_lower_arb384':A(lower).str(100,more=True),
    'passes_matched_m0':lower>m0,'passes_bare64_threshold':lower>bare,
    'radius_over_old_Phi_radius_exact':str(r/oldr),
    'time_bands':[{'left':str(t),'right':str(u),'mass_lower':str(lo),'mass_upper':str(hi),
      'J2_upper_term_arb384':value.str(95,more=True)}for (t,u,lo,hi),value in zip(bands,contributions)],
    'infinite_time_tail_added':False}
   rows.append(row);print('a',a,'r',r,'lower',row['K_cell_lower_arb384'][:55],
    'PASS_m0',row['passes_matched_m0'],'PASS64',row['passes_bare64_threshold'],flush=True)
 assert len(rows)==8
 toy=toy_checks();T48=Q(6,125);tube48=Q(1,5000000000)*(100000-40000)
 assert tube48==Q(3,250000)>T48/4096
 out={'status':'COMPLETE_FIXED_EIGHT_CASE_LOCAL_CALIBRATION_ONLY','pid':os.getpid(),'arb_bits':384,
  'kernel_calls':0,'new_gradient_calls':0,'MC_samples':0,'matlab_calls':0,'fixed_cases':8,
  'data_potential_family':'H','center':list(map(str,c)),'matched_m0':str(m0),'bare64_threshold':str(bare),
  'old_Phi_radius_exact':str(oldr),'old_Phi_radius_target':str(m0),'rows':rows,
  'pass_count_matched_m0':sum(row['passes_matched_m0']for row in rows),
  'pass_count_bare64':sum(row['passes_bare64_threshold']for row in rows),'simplex_toys':toy,
  'T48_tube_only':{'T':str(T48),'distance_upper':'1/100000','K_lower':str(tube48),'bare_threshold':str(T48/4096)},
  'global_cover_executed':False,'unresolved_root':True,'whole_cube_mT_certified':False,
  'T032_already_refuted_elsewhere':True,'new_C6_bound':False,
  'input_sha256':meta['input_sha256'],'source_sha256':meta['source_sha256'],'elapsed_seconds':perf_counter()-started}
 (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n');print(out['status'],'seconds',out['elapsed_seconds'])
if __name__=='__main__':main()
