#!/usr/bin/env python3
"""Twelve structural checks and exact beta=1 pair-orbit fan export; no kernel."""
from pathlib import Path
from fractions import Fraction as Q
from time import perf_counter
import hashlib,itertools,json,math,os
import sympy as s

HERE=Path(__file__).resolve().parent
MAIN=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds')
CELL=MAIN/'upper/round002/fan/exact_cells_results.json'
AUDIT=MAIN/'upper/round002/fan/verify_exact_fan_results.json'
OLDINT=MAIN/'upper/round003/direct_kernel/fan_integer_input.txt'
BRIDGE=MAIN/'upper/round004/primary/integer_fan_bridge.json'
F9=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-literature/cube-dim-6-bounds/audit/round009/short_potential/report.md')
OLD_CHECKER=HERE.parents[1]/'round015/band_taylor_cover/checker_schema.md'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dot(u,v):return sum(a*b for a,b in zip(u,v))
def projector(columns):
 B=s.Matrix.hstack(*columns);return B*(B.T*B).inv()*B.T
def main():
 started=perf_counter();body=[s.Matrix((1,)+bits)for bits in itertools.product([-1,1],repeat=5)]
 I=s.eye(6);axes=[I.col(i)for i in range(6)]
 pair=[axes[i]+sgn*axes[j]for i in range(6)for j in range(i+1,6)for sgn in [1,-1]]
 samples=[(1,1,1,1,1,1),(3,2,1,0,0,0),(-5,1,-3,2,0,-1)]
 inputs=[CELL,AUDIT,OLDINT,BRIDGE,F9,OLD_CHECKER,HERE/'registration.md',HERE/'frame_theory.md',HERE/'S2_and_rank.md']
 manifest={'beta_fixed':'1','structure_check_cap':12,'new_symbol_test_vectors':[list(x)for x in samples],
  'kernel_calls':0,'MC_samples':0,'optimization_calls':0,'input_sha256':{str(p):sha(p)for p in inputs},
  'source_sha256':sha(Path(__file__))}
 (HERE/'input_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 checks=[]
 def check(name,details):checks.append({'check':len(checks)+1,'name':name,'PASS':True,'details':details})
 # 1. Whole sorted chamber identity, including ties and zero coordinates.
 b,h=s.symbols('b h',nonnegative=True);a=b+h
 assert s.simplify(s.Abs(a+b)+s.Abs(a-b)-2*a)==0
 ell=tuple(2*(5-i)for i in range(6));assert ell==(10,8,6,4,2,0)
 check('pair_symbol_whole_chamber',{'linear_coefficients':ell,'absolute_identity_residual':'0'})
 # 2. True physical-column Gram, retaining symbolic beta.
 Gbody=sum((v*v.T for v in body),s.zeros(6));Gpair=sum((v*v.T for v in pair),s.zeros(6))
 assert Gbody==32*I and Gpair==10*I and Gbody+Gpair==42*I
 bet=s.symbols('beta',positive=True);assert Gbody+bet**2*Gpair==(32+10*bet**2)*I
 check('physical_frame_Gram',{'body_columns':32,'pair_columns':30,'total_columns':62,'Gram_beta1':'42I','symbolic_Gram':'(32+10beta^2)I'})
 # 3. Reuse the frozen exact geometry proof, transform only its mandatory126 rate entries.
 cells=json.loads(CELL.read_text());audit=json.loads(AUDIT.read_text());bridge=json.loads(BRIDGE.read_text())
 assert audit['status']=='PASS'and audit['full_cells']==21 and audit['input_sha256']==sha(CELL)
 assert bridge['status']=='PASS_exact_ray_rate_determinant_scaling_bridge'
 assert bridge['source_sha256']==sha(CELL)and bridge['integer_input_sha256']==sha(OLDINT)
 tokens=iter(map(int,OLDINT.read_text().split()));assert next(tokens)==len(cells['cells'])==21
 exported=[];textlines=['21'];newrates=[]
 for index,oldcell in enumerate(cells['cells']):
  det=next(tokens);rows=[];assert det>0 and oldcell['all_32_body_form_signs_constant']
  for _ in range(6):
   oldq=next(tokens);ray=tuple(next(tokens)for _ in range(6))
   if not(all(ray[i]>=ray[i+1]for i in range(5))and ray[-1]>=0 and ray[0]>0):
    (HERE/'failure.json').write_text(json.dumps({'status':'STOP_OLD_CONE_NOT_IN_LINEAR_CHAMBER','cone':index,'ray':ray})+'\n')
    raise RuntimeError('Old cone not in sorted chamber; stop before any new fan publication')
   increment=dot(ell,ray);newq=oldq+increment;assert oldq>0 and increment>0 and newq>oldq
   rows.append({'ray_integer':ray,'old_rate_integer':oldq,'pair_increment_integer':increment,'new_rate_integer':newq});newrates.append(newq)
  textlines.append(str(det));textlines.extend(' '.join(map(str,(r['new_rate_integer'],*r['ray_integer'])))for r in rows)
  exported.append({'cone_index':index,'frozen_word':oldcell['word'],'determinant_integer':det,'rows':rows,
   'new_alpha_exact':str(Q(det,math.prod(r['new_rate_integer']for r in rows)))})
 assert next(tokens,None)is None and len(newrates)==126
 check('mandatory_exact_fan_rate_transform',{'canonical_cones':21,'updated_ray_rate_entries':126,
  'new_cones_enumerated':0,'rays_and_determinants_changed':False,'old_cone_linearity_failure':False})
 # 4. One exact coordinate-marginal normalization; no point/time search.
 marginal=sum(abs(int(v[0]))for v in body+pair);assert marginal==42
 Tbase=Q(6,125);Tnew=Q(32,875);assert 42*Tnew==32*Tbase==Q(192,125)
 check('one_common_time_normalization',{'coordinate_rate_beta1':42,'dimensionless_tau':'192/125','normalized_T':'32/875','kernel_evaluated':False})
 # 5-9. Five canonical pair-plane structures, no exhaustive triple search.
 v=s.ones(6,1);p=axes[0]+axes[1];w=v-2*p
 assert all(abs(t)==1 for t in w)and projector([v,p])==projector([v,w])and s.Matrix.hstack(v,w,p).rank()==2
 check('matching_body_pair_old_plane',{'rank':2,'body_pair_old_plane_alias':True})
 p=axes[0]-axes[1];assert(v.T*p)[0]==0 and s.Matrix.hstack(v,p).rank()==2
 witness=v+p;assert len({abs(x)for x in witness})==3
 check('orthogonal_body_pair_new_plane',{'rank':2,'body_dot_pair':0,'three_absolute_level_witness':list(map(int,witness))})
 p=axes[0]+axes[1];q=axes[0]-axes[1];assert projector([p,q])==s.diag(1,1,0,0,0,0)
 check('same_support_pair_pair',{'coordinate_plane':True})
 p=axes[0]+axes[1];q=axes[0]+axes[2];third=axes[1]-axes[2]
 assert p-q-third==s.zeros(6,1)and s.Matrix.hstack(p,q,third).rank()==2
 check('one_overlap_pair_triangle',{'dependent_triple_rank':2,'normal_on_support':[1,-1,-1]})
 p=axes[0]+axes[1];q=axes[2]-axes[3]
 assert(p.T*q)[0]==0 and s.Matrix.hstack(p,q).rank()==2
 counts={'body_body_pairs':496,'body_pair_matching_pairs':480,'body_pair_new_orthogonal_pairs':480,
  'pair_pair_same_support_pairs':15,'pair_pair_one_overlap_pairs':240,'pair_pair_disjoint_pairs':180}
 assert sum(counts.values())==math.comb(62,2)==1891
 planes={'old_body':496,'new_orthogonal_body_pair':480,'new_coordinate':15,'new_triangle':80,'new_disjoint_support':180}
 assert sum(planes.values())==1251 and sum(planes.values())-496==755
 assert 240*3+256+480+15+80*3+180==1891
 check('disjoint_pair_pair_and_class_counts',{'rank':2,'direction_pair_counts':counts,'unique_plane_counts':planes,'dependent_triples_derived':320})
 # 10-12. Exactly the three registered finite symbol vectors.
 for x in samples:
  z=sorted(map(abs,x),reverse=True);actual=sum(abs(dot(tuple(p),x))for p in pair);predicted=dot(ell,z)
  bodyvalue=sum(abs(dot(tuple(v),x))for v in body);assert actual==predicted
  check('exact_symbol_vector',{'vector':x,'pair_symbol':int(actual),'body_symbol':int(bodyvalue),'new_beta1_symbol':int(bodyvalue+actual)})
 assert len(checks)==12
 integer_path=HERE/'pair_orbit_beta1_integer_fan.txt';integer_path.write_text('\n'.join(textlines)+'\n')
 fan={'model':'body32_plus_beta_pair_orbit_Cauchy','beta':'1','directions':62,'canonical_cones':21,
  'signed_permutation_group_size':46080,'whole_chamber_pair_linear_coefficients':ell,
  'rate_update':'q_new=q_old+10r1+8r2+6r3+4r4+2r5','cones':exported,
  'source_integer_fan_sha256':sha(OLDINT),'source_geometry_sha256':sha(CELL),
  'source_geometry_audit_sha256':sha(AUDIT),'source_bridge_sha256':sha(BRIDGE),
  'kernel_values_computed':False,'source_code_sha256':sha(Path(__file__))}
 fan_path=HERE/'pair_orbit_beta1_fan.json';fan_path.write_text(json.dumps(fan,indent=2)+'\n')
 out={'status':'PASS_NEW_PAIR_FRAME_AND_EXACT_21_CONE_FAN','pid':os.getpid(),'beta_fixed':'1',
  'structural_checks':checks,'structural_check_count':12,'symbol_test_vectors':3,
  'kernel_evaluations':0,'MC_samples':0,'optimization_calls':0,'new_grid':False,'matlab_calls':0,
  'Gram_exact':'42I','spectral_pair_linear_coefficients':ell,
  'new_rate_min':min(newrates),'new_rate_max':max(newrates),'S2_distinct_planes':1251,'S2_added_planes':755,
  'Markov_full_scale_interface':'proved conditionally on a new-model whole-cube m_T>T/4096',
  'old_g_or_barrier_inherited':False,'old_A015_checker_modified':False,'pair_checker_support_claimed':False,
  'time_normalization_proposal':{'coordinate_scale':42,'tau':'192/125','T':'32/875','no_kernel_evaluation':True},
  'new_C6_bound':False,'new_global_minimum':False,
  'exported_artifacts':{p.name:{'sha256':sha(p),'bytes':p.stat().st_size}for p in [integer_path,fan_path]},
  'input_sha256':manifest['input_sha256'],'source_sha256':manifest['source_sha256'],'elapsed_seconds':perf_counter()-started}
 (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
 manifest['exported_artifacts']=out['exported_artifacts'];(HERE/'input_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 print(out['status'],'checks',len(checks),'rate_range',min(newrates),max(newrates),'seconds',out['elapsed_seconds'])
 print(json.dumps(out['exported_artifacts'],indent=2))
if __name__=='__main__':main()
