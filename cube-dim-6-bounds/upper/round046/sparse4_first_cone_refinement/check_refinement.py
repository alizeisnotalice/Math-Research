#!/usr/bin/env python3
"""Exact first-parent arrangement/pulling certificate, no float hull or kernel."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations,product
from functools import lru_cache
import json,hashlib,math,os,sys,time
import flint
from flint import fmpq_mat,fmpq
HERE=Path(__file__).resolve().parent
START=time.time();DEADLINE=1792964915 # replaced below by task-start UTC +1200
# Explicit UTC task deadline; no computation beyond this parent cone.
import datetime
DEADLINE=datetime.datetime(2026,10,25,21,48,35,tzinfo=datetime.timezone.utc).timestamp()
MAX_SIMPLICES=100000
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
inp=json.loads((HERE/'input_manifest.json').read_text())
for p,h in inp['source_sha256'].items():assert sha(p)==h,p
assert not (HERE/'local_fan.json').exists() and not (HERE/'results.json').exists()
def stopcheck():
 if time.time()>=DEADLINE:raise RuntimeError('registered task wall cap reached')
def mat(rows):return fmpq_mat([[fmpq(z.numerator,z.denominator) if isinstance(z,F) else z for z in row] for row in rows])
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def detcols(cols):return F(str(mat(list(map(list,zip(*cols)))).det()))
def rankrows(rows):
 basis={}
 for z in rows:
  z=list(z)
  for j,b in sorted(basis.items()):
   if z[j]:
    c=z[j];z=[x-c*y for x,y in zip(z,b)]
  j=next((j for j,a in enumerate(z) if a),None)
  if j is not None:
   c=z[j];basis[j]=[a/c for a in z]
 return len(basis)
b=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-upper-astra/cube-dim-6-bounds/upper/round044/sparse4_frame_preflight')
f=json.loads((b/'frame.json').read_text());cols=[tuple(x) for x in f['column_numerators']]
assert len(cols)==152 and f['common_denominator']==80 and f['beta']=='1/80'
old=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-upper-astra/cube-dim-6-bounds/upper/round003/direct_kernel/fan_integer_input.txt')
lines=old.read_text().splitlines();assert lines[0]=='21'
rays=[tuple(map(int,l.split()[1:])) for l in lines[2:8]]
assert rays==[(2,1,1,1,1,0),(3,3,2,2,2,0),(2,2,2,1,1,0),(1,1,1,1,0,0),(1,1,1,1,1,0),(1,1,1,1,1,1)]
parent_det=detcols(rays);assert parent_det==1
new_parent_dots=[tuple(dot(w,r) for r in rays) for w in cols[32:]]
cross=[i+1 for i,x in enumerate(new_parent_dots) if min(x)<0<max(x)]
assert cross==[2,10,26,50,82]
H=[tuple(int(i==j) for i in range(6)) for j in range(6)]+[new_parent_dots[i-1] for i in cross]
print(json.dumps({'phase':'first_cone_structure_pass','crossing_indices':cross,'cuts':[list(x) for x in H[6:]],'parent_new_form_dots':720}),flush=True)
# The active boundary matrices do not depend on the cell signs. Cache once.
vertices={};active_count=0;singular=0;outside=0
for active in combinations(range(11),5):
 stopcheck();active_count+=1
 M=mat([[1]*6]+[H[j] for j in active])
 if not M.det():singular+=1;continue
 y=M.solve(mat([[1]]+[[0]]*5))
 u=tuple(F(str(y[i,0])) for i in range(6))
 if min(u)<0:outside+=1;continue
 assert sum(u)==1 and all(dot(H[j],u)==0 for j in active)
 vertices.setdefault(u,active)
assert active_count==462
U=sorted(vertices);N=len(U);uid={u:i for i,u in enumerate(U)}
zero=[tuple(dot(h,u)==0 for h in H) for u in U]
@lru_cache(None)
def affrank(S):
 if not S:return -1
 o=U[S[0]]
 return rankrows([[a-b for a,b in zip(U[i],o)] for i in S[1:]])
patterns=[]
for p,sign in enumerate(product((-1,1),repeat=5)):
 stopcheck()
 S=tuple(i for i,u in enumerate(U) if all(sign[j]*dot(H[6+j],u)>=0 for j in range(5)))
 patterns.append({'id':p,'signs':list(sign),'vertices':list(S),'affine_dimension':affrank(S)})
full=[p for p in patterns if p['affine_dimension']==5]
assert full
face_records={};pull_cache={}
def facets(S):
 d=affrank(S);out=set()
 for j in range(11):
  T=tuple(i for i in S if zero[i][j])
  if T and len(T)<len(S) and affrank(T)==d-1:out.add(T)
 return sorted(out)
def pulling(S):
 stopcheck()
 if S in pull_cache:return pull_cache[S]
 d=affrank(S);assert d>=0
 p=min(S)
 if d==0:
  assert len(S)==1
  out=[(p,)];fs=[]
 else:
  fs=facets(S);assert fs
  opposite=[T for T in fs if p not in T];assert opposite
  out=[]
  for T in opposite:
   for simplex in pulling(T):
    q=(p,)+simplex
    assert len(q)==d+1 and affrank(tuple(sorted(q)))==d
    out.append(q)
    if len(out)>MAX_SIMPLICES:raise RuntimeError('registered simplex output cap reached')
  assert len({tuple(sorted(q)) for q in out})==len(out)
 face_records[S]={'vertices':list(S),'dimension':d,'apex':p,'facets':[list(T) for T in fs],'simplices':[list(q) for q in out]}
 pull_cache[S]=out
 return out
# Primitive physical rays and all true-frame rates. Never import old q.
vertex_records=[];dots=[];scales=[];primitive=[]
for i,u in enumerate(U):
 x=tuple(sum(F(rays[j][a])*u[j] for j in range(6)) for a in range(6))
 D=math.lcm(*(z.denominator for z in x));ints=[int(z*D) for z in x]
 g=math.gcd(*ints);r=tuple(z//g for z in ints);scale=F(g,D)
 assert g>0 and all(xj==scale*rj for xj,rj in zip(x,r))
 assert all(r[j]>=r[j+1] for j in range(5)) and r[5]>=0
 ds=tuple(dot(a,r) for a in cols);q=F(sum(abs(z) for z in ds),80)
 assert q>0
 dots.append(ds);scales.append(scale);primitive.append(r)
 vertex_records.append({'id':i,'u':[str(z) for z in u],'active_witness':list(vertices[u]),'primitive_physical_ray':list(r),'slice_scale':str(scale),'q_true_152':str(q)})
all_simp=[];seen=set();total=F(0)
for p in full:
 S=tuple(p['vertices']);ts=pulling(S)
 slope=[];form_signs=[]
 for j in range(152):
  vs=[dots[i][j] for i in S];assert not(min(vs)<0<max(vs))
  form_signs.append(1 if max(vs)>0 else (-1 if min(vs)<0 else 0))
 slope=[F(sum(form_signs[j]*cols[j][a] for j in range(152)),80) for a in range(6)]
 p['linear_symbol_coefficient']=[str(z) for z in slope]
 p['simplex_indices']=[];pvol=F(0)
 for simplex in ts:
  stopcheck()
  key=tuple(sorted(simplex));assert key not in seen;seen.add(key)
  v=abs(detcols([U[i] for i in simplex]));assert v>0
  dr=detcols([primitive[i] for i in simplex])
  assert abs(dr)*math.prod(scales[i] for i in simplex)==v
  for j in range(152):
   ds=[dots[i][j] for i in simplex];assert not(min(ds)<0<max(ds))
  for i in simplex:assert dot(slope,primitive[i])==F(vertex_records[i]['q_true_152'])
  for i in simplex:
   assert min(U[i])>=0 and sum(U[i])==1
   assert all(p['signs'][j]*dot(H[j+6],U[i])>=0 for j in range(5))
  p['simplex_indices'].append(len(all_simp))
  all_simp.append({'id':len(all_simp),'pattern':p['id'],'vertices':list(simplex),'relative_slice_volume':str(v),'primitive_physical_determinant':str(dr)})
  pvol+=v
  if len(all_simp)>MAX_SIMPLICES:raise RuntimeError('registered global simplex cap reached')
 p['relative_slice_volume']=str(pvol);total+=pvol
assert total==1
# Every lower-dimensional closed pattern lies in one full cell closure.
lower_maps=[]
for p in patterns:
 if 0<=p['affine_dimension']<5:
  S=set(p['vertices']);witness=next((q['id'] for q in full if S.issubset(q['vertices'])),None)
  assert witness is not None
  p['containing_full_pattern']=witness;lower_maps.append([p['id'],witness])
# Exact common-face compatibility; compare induced maximal simplices.
shared_checked=0
for p,q in combinations(full,2):
 stopcheck()
 S=tuple(sorted(set(p['vertices'])&set(q['vertices'])))
 if not S:continue
 d=affrank(S);assert d<5;shared_checked+=1
 SS=set(S)
 def induced(cell):
  ans=set()
  for j in cell['simplex_indices']:
   T=tuple(sorted(SS&set(all_simp[j]['vertices'])))
   if len(T)==d+1:assert affrank(T)==d;ans.add(T)
  return ans
 assert induced(p)==induced(q) and induced(p)
# Stronger per-cell volume receipts and structural summary.
fan={'scope':'ONLY_OLD_PARENT_CONE_1','frame_sha256':sha(b/'frame.json'),'beta':'1/80','parent_rays':[list(r) for r in rays],'parent_determinant':'1','normalization':'sum u=1; normalized slice volume(parent)=1; u-coordinate Euclidean parent volume=sqrt(6)/120; physical xi-slice Euclidean parent volume=sqrt(3)/120','crossing_new_form_indices_1_based':cross,'boundary_forms_in_u':[list(x) for x in H],'vertices':vertex_records,'closed_patterns':patterns,'pulling_face_certificates':list(face_records.values()),'simplices':all_simp}
(HERE/'local_fan.json').write_text(json.dumps(fan,separators=(',',':'))+'\n')
results={'status':'PASS_COMPLETE_EXACT_FIRST_PARENT_REFINEMENT_ONLY','pid':os.getpid(),'python':sys.version,'python_flint':flint.__version__,'active_systems_solved_or_classified':active_count,'singular_active_systems':singular,'outside_parent_active_solutions':outside,'registered_active_budget':14784,'sign_patterns':32,'full_dimensional_patterns':len(full),'lower_dimensional_nonempty_patterns':len(lower_maps),'empty_patterns':sum(p['affine_dimension']<0 for p in patterns),'unique_slice_vertices':N,'unique_primitive_child_rays':len(set(primitive)),'output_simplices':len(all_simp),'pulling_faces_recorded':len(face_records),'shared_faces_independently_compared':shared_checked,'parent_new_form_dots':720,'child_all_frame_ray_dots':N*152,'new_form_crossings':cross,'total_relative_volume_exact':str(total),'all152_forms_constant_on_every_simplex':True,'all_true_rates_verified':True,'coverage_proof':'closed halfspace partition + exact complete vertex/facet enumeration + recursive pulling; volume identity is redundant check','boundary_patterns_retained':True,'local_fan_sha256':sha(HERE/'local_fan.json'),'source_sha256':sha(Path(__file__)),'registration_sha256':sha(HERE/'registration.md'),'input_manifest_sha256':sha(HERE/'input_manifest.json'),'elapsed_seconds':time.time()-START,'full_fan_generated':False,'kernel_native_G_H_K_p_calls':0,'MATLAB_calls':0,'new_C6_bound':False}
(HERE/'results.json').write_text(json.dumps(results,indent=2)+'\n')
print(json.dumps(results,indent=2),flush=True)
