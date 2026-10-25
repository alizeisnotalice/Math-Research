#!/usr/bin/env python3
"""Read-only arithmetic replay of the saved local certificate; no enumeration/kernel."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json,math,os,time
import sympy as s
r=Path(__file__).resolve().parent;start=time.time()
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
inp=json.loads((r/'input_manifest.json').read_text())
for p,h in inp['source_sha256'].items():assert sha(p)==h,p
fan=json.loads((r/'local_fan.json').read_text())
b=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-upper-astra/cube-dim-6-bounds/upper/round044/sparse4_frame_preflight')
f=json.loads((b/'frame.json').read_text());A=f['column_numerators'];den=f['common_denominator']
assert fan['frame_sha256']==sha(b/'frame.json') and den==80
U=[tuple(F(x) for x in v['u']) for v in fan['vertices']]
H=fan['boundary_forms_in_u'];R=s.Matrix.hstack(*map(s.Matrix,fan['parent_rays']))
assert R.det()==1
normal=R.T.inv()*s.ones(6,1)
assert normal==s.Matrix([0,-1,1,1,0,0]) and (normal.T*normal)[0]==3
assert all(h[4:]==[0,0] for h in H[6:])
dot=lambda a,b:sum(x*y for x,y in zip(a,b))
def rank(S):
 if not S:return -1
 o=U[S[0]]
 return s.Matrix([[a-b for a,b in zip(U[i],o)] for i in S[1:]]).rank() if len(S)>1 else 0
for v,u in zip(fan['vertices'],U):
 assert min(u)>=0 and sum(u)==1
 active=v['active_witness'];M=s.Matrix([[1]*6]+[H[j] for j in active]);assert M.det()!=0
 assert all(dot(H[j],u)==0 for j in active)
 pr=v['primitive_physical_ray'];scale=F(v['slice_scale'])
 x=R*s.Matrix(u);assert all(x[j]==scale*pr[j] for j in range(6))
 assert math.gcd(*pr)==1 and pr[5]>=0 and all(pr[i]>=pr[i+1] for i in range(5))
 assert F(v['q_true_152'])==F(sum(abs(dot(a,pr)) for a in A),den)
faces={tuple(d['vertices']):d for d in fan['pulling_face_certificates']}
for S,d in faces.items():
 dim=rank(S);assert dim==d['dimension'];p=min(S);assert p==d['apex']
 expected_facets=set()
 if dim:
  for h in H:
   T=tuple(i for i in S if dot(h,U[i])==0)
   if T and len(T)<len(S) and rank(T)==dim-1:expected_facets.add(T)
 assert expected_facets==set(map(tuple,d['facets']))
 if dim==0:expected={(p,)}
 else:
  expected=set()
  for T in expected_facets:
   if p not in T:
    assert T in faces
    for q in faces[T]['simplices']:expected.add(tuple(sorted([p]+q)))
 actual={tuple(sorted(t)) for t in d['simplices']}
 assert actual==expected and len(actual)==len(d['simplices'])
 for t in d['simplices']:assert len(t)==dim+1 and rank(tuple(t))==dim
patterns=fan['closed_patterns'];full=[p for p in patterns if p['affine_dimension']==5]
for p in patterns:
 S=tuple(p['vertices']);assert rank(S)==p['affine_dimension']
 expected=tuple(i for i,u in enumerate(U) if all(z*dot(h,u)>=0 for z,h in zip(p['signs'],H[6:])))
 assert S==expected
 if 0<=p['affine_dimension']<5:
  q=patterns[p['containing_full_pattern']];assert q['affine_dimension']==5 and set(S)<=set(q['vertices'])
vol=F(0);simpkeys=set()
for p in full:
 S=tuple(p['vertices']);record=faces[S]
 expected={tuple(sorted(t)) for t in record['simplices']}
 saved={tuple(sorted(fan['simplices'][j]['vertices'])) for j in p['simplex_indices']}
 assert expected==saved
 pv=F(0)
 for j in p['simplex_indices']:
  t=fan['simplices'][j];T=t['vertices'];key=tuple(sorted(T));assert key not in simpkeys;simpkeys.add(key)
  d=abs(s.Matrix.hstack(*[s.Matrix(U[i]) for i in T]).det());assert d==F(t['relative_slice_volume']) and d>0
  pr=[fan['vertices'][i]['primitive_physical_ray'] for i in T]
  pd=s.Matrix.hstack(*map(s.Matrix,pr)).det();assert pd==F(t['primitive_physical_determinant'])
  assert abs(pd)*math.prod(F(fan['vertices'][i]['slice_scale']) for i in T)==d
  for a in A:
   values=[dot(a,x) for x in pr];assert not(min(values)<0<max(values))
  slope=[F(x) for x in p['linear_symbol_coefficient']]
  for i,x in zip(T,pr):assert dot(slope,x)==F(fan['vertices'][i]['q_true_152'])
  pv+=F(d)
 assert pv==F(p['relative_slice_volume']);vol+=pv
assert vol==1 and len(simpkeys)==len(fan['simplices'])==14
out={'status':'PASS_SAVED_FIRST_PARENT_CERTIFICATE','pid':os.getpid(),'elapsed_seconds':time.time()-start,'active_set_reenumerations':0,'new_sign_patterns':0,'full_patterns_rechecked':len(full),'pulling_face_nodes_rechecked':len(faces),'simplices_rechecked':len(simpkeys),'relative_slice_volume_exact':str(vol),'u_Euclidean_parent_volume':'sqrt(6)/120','physical_xi_slice_normal':[0,-1,1,1,0,0],'physical_xi_Euclidean_parent_volume':'sqrt(3)/120','cuts_ignore_u5_u6_join_structure':True,'fan_sha256':sha(r/'local_fan.json'),'checker_sha256':sha(Path(__file__)),'kernel_calls':0,'full_fan':False}
(r/'checked_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
