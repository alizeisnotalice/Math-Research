"""Coordinate selection only: no kernel imports, values, optimization or tests."""
from pathlib import Path
from fractions import Fraction as Q
from hashlib import sha256
import json
ROOT=Path(__file__).resolve().parent
source=ROOT.parent/'round03/first_checker_summary.json'
d=json.loads(source.read_text())
points=[]
for i,box in enumerate(d['first_K_fail_examples'][:5]):
 u0,u1,v0,v1=map(Q,box[:4]);u=(u0+u1)/2;v=(v0+v1)/2
 selection='rational_box_center'
 if u+v>1:
  v=1-u;selection='u_midpoint_v_equal_1_minus_u'
 assert u0<=u<=u1 and v0<=v<=v1 and 0<v<=u and u+v<=1
 points.append({'id':f'K{i+1}','quantity':'K','threshold':'1','u':str(u),'v':str(v),
                'source_kind':'first_K_fail_examples','source_index':i,'box':box[:4],
                'failed_box_lower_bound':box[4],'selection':selection,'MATLAB_subset':True})
for i,box in enumerate(d['first_A_fail_examples'][:5]):
 u0,u1,v0,v1=map(Q,box[:4]);u=(u0+u1)/2;v=(v0+v1)/2
 assert 0<v<=u and u+v!=Q(33,20)
 points.append({'id':f'A{i+1}','quantity':'AK','threshold':'0','u':str(u),'v':str(v),
                'source_kind':'first_A_fail_examples','source_index':i,'box':box[:4],
                'failed_box_lower_bound_display':box[4],'selection':'rational_box_center','MATLAB_subset':False})
for j,i in enumerate((0,4)):
 box=d['first_K_fail_examples'][i];u0,u1,v0,v1=map(Q,box[:4])
 left,right=max(u0,1-v1),min(u1,1-v0)
 assert left<=right
 u=(left+right)/2;v=1-u
 assert u0<=u<=u1 and v0<=v<=v1 and 0<v<=u
 points.append({'id':f'D{j+1}','quantity':'K','threshold':'1','u':str(u),'v':str(v),
                'source_kind':'first_K_fail_examples','source_index':i,'box':box[:4],
                'failed_box_lower_bound':box[4],'selection':'midpoint_of_box_intersection_with_u_plus_v_equal_1',
                'MATLAB_subset':False})
assert len(points)==12 and len({(p['quantity'],p['u'],p['v']) for p in points})==12
assert all(0<Q(p['v'])<=Q(p['u']) and 0<Q(p['u'])+Q(p['v'])<Q(33,20) for p in points)
kernels=[]
for label,path in [('K1',ROOT.parent/'round01/scaled_certificate.json'),('K3',ROOT.parent/'round03/candidate.json')]:
 kernels.append({'id':label,'file':str(path.relative_to(ROOT.parent)), 'sha256':sha256(path.read_bytes()).hexdigest()})
output={'status':'coordinates_only_no_kernel_evaluation','source_summary_sha256':sha256(source.read_bytes()).hexdigest(),
        'selection':'first five K centers per MATLAB plan, first five A centers, and fixed first/fifth K diamond intersections',
        'points':points,'kernels':kernels,'point_count':12,'K_points':7,'AK_points':5,
        'implementations':['original_hlcert_adapter','independent_hl_core_point_methods'],
        'maximum_calls':48,'call_unit':'one selected scalar K or AK at one point for one kernel in one implementation'}
target=ROOT/'points.json'
if target.exists():raise SystemExit('points.json already exists; do not silently select a new list')
target.write_text(json.dumps(output,indent=2)+'\n')
print('Prepared 12 exact rational points, 0 kernel evaluations.')
