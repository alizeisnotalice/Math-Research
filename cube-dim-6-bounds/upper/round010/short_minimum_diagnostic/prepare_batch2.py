#!/usr/bin/env python3
from pathlib import Path
from fractions import Fraction as F
from math import lcm,gcd
import hashlib,json,time
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parents[2]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
gr=json.loads((HERE/'gradient_results.json').read_text());p1=json.loads((HERE/'point1_results.json').read_text());o=json.loads((HERE/'other_time_results.json').read_text())
assert gr['exit_code']==p1['exit_code']==o['exit_code']==0
assert gr['new_full_passes']+p1['new_full_passes_this_point']+o['new_H_passes']==5
assert p1['H_rows'][0]['refutes_strict_64_target'] and not o['refutes_strict_64_target']
receipt={'status':'BATCH1_COMPLETED_AND_STOPPED_BEFORE_BATCH2_AUTHORIZATION','new_full_passes':5,'new_spatial_points':1,'gradients':{'pid':gr['pid'],'seconds':gr['elapsed_seconds'],'exit_code':0},'point1_G_H064':{'pid':p1['pid'],'seconds':p1['elapsed_seconds'],'exit_code':0},'same_point_H032':{'pid':o['pid'],'seconds':o['elapsed_seconds'],'exit_code':0},'batch_wall_seconds':o['batch_elapsed_seconds'],'conclusion':'T=.064 refuted for <64; T=.032 refuted for 62.5 only. All original batch handles stopped.','source_results':[{'file':f,'sha256':sha(HERE/f)}for f in ['gradient_results.json','point1_results.json','other_time_results.json']]}
(HERE/'batch1_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
c=[F(v,gr['denominator'])for v in gr['center_numerators']];x=[F(v,p1['denominator'])for v in p1['point_numerators']]
y=[F(1),x[1]]+[x[i]+(x[i]-c[i])/50 for i in range(2,6)]
assert all(y[i]>y[i+1]>0 for i in range(5)) and y[1]<1
D=lcm(*(v.denominator for v in y));ns=[int(v*D)for v in y]
assert all(F(n,D)==v for n,v in zip(ns,y))
matches=[];token=str(ns[2])
for root in [PROJECT/'upper',Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-theory/cube-dim-6-bounds/upper'),Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper')]:
 for p in root.rglob('*.json'):
  if HERE in p.parents:continue
  if token in p.read_text():matches.append({'path':str(p),'sha256':sha(p)})
assert not matches,'Possible duplicate candidate'
out={'status':'BATCH2_POINT_REGISTERED_BEFORE_KERNEL','point_index':2,'T':'4/125','candidate_numerators':ns,'denominator':D,'coordinate_fractions':list(map(str,y)),'approximate_coordinates_for_display':[float(v)for v in y],'rule':'x1=1; x2 fixed at point1; for i=3..6 y_i=x_i+(x_i-c_i)/50 exactly','extension_fraction':'1/50','strict_order_and_cube_verified':True,'duplicate_scan_token':token,'duplicate_matches':matches,'max_new_passes':2,'only_new_passes':['exact G','H_4/125'],'max_seconds':1200,'batch2_start_epoch':time.time(),'first_batch_kept_completed':True,'inputs':[{'path':str(p.resolve()),'sha256':sha(p)}for p in [HERE/'batch1_receipt.json',HERE/'gradient_results.json',HERE/'point1_results.json',HERE/'other_time_results.json',Path(__file__)]]}
(HERE/'point2.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
