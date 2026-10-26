#!/usr/bin/env python3
"""Preparation only: exact generic-anchor/cache lookups; no H or log calls."""
from pathlib import Path
from fractions import Fraction as F
import gzip,hashlib,importlib.util,json,os,time
SOURCE=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_g105086/final/state.json')
SOURCE_SHA='8647a4543ef6a869e747ec06d93f5137cf11feb7b3a64a002995d25a681c4f9a'
MANIFEST=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-upper/cube-dim-5-bounds/upper_g105086/manifest.json')
MANIFEST_SHA='10461718679541d040e06241896c6f01da54165355f57ad6441084d0ddd9cbab'
API=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-upper/cube-dim-5-bounds/upper_cover046/arb_cover.py')
API_SHA='df3fb4165aea29e4aa3d1fb7ecacdb9cc5ad4a531bc03e3cc43a0820671fe028'
TABLE=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-upper/cube-dim-5-bounds/astra_vertex017/boundary_triangulation.json.gz')
TABLE_SHA='8dd2e9895d7465ac599d43feb1b5099506f4e43d1248215eba0a77e11f4d13f9'
OUT=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_g105089')
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  while b:=f.read(8*1024*1024):h.update(b)
 return h.hexdigest()
def write(p,x):p.write_text(json.dumps(x,sort_keys=True,indent=2)+'\n')
def main():
 began=time.monotonic();OUT.mkdir(exist_ok=True);assert sha(SOURCE)==SOURCE_SHA and sha(MANIFEST)==MANIFEST_SHA and sha(API)==API_SHA and sha(TABLE)==TABLE_SHA
 write(OUT/'preflight_registration.json',{'pid':os.getpid(),'code':str(Path(__file__).resolve()),'code_sha256':sha(__file__),
  'source':str(SOURCE),'source_sha256':SOURCE_SHA,'source086_manifest_sha256':MANIFEST_SHA,'target':'21/200000','epsilon':'1/64',
  'finite_cases':{'pending_nodes':1856,'immediate_children':3712},'purpose':'prepare fixed generic anchors and cache-hit census only;do not qualify children or launch kernels',
  'new_H_physical_gradient_kernel_log_calls':0,'production_authorized':False})
 state=json.load(open(SOURCE));assert [len(state[k]) for k in ('passed','splits','queue','cache')]==[96104,45600,1856,180040]
 assert state['epsilon']==state['cache_model_epsilon']=='1/64' and state['target']=='21/200000'
 spec=importlib.util.spec_from_file_location('preflight_locked046',API);api=importlib.util.module_from_spec(spec);spec.loader.exec_module(api)
 kernel=api.FirstJetKernel(json.loads(gzip.decompress(TABLE.read_bytes())),F(1,64));pending={tuple(r['node']):r for r in state['pending_certificates086']}
 rows=[];childrows=[]
 def inspect(node):
  point,candidate,tests=api.choose_anchor(api.box_of(node,state['roots']),kernel)
  assert point is not None;anchor=api.pointkey(point);cached=anchor in state['cache']
  return {'node':node,'anchor':anchor,'point':list(map(str,point)),'candidate':candidate,'geometry_tests':tests,'cached':cached,
   'H_interval_if_cached':state['cache'][anchor]['values'][0] if cached else None}
 for node in state['queue']:
  row=inspect(node);old=pending[tuple(node)];row['existing_g105_failed_certificate_matches_anchor']=old['certificate']['anchor']==row['anchor']
  assert old['certificate']['target']=='21/200000' and old['certificate']['passed'] is False
  rows.append(row)
  for bit in '01':childrows.append(inspect([node[0],node[1]+bit]))
 assert len({tuple(r['node']) for r in childrows})==3712
 write(OUT/'pending_anchor_preflight.json',rows);write(OUT/'child_anchor_preflight.json',childrows)
 summary={'status':'preparation_only_no_H','pid':os.getpid(),'source_sha256':SOURCE_SHA,'target':'21/200000','epsilon':'1/64',
  'pending':len(rows),'pending_cached':sum(r['cached'] for r in rows),'pending_matching_existing_failed_certificate':sum(r['existing_g105_failed_certificate_matches_anchor'] for r in rows),
  'immediate_children':len(childrows),'child_cached':sum(r['cached'] for r in childrows),'child_not_cached':sum(not r['cached'] for r in childrows),
  'all_candidates_pass_280_generic_guard':True,'new_H_physical_gradient_kernel_log_calls':0,'wall_seconds':time.monotonic()-began,
  'production_authorized':False,'pending_preflight_sha256':sha(OUT/'pending_anchor_preflight.json'),'child_preflight_sha256':sha(OUT/'child_anchor_preflight.json')}
 write(OUT/'preflight_summary.json',summary);assert sha(SOURCE)==SOURCE_SHA;print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
