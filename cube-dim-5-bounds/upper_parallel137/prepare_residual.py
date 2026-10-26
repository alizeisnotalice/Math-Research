"""Concrete residual EPS128 point/time registration only; never evaluate H."""
from pathlib import Path
from fractions import Fraction as F
from math import lcm
from collections import Counter
import gzip,json,os,time,resource
from common import sha,digest,write
H=Path(__file__).resolve().parent;D=Path('/Users/zhengzhihao/Desktop/T');S=D/'cube-finite-dim-scratch';W=S/'upper_parallel137';O=W/'residual_preflight'
K=H.parent/'upper_mixed098/epsilon_kernel.py';KSHA='128d98c94e55632336edafe024aa3a27c1486898670a52be89cb39cc22ed6cc0';V=H.parent/'upper_g110107/worker_v2.py';VSHA='62375ec07bf0a030c397c4a925b0809ab1f625edf75a19f434e4eb7f8a492e8b'
TABLE=D/'Math-Research-cube5-astra-upper/cube-dim-5-bounds/astra_vertex017/boundary_triangulation.json.gz';TSHA='8dd2e9895d7465ac599d43feb1b5099506f4e43d1248215eba0a77e11f4d13f9'
def pk(p):return ','.join(map(str,map(F,p)))
def main():
 began=time.monotonic();assert not O.exists();O.mkdir();assert sha(K)==KSHA and sha(V)==VSHA and sha(TABLE)==TSHA
 result_manifest=H/'result_manifest.json';assert sha(result_manifest)=='16d98997cd04258b28bd997d5489c6a7ad4b52a359f56403e9943017bd11fd2e';result=json.load(open(result_manifest))
 for p,h in result['result_external_files'].items():assert sha(p)==h
 plan123=json.load(open(S/'upper_continue123/finite_plan.json'));refs=plan123['I_scalar_planned_registry_refs'];planned={};locks={str(K):KSHA,str(V):VSHA,str(TABLE):TSHA,str(result_manifest):sha(result_manifest),str(S/'upper_continue123/finite_plan.json'):sha(S/'upper_continue123/finite_plan.json')}
 for name,ref in refs.items():
  assert sha(ref['jobs'])==ref['jobs_sha256'] and sha(ref['registration'])==ref['registration_sha256'];reg=json.load(open(ref['registration']));times=list(map(F,reg['time_endpoints']));assert F(1,128) in times
  jobs=json.loads(gzip.decompress(Path(ref['jobs']).read_bytes()));planned[name]={pk(j['point']) for j in jobs};locks[ref['jobs']]=ref['jobs_sha256'];locks[ref['registration']]=ref['registration_sha256']
 for name,path in [('old107_all4973planned',S/'upper_g110107/planned_requests.json'),('old123_all6237planned',S/'upper_continue123/planned_tasks.json')]:
  planned[name]={pk(t['point']) for t in json.load(open(path))};locks[str(path)]=sha(path)
 sourcepool=S/'upper_crossanchor131/source_pool.jsonl';assert sha(sourcepool)=='56c0b6c376283f518f8ba848f9b84bec80215b7c57e7f019e378d552a1d01cd8';locks[str(sourcepool)]=sha(sourcepool)
 full={}
 for row in map(json.loads,open(sourcepool)):
  if row['core']['epsilon']=='1/128':
   pointkey=pk(row['core']['point']);assert pointkey not in full;full[pointkey]={'source_id':row['source_id'],'source_index':row['source_index'],'source_pool_object_sha256':digest(row),'record_sha256':row['raw_record_sha256'],'point':row['core']['point'],'curvature_rule':row['core']['curvature_rule'],'epsilon':'1/128','target':'11/100000'}
 assert len(full)==10527
 # Every tested source/region is a complete whole-box test, never corner switching.
 failed=set()
 for row in map(json.loads,open(S/'upper_curvature128/evaluation001/candidate_results.jsonl')):
  if not row['passed']:failed.add((row['child_node'][0],row['child_node'][1],row['point_key']))
 for path in (W/'run001').glob('shard*/results.jsonl'):
  for row in map(json.loads,open(path)):
   if row['kind']=='evaluated' and not row['certificate']['passed'] and row['certificate']['epsilon']=='1/128':failed.add((*row['certificate']['child_node'],pk(row['certificate']['point'])))
 write(O/'registration.json',{'pid':os.getpid(),'code_sha256':sha(__file__),'source_locks':locks,'residual_exact':4274,'target':'11/100000','epsilon':'1/128','proposed_worker_count':15,'finite_batches':[4096,178],'wall_deadline':None,'new_H_gradient_kernel_log_exp_calls':0,'actual_corner_checks':0,'dispatch_authorized':False,'selection':'exact residual gap order;midpoint then046 eta=minside/4096 beta1/8,1/16,1/32;completed typedfulljet reuse else avoid all frozen point/time registries','no_new_child_dispatch':True})
 import importlib.util
 spec=importlib.util.spec_from_file_location('readonly098_geometry_only',K);api=importlib.util.module_from_spec(spec);spec.loader.exec_module(api);geometry=api.FirstJetKernel(json.loads(gzip.decompress(TABLE.read_bytes())),F(1,128))
 residual=list(map(json.loads,open(W/'remaining_gap_order.jsonl')));assert len(residual)==4274;chosen=set();tasks=[];blocked=[];counts=Counter();collisions=Counter()
 for ordinal,row in enumerate(residual):
  lo=list(map(F,row['box_lower']));hi=list(map(F,row['box_upper']));mid=[(a+b)/2 for a,b in zip(lo,hi)];eta=min(b-a for a,b in zip(lo,hi))/4096;points=[[F(1),*mid]]+[[F(1),*[x+eta*beta**j for j,x in enumerate(mid)]] for beta in (F(1,8),F(1,16),F(1,32))];tests=[];selection=None;node=row['node']
  ancestor=row['prior_direct_parent_newq_failure']['H64_lower_source_node'];assert ancestor[0]==node[0] and node[1].startswith(ancestor[1]);assert F(row['prior_direct_parent_newq_failure']['H64_lower'])==F(21,200000)
  for j,point in enumerate(points):
   p=pk(point);g=geometry.eligible(point);assert all(a<x<b for a,x,b in zip(lo,point[1:],hi));cached=p in full;known=(*node,p) in failed;collision=[name for name,ps in planned.items() if p in ps]
   reason='generic280_guard' if not g['eligible'] else 'already_failed_same_region_EPS128_point' if known else 'duplicate_new_residual_point' if p in chosen else 'frozen_scalar_or_planned_H128_object' if collision and not cached else None
   tests.append({'candidate':j,'point':list(map(str,point)),'generic':g,'inside_exact_child':True,'completed_typed_EPS128_fulljet':cached,'known_same_region_failed_point':known,'registered_point_time_collisions':collision,'rejection':reason})
   if g['eligible']:
    for name in collision:collisions[name]+=1
   if reason is None:selection=(j,point,p,cached);break
  if selection is None:blocked.append({'ordinal':ordinal,'node':node,'candidate_tests':tests});continue
  j,point,p,cached=selection;chosen.add(p);counts['candidate'+str(j)]+=1;counts['typed_fulljet_hit' if cached else 'fresh_H128_request']+=1
  batch=ordinal//4096+1;local=ordinal%4096;tasks.append({'ordinal':ordinal,'batch':batch,'batch_ordinal':local,'node':node,'source137_child_index':row['child_index'],'source137_gap_rank':row['gap_rank'],'source137_failed_margin':row['best_crossanchor_failed_certificate']['minimum_log_margin'],'box_lower':row['box_lower'],'box_upper':row['box_upper'],'physical_volume':row['physical_volume'],'epsilon':'1/128','target':'11/100000','point':list(map(str,point)),'pair_cache_key':['1/128',p],'candidate':j,'original_H64_candidate':-1,'original_H64_point':None,'same_point_H64_guard':None,'candidate_tests':tests,'geometry_tests':[{'candidate':t['candidate'],**t['generic']} for t in tests],'cached':cached,'cache_reuse_source':full.get(p),'H64_g105_ancestor':{'node':ancestor,'lower':'21/200000','source':'locked128/130 H64_lower_source_node;selectedpoint inside descendant child','newpoint_inside_child':True,'child_is_descendant':True},'future_curvature':'129 anchor-fixed Jensen q=1024-768*(21/200000)/newH128_upper after fulljet;no mixed gradient','point_checkpoint':str(O/f'future_batch{batch:04d}'/f'point{local:04d}.json'),'stop_file':str(O/'STOP'),'selection_reason':'newchild center/generic anchor;no repeated failed parent or same-region EPS128 anchor','scope':'preflight only;not dispatch;no newborn child scheduling'})
 assert not blocked and len(tasks)==4274 and len(chosen)==4274
 write(O/'tasks.json',tasks);write(O/'blocked.json',blocked);write(O/'completed_fulljet_reuse_registry.json',full);batches=[]
 for i,start in enumerate((0,4096),1):
  part=tasks[start:start+4096];path=O/f'input_batch{i:04d}.json';write(path,part);batches.append({'batch':i,'count':len(part),'path':str(path),'sha256':sha(path),'maximum_new_H128':sum(not t['cached'] for t in part),'fulljet_hits':sum(t['cached'] for t in part),'worker_count':15,'fixed_input_only':True,'no_automatic_geometric_children':True})
 summary={'status':'concrete4274_residual_EPS128_inputs_prepared_NOT_dispatched','pid':os.getpid(),'actual_preflight_workers':1,'proposed_production_workers':15,'exact_residual':4274,'finite_batches':batches,'counts':dict(counts),'completed_fulljet_keys':10527,'planned_registry_counts':{k:len(v) for k,v in planned.items()},'eligible_candidate_collision_counts':dict(collisions),'same_region_failed_keys':len(failed),'new_selected_points_in_scalar_or_planned_registry':sum(not t['cached'] and any(t['pair_cache_key'][1] in ps for ps in planned.values()) for t in tasks),'all_selected_inside_actual_child_and_generic280':True,'all_selected_unique':True,'actual_corner_checks':0,'new_H_gradient_kernel_log_exp_calls':0,'dispatch_authorized':False,'wall_seconds':time.monotonic()-began,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 for p,h in locks.items():assert sha(p)==h
 write(O/'summary.json',summary);print(json.dumps(summary),flush=True)
if __name__=='__main__':main()
