"""Author structural/source/durable-prefix gates, no repeated corner arithmetic."""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict,Counter
import json,os,time,subprocess
from common import sha,digest,write
from shard_worker import recover
H=Path(__file__).resolve().parent;S=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch');W=S/'upper_parallel137';O=W/'run001'
def nodekey(n):return n[0],n[1]
def physical(n,roots):
 lo=[F(v,32) for v in roots[n[0]]];hi=[v+F(1,32) for v in lo]
 for j,b in enumerate(n[1]):
  k=j%4;m=(lo[k]+hi[k])/2
  if b=='0':hi[k]=m
  else:assert b=='1';lo[k]=m
 return lo,hi
def main():
 began=time.monotonic();assert not (W/'author_validation.json').exists();m=json.load(open(H/'prepared_manifest.json'));assert sha(H/'prepared_manifest.json')=='3b87ed39e6d1d0f7b7e9c1ad880d0e60708a8126ab4ab14a779da59524eab6c3'
 for n,h in m['repository_files'].items():assert sha(H/n)==h
 for p,h in {**m['source_locks'],**m['prepared_external_files']}.items():assert sha(p)==h
 write(W/'author_validation_registration.json',{'pid':os.getpid(),'code_sha256':sha(__file__),'prepared_manifest_sha256':sha(H/'prepared_manifest.json'),'repeat_corner_arithmetic':False,'new_H_gradient_kernel_log_exp_calls':0})
 plan=json.load(open(W/'shard_plan.json'));summary=json.load(open(O/'master_attempt0001_summary.json'));master=json.load(open(O/'master_attempt0001_state.json'));assert summary['status']=='complete_fixed89426_objects' and summary['master_fault'] is None and master['faults'] is None
 children=json.load(open(S/'upper_crossanchor131/child_inputs.json'));qual={};best={};counts=Counter();allpairs=set()
 for spec in plan['shards']:
  data=json.load(open(spec['input_path']));out=O/f'shard{spec["shard_id"]:02d}';cp=json.load(open(out/'checkpoint.json'));rec=recover(out/'results.jsonl',out/'checkpoint.json',data['pairs'],cp['fingerprint'],True);assert rec=={k:v for k,v in cp.items() if k!='status'} and cp['status']=='complete';assert cp['next_pair']==spec['pair_count']
  assert sha(out/'checkpoint.json')==summary['worker_summaries'][spec['shard_id']]['checkpoint_sha256'] and sha(out/'results.jsonl')==summary['worker_summaries'][spec['shard_id']]['results_sha256']
  for raw in map(json.loads,open(out/'results.jsonl')):
   i=raw['global_pair_index'];assert i not in allpairs;allpairs.add(i);counts[raw['kind']]+=1
   if raw['kind']=='evaluated':
    cert=raw['certificate'];source=data['sources'][raw['source_id']];core=source['core'];assert cert['source_pool_object_sha256']==digest(source)
    for key in ('epsilon','point','target','Delta_lower','loggradient_directed','curvature_rule','penalty_coefficient'):assert cert[key]==core[key]
    assert cert['checked_points']==16 and cert['passed']==(F(cert['minimum_log_margin'])>0)
    k=raw['child_index'];b=best.get(k)
    if b is None or F(cert['minimum_log_margin'])>F(b['minimum_log_margin']):best[k]=cert
  for k,c in cp['qualified_children'].items():assert k not in qual;qual[k]=c
 assert allpairs==set(range(89426)) and qual==master['qualified_children'] and len(qual)==1153
 assert counts=={'evaluated':71812,'skipped':17614}
 pending=[r for r in children if str(r['child_index']) not in qual];assert [r['node'] for r in pending]==master['pending_children'] and len(pending)==4274
 seed=json.load(open(S/'upper_continue123/seed.json'));state=json.load(open(S/'upper_continue123/batch0002/state.json'));roots=seed['roots'];assert roots==[[a,b,c,d] for a in range(32) for b in range(a+1) for c in range(b+1) for d in range(c+1)]
 assert state['source_seed_sha256']==sha(S/'upper_continue123/seed.json') if 'source_seed_sha256' in state else True
 sourcepassed=[r['node'] for r in seed['passed']]+[r['node'] for r in state['new_passed']];splits=seed['splits']+state['new_splits'];oldqueue=state['queue'];results=list(map(json.loads,open(S/'upper_curvature128/evaluation001/candidate_results.jsonl')));assert {nodekey(r['child_node']) for r in results}=={nodekey(n) for n in oldqueue} and len(results)==len(oldqueue)==11330
 assert {nodekey(r['node']) for r in children}=={nodekey(r['child_node']) for r in results if not r['passed']}
 passed=sourcepassed+[r['child_node'] for r in results if r['passed']]+[c['child_node'] for c in qual.values()];queue=[r['node'] for r in pending];by=defaultdict(list);leaves=set();inner=set()
 for n in passed+queue:
  k=nodekey(n);assert k not in leaves;leaves.add(k);by[n[0]].append(n[1])
 for r in splits:
  k=nodekey(r['node']);assert k not in inner;inner.add(k);root,path=k;assert r['axis']==len(path)%4 and r['children']==[[root,path+'0'],[root,path+'1']]
 assert len(passed)==105635 and len(queue)==4274 and len(inner)==57549 and len(leaves)==52360+len(inner) and not leaves&inner and len(by)==52360
 for root in range(52360):
  paths=sorted(by[root]);assert all(not b.startswith(a) for a,b in zip(paths,paths[1:]));depth=max(map(len,paths));assert sum(1<<(depth-len(p)) for p in paths)==1<<depth
 for r,p in leaves|inner:
  if p:assert (r,p[:-1]) in inner
 for r in children:
  lo,hi=physical(r['node'],roots);assert list(map(str,lo))==r['box_lower'] and list(map(str,hi))==r['box_upper']
 ordered=[]
 for child in pending:
  c=best[child['child_index']];assert not c['passed'];volume=F(1,32**4*2**len(child['node'][1]));ordered.append({'child_index':child['child_index'],'node':child['node'],'box_lower':child['box_lower'],'box_upper':child['box_upper'],'physical_volume':str(volume),'best_crossanchor_failed_certificate':c,'prior_direct_parent_newq_failure':child['prior_direct_parent_newq_failure']})
 ordered.sort(key=lambda r:(F(r['best_crossanchor_failed_certificate']['minimum_log_margin']),-F(r['physical_volume']),*r['node']))
 with open(W/'remaining_gap_order.jsonl','x') as f:
  for i,r in enumerate(ordered):r['gap_rank']=i;f.write(json.dumps(r,sort_keys=True,separators=(',',':'))+'\n')
 pids=[summary['master_pid']]+[r['pid'] for r in summary['process_evidence']];alive=[]
 for pid in pids:
  try:os.kill(pid,0);alive.append(pid)
  except ProcessLookupError:pass
 assert not alive and all(r['exitcode']==0 and not r['alive'] for r in summary['process_evidence'])
 result={'status':'PASS_full_sources_prefix_objects_closedtree_and_process_exit','pid':os.getpid(),'actual_master_pid':summary['master_pid'],'actual_worker_pids':pids[1:],'all16_production_PIDs_OS_absent':True,'production_wall_seconds':summary['wall_seconds'],'fixed_pairs':89426,'single_source_tests':71812,'exact_corners':1148992,'skip_receipts':17614,'new_passes':1153,'passed':105635,'pending':4274,'splits':57549,'active_leaves':109909,'ordered_roots':52360,'exact_root_Kraft_and_closed_box_mapper':True,'source_cache_and_code_hashes_unchanged':True,'faults':0,'new_H_gradient_kernel_log_exp_calls':0,'repeated_corner_checks':0,'remaining_gap_order_sha256':sha(W/'remaining_gap_order.jsonl'),'master_state_sha256':sha(O/'master_attempt0001_state.json'),'master_summary_sha256':sha(O/'master_attempt0001_summary.json'),'globalg110_qualified':False,'S_credit':0,'wall_seconds':time.monotonic()-began};write(W/'author_validation.json',result);print(json.dumps(result),flush=True)
if __name__=='__main__':main()
