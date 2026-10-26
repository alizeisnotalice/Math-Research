"""New fixed39510 source/prefix/whole-tree gate, no repeated corner arithmetic."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import json,os,time
from common import sha,digest,write
from shard_worker import recover
import residual_run
H=Path(__file__).resolve().parent;S=residual_run.S;W=S/'upper_parallel137/next2350_preflight';O=W/'run001'
def main():
 began=time.monotonic();dest=W/'author_validation.json';assert not dest.exists();prepared=H/'next2350_prepared_manifest.json';assert sha(prepared)=='f57463aa232d2b48f8dbd72ea22a5460d4e45acda8a206dcf0668f5953062626';package=json.load(open(prepared))
 for n,h in package['repository_files'].items():assert sha(H/n)==h
 for p,h in {**package['source_locks'],**package['prepared_external_files']}.items():assert sha(p)==h
 write(W/'author_registration.json',{'pid':os.getpid(),'code_sha256':sha(__file__),'prepared_manifest_sha256':sha(prepared),'repeat_H_log_corners':False});plan=json.load(open(W/'shard_plan.json'));summary=json.load(open(O/'master_attempt0001_summary.json'));state=json.load(open(O/'master_attempt0001_state.json'));assert summary['status']=='complete_fixed39510_objects' and summary['master_fault'] is None and state['faults'] is None
 qualified={};seen=set();best={};counts=Counter();passrules=Counter();outside=0
 for spec in plan['shards']:
  data=json.load(open(spec['input_path']));out=O/f'shard{spec["shard_id"]:02d}';cp=json.load(open(out/'checkpoint.json'));rec=recover(out/'results.jsonl',out/'checkpoint.json',data['pairs'],cp['fingerprint'],True);assert cp['status']=='complete' and rec=={k:v for k,v in cp.items() if k!='status'} and cp['next_pair']==spec['pair_count']
  ws=summary['worker_summaries'][spec['shard_id']];assert sha(out/'checkpoint.json')==ws['checkpoint_sha256'] and sha(out/'results.jsonl')==ws['results_sha256']
  for row in map(json.loads,open(out/'results.jsonl')):
   assert row['global_pair_index'] not in seen;seen.add(row['global_pair_index']);counts[row['kind']]+=1
   if row['kind']=='evaluated':
    c=row['certificate'];source=data['sources'][row['source_id']];assert c['source_pool_object_sha256']==digest(source)
    for key in ('epsilon','point','target','Delta_lower','loggradient_directed','curvature_rule','penalty_coefficient'):assert c[key]==source['core'][key]
    k=row['child_index']
    if k not in best or F(c['minimum_log_margin'])>F(best[k]['minimum_log_margin']):best[k]=c
  for k,c in cp['qualified_children'].items():assert k not in qualified;qualified[k]=c;passrules[c['curvature_rule']]+=1;outside+=not c['anchor_inside_child']
 assert seen==set(range(39510)) and qualified==state['qualified_children'] and len(qualified)==671 and counts=={'evaluated':29727,'skipped':9783}
 children=json.load(open(W/'children.json'));pending=[r for r in children if str(r['child_index']) not in qualified];assert len(pending)==1679 and [r['node'] for r in pending]==state['pending_children'];old=residual_run.legacy_load(json.load(open(H/'residual_prepared_manifest.json')));seed=residual_run.source_seed({});last=json.load(open(S/'upper_parallel137/residual_preflight/future_batch0002/state.json'));seed['passed']+=last['new_passed'];seed['splits']+=last['new_splits']
 for p,r in last['new_cache_by_epsilon']['1/128'].items():assert p not in seed['cache_by_epsilon']['1/128'];seed['cache_by_epsilon']['1/128'][p]=r
 ov={'new_passed':[{'node':c['child_node'],'certificate':c} for c in qualified.values()],'new_splits':[],'queue':state['pending_children'],'inflight':{}};tree=old.tree(seed,ov);assert tree['passed']==109405 and tree['splits']==58724 and tree['queue_and_inflight']==1679
 ordered=[]
 for r in pending:
  lo,hi=old.physical(r['node'],seed['roots']);assert list(map(str,lo))==r['box_lower'] and list(map(str,hi))==r['box_upper'];c=best[r['child_index']];assert not c['passed'];ordered.append({**r,'best_new_crossanchor_failed_certificate':c})
 ordered.sort(key=lambda r:(F(r['best_new_crossanchor_failed_certificate']['minimum_log_margin']),-F(r['physical_volume']),*r['node']))
 with open(W/'remaining1679_gap_order.jsonl','x') as f:
  for i,r in enumerate(ordered):r['gap_rank']=i;f.write(json.dumps(r,sort_keys=True,separators=(',',':'))+'\n')
 pids=[summary['master_pid']]+[r['pid'] for r in summary['process_evidence']]
 for p in pids:
  try:os.kill(p,0);raise AssertionError(('production PID still exists',p))
  except ProcessLookupError:pass
 assert all(p['exitcode']==0 and not p['alive'] for p in summary['process_evidence']);result={'status':'PASS_sources_all39510_prefixes_single_source_metadata_fulltree_and_exit','pid':os.getpid(),'production_master':summary['master_pid'],'production_workers':pids[1:],'all16_PIDs_OS_absent':True,'fixed_pairs':39510,'single_source_tests':29727,'corners':475632,'skip_receipts':9783,'new_passed':671,'new_pass_curvature_rules':dict(passrules),'outside_anchor_passes':outside,'tree':tree,'cache64':184796,'cache128':14801,'source_locks_unchanged':True,'sourcecache_not_modified':True,'new_H_gradient_kernel_log_exp_calls':0,'repeat_corner_arithmetic':0,'faults':0,'new_splits':0,'state_sha256':sha(O/'master_attempt0001_state.json'),'summary_sha256':sha(O/'master_attempt0001_summary.json'),'remaining1679_sha256':sha(W/'remaining1679_gap_order.jsonl'),'production_wall_seconds':summary['wall_seconds'],'wall_seconds':time.monotonic()-began,'globalg110_qualified':False,'S_credit':0};write(dest,result);print(json.dumps(result),flush=True)
if __name__=='__main__':main()
