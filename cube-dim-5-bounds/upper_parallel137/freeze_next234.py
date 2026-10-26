"""Lock exact234 frontier/source union and3904 fixed pair consumer, no core call."""
from pathlib import Path
import ast,difflib,json,os
from common import sha,digest,write,CORE_PATH,CORE_SHA
H=Path(__file__).resolve().parent;S=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch');W=S/'upper_parallel137/next234_preflight'
def main():
 dest=H/'next234_prepared_manifest.json';assert not dest.exists() and not (W/'run001').exists();refs=json.load(open(W/'source_pool_refs.json'))
 for key in ('old_pool','old137_new4274_pool','new_pool'):assert sha(refs[key]['path'])==refs[key]['sha256']
 pairs=list(map(json.loads,open(W/'pairs.jsonl')));children=json.load(open(W/'children.json'));plan=json.load(open(W/'shard_plan.json'));assert len(pairs)==3904 and len(children)==234 and plan['total_pairs']==3904 and plan['total_children']==234;seen=set();owned=set()
 for spec in plan['shards']:
  assert sha(spec['input_path'])==spec['input_sha256'];data=json.load(open(spec['input_path']))
  for p in data['pairs']:
   i=p['pair_index'];assert i not in seen and p==pairs[i] and p['child_index']%15==spec['shard_id'];seen.add(i);source=data['sources'][p['source_id']];assert digest(source['core'])==p['source_id'] and source['source_index']==p['source_index'] and source['core']['curvature_rule']==p['curvature_rule'] and source['core']['target']=='11/100000'
  for k,c in data['children'].items():assert k not in owned and c==children[int(k)];owned.add(k)
 assert seen==set(range(3904)) and owned==set(map(str,range(234))) and sha(CORE_PATH)==CORE_SHA
 diff=H/'next234_vs_previous_consumer.diff';assert not diff.exists();diff.write_text(''.join(difflib.unified_diff((H/'next_consume.py').read_text().splitlines(True),(H/'next234_consume.py').read_text().splitlines(True),fromfile='locked_next_consume.py',tofile='next234_consume.py')))
 names=['common.py','shard_worker.py','prepare_next234.py','prepare_next234_shards.py','next234_consume.py','next234_interface.md','next234_vs_previous_consumer.diff','freeze_next234.py'];repo={n:sha(H/n) for n in names}
 for n in names:
  if n.endswith('.py'):ast.parse((H/n).read_text())
 result=H/'H1679_results_manifest.json';m=json.load(open(result));locks={**m['source_locks'],**m['result_external_files'],str(result):sha(result),str(CORE_PATH):CORE_SHA}
 for key in ('old_pool','old137_new4274_pool'):locks[refs[key]['path']]=refs[key]['sha256']
 for p,h in locks.items():assert sha(p)==h
 external={str(p):sha(p) for p in W.iterdir() if p.is_file()};write(dest,{'schema':'next234_fixed3904_candidate_prepared_v1','status':'prepared_only_wait_primary_zeroH_dispatch_review','freeze_pid':os.getpid(),'repository_files':repo,'source_locks':locks,'prepared_external_files':external,'fixed_pairs':3904,'fixed_children':234,'source_library':117284,'selected_distinct_sources':2346,'worker_count':15,'core_sha256':CORE_SHA,'source_passed':110967,'source_splits':58841,'roots':52360,'cache64':184796,'cache128':16480,'numerical_dependency':'143 new1679 plus138 old4274 and primary acceptance','integer_dependency':'prior141 and separate newpair review','H_gradient_kernel_log_exp_corner_calls':0,'new_splits':0,'backup_H_points':234,'backup_H_dispatch_authorized':False,'production_authorized':False,'candidate_only':True,'globalg110_qualified':False,'S_credit':0,'wall_deadline':None,'memory_estimate_GiB':[0.3,1.0],'memory_estimate_not_cap':True});print(json.dumps({'manifest':str(dest),'sha256':sha(dest),'consumer_sha256':sha(H/'next234_consume.py'),'pairs_sha256':sha(W/'pairs.jsonl'),'source_locks':len(locks),'external_files':len(external),'actual_core_calls':0}),flush=True)
if __name__=='__main__':main()
