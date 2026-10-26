"""Lock concrete new-pair inputs and exact dispatcher-only delta, no cert call."""
from pathlib import Path
import ast,difflib,json,os
from common import sha,digest,write,CORE_PATH,CORE_SHA
H=Path(__file__).resolve().parent;S=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch');W=S/'upper_parallel137/next2350_preflight';OLD=H.parent/'upper_crossanchor131'
def main():
 dest=H/'next2350_prepared_manifest.json';assert not dest.exists() and not (W/'run001').exists();prepared=json.load(open(OLD/'prepared_manifest.json'));oldpool=S/'upper_crossanchor131/source_pool.jsonl';refs=json.load(open(W/'source_pool_refs.json'));assert sha(oldpool)==refs['old_pool']['sha256'] and sha(refs['new_pool']['path'])==refs['new_pool']['sha256']
 plan=json.load(open(W/'shard_plan.json'));pairs=list(map(json.loads,open(W/'pairs.jsonl')));children=json.load(open(W/'children.json'));assert len(pairs)==plan['total_pairs']==39510 and len(children)==2350;allpairs={};allchildren={}
 for spec in plan['shards']:
  assert sha(spec['input_path'])==spec['input_sha256'];data=json.load(open(spec['input_path']));assert data['core_sha256']==CORE_SHA
  for p in data['pairs']:
   i=p['pair_index'];assert i not in allpairs and p==pairs[i] and p['child_index']%15==spec['shard_id'];allpairs[i]=p;source=data['sources'][p['source_id']];assert source['source_id']==digest(source['core']) and source['source_index']==p['source_index'] and source['core']['target']=='11/100000';assert source['core']['curvature_rule']==p['curvature_rule']
  for k,c in data['children'].items():assert k not in allchildren and c==children[int(k)];allchildren[k]=c
 assert set(allpairs)==set(range(39510)) and set(allchildren)==set(map(str,range(2350))) and sha(CORE_PATH)==CORE_SHA
 delta=H/'next_consumer_vs137.diff';assert not delta.exists();delta.write_text(''.join(difflib.unified_diff((H/'run.py').read_text().splitlines(True),(H/'next_consume.py').read_text().splitlines(True),fromfile='locked137/run.py',tofile='next_consume.py')))
 names=['common.py','shard_worker.py','prepare_next2350.py','prepare_next_shards.py','next_consume.py','next2350_interface.md','next_consumer_vs137.diff','freeze_next2350.py'];repo={n:sha(H/n) for n in names}
 for n in names:
  if n.endswith('.py'):ast.parse((H/n).read_text())
 result=H/'residual_results_manifest.json';locks={**prepared['source_locks'],str(OLD/'prepared_manifest.json'):sha(OLD/'prepared_manifest.json'),str(oldpool):sha(oldpool),str(CORE_PATH):CORE_SHA,str(result):sha(result)}
 for batch in (1,2):
  b=S/f'upper_parallel137/residual_preflight/future_batch{batch:04d}'
  for name in ['state.json','point_manifest.json','manifest.json','author_validation.json']:locks[str(b/name)]=sha(b/name)
 for p,h in locks.items():assert sha(p)==h
 external={str(p):sha(p) for p in sorted(W.iterdir()) if p.is_file()};write(dest,{'schema':'next2350_fixed39510_candidate_prepared_v1','status':'prepared_only_pending_primary_corner_dispatch','freeze_pid':os.getpid(),'repository_files':repo,'source_locks':locks,'prepared_external_files':external,'fixed_pairs':39510,'fixed_children':2350,'source_library':115605,'selected_distinct_sources':8296,'worker_count':15,'core_sha256':CORE_SHA,'source_passed':108734,'source_splits':58724,'root_count':52360,'new4274_numerical_qualification_dependency':'audit138 plus primary acceptance;140 is unrelated J audit','existing_137_zeroH_integer_qualification_dependency':'independent137 audit plus primary acceptance','source_math_core_unchanged':True,'actual_H_gradient_kernel_log_exp_corner_calls':0,'new_splits':0,'production_authorized':False,'candidate_only':True,'globalg110_qualified':False,'S_credit':0,'memory_estimate_GiB':[0.8,2.5],'memory_estimate_not_cap':True,'wall_deadline':None,'automatic_backup_H_dispatch':False});print(json.dumps({'manifest':str(dest),'sha256':sha(dest),'consumer_sha256':sha(H/'next_consume.py'),'pair_file_sha256':sha(W/'pairs.jsonl'),'source_locks':len(locks),'external_files':len(external),'core_calls':0}),flush=True)
if __name__=='__main__':main()
