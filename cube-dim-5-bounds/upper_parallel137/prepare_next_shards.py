"""Freeze child-owned39510 pair shards and unchanged-core finite consumer."""
from pathlib import Path
import ast,json,os,time
from common import sha,digest,write,CORE_SHA
H=Path(__file__).resolve().parent;S=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch');W=S/'upper_parallel137/next2350_preflight';OLD=S/'upper_crossanchor131'
def main():
 began=time.monotonic();assert not (W/'shard_plan.json').exists();pool=W/'new_source_pool.jsonl';before=sha(pool);new=[];beforecore=[]
 for r in map(json.loads,open(pool)):beforecore.append(digest(r['core']));r['qualification']=r['qualification'].replace('audit pending','audit138 pending');new.append(r)
 # The creator's old source_refs used an incorrect140 label; math fields stay identical.
 refs=json.load(open(W/'source_pool_refs.json'));oldlabel=refs['new_source_independent_numeric_qualification'];refs['new_source_independent_numeric_qualification']='pending audit138/primary;preparation only'
 with open(pool.with_suffix('.tmp'),'x') as f:
  for r in new:f.write(json.dumps(r,sort_keys=True,separators=(',',':'))+'\n')
 os.replace(pool.with_suffix('.tmp'),pool);assert beforecore==[digest(r['core']) for r in new];refs['new_pool']['sha256']=sha(pool);write(W/'source_pool_refs.json',refs);write(W/'qualification_metadata_correction.json',{'old_newpool_sha256':before,'new_newpool_sha256':sha(pool),'creator_source_refs_label':oldlabel,'correct_numeric_qualification':'audit138 new4274 H128;140 belongs to unrelated J audit','all_4274_core_sha256_sequences_equal':True,'core_sequence_sha256':digest(beforecore),'new_H_gradient_kernel_log_exp_corner_calls':0})
 pairs=list(map(json.loads,open(W/'pairs.jsonl')));children=json.load(open(W/'children.json'));assert len(pairs)==39510 and len(children)==2350;need={r['source_id'] for r in pairs};sources={}
 for path in [OLD/'source_pool.jsonl',pool]:
  for r in map(json.loads,open(path)):
   if r['source_id'] in need:assert r['source_id'] not in sources;sources[r['source_id']]=r
 assert len(sources)==len(need)==8296;per=[[] for _ in range(15)]
 for i,p in enumerate(pairs):assert p['pair_index']==i;per[p['child_index']%15].append(p)
 plan=[];seen=set();owner={}
 for shard,pairlist in enumerate(per):
  ids=sorted({p['source_id'] for p in pairlist});ci=sorted({p['child_index'] for p in pairlist});obj={'shard_id':shard,'core_sha256':CORE_SHA,'source_pool_sha256':sha(W/'source_pool_refs.json'),'child_inputs_sha256':sha(W/'children.json'),'candidate_pairs_sha256':sha(W/'pairs.jsonl'),'source_pool_scope':'immutable131 plus new4274 jets;new numerical qualification pending138','pairs':pairlist,'children':{str(i):children[i] for i in ci},'sources':{sid:sources[sid] for sid in ids}};path=W/f'shard{shard:02d}_input.json';write(path,obj)
  for p in pairlist:assert p['pair_index'] not in seen;seen.add(p['pair_index']);assert owner.setdefault(p['child_index'],shard)==shard
  plan.append({'shard_id':shard,'input_path':str(path),'input_sha256':sha(path),'pair_count':len(pairlist),'child_count':len(ci),'source_count':len(ids),'input_bytes':path.stat().st_size})
 assert seen==set(range(39510)) and set(owner)==set(range(2350));write(W/'shard_plan.json',{'worker_count':15,'shards':plan,'total_pairs':39510,'total_children':2350,'unique_needed_sources':8296,'whole_child_unique_owner':True,'all_global_pairs_exactly_once':True,'core_sha256':CORE_SHA,'source137_newH_results_manifest_sha256':'86a877cab4efcd4f8b7c3d5a22196d6593bb81c01c0d4179129d0077810e7ad9','candidate_only_pending_numeric138':True})
 code=(H/'run.py').read_text()
 for a,b in [('89426','39510'),('5427','2350'),('104482','108734'),('57549','58724'),("scratch/upper_parallel137')","scratch/upper_parallel137/next2350_preflight')"),("'prepared_manifest.json'","'next2350_prepared_manifest.json'"),('source128_base','source137_newH_base')]:code=code.replace(a,b)
 assert "scratch/upper_parallel137/next2350_preflight')" in code and '89426' not in code and '104482' not in code;ast.parse(code);path=H/'next_consume.py';assert not path.exists();path.write_text(code)
 result={'status':'PASS_15_fixed39510_shards_AND_candidate_consumer_prepared_ONLY','pid':os.getpid(),'pairs':39510,'children':2350,'needed_sources':8296,'workers':15,'pair_counts':[r['pair_count'] for r in plan],'child_counts':[r['child_count'] for r in plan],'source_counts':[r['source_count'] for r in plan],'source_copies':sum(r['source_count'] for r in plan),'total_shard_bytes':sum(r['input_bytes'] for r in plan),'core_sha256':CORE_SHA,'worker_sha256':sha(H/'shard_worker.py'),'consumer_sha256':sha(path),'actual_H_gradient_kernel_log_exp_corner_calls':0,'dispatch_authorized':False,'numeric_qualification_dependency':'audit138 plus primary acceptance','wall_seconds':time.monotonic()-began};write(W/'shard_preparation_receipt.json',result);print(json.dumps(result),flush=True)
if __name__=='__main__':main()
