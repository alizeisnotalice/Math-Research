#!/usr/bin/env python3
"""Freeze15 child-owned shards, only selected10549 sources, zero corner tests."""
from pathlib import Path
from collections import Counter
import json,os,time,resource
from common import sha,digest,write,CORE_PATH,CORE_SHA
H=Path(__file__).resolve().parent;W=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_parallel137');OLD=H.parent/'upper_crossanchor131';S=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_crossanchor131')
OLD_SHA='1520dc778a8d92631ed33179690605a9eefdca6ae3e7c209dd6c8ffbcc4638f3'
def main():
 began=time.monotonic();W.mkdir(parents=True,exist_ok=True);assert not (W/'preparation_registration.json').exists();assert sha(OLD/'prepared_manifest.json')==OLD_SHA and sha(CORE_PATH)==CORE_SHA
 m=json.load(open(OLD/'prepared_manifest.json'))
 for name,h in m['repository_files'].items():assert sha(OLD/name)==h
 for p,h in m['prepared_external_files'].items():assert sha(p)==h
 for p,h in m['source_locks'].items():assert sha(p)==h
 write(W/'preparation_registration.json',{'task':'upper_parallel137_prepare_only','pid':os.getpid(),'code_sha256':sha(__file__),'core_sha256':CORE_SHA,'source131_prepared_manifest_sha256':OLD_SHA,'fixed_global_pairs':89426,'fixed_children':5427,'worker_count':15,'child_owner':'child_index modulo15;whole candidate set stays with one worker','source_loading':'one selected10549 scan in preparation;runtime each shard loads only its subset','peak_memory_estimate_GiB':[0.8,2.5],'memory_estimate_is_not_cap':True,'wall_deadline':None,'actual_new_corner_checks':0,'new_H_gradient_kernel_log_exp_calls':0,'dispatch_authorized':False})
 children=json.load(open(S/'child_inputs.json'));pairs=list(map(json.loads,open(S/'candidate_pairs.jsonl')));assert len(children)==5427 and len(pairs)==89426
 needed={r['source_id'] for r in pairs};sources={}
 for r in map(json.loads,open(S/'source_pool.jsonl')):
  if r['source_id'] in needed:sources[r['source_id']]=r
 assert len(sources)==len(needed)==10549
 per=[[] for _ in range(15)];owners={};seen=set()
 for p in pairs:
  child=p['child_index'];owner=child%15;assert owners.setdefault(child,owner)==owner and p['pair_index'] not in seen;seen.add(p['pair_index']);per[owner].append(p)
 assert seen==set(range(89426)) and set(owners)==set(range(5427))
 plan=[];total_sources=0
 for shard,pairlist in enumerate(per):
  ids=sorted({p['source_id'] for p in pairlist});ci=sorted({p['child_index'] for p in pairlist});obj={'shard_id':shard,'child_owner_rule':'child_index modulo15','source131_prepared_manifest_sha256':OLD_SHA,'source_pool_sha256':m['prepared_external_files'][str(S/'source_pool.jsonl')],'child_inputs_sha256':m['prepared_external_files'][str(S/'child_inputs.json')],'candidate_pairs_sha256':m['prepared_external_files'][str(S/'candidate_pairs.jsonl')],'core_sha256':CORE_SHA,'pairs':pairlist,'children':{str(i):children[i] for i in ci},'sources':{sid:sources[sid] for sid in ids}}
  path=W/f'shard{shard:02d}_input.json';assert not path.exists();write(path,obj);total_sources+=len(ids)
  plan.append({'shard_id':shard,'input_path':str(path),'input_sha256':sha(path),'pair_count':len(pairlist),'child_count':len(ci),'source_count':len(ids),'input_bytes':path.stat().st_size,'global_pair_indices_sha256':digest([p['pair_index'] for p in pairlist]),'child_indices_sha256':digest(ci)})
 assert sum(r['pair_count'] for r in plan)==89426 and sum(r['child_count'] for r in plan)==5427
 write(W/'shard_plan.json',{'worker_count':15,'shards':plan,'total_pairs':89426,'total_children':5427,'unique_needed_sources':10549,'all_source_copies_across_shards':total_sources,'whole_child_unique_owner':True,'all_global_pairs_exactly_once':True,'source131_prepared_manifest_sha256':OLD_SHA,'core_sha256':CORE_SHA,'production_authorized':False})
 result={'status':'frozen15shards_prepared_NOT_corner_executed','pid':os.getpid(),'pairs':89426,'children':5427,'workers':15,'pair_counts':[r['pair_count'] for r in plan],'child_counts':[r['child_count'] for r in plan],'source_counts':[r['source_count'] for r in plan],'total_shard_source_copies':total_sources,'unique_sources':10549,'total_input_bytes':sum(r['input_bytes'] for r in plan),'peak_memory_estimate_GiB':[0.8,2.5],'actual_new_corner_checks':0,'new_H_gradient_kernel_log_exp_calls':0,'wall_seconds':time.monotonic()-began,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'dispatch_authorized':False}
 write(W/'preparation_summary.json',result);print(json.dumps(result),flush=True)
if __name__=='__main__':main()
