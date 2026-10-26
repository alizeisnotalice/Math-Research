"""Hash-lock prepared code/inputs without executing any certificate core."""
from pathlib import Path
import ast,json,os
from common import sha,write,CORE_SHA
H=Path(__file__).resolve().parent;W=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_parallel137');OLD=H.parent/'upper_crossanchor131'
def main():
 dest=H/'prepared_manifest.json';assert not dest.exists() and not (W/'run001').exists()
 old=json.load(open(OLD/'prepared_manifest.json'));assert sha(OLD/'prepared_manifest.json')=='1520dc778a8d92631ed33179690605a9eefdca6ae3e7c209dd6c8ffbcc4638f3'
 locks=dict(old['source_locks']);locks[str(OLD/'prepared_manifest.json')]=sha(OLD/'prepared_manifest.json')
 for name,h in old['repository_files'].items():locks[str(OLD/name)]=h
 for p,h in old['prepared_external_files'].items():locks[p]=h
 for p,h in locks.items():assert sha(p)==h
 repository={p.name:sha(p) for p in sorted(H.iterdir()) if p.is_file() and p.name!='prepared_manifest.json'}
 for p in H.glob('*.py'):ast.parse(p.read_text())
 plan=json.load(open(W/'shard_plan.json'));receipt=json.load(open(W/'prepared_check_receipt.json'));assert receipt['status']=='PASS_input_equality_and_synthetic_recovery_only' and receipt['actual_source_corner_calls']==0
 external={str(p):sha(p) for p in sorted(W.iterdir()) if p.is_file()};assert len(plan['shards'])==15 and sum(x['pair_count'] for x in plan['shards'])==89426
 for spec in plan['shards']:assert external[spec['input_path']]==spec['input_sha256']
 write(dest,{'schema':'upper_parallel137_prepared_v1','status':'prepared_only_waiting_primary_review_before_corner_dispatch','freeze_pid':os.getpid(),'target':'11/100000','worker_count':15,'fixed_pairs':89426,'fixed_children':5427,'unique_needed_sources':10549,'source_copies_across_shards':64568,'source131_prepared_manifest_sha256':sha(OLD/'prepared_manifest.json'),'exact_core_sha256':CORE_SHA,'repository_files':repository,'source_locks':locks,'prepared_external_files':external,'shard_plan_sha256':sha(W/'shard_plan.json'),'check_receipt_sha256':sha(W/'prepared_check_receipt.json'),'actual_source_corner_checks':0,'new_H_gradient_kernel_log_exp_calls':0,'new_splits':0,'wall_deadline':None,'memory_estimate_GiB':[0.8,2.5],'memory_estimate_not_cap':True,'production_authorized':False,'scope':'fixed existing jets single-source whole closedbox candidates; no globalg110 or S credit'})
 print(json.dumps({'manifest':str(dest),'sha256':sha(dest),'freeze_pid':os.getpid(),'repository_files':len(repository),'source_locks':len(locks),'prepared_external_files':len(external),'run_exists':(W/'run001').exists()}),flush=True)
if __name__=='__main__':main()
