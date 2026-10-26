"""Freeze exact finite prepared inputs and code only; no newmath call."""
from pathlib import Path
import ast,json,hashlib,os
H=Path(__file__).resolve().parent;W=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_retarget155')
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  while b:=f.read(8*1024*1024):h.update(b)
 return h.hexdigest()
def main():
 dest=H/'prepared_manifest.json';assert not dest.exists() and not (W/'run001').exists();summary=json.load(open(W/'summary.json'));assert summary['objects']==107530 and summary['leaves']==111201 and summary['actual_new_q_log_H_gradient_corner_calls']==0;plan=json.load(open(W/'finite_plan.json'));seen_obj=set();seen_leaf=set()
 for spec in plan['shards']:
  assert sha(spec['input_path'])==spec['input_sha256'];data=json.load(open(spec['input_path']));ids=set()
  for s in data['sources']:
   oid=s['canonical_id'];assert oid not in seen_obj and s['canonical_index']%15==spec['shard_id'] and s['new_q_status']=='NOT_COMPUTED';seen_obj.add(oid);ids.add(oid)
  for r in data['leaves']:assert r['leaf_index'] not in seen_leaf and r['canonical_id'] in ids;seen_leaf.add(r['leaf_index'])
 assert len(seen_obj)==107530 and seen_leaf==set(range(111201));locks=json.load(open(W/'source_locks.json'))
 for p,h in locks.items():assert sha(p)==h
 names=['prepare155.py','util155.py','core155.py','worker155.py','run155.py','mechanism155.md','freeze155.py'];files={n:sha(H/n) for n in names}
 for n in names:
  if n.endswith('.py'):ast.parse((H/n).read_text(),filename=n)
 external={str(p):sha(p) for p in W.iterdir() if p.is_file()};penalty='/Users/zhengzhihao/Desktop/T/Math-Research-cube5-curvature128-sol/cube-dim-5-bounds/upper_curvature153/penalty.py';logsource=H.parent/'upper_retarget149/shared_log.py';payload={'schema':'upper_retarget155_joint153_g115_prepared_v1','status':'PREPARED_ONLY_wait_primary_dispatch_review','freeze_pid':os.getpid(),'repository_files':files,'source_locks':locks,'prepared_external_files':external,'penalty_source':{'path':penalty,'sha256':sha(penalty)},'shared_log_source':{'path':str(logsource),'sha256':sha(logsource),'ratio':'23/22','planned_new_constants':1,'actual_evaluated':0},'source_objects':107530,'fixed_old_closed_leaves':111201,'worker_count':15,'ownership':'canonical_index modulo15;allaliasleaves remainwith one qobjectowner','oldq129_qualified_objects':8135,'old_H64g105_support_cells':3812,'newq_mass_policy':'only original H64g105;no newrawL and no globalH128 gstar as H64','new_q_rule':'153;allowbelow256;positiveand<=oldq only','single16corner_outputs':['margin110_new','margin115=margin110_new-Ulog'],'required_full110_invariant':'margin110_new>=saved_oldmargin>0 foreach111201leaf','full_domain_candidate':'g0*(1+min margin110_new);old completeforest retained eveniftarget115 pending','actual_new_q_log_H_gradient_corner_calls':0,'oldq_corner_baseline_runs':0,'new_exp_calls':0,'H_or_split_or_newcutoff_authorized':False,'execution_authorized':False,'formal_g_or_C5_updated':False,'wall_deadline':None,'memory_estimate_GiB':[1.5,4],'memory_estimate_not_cap':True}
 with open(dest,'x') as f:json.dump(payload,f,sort_keys=True,separators=(',',':'));f.write('\n')
 print(json.dumps({'manifest':str(dest),'sha256':sha(dest),'plan_sha256':sha(W/'finite_plan.json'),'source_locks':len(locks),'external_files':len(external),'q_log_corner_calls':0}),flush=True)
if __name__=='__main__':main()
