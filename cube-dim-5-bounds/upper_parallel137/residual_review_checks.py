"""New dispatcher source-codec and task guards only; no numerical kernel call."""
from pathlib import Path
from fractions import Fraction as F
import ast,json,os,time
from common import sha,digest,write
import residual_run
H=Path(__file__).resolve().parent;S=residual_run.S;W=residual_run.W;OLD=residual_run.LEGACY
def main():
 began=time.monotonic();assert not (W/'source_codec_check.json').exists();prepared=json.load(open(OLD/'production_prepared_manifest.json'));names=['run.py','binding.py','worker.py','cache_model.py','certificate.py','reuse.py'];locks={str(OLD/n):prepared['repository_files'][n] for n in names}
 assert locks[str(OLD/'worker.py')]=='62375ec07bf0a030c397c4a925b0809ab1f625edf75a19f434e4eb7f8a492e8b'
 old=residual_run.legacy_load({'legacy_code_locks':locks});seed=residual_run.source_seed({});tasks=json.load(open(W/'tasks.json'));residual=list(map(json.loads,open(S/'upper_parallel137/remaining_gap_order.jsonl')));assert len(tasks)==len(residual)==4274 and [t['node'] for t in tasks]==[r['node'] for r in residual]
 old.worker.RUN_TARGET=residual_run.TARGET
 for t,r in zip(tasks,residual):
  assert old.worker.validate_task(t)==F(11,100000);lo,hi=old.physical(t['node'],seed['roots']);assert list(map(str,lo))==t['box_lower']==r['box_lower'] and list(map(str,hi))==t['box_upper']==r['box_upper'];point=list(map(F,t['point']));assert all(a<x<b for a,x,b in zip(lo,point[1:],hi))
  ancestor=t['H64_g105_ancestor']['node'];assert ancestor==r['prior_direct_parent_newq_failure']['H64_lower_source_node'];alo,ahi=old.physical(ancestor,seed['roots']);assert all(a<=x<=b for a,x,b in zip(alo,point[1:],ahi)) and F(t['H64_g105_ancestor']['lower'])==F(r['prior_direct_parent_newq_failure']['H64_lower'])==F(21,200000)
  assert t['same_point_H64_guard'] is None and t['pair_cache_key'][1] not in seed['cache_by_epsilon']['1/128'] and not t['cached'];assert t['geometry_tests'][-1]['eligible'] and t['geometry_tests'][-1]['distinct_count']==280 and t['geometry_tests'][-1]['zero_count']==0
 ov={'new_passed':[],'new_splits':[],'queue':[t['node'] for t in tasks],'inflight':{}};tree=old.tree(seed,ov);assert tree['passed']==105635 and tree['queue_and_inflight']==4274 and tree['splits']==57549
 for p in H.glob('residual_*.py'):ast.parse(p.read_text(),filename=str(p))
 result={'status':'PASS_new_sourcecodec_fulltree_exactbox_ancestry_taskguard_only','pid':os.getpid(),'legacy_code_locks':locks,'source_seed123_sha256':sha(S/'upper_continue123/seed.json'),'source123_final_state_sha256':sha(S/'upper_continue123/batch0002/state.json'),'typedcache64':len(seed['cache_by_epsilon']['1/64']),'typedcache128':len(seed['cache_by_epsilon']['1/128']),'tree':tree,'all4274_fresh_generic_inside_child_and_H64g105_ancestor':True,'all_samepoint_H64_guards_null':True,'actual_H_gradient_kernel_evaluate_log_exp_calls':0,'actual_corner_certificate_calls':0,'production_authorized':False,'wall_seconds':time.monotonic()-began};write(W/'source_codec_check.json',result);print(json.dumps(result),flush=True)
if __name__=='__main__':main()
