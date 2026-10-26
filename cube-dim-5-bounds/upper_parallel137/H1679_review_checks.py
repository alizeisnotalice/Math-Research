"""Verify new source codec, closed tree and fresh-point guards without H/corners."""
from pathlib import Path
from fractions import Fraction as F
import ast,json,os,time
from common import sha,write
import H1679_run as driver
H=Path(__file__).resolve().parent;S=driver.S;W=driver.W
def main():
 began=time.monotonic();dest=W/'source_codec_check.json';assert not dest.exists();package=json.load(open(H/'residual_prepared_manifest.json'));old=driver.legacy_load(package);seed=driver.source_seed(package);old.worker.RUN_TARGET=F(11,100000);tasks=json.load(open(W/'tasks.json'));assert len(tasks)==1679 and {tuple(t['node']) for t in tasks}=={tuple(n) for n in seed['queue']};keys=set()
 for t in tasks:
  assert old.worker.validate_task(t)==F(11,100000) and t['same_point_H64_guard'] is None and not t['cached'];p=t['pair_cache_key'][1];assert p not in seed['cache_by_epsilon']['1/128'] and p not in keys;keys.add(p);lo,hi=old.physical(t['node'],seed['roots']);assert list(map(str,lo))==t['box_lower'] and list(map(str,hi))==t['box_upper'];point=list(map(F,t['point']));assert all(a<x<b for a,x,b in zip(lo,point[1:],hi));ancestor=t['H64_g105_ancestor']['node'];alo,ahi=old.physical(ancestor,seed['roots']);assert ancestor[0]==t['node'][0] and t['node'][1].startswith(ancestor[1]) and all(a<=x<=b for a,x,b in zip(alo,point[1:],ahi));assert F(t['H64_g105_ancestor']['lower'])==F(21,200000)
 tree=old.tree(seed,{'new_passed':[],'new_splits':[],'queue':[t['node'] for t in tasks],'inflight':{}});assert tree['passed']==109405 and tree['queue_and_inflight']==1679 and tree['splits']==58724
 for p in [H/'H1679_run.py',H/'prepare_H1679.py',H/'H1679_review_checks.py']:ast.parse(p.read_text())
 result={'status':'PASS_true1679_sourcecodec_freshpoints_fulltree_and_taskguards_ONLY','pid':os.getpid(),'tree':tree,'cache64':len(seed['cache_by_epsilon']['1/64']),'cache128':len(seed['cache_by_epsilon']['1/128']),'target':'11/100000','epsilon':'1/128','fresh_unique_anchors':len(keys),'samepoint_H64_guards':0,'source_passed_retained':109405,'source_splits_retained':58724,'actual_H_gradient_kernel_log_exp_corner_calls':0,'dispatch_authorized':False,'wall_seconds':time.monotonic()-began};write(dest,result);print(json.dumps(result),flush=True)
if __name__=='__main__':main()
