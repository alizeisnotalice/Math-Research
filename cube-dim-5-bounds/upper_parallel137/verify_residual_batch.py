"""Replay new production/source/point/typedcache/tree gates without H/log/corners."""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import argparse,json,os,time
from common import sha,digest,write
import residual_run
H=Path(__file__).resolve().parent;S=residual_run.S;W=residual_run.W
def main():
 began=time.monotonic();ap=argparse.ArgumentParser();ap.add_argument('--batch',type=int,choices=(1,2),required=True);a=ap.parse_args();bdir=W/f'future_batch{a.batch:04d}';dest=bdir/'author_validation.json';assert not dest.exists();package=json.load(open(H/'residual_prepared_manifest.json'));assert sha(H/'residual_prepared_manifest.json')=='99e932b0b86c9d89efbebd6a187cd26f977679254b4b244fba82b7b63cd56260'
 for n,h in package['repository_files'].items():assert sha(H/n)==h
 for p,h in {**package['source_locks'],**package['external_files'],**package['legacy_code_locks']}.items():assert sha(p)==h
 old=residual_run.legacy_load(package);seed=residual_run.source_seed(package);summary=json.load(open(bdir/'summary.json'));reg=json.load(open(bdir/'registration.json'));manifest=json.load(open(bdir/'manifest.json'));assert summary['status']=='finite_original_prefix_complete' and summary['faults']==summary['obstructions']==summary['inflight']==0
 for n,h in manifest['files'].items():assert sha(bdir/n)==h
 write(bdir/'author_registration.json',{'pid':os.getpid(),'code_sha256':sha(__file__),'batch':a.batch,'point_manifest_sha256':sha(bdir/'point_manifest.json'),'state_sha256':sha(bdir/'state.json'),'journal_sha256':sha(bdir/'journal.jsonl'),'new_H_gradient_kernel_log_exp_calls':0,'repeat_corner_certificates':0})
 tasks=json.load(open(W/f'input_batch{a.batch:04d}.json'));ov=json.load(open(bdir/'initial_overlay.json'));previous='0'*64;seq=0;received={};merged={};counts=Counter()
 assert reg['source_cache_canonical_sha256']==digest(seed['cache_by_epsilon']) and reg['source_passed_canonical_sha256']==digest(seed['passed'])
 for event in map(json.loads,open(bdir/'journal.jsonl')):
  seq+=1;h=event.pop('hash');assert event['seq']==seq and event['previous']==previous and digest(event)==h;event['hash']=h;previous=h;counts[event['kind']]+=1
  body=event['body']
  if event['kind']=='dispatch':assert body['task']==tasks[body['task']['batch_ordinal']]
  elif event['kind']=='received':assert body['ordinal'] not in received;received[body['ordinal']]=body['returned']
  elif event['kind']=='merge':assert body['ordinal'] not in merged;merged[body['ordinal']]=body['certificate']
  old.apply(ov,event,seed['cache_by_epsilon'])
 state=json.load(open(bdir/'state.json'));assert ov==state and len(received)==len(merged)==len(tasks) and counts=={'dispatch':len(tasks),'received':len(tasks),'merge':len(tasks),'finish':1};assert all(v['kind']=='complete' for v in received.values())
 points=json.load(open(bdir/'point_manifest.json'))['records'];assert len(points)==len(tasks);bits=Counter();segments=Counter();minhl=None;maxwidth=F(0);fresh={}
 for row in points:
  i=row['ordinal'];task=tasks[i];p=Path(row['point_checkpoint']);assert sha(p)==row['point_checkpoint_sha256'];cp=json.load(open(p));got=received[i];record=got['record'];c=merged[i]
  assert cp['record']==row['record']==record and cp['partial'] is None and cp['stage']=='complete' and cp['calls']==got['calls'] and cp['attempts']==got['attempts'];assert cp['point']==task['point'] and cp['epsilon']=='1/128' and cp['target']=='11/100000' and cp['kernel_sha256']==package['kernel_sha256']
  old.validate_record(record,F(11,100000));assert record['pair_cache_key']==task['pair_cache_key'] and record['point']==task['point'];key=record['pair_cache_key'][1];assert key not in seed['cache_by_epsilon']['1/128'] and key not in fresh;fresh[key]=record
  hl,hu=map(F,record['values'][0]);assert F(11,100000)<hl<=hu;minhl=hl if minhl is None else min(minhl,hl);maxwidth=max(maxwidth,hu-hl);bits[record['bits']]+=1;segments[cp['calls']]+=1
  assert c['target']=='11/100000' and c['leaf_epsilon']=='1/128' and c['point_key']==key and c['checked_points']==16 and c['H64_g105_ancestor_proof']==task['H64_g105_ancestor'];assert F(c['H128_upper'])==hu and F(c['H64_lower'])==F(21,200000) and F(c['penalty_coefficient'])==1024-768*F(21,200000)/hu
  assert c['passed']==(F(c['minimum_log_margin'])>0)
  for pair,(lo,hi) in zip(record['values'][1:],c['loggradient_directed']):
   u,v=map(F,pair);assert u<=v and F(lo)<=min(u/hl,u/hu,v/hl,v/hu)<=max(u/hl,u/hu,v/hl,v/hu)<=F(hi)
 assert all(state['new_cache_by_epsilon']['1/128'][k]==r for k,r in fresh.items());tree=old.tree(seed,state);assert tree==summary['tree'];assert digest(seed['cache_by_epsilon'])==reg['source_cache_canonical_sha256'] and digest(seed['passed'])==reg['source_passed_canonical_sha256']
 pids=[reg['master_pid']]+[p['pid'] for p in summary['pool_evidence']];alive=[]
 for pid in pids:
  try:os.kill(pid,0);alive.append(pid)
  except ProcessLookupError:pass
 assert not alive and len(summary['pool_evidence'])==15 and all(p['exitcode']==0 and not p['alive'] for p in summary['pool_evidence'])
 result={'status':'PASS_all_events_points_sourcecache_qbinding_fulltree_and_exit','pid':os.getpid(),'batch':a.batch,'production_master_pid':reg['master_pid'],'production_worker_pids':pids[1:],'all_production_PIDs_OS_absent':True,'new_H_count':len(fresh),'new_physical_gradient_count':4*len(fresh),'events':seq,'event_counts':dict(counts),'bits':dict(bits),'segments_per_point':dict(segments),'min_new_H_lower':str(minhl),'max_new_H_interval_width':str(maxwidth),'new_passed':summary['new_passed'],'new_splits':summary['new_splits'],'tree':tree,'faults':0,'obstructions':0,'inflight':0,'sourcecache_unchanged':True,'state_sha256':sha(bdir/'state.json'),'journal_sha256':sha(bdir/'journal.jsonl'),'point_manifest_sha256':sha(bdir/'point_manifest.json'),'new_H_gradient_kernel_log_exp_calls':0,'repeat_corner_certificates':0,'production_wall_seconds':summary['wall_seconds'],'wall_seconds':time.monotonic()-began,'globalg110_qualified':False,'S_credit':0};write(dest,result);print(json.dumps(result),flush=True)
if __name__=='__main__':main()
