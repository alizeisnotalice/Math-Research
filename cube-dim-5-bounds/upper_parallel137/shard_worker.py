"""One child-owned shard, durable hash-chain prefix, no automatic restart."""
from pathlib import Path
from fractions import Fraction as F
from datetime import datetime,timezone
import json,os,resource,traceback
from common import sha,digest,write,load_core,CORE_SHA
def recover(rows_path,checkpoint_path,pairs,fingerprint,resume):
 qualified={};cursor=evaluated=skipped=0;previous='0'*64;offset=0;cp=None;at_cp=None
 if checkpoint_path.exists():
  assert resume,'existing checkpoint requires explicit reviewed resume';cp=json.load(open(checkpoint_path));assert cp['fingerprint']==fingerprint
 if rows_path.exists():
  assert resume,'existing JSONL must never be overwritten'
  with open(rows_path,'rb') as f:
   for raw in f:
    assert raw.endswith(b'\n'),'incomplete tail retained:inspect/archive before reviewed recovery,do not blind truncate'
    row=json.loads(raw);h=row.pop('hash');assert row['seq']==cursor+1 and row['previous']==previous and digest(row)==h
    pair=pairs[cursor];assert row['global_pair_index']==pair['pair_index'] and row['child_index']==pair['child_index'] and row['source_id']==pair['source_id'];k=str(pair['child_index'])
    if row['kind']=='skipped':assert k in qualified and row['reason']=='child already certified by one full16corner source';skipped+=1
    else:
     assert row['kind']=='evaluated' and k not in qualified;cert=row['certificate'];assert cert['source_id']==pair['source_id'] and cert['child_node']==pair['child_node'] and cert['curvature_rule']==pair['curvature_rule'] and cert['checked_points']==16 and cert['passed']==(F(cert['minimum_log_margin'])>0);evaluated+=1
     if cert['passed']:qualified[k]=cert
    cursor+=1;previous=h;offset+=len(raw)
    if cp and cursor==cp['next_pair']:at_cp=(previous,offset,evaluated,skipped,dict(qualified))
 if cp:
  assert cp['next_pair']<=cursor
  expected=('0'*64,0,0,0,{}) if cp['next_pair']==0 else at_cp
  assert expected==(cp['last_event_hash'],cp['valid_prefix_bytes'],cp['evaluated'],cp['skipped'],cp['qualified_children'])
 return {'next_pair':cursor,'last_event_hash':previous,'valid_prefix_bytes':offset,'evaluated':evaluated,'skipped':skipped,'qualified_children':qualified,'fingerprint':fingerprint}
def execute(spec,out_dir,prepared_sha,resume,stop_event,start_event,ready_queue,result_queue):
 shard=spec['shard_id'];out=Path(out_dir)/f'shard{shard:02d}';out.mkdir(exist_ok=True);state=None;data=None;rows=None
 try:
  os.nice(5);assert sha(spec['input_path'])==spec['input_sha256'];data=json.load(open(spec['input_path']));assert data['shard_id']==shard and data['core_sha256']==CORE_SHA
  exact=load_core();pairs=data['pairs'];assert len(pairs)==spec['pair_count'] and len(data['children'])==spec['child_count'];assert all(p['child_index']%15==shard for p in pairs)
  fingerprint={'shard_id':shard,'input_sha256':spec['input_sha256'],'core_sha256':CORE_SHA,'worker_code_sha256':sha(__file__),'prepared_manifest_sha256':prepared_sha,'source_pool_sha256':data['source_pool_sha256'],'child_inputs_sha256':data['child_inputs_sha256'],'candidate_pairs_sha256':data['candidate_pairs_sha256'],'target':'11/100000'}
  original=out/'worker_registration.json'
  if original.exists():assert resume and json.load(open(original))['fingerprint']==fingerprint
  else:write(original,{'initial_pid':os.getpid(),'fingerprint':fingerprint,'pair_count':len(pairs),'child_count':len(data['children']),'source_count':len(data['sources'])})
  rows_path=out/'results.jsonl';checkpoint=out/'checkpoint.json';state=recover(rows_path,checkpoint,pairs,fingerprint,resume)
  attempt=1
  while (out/f'attempt{attempt:04d}_registration.json').exists():attempt+=1
  write(out/f'attempt{attempt:04d}_registration.json',{'pid':os.getpid(),'attempt':attempt,'registered_utc':datetime.now(timezone.utc).isoformat(),'resume_requested':resume,'start_pair':state['next_pair'],'source_count':len(data['sources']),'fingerprint':fingerprint,'nice':os.getpriority(os.PRIO_PROCESS,0)})
  ready_queue.put({'shard_id':shard,'pid':os.getpid(),'pair_count':len(pairs),'child_count':len(data['children']),'loaded_sources':len(data['sources']),'recovered_successful_prefix':state['next_pair'],'fingerprint':fingerprint,'nice':os.getpriority(os.PRIO_PROCESS,0)})
  start_event.wait();before_eval=state['evaluated'];before_skip=state['skipped']
  with open(rows_path,'ab') as rows:
   while state['next_pair']<len(pairs) and not stop_event.is_set():
    cursor=state['next_pair'];pair=pairs[cursor];k=str(pair['child_index']);row={'seq':cursor+1,'previous':state['last_event_hash'],'global_pair_index':pair['pair_index'],'child_index':pair['child_index'],'source_id':pair['source_id']}
    if k in state['qualified_children']:row.update({'kind':'skipped','reason':'child already certified by one full16corner source'});state['skipped']+=1
    else:
     source=data['sources'][pair['source_id']];child=data['children'][k];assert source['source_index']==pair['source_index'] and child['node']==pair['child_node'] and source['core']['curvature_rule']==pair['curvature_rule']
     cert=exact(source,child);row.update({'kind':'evaluated','certificate':cert});state['evaluated']+=1
     if cert['passed']:state['qualified_children'][k]=cert
    row['hash']=digest(row);encoded=(json.dumps(row,sort_keys=True,separators=(',',':'))+'\n').encode();rows.write(encoded);rows.flush();os.fsync(rows.fileno())
    state['next_pair']+=1;state['last_event_hash']=row['hash'];state['valid_prefix_bytes']+=len(encoded)
    if state['next_pair']%128==0:write(checkpoint,state)
  state['status']='complete' if state['next_pair']==len(pairs) else 'stopped_prefix';write(checkpoint,state)
  summary={'status':state['status'],'pid':os.getpid(),'shard_id':shard,'pair_count':len(pairs),'acknowledged_pairs':state['next_pair'],'evaluated_total':state['evaluated'],'skipped_total':state['skipped'],'new_evaluated':state['evaluated']-before_eval,'new_skipped':state['skipped']-before_skip,'qualified_child_count':len(state['qualified_children']),'child_count':len(data['children']),'loaded_sources':len(data['sources']),'checkpoint_sha256':sha(checkpoint),'results_sha256':sha(rows_path),'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'new_H_gradient_kernel_log_exp_calls':0};write(out/f'attempt{attempt:04d}_summary.json',summary);result_queue.put(summary)
 except BaseException as exc:
  stop_event.set();start_event.set();fault={'status':'worker_fault','pid':os.getpid(),'shard_id':shard,'exception':repr(exc),'traceback':traceback.format_exc(),'prefix_next_pair':state['next_pair'] if state else None,'action':'preserve JSONL/checkpoint;no automatic restart or truncation','new_H_gradient_kernel_log_exp_calls':0}
  write(out/f'fault_pid{os.getpid()}.json',fault)
  if state is not None:state['status']='fault_prefix_preserved';write(out/'fault_prefix_state.json',state)
  ready_queue.put({'shard_id':shard,'pid':os.getpid(),'status':'worker_fault'});result_queue.put(fault)
