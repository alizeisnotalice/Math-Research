"""One canonical-owned finite shard. Durable q/leaf prefixes; no automatic retry."""
from pathlib import Path
from fractions import Fraction as F
from datetime import datetime,timezone
import json,os,traceback,resource
from util155 import sha,digest,atomic
from core155 import load_penalty_after_authorization,new_parameters_after_authorization,one_new16corner_calculation
def append(f,row,previous,seq):
 event={'seq':seq+1,'previous':previous,'body':row};event['hash']=digest(event);b=(json.dumps(event,sort_keys=True,separators=(',',':'))+'\n').encode();f.write(b);f.flush();os.fsync(f.fileno());return event['hash'],seq+1
def scan(path,items,kind,resume):
 rows=[];previous='0'*64
 if path.exists():
  assert resume,'existing prefix needs explicit fingerprinted resume'
  for b in open(path,'rb'):
   assert b.endswith(b'\n'),'truncated tail retained,inspect before recovery';event=json.loads(b);h=event.pop('hash');assert event['seq']==len(rows)+1 and event['previous']==previous and digest(event)==h;row=event['body'];item=items[len(rows)]
   assert row['canonical_id']==item['canonical_id']
   if kind=='leaf':assert row['leaf_index']==item['leaf_index'] and row['node']==item['node'] and row['margin110_new_ge_old_positive'] and F(row['margin110_new'])>=F(item['saved_g110_margin'])>0
   rows.append(row);previous=h
 return rows,previous
def execute(spec,outroot,reviewed_sha,resume,start,stop,ready,result):
 out=Path(outroot)/f'shard{spec["shard_id"]:02d}';out.mkdir(exist_ok=True);state=None
 try:
  os.nice(5);assert sha(spec['input_path'])==spec['input_sha256'];data=json.load(open(spec['input_path']));fingerprint={'shard_id':spec['shard_id'],'input_sha256':spec['input_sha256'],'reviewed_manifest_sha256':reviewed_sha,'worker_sha256':sha(__file__)};reg=out/'registration.json'
  if reg.exists():assert resume and json.load(open(reg))['fingerprint']==fingerprint
  else:atomic(reg,{'pid':os.getpid(),'registered_utc':datetime.now(timezone.utc).isoformat(),'fingerprint':fingerprint,'objects':len(data['sources']),'leaves':len(data['leaves']),'nice':os.getpriority(os.PRIO_PROCESS,0)},True)
  ready.put({'shard_id':spec['shard_id'],'pid':os.getpid(),'fingerprint':fingerprint,'objects':len(data['sources']),'leaves':len(data['leaves'])});start.wait()
  if stop.is_set():result.put({'status':'stopped_before_newmath','shard_id':spec['shard_id'],'pid':os.getpid()});return
  mathstart=json.load(open(Path(outroot)/'math_start.json'));assert mathstart['reviewed_manifest_sha256']==reviewed_sha and sha(mathstart['common_log_path'])==mathstart['common_log_sha256'];common=json.load(open(mathstart['common_log_path']));assert common['ratio']=='23/22';Ulog=F(common['upper']);assert Ulog>0
  api=load_penalty_after_authorization();qrows,qhash=scan(out/'q_updates.jsonl',data['sources'],'q',resume);leafrows,lhash=scan(out/'leaf_results.jsonl',data['leaves'],'leaf',resume);parameters={r['canonical_id']:r for r in qrows};qbegin=len(qrows);lbegin=len(leafrows);minimum=None;passed115=0
  assert not leafrows or len(qrows)==len(data['sources'])
  if (out/'checkpoint.json').exists():
   assert resume;cp=json.load(open(out/'checkpoint.json'));assert cp['fingerprint']==fingerprint and cp['common_log_sha256']==mathstart['common_log_sha256'] and cp['next_q']<=len(qrows) and cp['next_leaf']<=len(leafrows)
   for filename,count,key in [('q_updates.jsonl',cp['next_q'],'q_hash'),('leaf_results.jsonl',cp['next_leaf'],'leaf_hash')]:
    expected='0'*64
    for i,b in enumerate(open(out/filename,'rb')):
     if i>=count:break
     expected=json.loads(b)['hash']
    assert cp[key]==expected
  for r in leafrows:
   assert F(r['margin115'])==F(r['margin110_new'])-Ulog;key=(F(r['margin110_new']),r['leaf_index']);minimum=key if minimum is None else min(minimum,key);passed115+=r['passed115']
  state={'fingerprint':fingerprint,'common_log_sha256':mathstart['common_log_sha256'],'next_q':len(qrows),'next_leaf':len(leafrows),'q_hash':qhash,'leaf_hash':lhash,'passed115':passed115,'minimum110':None if minimum is None else [str(minimum[0]),minimum[1]]}
  with open(out/'q_updates.jsonl','ab') as qf:
   for source in data['sources'][len(qrows):]:
    if stop.is_set():break
    row=new_parameters_after_authorization(source,Ulog,api);qhash,n=append(qf,row,qhash,state['next_q']);parameters[row['canonical_id']]=row;state.update({'next_q':n,'q_hash':qhash})
    if n%128==0:atomic(out/'checkpoint.json',state)
  if state['next_q']==len(data['sources']):
   with open(out/'leaf_results.jsonl','ab') as lf:
    for leaf in data['leaves'][len(leafrows):]:
     if stop.is_set():break
     row=one_new16corner_calculation(parameters[leaf['canonical_id']],leaf,Ulog);lhash,n=append(lf,row,lhash,state['next_leaf']);passed115+=row['passed115'];key=(F(row['margin110_new']),row['leaf_index']);minimum=key if minimum is None else min(minimum,key);state.update({'next_leaf':n,'leaf_hash':lhash,'passed115':passed115,'minimum110':[str(minimum[0]),minimum[1]]})
     if n%128==0:atomic(out/'checkpoint.json',state)
  state['status']='complete' if state['next_q']==len(data['sources']) and state['next_leaf']==len(data['leaves']) else 'stopped_successful_prefix';atomic(out/'checkpoint.json',state);summary={**state,'shard_id':spec['shard_id'],'pid':os.getpid(),'new_q_objects':state['next_q']-qbegin,'new_leaf16corner_cases':state['next_leaf']-lbegin,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'new_H_gradient_or_split_calls':0};atomic(out/'summary.json',summary);result.put(summary)
 except BaseException as exc:
  stop.set();start.set();fault={'status':'new155_mechanism_fault_prefix_preserved','shard_id':spec['shard_id'],'pid':os.getpid(),'exception':repr(exc),'traceback':traceback.format_exc(),'action':'no blind restart;old complete110 forest unchanged'};atomic(out/f'fault_pid{os.getpid()}.json',fault)
  if state is not None:atomic(out/'fault_prefix.json',state)
  ready.put({'shard_id':spec['shard_id'],'pid':os.getpid(),'status':'fault'});result.put(fault)
