"""Future two finite phases, <=15 live children; durable source/pair prefixes."""
from pathlib import Path
import json,os,resource,time,traceback
from fractions import Fraction as F
from util173 import sha,digest,read,atomic,append
from core173 import band_source_certificate,compile_curve,closedbox_certificate,initial_best,update_best
def execute(phase,spec,outroot,manifest_sha,runtime,ready,result,start,stop,resume):
 sid=spec['shard_id'];out=Path(outroot)/phase/f'shard{sid:02d}';out.mkdir(parents=True,exist_ok=True);state=None;current=None;began=time.monotonic()
 try:
  os.nice(5);inputpath=spec['path'] if phase=='mass' else spec['input_path'];inputsha=spec['sha256'];assert sha(inputpath)==inputsha;data=read(inputpath);assert sha(runtime['source_registry_path'])==runtime['source_registry_sha256'];sources=read(runtime['source_registry_path']);assert sha(runtime['old_common_log_path'])==runtime['old_common_log_sha256'];C=F(read(runtime['old_common_log_path'])['upper']);fp={'phase':phase,'shard_id':sid,'input_sha256':inputsha,'prepared_manifest_sha256':manifest_sha,'worker_sha256':sha(__file__),'source_registry_sha256':runtime['source_registry_sha256'],'common_log_sha256':runtime['old_common_log_sha256']};best={};ginput=qregistry=None
  if phase=='mass':
   gin=runtime['g_input_by_shard'][sid];assert sha(gin['g_input_path'])==gin['g_input_sha256'];ginput=read(gin['g_input_path']);fp['g_input_sha256']=gin['g_input_sha256'];jobs=data['source_ids']
  else:
   assert phase=='corners' and sha(runtime['q_registry_path'])==runtime['q_registry_sha256'];qregistry=read(runtime['q_registry_path']);fp['q_registry_sha256']=runtime['q_registry_sha256'];jobs=data['pairs'];best={k:initial_best(r) for k,r in data['receivers'].items()}
   for b in best.values():assert F(b['best_margin110'])-C==F(b['best_margin115'])
  reg=out/'registration.json'
  if reg.exists():assert resume and read(reg)['fingerprint']==fp
  else:atomic(reg,{'pid':os.getpid(),'fingerprint':fp,'input_count':len(jobs),'phase':phase})
  journal=out/'events.jsonl';n=0;previous='0'*64;counts={'mass':0,'evaluated':0,'skip_after_firstPASS':0,'skip_oldpair_q_unchanged':0};certificates={};cp=read(out/'checkpoint.json') if (out/'checkpoint.json').exists() else None;atcp=None
  if cp:assert resume and cp['fingerprint']==fp
  if journal.exists():
   assert resume
   for raw in open(journal,'rb'):
    assert raw.endswith(b'\n'),'truncated tail preserved for review';row=json.loads(raw);h=row.pop('hash');assert row['seq']==n+1 and row['previous']==previous and digest(row)==h
    if phase=='mass':assert row['canonical_id']==jobs[n] and row['kind']=='mass_certificate';assert row['certificate']['canonical_id']==jobs[n];certificates[jobs[n]]={'certificate':row['certificate'],'row_hash':h};counts['mass']+=1
    else:
     pair=jobs[n];assert all(row[k]==pair[k] for k in ('pair_index','receiver_index','canonical_id'));b=best[str(pair['receiver_index'])];q=qregistry[pair['canonical_id']]['certificate'];assert F(q['q_new'])<=F(q['oldqualified155q'])
     if row['kind']=='skip_after_firstPASS':assert b['first_pass115'] is not None;counts['skip_after_firstPASS']+=1
     elif row['kind']=='skip_oldpair_q_unchanged':assert b['first_pass115'] is None and pair['prior_test_same_old155_q'] is not None and F(q['q_new'])==F(q['oldqualified155q']);counts['skip_oldpair_q_unchanged']+=1
     else:
      assert row['kind']=='evaluated' and b['first_pass115'] is None;c=row['certificate'];assert c['canonical_id']==pair['canonical_id'] and c['node']==pair['node'] and c['source_mass_row_hash']==qregistry[pair['canonical_id']]['row_hash'] and F(c['margin110'])==F(c['margin115'])+C;assert c['passed115']==(F(c['margin115'])>0);update_best(b,c);counts['evaluated']+=1
    n+=1;previous=h
    if cp and n==cp['next_index']:atcp={'last_hash':previous,'counts':dict(counts),'best_sha256':digest(best)}
  if cp:
   empty={'last_hash':'0'*64,'counts':{'mass':0,'evaluated':0,'skip_after_firstPASS':0,'skip_oldpair_q_unchanged':0},'best_sha256':digest({k:initial_best(r) for k,r in data.get('receivers',{}).items()})};assert cp['next_index']<=n and (atcp if cp['next_index'] else empty)=={k:cp[k] for k in ('last_hash','counts','best_sha256')}
  state={'fingerprint':fp,'next_index':n,'last_hash':previous,'counts':counts,'best_sha256':digest(best)};ready.put({'phase':phase,'shard_id':sid,'pid':os.getpid(),'recovered_prefix':n,'input_count':len(jobs)});start.wait();compiled={};before=n
  with open(journal,'ab') as f:
   while state['next_index']<len(jobs) and not stop.is_set():
    n=state['next_index'];row={'seq':n+1,'previous':state['last_hash']}
    if phase=='mass':
     cid=jobs[n];current=cid;cert=band_source_certificate(sources[cid],ginput[cid],runtime['time_endpoints'],C);row.update({'kind':'mass_certificate','canonical_id':cid,'certificate':cert});counts['mass']+=1
    else:
     pair=jobs[n];current=pair['pair_index'];k=str(pair['receiver_index']);b=best[k];cid=pair['canonical_id'];q=qregistry[cid];assert 0<F(q['certificate']['q_new'])<=F(q['certificate']['oldqualified155q']);row.update({a:pair[a] for a in ('pair_index','receiver_index','canonical_id')})
     if b['first_pass115'] is not None:row['kind']='skip_after_firstPASS';counts['skip_after_firstPASS']+=1
     elif pair['prior_test_same_old155_q'] is not None and not q['certificate']['q_strictly_decreased']:assert F(q['certificate']['q_new'])==F(q['certificate']['oldqualified155q']);row['kind']='skip_oldpair_q_unchanged';counts['skip_oldpair_q_unchanged']+=1
     else:
      if cid not in compiled:compiled[cid]=compile_curve(sources[cid],q['certificate'],q['row_hash'],C)
      cert=closedbox_certificate(compiled[cid],data['receivers'][k]);update_best(b,cert);row.update({'kind':'evaluated','certificate':cert});counts['evaluated']+=1
    h=append(f,row)
    if phase=='mass':certificates[cid]={'certificate':cert,'row_hash':h}
    state.update({'next_index':n+1,'last_hash':h})
    if phase=='mass' or (n+1)%128==0:state['best_sha256']=digest(best);atomic(out/'checkpoint.json',state)
  state['best_sha256']=digest(best);state['status']='complete' if state['next_index']==len(jobs) else 'stopped_successful_prefix';atomic(out/'checkpoint.json',state)
  if phase=='mass':atomic(out/'source_certificates.json',certificates)
  else:atomic(out/'receiver_best.json',best)
  u=resource.getrusage(resource.RUSAGE_SELF);summary={**state,'shard_id':sid,'phase':phase,'pid':os.getpid(),'new_cases_this_attempt':state['next_index']-before,'wall_seconds':time.monotonic()-began,'cpu_user_seconds':u.ru_utime,'cpu_system_seconds':u.ru_stime,'peak_rss_bytes':u.ru_maxrss,'new_H_log_gradient_mean_calls':0};atomic(out/'summary.json',summary);result.put(summary)
 except BaseException as exc:
  stop.set();start.set();fault={'status':'fault_successful_prefix_preserved','phase':phase,'shard_id':sid,'pid':os.getpid(),'current_object':current,'exception':repr(exc),'details':getattr(exc,'details',None),'traceback':traceback.format_exc()};atomic(out/f'fault{os.getpid()}.json',fault)
  if state is not None:atomic(out/'fault_prefix.json',state)
  ready.put({'phase':phase,'shard_id':sid,'pid':os.getpid(),'status':'fault'});result.put(fault)
