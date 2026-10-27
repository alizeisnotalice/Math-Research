"""Saved-only preparation: exact terminal metadata/unknown/budget checks and new-certificate jobs; no Arb/native."""
import sys,os,json,hashlib,collections
from pathlib import Path
sys.set_int_max_str_digits(0);sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent;A=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-numerics/cube-dim-6-bounds/upper/round129/recovered_continuation');REC=A.parents[2]/'audit/round127/crash_recovery'
rd=lambda p:json.loads(Path(p).read_text());sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(n,x):(P/n).write_text(json.dumps(x,indent=2)+'\n')
if(P/'prepared_new_jobs.json').exists():raise SystemExit('preparation already frozen; inspect before any revision')
run=rd(A/'run_receipt_v129.json');final=rd(A/'final_runtime_receipt_v129.json');exit=rd(A/'detached_exit_receipt.json');idx=rd(A/'semantic_runtime_index_v129.json');certidx=rd(A/'certificate_index_v129.json');oldcache=rd(REC/'combined_cache.json');unknown=rd(REC/'unknown_no_reissue.json');budgets=rd(REC/'root_recovery_and_budget.json')
assert exit['exit_code']==0 and run['failure']is None and final['failure']is None and final['live_owned_processes']=={} and run['status']==final['status']=='FINITE_FRONTIER_RUN_COMPLETE'
assert sha(A/'semantic_runtime_index_v129.json')==run['semantic_runtime_index_sha256']==final['semantic_runtime_index_sha256']and sha(A/'reservations_v129.jsonl')==run['reservation_journal_sha256']==final['reservation_journal_sha256']
unknownkeys={u['semantic_key']for u in unknown};states={};unknown_states={};leaves=[]
for key,entry in idx['entries'].items():
 assert sha(entry['state_path'])==entry['state_sha256'];s=rd(entry['state_path']);assert s['semantic_key']==key and s['status']==entry['status']
 if s['status']=='UNKNOWN_NO_REISSUE':
  assert key in unknownkeys and s['levels']==[]and entry['level_count']==0;unknown_states[key]=s
 else:
  assert s['status']=='AVAILABLE'and s['artifact_key']=='C129.'+key and key not in oldcache['entries']and key not in unknownkeys
  assert len(s['levels'])==entry['level_count']==1 and s['levels'][0]['status']=='ACCEPTED';states[key]=s
assert len(states)==18257 and len(unknown_states)==15 and set(unknown_states)==unknownkeys
for u in unknown:
 s=unknown_states[u['semantic_key']];assert s['expected']==u['expected']and s['owner_root']==u['owner_root']and s['owner_node_id']==u['owner_node_id']and s['prior_C127_level']==u['level']
journal=[json.loads(line)for line in(A/'reservations_v129.jsonl').read_text().splitlines()];byroot=collections.Counter();seen=set()
for ordinal,row in enumerate(journal,3128):
 k=row['semantic_key'];assert row['global_ordinal']==ordinal and k in states and k not in seen;seen.add(k);byroot[row['owner_root']]+=1
 assert row['root_ordinal']==budgets[row['owner_root']]['native_used_127']+byroot[row['owner_root']]and row['precision_bits']==states[k]['accepted_precision_bits']
assert len(journal)==18257 and seen==set(states)and run['total_native_attempt_reservations']==21384==3127+len(journal)
jobs=[];new_cert_sources=[];old_cert_sources=[];jobdir=P/'jobs';jobdir.mkdir(exist_ok=True)
assert not list(jobdir.iterdir()),'job directory must be empty before first frozen preparation'
def new_source_job(k):
 s=states[k];level=s['levels'][0];assert sha(s['native_path'])==s['native_sha256']and sha(s['receipt_path'])==level['receipt_sha256'];receipt=rd(s['receipt_path'])
 assert receipt['exit_code']==0 and receipt['raw_sha256']==s['native_sha256']and receipt['raw_path_after_move']==s['native_path']
 intent=Path(s['receipt_path']).parent/'intent.json';process=Path(s['receipt_path']).parent/'original_process.json';assert sha(intent)==receipt['intent_sha256']and sha(process)==receipt['process_sha256']
 return {'source_kind':'NEW129','semantic_key':k,'state_ref':{'path':idx['entries'][k]['state_path'],'sha256':idx['entries'][k]['state_sha256']},'jet_sha256':sha(s['jet_path']),'receipt_sha256':sha(s['receipt_path']),'intent_path':str(intent),'intent_sha256':sha(intent),'process_path':str(process),'process_sha256':sha(process)}
for e in sorted(certidx['entries'],key=lambda e:e['key']):
 assert sha(e['path'])==e['sha256']
 if e['native_artifact_key'].startswith('C129.'):
  k=e['native_artifact_key'][5:];s=states[k];assert e['expected']['center']==s['expected']['center']and e['expected']['window']==s['expected']['window'];job=new_source_job(k);new_cert_sources.append(k)
 else:
  k=oldcache['artifact_to_semantic'][e['native_artifact_key']];m=oldcache['entries'][k];assert m['artifact_key']==e['native_artifact_key']and e['expected']['center']==m['expected']['center']and e['expected']['window']==m['expected']['window']
  if 'jet_path'in m:jr={'path':m['jet_path'],'sha256':m['jet_sha256'],'selector':[]}
  else:jr={**m['jet_ref'],'sha256':sha(m['jet_ref']['path'])}
  assert sha(jr['path'])==jr['sha256'];job={'source_kind':'FROZEN_BASELINE','semantic_key':k,'frozen_source_metadata':m,'jet_ref':jr};old_cert_sources.append(k)
 job['certificate']=e;jobs.append(job)
assert len(new_cert_sources)==len(set(new_cert_sources))==18253 and len(old_cert_sources)==4 and len(jobs)==certidx['new_certificates']==18257
native_only=sorted(set(states)-set(new_cert_sources));assert len(native_only)==4
for k in native_only:job=new_source_job(k);job['certificate']=None;jobs.append(job)
shards=[]
for number in range(15):
 rows=jobs[number::15];q=jobdir/f'new_certificates_{number:02d}.jsonl';q.write_text(''.join(json.dumps(x,separators=(',',':'))+'\n'for x in rows));nc=sum(row['certificate']is not None for row in rows);shards.append({'id':number,'path':str(q),'sha256':sha(q),'jobs':len(rows),'certificates':nc,'corner_exp_max':32*nc})
sourcefiles=[A/n for n in ['registration_129.json','run_prepared_129.py','parent_release_v129.json','release_requirements.json','prepared_packet.json','continuation_seed129.json','unknown_geometry_receipt.json','remaining_budgets129.json','forest_v129.json','certificate_index_v129.json','semantic_runtime_index_v129.json','reservations_v129.jsonl','run_receipt_v129.json','final_runtime_receipt_v129.json','detached_exit_receipt.json']]+[REC/n for n in ['recovery_receipt.json','combined_cache.json','unknown_no_reissue.json','root_recovery_and_budget.json','recovered_forest.json','recovered_certificate_index.json']]
save('terminal_source_pins.json',{'files':[{'path':str(q),'sha256':sha(q),'bytes':q.stat().st_size}for q in sourcefiles]})
out={'status':'PREPARED_NEW_CERTIFICATE_JOBS_NO_ARB_PENDING_MAIN_LAUNCH','new_certificates':18257,'new_native_sources_with_certificates':18253,'frozen_old_jets_used_for_new_geometry_only':4,'new_native_only_jobs':4,'total_jobs':18261,'new_Arb768_corner_exp_max':584224,'shards':shards,'workers_max':15,'old15193_revalidated':False,'old_baseline_Arb_exp_calls':0,'native_calls':0,'unknown15_preserved_no_new_levels':True,'new_journal_rows':18257,'cumulative127_129_reservations':21384,'per_root_new_reservations':dict(byroot),'accepted_precision_levels':dict(collections.Counter(s['accepted_precision_bits']for s in states.values())),'terminal_source_pins_sha256':sha(P/'terminal_source_pins.json')}
save('prepared_new_jobs.json',out);print(json.dumps({k:v for k,v in out.items()if k not in ['shards','per_root_new_reservations']},indent=2))
