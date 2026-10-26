"""One finite2048-parent batch. Newchildren only; complete-parent prefix is always fullF."""
from pathlib import Path
from fractions import Fraction as F
import gzip,json,hashlib,multiprocessing as mp,os,time,sys
from flint import arb,ctx
import kernel as K
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
H=Path(__file__).resolve().parent;S=K.S;N=1<<160

def budget(low,up,frontcount,reg):
 lo=F(low,N);uh=F(up,N);E=F(reg['E076_upper']);g=F(reg['valid_g']);B=K.ep(arb(K.rq(uh))**(arb(4)/5)+arb(K.rq(E)),True);I=K.ep(arb(K.rq(B))**(arb(5)/4),True);cost=K.ep((5*arb(K.rq(B))/arb(K.rq(g)))**(arb(5)/4),True);margin=32**4*g**5-5**5*B**5;return {'LH_truncated_integral':str(lo),'UH_truncated_integral':str(uh),'UH_display':float(uh),'B_upper':str(B),'I_upper':str(I),'I_upper_display':float(I),'valid_g':str(g),'cost_upper':str(cost),'cost_upper_display':float(cost),'strict32':margin>0,'exact_32_gate_margin':str(margin),'source_E076_upper':str(E),'I_strictly_below_old_cap':I<F(reg['old_complete_I_cap']),'old_complete_I_cap':reg['old_complete_I_cap'],'physical_frontier':frontcount,'no_tube_added':True,'formal_C5_changed_by_author':False,'qualification_status':'candidate_continuation_pending095_102_and_main'}

def install(state,parent,children,rows,jobmap):
 nodes=state['nodes'];p=nodes[parent];assert p['children'] is None;root=parent
 while nodes[root]['parent'] is not None:root=nodes[root]['parent']
 oldL,oldU=nodes[root]['effL'],nodes[root]['effU']
 for k in children:
  assert k not in nodes;r=rows[k];j=jobmap[k];b=r['bounds'];lo=j['weight']*b['integral_lower_i'];up=j['weight']*b['integral_upper_i'];nodes[k]={'key':k,'level':j['level'],'indices':j['indices'],'weight':j['weight'],'selfL':lo,'selfU':up,'effL':lo,'effU':up,'children':None,'parent':parent,'source':'103'}
 assert sum(nodes[k]['weight'] for k in children)==16*p['weight'];rawL=sum(nodes[k]['selfL'] for k in children);rawU=sum(nodes[k]['selfU'] for k in children);p['children']=children;k=parent
 while k is not None:
  n=nodes[k];n['effL']=max(n['selfL'],sum(nodes[c]['effL'] for c in n['children'])) if n['children'] else n['selfL'];n['effU']=min(n['selfU'],sum(nodes[c]['effU'] for c in n['children'])) if n['children'] else n['selfU'];assert n['effL']<=n['effU'],('hierarchy_not_separated',k,n['effL'],n['effU']);k=n['parent']
 return nodes[root]['effL']-oldL,nodes[root]['effU']-oldU,rawL,rawU

def main():
 began=time.monotonic();ctx.prec=384;import argparse
 ap=argparse.ArgumentParser();ap.add_argument('--resume',action='store_true');args=ap.parse_args();reg=json.loads((H/'registration.json').read_text());assert K.sha(__file__)==reg['run_sha256'] and K.sha(H/'kernel.py')==reg['kernel_sha256'];assert K.sha(H/'jobs.json.gz')==reg['jobs_sha256'] and K.sha(H/'ratio_catalog.json')==reg['catalog_sha256'];fp=hashlib.sha256((K.sha(__file__)+K.sha(H/'kernel.py')+K.sha(H/'registration.json')).encode()).hexdigest();jobs=json.loads(gzip.decompress((H/'jobs.json.gz').read_bytes()));jobmap={j['key']:j for j in jobs};pre=json.loads((H/'preflight_registration.json').read_text());edges=pre['children_by_parent'];parents=[r['key'] for r in pre['parents']];sp=S/'state.json.gz';initial=Path(reg['initial_state_path']);rp=H/'run_receipt.json';completed_file=S/'completed_child_results.json.gz'
 if sp.exists():
  assert args.resume;state=json.loads(gzip.decompress(sp.read_bytes()));assert state['adaptive103']['fingerprint']==fp and not state['adaptive103']['complete'];results={r['key']:r for r in json.loads(gzip.decompress(completed_file.read_bytes()))} if completed_file.exists() else {}
 else:
  assert not args.resume and not rp.exists();assert K.sha(initial)==reg['initial_state_sha256'];state=json.loads(gzip.decompress(initial.read_bytes()));state['adaptive103']={'fingerprint':fp,'parents_complete':[],'complete':False,'globalL_units':sum(state['nodes'][k]['effL'] for k in state['roots']),'globalU_units':sum(state['nodes'][k]['effU'] for k in state['roots']),'event_hash':'0'*64,'event_sequence':0};results={};K.gz(sp,state)
 meta=state['adaptive103'];front=set(state['frontier']);before=budget(meta['globalL_units'],meta['globalU_units'],len(front),reg);applied=set(meta['parents_complete']);journal=S/'parent_events.jsonl';receipt={'fingerprint':fp,'pid':os.getpid(),'began_unix':time.time(),'registered_parents':2048,'registered_children':len(jobs),'before':before,'status':'running','source_registration_sha256':K.sha(H/'registration.json')};K.atomic(rp,receipt)
 def checkpoint():
  state['frontier']=sorted(front);K.gz(completed_file,list(results.values()));K.gz(sp,state);current=budget(meta['globalL_units'],meta['globalU_units'],len(front),reg);K.atomic(H/'progress.json',{'status':'running','master_pid':os.getpid(),'completed_children':len(results),'complete_parents':len(applied),'total_parents':2048,'total_children':len(jobs),'state_sha256':K.sha(sp),'budget':current})
 def try_parent(pk):
  if pk in applied or not all(c in results for c in edges[pk]):return
  dL,dU,rawL,rawU=install(state,pk,edges[pk],results,jobmap);assert dL>=0 and dU<=0;meta['globalL_units']+=dL;meta['globalU_units']+=dU;applied.add(pk);meta['parents_complete'].append(pk);front.remove(pk);front.update(edges[pk]);now=budget(meta['globalL_units'],meta['globalU_units'],len(front),reg);event={'seq':meta['event_sequence']+1,'parent':pk,'children':edges[pk],'child_raw_L_units':rawL,'child_raw_U_units':rawU,'root_delta_L_units':dL,'root_delta_U_units':dU,'global_L_units':meta['globalL_units'],'global_U_units':meta['globalU_units'],'previous_hash':meta['event_hash']};eh=hashlib.sha256(K.serial(event).encode()).hexdigest();event['event_hash']=eh
  with journal.open('a') as f:f.write(K.serial(event)+'\n');f.flush()
  meta['event_hash']=eh;meta['event_sequence']=event['seq']
  if now['strict32'] and not (H/'first_strict32.json').exists():
   checkpoint();snap=S/'first_strict32_state.json.gz';K.gz(snap,state);K.atomic(H/'first_strict32.json',{'status':'numerical_gate_pending_main_and_independent_qualification','complete_parents':len(applied),'completed_children':len(results),'budget':now,'state_sha256':K.sha(snap),'state_path':str(snap),'event_hash':eh});print('FIRST_STRICT32',K.serial(now),flush=True)
  if len(applied)%32==0:checkpoint();print('PARENTS_COMPLETE',len(applied),'children',len(results),'wall',round(time.monotonic()-began,3),'I',now['I_upper_display'],'cost',now['cost_upper_display'],flush=True)
 # Existing successful outputs are reused, never rerun; partially completed point cp resumes insidekernel.
 for j in jobs:
  if j['key'] in results:continue
  p=K.path(j['key'])
  if p.exists():
   ev=json.loads(gzip.decompress(p.read_bytes()));assert ev['fingerprint']==fp and ev['point']==j['point'];results[j['key']]={'key':j['key'],'status':'reused_completed','bounds':ev['bounds'],'evidence_sha256':K.sha(p),'runtime':ev['runtime']}
 for pk in parents:try_parent(pk)
 todo=[j for j in jobs if j['key'] not in results];table=json.loads(gzip.decompress(Path(reg['table_path']).read_bytes()));assert K.sha(reg['table_path'])==reg['table_sha256'];catalog=json.loads((H/'ratio_catalog.json').read_text());pool=mp.get_context('spawn').Pool(12,initializer=K.init,initargs=(table,reg['time_endpoints'],catalog,fp));receipt['worker_pids_dispatch']=[p.pid for p in pool._pool];K.atomic(rp,receipt);fault=None
 try:
  for r in pool.imap_unordered(K.worker,todo):
   if r['status']=='failed':fault=r;pool.terminate();break
   results[r['key']]=r;try_parent(jobmap[r['key']]['parent'])
   if len(results)%128==0:checkpoint();print('CHILDREN_COMPLETE',len(results),'parents',len(applied),'wall',round(time.monotonic()-began,3),flush=True)
  if fault is None:pool.close()
  pool.join();receipt['workers_joined']=True
 except BaseException as e:
  pool.terminate();pool.join();receipt.update({'workers_joined':True,'exception':repr(e)})
 checkpoint();complete=fault is None and len(results)==len(jobs) and len(applied)==2048;meta['complete']=complete;state['frontier']=sorted(front);K.gz(sp,state);after=budget(meta['globalL_units'],meta['globalU_units'],len(front),reg);receipt.update({'status':'complete' if complete else 'stopped_on_real_fault','fault':fault,'after':after,'complete_parents':len(applied),'completed_children':len(results),'worker_pids':sorted({r['runtime']['pid'] for r in results.values()}),'wall_seconds':time.monotonic()-began,'new_H_calls':sum(r['runtime']['new_H_calls'] for r in results.values()),'reused_H_values':sum(r['runtime']['reused_H_values'] for r in results.values()),'new_DD_cones_this_execution':sum(r['runtime']['new_DD_cones'] for r in results.values()),'event_hash':meta['event_hash'],'state_sha256':K.sha(sp),'first_strict32_saved':(H/'first_strict32.json').exists()});K.atomic(rp,receipt)
 if complete:
  compact=[{'key':k,'level':jobmap[k]['level'],'indices':jobmap[k]['indices'],'weight':jobmap[k]['weight'],'parent':jobmap[k]['parent'],'bounds':r['bounds'],'evidence_sha256':r['evidence_sha256']} for k,r in sorted(results.items())];K.gz(H/'children_certificates.json.gz',compact);K.gz(H/'physical_frontier.json.gz',state['frontier']);K.gz(H/'numeric_hierarchy.json.gz',{'old_source095_manifest_sha256':reg['source095_manifest_sha256'],'expanded_parent_nodes':{k:state['nodes'][k] for k in parents},'budget':after});receipt['certificates_sha256']=K.sha(H/'children_certificates.json.gz');K.atomic(rp,receipt)
 K.atomic(H/'progress.json',{'status':receipt['status'],'complete_parents':len(applied),'completed_children':len(results),'budget':after,'state_sha256':K.sha(sp)});print(K.serial(receipt),flush=True)
if __name__=='__main__':main()
