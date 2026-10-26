#!/usr/bin/env python3
"""Explicitly reviewed15-process finite dispatcher. No automatic next experiment."""
from pathlib import Path
from datetime import datetime,timezone
import argparse,json,multiprocessing as mp,os,time,queue,traceback
from common import sha,digest,write,CORE_SHA
from shard_worker import execute
H=Path(__file__).resolve().parent;W=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_parallel137/next2350_preflight')
def main():
 began=time.monotonic();ap=argparse.ArgumentParser();ap.add_argument('--dispatch-reviewed',action='store_true');ap.add_argument('--reviewed-manifest-sha',required=True);ap.add_argument('--resume-reviewed',action='store_true');args=ap.parse_args();assert args.dispatch_reviewed,'main must review this dispatcher before actual corner execution'
 manifest=H/'next2350_prepared_manifest.json';assert sha(manifest)==args.reviewed_manifest_sha;m=json.load(open(manifest))
 for name,h in m['repository_files'].items():assert sha(H/name)==h
 for p,h in m['source_locks'].items():assert sha(p)==h
 for p,h in m['prepared_external_files'].items():assert sha(p)==h
 plan=json.load(open(W/'shard_plan.json'));assert plan['worker_count']==15 and plan['total_pairs']==39510 and plan['total_children']==2350 and plan['core_sha256']==CORE_SHA
 out=W/'run001'
 if out.exists():assert args.resume_reviewed,'existing run requires fingerprint-verified explicit resume'
 else:assert not args.resume_reviewed;out.mkdir()
 attempt=1
 while (out/f'master_attempt{attempt:04d}_registration.json').exists():attempt+=1
 write(out/f'master_attempt{attempt:04d}_registration.json',{'master_pid':os.getpid(),'registered_utc':datetime.now(timezone.utc).isoformat(),'attempt':attempt,'prepared_manifest_sha256':args.reviewed_manifest_sha,'core_sha256':CORE_SHA,'worker_count':15,'fixed_shards':15,'fixed_pairs':39510,'fixed_children':2350,'resume_reviewed':args.resume_reviewed,'shard_plan_sha256':sha(W/'shard_plan.json'),'source_loading':'each shard loads only frozenneeded subset,never full111331','peak_memory_estimate_GiB':[0.8,2.5],'memory_estimate_not_cap':True,'wall_deadline':None,'new_H_gradient_kernel_log_exp_calls':0})
 ctx=mp.get_context('spawn');stop=ctx.Event();start=ctx.Event();ready=ctx.Queue();results=ctx.Queue();processes=[];readies={};summaries={};master_fault=None
 try:
  for spec in plan['shards']:
   proc=ctx.Process(target=execute,args=(spec,str(out),args.reviewed_manifest_sha,args.resume_reviewed,stop,start,ready,results),name=f'upper137-shard{spec["shard_id"]:02d}');proc.start();processes.append(proc)
  # No arithmetic begins until all actualworker PIDs have been persisted.
  while len(readies)<15 and not stop.is_set():
   try:r=ready.get(timeout=.5);assert r['shard_id'] not in readies;readies[r['shard_id']]=r
   except queue.Empty:pass
   for proc in processes:
    if proc.exitcode is not None and len(readies)<15:raise RuntimeError(('worker exited before start gate',proc.pid,proc.exitcode))
  if stop.is_set():raise RuntimeError('worker initialization/recovery fault,not starting corner work')
  assert len({r['pid'] for r in readies.values()})==15 and all(r['fingerprint']['core_sha256']==CORE_SHA for r in readies.values())
  write(out/f'master_attempt{attempt:04d}_pool_registration.json',{'master_pid':os.getpid(),'workers':[readies[i] for i in range(15)],'persisted_before_any_new_corner':True,'core_sha256':CORE_SHA})
  print(json.dumps({'event':'15workers_registered_before_corner_start','master_pid':os.getpid(),'worker_pids':[readies[i]['pid'] for i in range(15)],'shards':15,'pairs':39510,'prepared_manifest_sha256':args.reviewed_manifest_sha}),flush=True);start.set()
  while len(summaries)<15:
   try:r=results.get(timeout=.5);assert r['shard_id'] not in summaries;summaries[r['shard_id']]=r
   except queue.Empty:pass
   for i,proc in enumerate(processes):
    if proc.exitcode not in (None,0) and i not in summaries:stop.set();summaries[i]={'status':'process_fault','shard_id':i,'pid':proc.pid,'exitcode':proc.exitcode}
   if any(r['status'] not in ('complete','stopped_prefix') for r in summaries.values()):stop.set()
   # Each surviving finiteworker observes stop at the next pair boundary.
   if all(p.exitcode is not None for p in processes) and len(summaries)<15:
    while True:
     try:r=results.get_nowait();summaries.setdefault(r['shard_id'],r)
     except queue.Empty:break
    for i,p in enumerate(processes):summaries.setdefault(i,{'status':'missing_result_fault','shard_id':i,'pid':p.pid,'exitcode':p.exitcode})
 except BaseException as exc:
  master_fault={'exception':repr(exc),'traceback':traceback.format_exc()};stop.set();start.set();write(out/f'master_attempt{attempt:04d}_fault.json',master_fault)
 finally:
  if stop.is_set():start.set()
  for proc in processes:proc.join()
  ready.close();results.close();ready.join_thread();results.join_thread()
 # Preserve all valid child proofs from checkpoints, even on a stopped run.
 qualified={};ack=evaluated=skipped=0;complete=True
 for spec in plan['shards']:
  cp=out/f'shard{spec["shard_id"]:02d}'/'checkpoint.json'
  if not cp.exists():complete=False;continue
  state=json.load(open(cp));assert state['fingerprint']['input_sha256']==spec['input_sha256'] and state['fingerprint']['prepared_manifest_sha256']==args.reviewed_manifest_sha
  ack+=state['next_pair'];evaluated+=state['evaluated'];skipped+=state['skipped'];complete&=state.get('status')=='complete' and state['next_pair']==spec['pair_count']
  for k,cert in state['qualified_children'].items():assert int(k)%15==spec['shard_id'] and k not in qualified;qualified[k]=cert
 allchildren=[]
 for spec in plan['shards']:
  data=json.load(open(spec['input_path']));allchildren.extend((int(k),v['node']) for k,v in data['children'].items())
 allchildren.sort();assert [i for i,n in allchildren]==list(range(2350));pending=[node for i,node in allchildren if str(i) not in qualified]
 if complete and not master_fault and all(r.get('status')=='complete' for r in summaries.values()):assert ack==evaluated+skipped==39510;status='complete_fixed39510_objects'
 else:status='fault_or_stopped_prefix_preserved'
 state={'status':status,'prepared_manifest_sha256':args.reviewed_manifest_sha,'qualified_children':qualified,'pending_children':pending,'source137_newH_base_passed':108734,'source137_newH_base_pending':2350,'candidate_total_passed':108734+len(qualified),'source_splits_unchanged':58724,'new_splits':[],'faults':master_fault,'S_credit':0,'globalg110_qualified':False};write(out/f'master_attempt{attempt:04d}_state.json',state)
 evidence=[{'pid':p.pid,'exitcode':p.exitcode,'alive':p.is_alive()} for p in processes]
 summary={'status':status,'master_pid':os.getpid(),'actual_worker_count':len(processes),'process_evidence':evidence,'acknowledged_pairs':ack,'exact_single_source_tests':evaluated,'exact_corners':16*evaluated,'skip_receipts':skipped,'new_child_passes':len(qualified),'pending':len(pending),'candidate_total_passed':108734+len(qualified),'new_H_gradient_kernel_log_exp_calls':0,'new_splits':0,'core_sha256':CORE_SHA,'master_fault':master_fault,'worker_summaries':[summaries.get(i) for i in range(15)],'wall_seconds':time.monotonic()-began,'no_automatic_residual_H':True};write(out/f'master_attempt{attempt:04d}_summary.json',summary);print(json.dumps(summary),flush=True)
if __name__=='__main__':main()
