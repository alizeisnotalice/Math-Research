#!/usr/bin/env python3
"""Future finite111201-case153+g115 experiment. Not run in preparation."""
from pathlib import Path
from fractions import Fraction as F
from datetime import datetime,timezone
import argparse,importlib.util,json,multiprocessing as mp,os,queue,time,traceback
from util155 import sha,digest,atomic
from worker155 import execute
H=Path(__file__).resolve().parent;W=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_retarget155')
def main():
 began=time.monotonic();ap=argparse.ArgumentParser();ap.add_argument('--dispatch-reviewed',action='store_true');ap.add_argument('--reviewed-manifest-sha',required=True);ap.add_argument('--resume-reviewed',action='store_true');args=ap.parse_args();assert args.dispatch_reviewed
 manifest=H/'prepared_manifest.json';assert sha(manifest)==args.reviewed_manifest_sha;package=json.load(open(manifest))
 for name,h in package['repository_files'].items():assert sha(H/name)==h
 for p,h in {**package['source_locks'],**package['prepared_external_files']}.items():assert sha(p)==h
 plan=json.load(open(W/'finite_plan.json'));assert plan['source_objects']==107530 and plan['source_leaves']==111201 and plan['worker_count']==15
 out=W/'run001'
 if out.exists():assert args.resume_reviewed,'existing output requires explicit fingerprint-reviewed resume'
 else:assert not args.resume_reviewed;out.mkdir()
 assert not (out/'full_domain_candidate.json').exists(),'completed run is immutable,not an unfinished prefix to resume'
 attempt=1
 while (out/f'actual_attempt{attempt:04d}_registration.json').exists():attempt+=1
 atomic(out/f'actual_attempt{attempt:04d}_registration.json',{'master_pid':os.getpid(),'registered_utc':datetime.now(timezone.utc).isoformat(),'reviewed_manifest_sha256':args.reviewed_manifest_sha,'fixed_unique_q_objects':107530,'fixed_closed_leaves':111201,'workers':15,'new_log_math_objects':1,'no_H_gradient_or_split':True,'wall_deadline':None,'source_owner':'canonical_index modulo15;one newq per uniqueobject'},True)
 ctx=mp.get_context('spawn');start=ctx.Event();stop=ctx.Event();ready=ctx.Queue();results=ctx.Queue();processes=[];infos={};summaries={};fault=None
 try:
  for spec in plan['shards']:
   p=ctx.Process(target=execute,args=(spec,str(out),args.reviewed_manifest_sha,args.resume_reviewed,start,stop,ready,results),name=f'155-canonical{spec["shard_id"]:02d}');p.start();processes.append(p)
  while len(infos)<15 and not stop.is_set():
   try:r=ready.get(timeout=.1);assert r['shard_id'] not in infos;infos[r['shard_id']]=r
   except queue.Empty:pass
   if any(p.exitcode is not None for p in processes):raise RuntimeError('worker exited before math start')
  assert not stop.is_set() and len(infos)==15 and len({r['pid'] for r in infos.values()})==15
  atomic(out/f'actual_attempt{attempt:04d}_pool_registration.json',{'master_pid':os.getpid(),'workers':[infos[i] for i in range(15)],'before_any_new_log_q_or_corner':True},True);print(json.dumps({'event':'15workers_registered_before_newmath','master_pid':os.getpid(),'worker_pids':[infos[i]['pid'] for i in range(15)],'objects':107530,'leaves':111201}),flush=True)
  logpath=out/'common_log_certificate.json'
  if logpath.exists():assert args.resume_reviewed;common=json.load(open(logpath));assert common['ratio']=='23/22'
  else:
   assert not args.resume_reviewed or not (out/'math_start.json').exists();path=Path(package['shared_log_source']['path']);assert sha(path)==package['shared_log_source']['sha256'];spec=importlib.util.spec_from_file_location('verified_one_shared_log155',path);logmodule=importlib.util.module_from_spec(spec);spec.loader.exec_module(logmodule);common=logmodule.certify_upper_once(64);atomic(logpath,common,True)
  mathstart={'reviewed_manifest_sha256':args.reviewed_manifest_sha,'common_log_path':str(logpath),'common_log_sha256':sha(logpath),'single_log_math_object':True,'new_penalty_code_sha256':package['penalty_source']['sha256']}
  if (out/'math_start.json').exists():assert args.resume_reviewed and json.load(open(out/'math_start.json'))==mathstart
  else:atomic(out/'math_start.json',mathstart,True)
  start.set()
  while len(summaries)<15:
   try:r=results.get(timeout=.1);assert r['shard_id'] not in summaries;summaries[r['shard_id']]=r
   except queue.Empty:pass
   if any(r['status'] not in ('complete','stopped_successful_prefix') for r in summaries.values()):stop.set()
   for i,p in enumerate(processes):
    if p.exitcode not in (None,0) and i not in summaries:stop.set();summaries[i]={'status':'worker_process_fault','shard_id':i,'pid':p.pid,'exitcode':p.exitcode}
   if all(p.exitcode is not None for p in processes) and len(summaries)<15:
    while True:
     try:r=results.get_nowait();summaries.setdefault(r['shard_id'],r)
     except queue.Empty:break
    for i,p in enumerate(processes):summaries.setdefault(i,{'status':'missing_result_fault','shard_id':i,'pid':p.pid})
 except BaseException as exc:
  fault={'exception':repr(exc),'traceback':traceback.format_exc(),'scope':'newmechanism fault;old completeg110 unchanged'};atomic(out/f'actual_attempt{attempt:04d}_fault.json',fault);stop.set();start.set()
 finally:
  if stop.is_set():start.set()
  for p in processes:p.join()
  ready.close();results.close();ready.join_thread();results.join_thread()
 complete=not fault and len(summaries)==15 and all(r['status']=='complete' for r in summaries.values());evidence=[{'pid':p.pid,'exitcode':p.exitcode,'alive':p.is_alive()} for p in processes];pending115=[];minimum=None;seen=set();leaf_count=q_count=passed115=0
 if complete:
  assert all(r['exitcode']==0 and not r['alive'] for r in evidence);Ulog=F(json.load(open(out/'common_log_certificate.json'))['upper'])
  for i in range(15):
   r=summaries[i];leaf_count+=r['next_leaf'];q_count+=r['next_q'];passed115+=r['passed115'];key=(F(r['minimum110'][0]),r['minimum110'][1]);minimum=key if minimum is None else min(minimum,key)
   for line in open(out/f'shard{i:02d}'/'leaf_results.jsonl'):
    row=json.loads(line)['body'];assert row['leaf_index'] not in seen;seen.add(row['leaf_index']);assert F(row['margin110_new'])>=F(row['saved_old_margin110'])>0 and F(row['margin115'])==F(row['margin110_new'])-Ulog
    if not row['passed115']:pending115.append({'leaf_index':row['leaf_index'],'node':row['node'],'canonical_id':row['canonical_id'],'margin115':row['margin115'],'margin110_new':row['margin110_new']})
  assert seen==set(range(111201)) and q_count==107530 and leaf_count==111201 and passed115+len(pending115)==111201;pending115.sort(key=lambda r:r['leaf_index']);floor=F(11,100000)*(1+minimum[0]);atomic(out/'target115_pending.json',pending115,True);full={'complete_old110_domain_preserved':True,'full110_leaf_count':111201,'source_forest_unchanged':True,'minimum_margin110_new':str(minimum[0]),'minimum_witness_leaf_index':minimum[1],'minimum_margin115':str(minimum[0]-Ulog),'candidate_floor':str(floor),'candidate_floor_derivation':'g0*(1+m_min),exp(m)>=1+m;no exp evaluation','candidate_floor_domain':'H128 on F;fullG extension onlyvia029 afterindependent/primary acceptance','target115_passed':passed115,'target115_pending':len(pending115),'old_complete_domain_not_pending':True,'acceptance_status':'candidate only'};atomic(out/'full_domain_candidate.json',full,True)
  for p,h in package['source_locks'].items():assert sha(p)==h
 summary={'status':'complete_fixed153_jointretarget_candidate' if complete else 'fault_or_stopped_prefix_preserved','master_pid':os.getpid(),'worker_evidence':evidence,'q_objects':q_count,'leaf_cases':leaf_count,'target115_passed':passed115,'target115_pending':len(pending115) if complete else None,'old_qualified110_domain_unchanged':True,'new110_margin_updates_complete':complete,'fault':fault,'new_H_gradient_or_split_calls':0,'new_log_math_objects':1 if (out/'common_log_certificate.json').exists() else 0,'new_exp_calls':0,'oldq_corner_baseline_runs':0,'formal_g_or_C5_updated':False,'wall_seconds':time.monotonic()-began};atomic(out/f'actual_attempt{attempt:04d}_summary.json',summary,True);print(json.dumps(summary),flush=True)
if __name__=='__main__':main()
