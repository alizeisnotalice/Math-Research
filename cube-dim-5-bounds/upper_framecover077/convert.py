#!/usr/bin/env python3
"""Frozen full-frame Q=I4/16 conversion; one finite process per block.

No previous implementation is imported. All numerical work uses Fraction.
The immutable source contains all retired records, retained by source index.
"""
import argparse
from collections import defaultdict
from datetime import datetime, timezone
from fractions import Fraction as F
import hashlib
from itertools import product
import json
import os
from pathlib import Path
import resource
import signal
import sys
import time

SOURCE=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_cover066/final/state.json')
SOURCE_HASH='42c67156078cb0eb134910a5248821eb3c826c96981afb3d1d5b2fe3242a871b'
LOCKS={
 '/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-lower/cube-dim-5-bounds/audit_frame075/report.md':'3308d1e5c660cd11d73d47a013acf794bd1f596afdb5a113573a5bfd7481f225',
 '/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-lower/cube-dim-5-bounds/audit_frame075/manifest.json':'1d5e270eb1664d714ec80f0b371420d40028b24c3943196b9ab41a0b3de1f9d1',
 '/Users/zhengzhihao/Desktop/T/Math-Research-cube5-lower-audit-cont/cube-dim-5-bounds/audit_conversion069/final/report.md':'bda9cc665e71c5004518b5881f5c625010a6d988a68a760c8666af38d183675f',
}
ENC=json.JSONEncoder(sort_keys=True,separators=(',',':'))
def sha(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  while chunk:=f.read(8*1024*1024):h.update(chunk)
 return h.hexdigest()
def digest(obj):
 h=hashlib.sha256()
 for chunk in ENC.iterencode(obj):h.update(chunk.encode())
 return h.hexdigest()
def write(path,obj):
 temp=path.with_suffix(path.suffix+'.tmp')
 with open(temp,'w') as f:
  for chunk in ENC.iterencode(obj):f.write(chunk)
  f.write('\n');f.flush();os.fsync(f.fileno())
 os.replace(temp,path)
def key(node):return (node[0],node[1])
def physical(node,roots):
 lo=[F(x,32) for x in roots[node[0]]];hi=[x+F(1,32) for x in lo]
 for b in node[1]:
  widths=[y-x for x,y in zip(lo,hi)];j=max(range(4),key=lambda j:widths[j]);m=(lo[j]+hi[j])/2
  assert b in '01'
  if b=='0':hi[j]=m
  else:lo[j]=m
 return lo,hi
def partition(roots,passed,queue,splits):
 leaves=set();byroot=defaultdict(list)
 for n in [r['node'] for r in passed]+queue:
  k=key(n);assert k not in leaves;leaves.add(k);byroot[k[0]].append(k[1])
 assert len(byroot)==len(roots)
 for r in range(len(roots)):
  paths=sorted(byroot[r]);assert all(not b.startswith(a) for a,b in zip(paths,paths[1:]))
  depth=max(map(len,paths));assert sum(1<<(depth-len(p)) for p in paths)==1<<depth
 # Every retained split is an internal tree node, and its recorded split axis
 # equals deterministic longest-side splitting (ties select first axis).
 splitkeys=set()
 for s in splits:
  r,p=s['node'];assert (r,p) not in splitkeys;splitkeys.add((r,p))
  assert s['children']==[[r,p+'0'],[r,p+'1']]
  assert s['axis']==len(p)%4
 for r,p in leaves|splitkeys:
  if p:assert (r,p[:-1]) in splitkeys
 assert len(leaves)==len(roots)+len(splits)
 return {'antichain':True,'roots_with_Kraft_sum_one':len(roots),'leaf_count':len(leaves),'internal_count':len(splitkeys)}
def qualify(split,state,index):
 c=split['failed_certificate'];assert c and not c['passed']
 lo,hi=physical(split['node'],state['roots']);rec=state['cache'][c['anchor']]
 point=list(map(F,rec['point']));assert point[0]==1 and ','.join(map(str,point))==c['anchor']
 assert all(a<=x<=b for a,x,b in zip(lo,point[1:],hi))
 values=[tuple(map(F,v)) for v in rec['values']];assert len(values)==5
 Hlo,Hhi=values[0];assert F(1,10000)<Hlo<=Hhi
 gd=[tuple(map(F,v)) for v in c['loggradient_directed']];assert len(gd)==4
 for (a,b),(gl,gu) in zip(values[1:],gd):
  assert a<=b
  ratios=[a/Hlo,a/Hhi,b/Hlo,b/Hhi];assert gl<=min(ratios)<=max(ratios)<=gu
 Delta=F(c['Delta_lower']);assert Delta>0
 lower=[a-x for a,x in zip(lo,point[1:])];upper=[b-x for b,x in zip(hi,point[1:])]
 linear=lambda v:sum(min(gl*x,gu*x) for x,(gl,gu) in zip(v,gd))
 old=[]
 for bits in product((0,1),repeat=4):
  v=tuple((upper[j] if bits[j] else lower[j]) for j in range(4))
  q=sum(x*x for x in v)/4-sum(v)**2/36
  old.append((Delta+linear(v)-4096*q,bits))
 om,ob=min(old);assert om==F(c['minimum_log_margin'])<=0 and list(ob)==c['worst_corner_bits']
 rows=[]
 for bits in product((0,1),repeat=4):
  v=tuple((upper[j] if bits[j] else lower[j]) for j in range(4))
  q=sum(x*x for x in v)/16
  rows.append((Delta+linear(v)-4096*q,bits))
 nm,nv=min(rows);assert nm>=om
 cert={'mechanism':'full_frame_16','passed':nm>0,'minimum_log_margin':str(nm),
  'worst_corner_bits':list(nv),'checked_points':len(rows),
  'quadratic':'Q_full=I4/16',
  'Delta_lower':c['Delta_lower'],'loggradient_directed':c['loggradient_directed'],'anchor':c['anchor']}
 return {'source_split_index':index,'node':split['node'],'source_failed_certificate_sha256':digest(c),
  'box_lower':list(map(str,lo)),'box_upper':list(map(str,hi)),'anchor_inside':True,
  'old_minimum_log_margin':str(om),'certificate':cert}

def run():
 began=time.monotonic()
 ap=argparse.ArgumentParser();ap.add_argument('--work-dir',type=Path,required=True)
 ap.add_argument('--resume',action='store_true');ap.add_argument('--reconstruct',action='store_true')
 args=ap.parse_args();work=args.work_dir.resolve();work.mkdir(parents=True,exist_ok=True)
 rows=[];cursor=0;prefix=None;attempt_dir=None
 def deadline(signum,frame):
  if prefix is not None:
   write(prefix,{'stage':'wall_guard_recoverable','registration_sha256':sha(prefix.parent/'registration.json'),
    'pid':os.getpid(),'cursor':cursor,'rows':rows,'total_wall_seconds':time.monotonic()-began})
  write(work/f'wall_guard_{os.getpid()}.json',{'pid':os.getpid(),'stage':'590_second_total_wall_guard',
   'wall_seconds':time.monotonic()-began,'exit_code':75,'new_kernel_chain_started':False})
  raise SystemExit(75)
 signal.signal(signal.SIGALRM,deadline);signal.alarm(590)
 codehash=sha(__file__);runreg=work/'registration.json'
 assert sha(SOURCE)==SOURCE_HASH
 for p,h in LOCKS.items():assert sha(p)==h
 mathroot=Path(next(p for p in LOCKS if p.endswith('audit_frame075/manifest.json'))).parent
 for f in json.load(open(mathroot/'manifest.json'))['files']:assert sha(mathroot/f['path'])==f['sha256']
 if not runreg.exists():
  assert not args.resume and not args.reconstruct
  write(runreg,{'task':'upper_framecover077','source':str(SOURCE),'source_sha256':SOURCE_HASH,
   'source_counts':{'passed':94644,'splits':56048,'queue':13764,'cache':178722,'roots':52360},
   'code':str(Path(__file__).resolve()),'code_sha256':codehash,'locks':LOCKS,
   'registered_utc':datetime.now(timezone.utc).isoformat(),'epsilon':'1/64','target':'1/10000',
   'mathematics':'075 full frame Q_full=I4/16; globally concave interval-affine minimum minus positive quadratic; exactly 16 corners',
   'resources':{'processes':1,'threads':1,'MATLAB_calls':0,'new_H_gradient_kernel_log_exp_Gamma_calls':0},
   'per_process_limits':{'tested_parents':4096,'total_wall_seconds':600,'normal_checkpoint_seconds':540,'alarm_seconds':590},
   'selection':'live untested source split with current queue descendant; (path depth,-current queue descendant count,root index,path); exact prefix registered before computation',
   'schema':{'registration':'source/code/math hashes; each block parent checkpoint and logical frontier hashes; exact selected source index/node/queue-count list',
    'certificate':'source_split_index,node,source_failed_certificate_sha256,physical box; certificate full_frame_16,Delta,gradient,anchor,min margin,worst_corner_bits,checked_points=16',
    'retired':'per-block source passed/queue/split index,node,unique accepted ancestor; source immutable records retain full provenance',
    'checkpoint':'immutable chain of parent checkpoint SHA, registration/certificates/retired/summary SHA, frontier hash',
    'logical_frontier':'digest({source_sha256,active_source_indices,accepted_certificates,tested_source_split_indices})',
    'final_state':'only passed/queue/splits transformed; task077/status077/fingerprint077/conversion077 added; all other source values unchanged'},
   'next_process_command':f'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 {sys.executable} -B {Path(__file__).resolve()} --work-dir {work}',
   'reconstruct_command':f'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 {sys.executable} -B {Path(__file__).resolve()} --work-dir {work} --reconstruct',
   'new_kernel_chain_authorized':False,'one_block_per_process':True,'STOP_file':str(work/'STOP')})
 reg=json.load(open(runreg));assert reg['code_sha256']==codehash and reg['source_sha256']==SOURCE_HASH
 state=json.load(open(SOURCE));assert not state['inflight'] and not state['resume_jobs'] and not state.get('active')
 assert state['epsilon']=='1/64' and state['target']=='1/10000'
 assert [len(state[k]) for k in ('passed','splits','queue','cache','roots')]==[94644,56048,13764,178722,52360]
 arrays={k:state[k] for k in ('passed','queue','splits')};active={k:set(range(len(v))) for k,v in arrays.items()}
 newleaves=[];tested=set();allaccepted=defaultdict(list)
 cachehash=digest(state['cache']);rootshash=digest(state['roots'])
 source_other_hashes={k:digest(v) for k,v in state.items() if k not in ('passed','splits','queue','cache','roots')}
 partition(state['roots'],arrays['passed'],arrays['queue'],arrays['splits'])
 def frontier_hash():
  return digest({'source_sha256':SOURCE_HASH,'active_source_indices':{k:sorted(v) for k,v in active.items()},
   'accepted_certificates':newleaves,'tested_source_split_indices':sorted(tested)})
 def selection():
  counts=defaultdict(int)
  for i in active['queue']:
   r,p=arrays['queue'][i]
   for depth in range(len(p)):counts[(r,p[:depth])]+=1
  choices=sorted((i for i in active['splits'] if i not in tested and counts[key(arrays['splits'][i]['node'])]>0),
   key=lambda i:(len(arrays['splits'][i]['node'][1]),-counts[key(arrays['splits'][i]['node'])],*arrays['splits'][i]['node']))
  return choices,counts
 def current_lists():return {k:[v[i] for i in sorted(active[k])] for k,v in arrays.items()}
 def apply_rows(got,expected_retired=None):
  accepted=defaultdict(list)
  for row in got:
   idx=row['source_split_index'];assert idx not in tested and idx in active['splits'];tested.add(idx)
   assert arrays['splits'][idx]['node']==row['node']
   assert digest(arrays['splits'][idx]['failed_certificate'])==row['source_failed_certificate_sha256']
   if row['certificate']['passed']:
    r,p=row['node'];assert not any(p.startswith(q) or q.startswith(p) for q in allaccepted[r])
    accepted[r].append(p);allaccepted[r].append(p);newleaves.append({'node':row['node'],'certificate':row['certificate']})
  retired={k:[] for k in arrays}
  for name,items in arrays.items():
   remove=[]
   for idx in sorted(active[name]):
    n=items[idx] if name=='queue' else items[idx]['node'];r,p=n
    matches=[q for q in accepted[r] if p.startswith(q) and (name=='splits' or len(p)>len(q))]
    assert len(matches)<=1
    if matches:remove.append(idx);retired[name].append({'source_index':idx,'node':n,'unique_new_parent':[r,matches[0]]})
   active[name].difference_update(remove)
  if expected_retired is not None:assert retired==expected_retired
  return retired
 prevcp=None;block=0
 while (work/f'batch{block+1:04d}'/'checkpoint.json').exists():
  bdir=work/f'batch{block+1:04d}';cp=json.load(open(bdir/'checkpoint.json'))
  assert cp['stage']=='complete' and cp['parent_checkpoint_sha256']==prevcp
  br=json.load(open(bdir/'registration.json'));assert br['parent_frontier_sha256']==frontier_hash()
  assert br['parent_checkpoint_sha256']==prevcp and cp['registration_sha256']==sha(bdir/'registration.json')
  assert cp['certificates_sha256']==sha(bdir/'certificates.json') and cp['retired_sha256']==sha(bdir/'retired.json')
  choices,counts=selection()
  expected=[{'source_split_index':i,'node':arrays['splits'][i]['node'],'queue_descendants':counts[key(arrays['splits'][i]['node'])]} for i in choices[:4096]]
  assert br['selected']==expected
  apply_rows(json.load(open(bdir/'certificates.json')),json.load(open(bdir/'retired.json')))
  assert cp['frontier_sha256']==frontier_hash()
  prevcp=sha(bdir/'checkpoint.json');block+=1
 def audit():
  lists=current_lists();a=partition(state['roots'],lists['passed']+newleaves,lists['queue'],lists['splits'])
  assert digest(state['cache'])==cachehash and digest(state['roots'])==rootshash
  assert {k:digest(state[k]) for k in source_other_hashes}==source_other_hashes
  return a
 choices,counts=selection()
 if not args.reconstruct and choices and active['queue']:
  if (work/'STOP').exists():raise SystemExit('STOP before dispatch; no source mutated')
  bdir=work/f'batch{block+1:04d}';bdir.mkdir(exist_ok=True);brpath=bdir/'registration.json'
  selected=choices[:4096]
  br={'task':'upper_framecover077','block':block+1,'source_state':str(SOURCE),'source_state_sha256':SOURCE_HASH,
   'run_registration_sha256':sha(runreg),'code_sha256':codehash,'parent_checkpoint_sha256':prevcp,
   'parent_frontier_sha256':frontier_hash(),'max_tested_per_process':4096,'total_wall_seconds':600,
   'normal_checkpoint_seconds':540,'alarm_seconds':590,
   'selected':[{'source_split_index':i,'node':arrays['splits'][i]['node'],'queue_descendants':counts[key(arrays['splits'][i]['node'])]} for i in selected]}
  if brpath.exists():assert args.resume and json.load(open(brpath))==br
  else:write(brpath,br)
  attempt_dir=bdir/'attempts';attempt_dir.mkdir(exist_ok=True)
  ordinal=len(list(attempt_dir.glob('*_registration.json')))+1;ar=attempt_dir/f'{ordinal:04d}_registration.json'
  prefix=bdir/'prefix_checkpoint.json'
  write(ar,{'pid':os.getpid(),'registered_utc':datetime.now(timezone.utc).isoformat(),
   'block_registration_sha256':sha(brpath),'parent_prefix_sha256':sha(prefix) if prefix.exists() else None,
   'total_wall_seconds':600,'normal_checkpoint_seconds':540,'alarm_seconds':590,
   'wall_elapsed_at_attempt_registration':time.monotonic()-began})
  accepted=defaultdict(list)
  if prefix.exists():
   pref=json.load(open(prefix));assert args.resume and pref['registration_sha256']==sha(brpath)
   rows=pref['rows'];cursor=pref['cursor']
   for row in rows:
    if row['certificate']['passed']:accepted[row['node'][0]].append(row['node'][1])
  def saveprefix(stage):write(prefix,{'stage':stage,'registration_sha256':sha(brpath),'pid':os.getpid(),'cursor':cursor,'rows':rows,'total_wall_seconds':time.monotonic()-began})
  saveprefix('ready');qstart=time.monotonic()
  try:
   while cursor<len(selected) and len(rows)<4096 and time.monotonic()-began<540:
    if (work/'STOP').exists():saveprefix('stopped_before_new_dispatch');return
    idx=selected[cursor];cursor+=1;n=arrays['splits'][idx]['node']
    if any(n[1].startswith(p) for p in accepted[n[0]]):continue
    row=qualify(arrays['splits'][idx],state,idx);rows.append(row)
    if row['certificate']['passed']:accepted[n[0]].append(n[1])
    if len(rows)%16==0:saveprefix('qualification_prefix')
   qseconds=time.monotonic()-qstart;saveprefix('qualification_complete_or_540_wall_checkpoint')
   retired=apply_rows(rows);write(bdir/'certificates.json',rows);write(bdir/'retired.json',retired)
   a=audit();fh=frontier_hash();remaining,_=selection()
   summary={'task':'upper_framecover077','block':block+1,'pid':os.getpid(),'tested':len(rows),
    'accepted':sum(r['certificate']['passed'] for r in rows),'failed_uncertified':sum(not r['certificate']['passed'] for r in rows),
    'points':sum(r['certificate']['checked_points'] for r in rows),'selected':len(selected),'cursor':cursor,
    'retired':{k:len(v) for k,v in retired.items()},'passed':len(active['passed'])+len(newleaves),
    'split':len(active['splits']),'queue':len(active['queue']),'qualification_seconds':qseconds,
    'cumulative_tested':len(tested),'cumulative_accepted':len(newleaves),'remaining_untested_eligible':len(remaining),
    'cache_entries':len(state['cache']),'cache_canonical_sha256':cachehash,'frontier_sha256':fh,
    'kernel_calls':0,'audit':a,'wall_seconds_to_block':time.monotonic()-began,
    'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
   write(bdir/'summary.json',summary)
   cp={'stage':'complete','block':block+1,'pid':os.getpid(),'parent_checkpoint_sha256':prevcp,
    'registration_sha256':sha(brpath),'certificates_sha256':sha(bdir/'certificates.json'),
    'retired_sha256':sha(bdir/'retired.json'),'summary_sha256':sha(bdir/'summary.json'),
    'frontier_sha256':fh,'cumulative_tested':len(tested),'cumulative_accepted':len(newleaves)}
   write(bdir/'checkpoint.json',cp);prevcp=sha(bdir/'checkpoint.json');block+=1
   write(work/'CURRENT.json',{'block':block,'checkpoint_path':str(bdir/'checkpoint.json'),'checkpoint_sha256':prevcp,'pid':os.getpid()})
   print(json.dumps({'event':'block_complete','directory':str(bdir),'checkpoint_sha256':prevcp,**summary}),flush=True)
   write(attempt_dir/f'{ordinal:04d}_receipt.json',{'pid':os.getpid(),'block':block,'block_checkpoint_sha256':prevcp,
    'block_wall_seconds':time.monotonic()-began,'process_total_wall_limit':600,'stage':'completed_block',
    'no_other_block_dispatched_in_this_process':True})
  except BaseException as exc:
   saveprefix('actual_fault_recoverable');write(work/f'fault_{os.getpid()}.json',{'pid':os.getpid(),'block':block+1,'error':repr(exc),'prefix':str(prefix)});raise
 choices,_=selection()
 if args.reconstruct or not choices or not active['queue']:
  audit();lists=current_lists()
  for k in arrays:state[k]=lists[k]+(newleaves if k=='passed' else [])
  state.update({'task077':'upper_framecover077','status077':'full_frame_zero_kernel_complete_waiting_main_review',
   'fingerprint077':digest({'registration':sha(runreg),'last_checkpoint':prevcp}),
   'conversion077':{'source_state':str(SOURCE),'source_state_sha256':SOURCE_HASH,'source_fingerprint':state['fingerprint'],
    'run_registration_sha256':sha(runreg),'last_block':block,'last_checkpoint_sha256':prevcp,
    'frontier_sha256':frontier_hash(),'tested':len(tested),'accepted':len(newleaves),
    'cache_canonical_sha256':cachehash,'new_kernel_chain_started':False,'math':'075 Q_full=I4/16',
    'closed_root_cover_candidate':not active['queue']}})
  final=work/'final';final.mkdir(exist_ok=True);write(final/'state.json',state)
  assert sha(SOURCE)==SOURCE_HASH
  for p,h in LOCKS.items():assert sha(p)==h
  summary={'task':'upper_framecover077','status':'complete_waiting_main_review','pid':os.getpid(),'blocks':block,
   'tested':len(tested),'accepted':len(newleaves),'failed_uncertified':len(tested)-len(newleaves),
   'passed':len(state['passed']),'split':len(state['splits']),'queue':len(state['queue']),
   'remaining_untested_eligible':len(choices),'closed_root_cover_candidate':not active['queue'],
   'state_path':str(final/'state.json'),'state_sha256':sha(final/'state.json'),
   'source_sha256_after':sha(SOURCE),'last_checkpoint_sha256':prevcp,'cache_canonical_sha256':cachehash,
   'kernel_calls':0,'wall_seconds':time.monotonic()-began,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
   'exit_code':0,'global_g_or_C5_claim':False}
  write(final/'summary.json',summary);print(json.dumps({'event':'final_complete',**summary}),flush=True)
 receipt={'pid':os.getpid(),'exit_code':0,'total_wall_seconds':time.monotonic()-began,'total_wall_limit':600,
  'source_sha256_after':sha(SOURCE),'last_complete_block':block,'checkpoint_sha256':prevcp,
  'new_kernel_calls':0,'one_block_per_process':True}
 assert receipt['total_wall_seconds']<600
 if attempt_dir is not None:write(attempt_dir/f'{ordinal:04d}_process_receipt.json',receipt)
 else:write(work/f'reconstruction_{os.getpid()}_receipt.json',receipt)
 signal.alarm(0);print(json.dumps({'event':'process_exit',**receipt}),flush=True)

if __name__=='__main__':run()

