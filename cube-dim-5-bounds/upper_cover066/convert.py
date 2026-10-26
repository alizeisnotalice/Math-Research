#!/usr/bin/env python3
"""Frozen queue-targeted zero-kernel conversion, chained sparse deltas.

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
import sys
import time

SOURCE=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_cover061/batch0001/state.json')
SOURCE_HASH='d0d62d99d57b3f15d3acb1f4a635b188e196b0de6e155b461460ccbd64a72d31'
LOCKS={
 '/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-I/cube-dim-5-bounds/upper_multibasis056/manifest.json':'0f3261b09dabc9557903547c0670125392e2bf106a8f4e7c3d952005b3978ba4',
 '/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-I/cube-dim-5-bounds/upper_multibasis056/report.md':'7aee0b3a6f685f5e722f2abd95cbf303357ef7ef4d4dfe8659c0b34e487c5303',
 '/Users/zhengzhihao/Desktop/T/Math-Research-cube5-gprobe/cube-dim-5-bounds/audit_multibasis060/report.md':'b2d7797bb487472383797399a75c4ee9856d0a0141f2564dcea28898a41dc081',
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
 sets=[sorted(set([a,F(0),b])) for a,b in zip(lower,upper)];rows=[]
 for v in product(*sets):
  norm=sum(x*x for x in v);q=norm/4-sum(abs(x) for x in v)**2/36
  assert F(5,36)*norm<=q<=F(2,9)*norm
  rows.append((Delta+linear(v)-4096*q,v))
 nm,nv=min(rows);assert nm>=om
 cert={'mechanism':'signed_basis_81','passed':nm>0,'minimum_log_margin':str(nm),
  'worst_offset':list(map(str,nv)),'checked_points':len(rows),
  'quadratic':'psi(v)=sum(v^2)/4-sum(abs(v))^2/36',
  'Delta_lower':c['Delta_lower'],'loggradient_directed':c['loggradient_directed'],'anchor':c['anchor']}
 return {'source_split_index':index,'node':split['node'],'source_failed_certificate_sha256':digest(c),
  'box_lower':list(map(str,lo)),'box_upper':list(map(str,hi)),'anchor_inside':True,
  'old_minimum_log_margin':str(om),'certificate':cert}

def run():
 ap=argparse.ArgumentParser();ap.add_argument('--work-dir',type=Path,required=True)
 ap.add_argument('--resume',action='store_true');ap.add_argument('--reconstruct',action='store_true')
 args=ap.parse_args();work=args.work_dir.resolve();work.mkdir(parents=True,exist_ok=True)
 began=time.monotonic();codehash=sha(__file__);runreg=work/'registration.json'
 assert sha(SOURCE)==SOURCE_HASH
 for p,h in LOCKS.items():assert sha(p)==h
 if not args.resume and not args.reconstruct:
  assert not runreg.exists()
  write(runreg,{'task':'upper_cover066','source':str(SOURCE),'source_sha256':SOURCE_HASH,
   'code':str(Path(__file__).resolve()),'code_sha256':codehash,'locks':LOCKS,
   'registered_utc':datetime.now(timezone.utc).isoformat(),'epsilon':'1/64','target':'1/10000',
   'resources':{'processes':1,'threads':1,'MATLAB':0,'new_kernel_calls':0},
   'block_limits':{'selected_eligible_source_parents':2048,'tested_parents':2048,'qualification_seconds':120},
   'selection':'only live untested source splits with a remaining queue descendant; (path depth,-initial queue descendant count,root index,path); exact candidate prefix registered before qualification',
   'schema':{'block_registration':'parent checkpoint file SHA and logical frontier SHA, selected source split indices+nodes+queue counts',
    'certificates':'one full 81-point margin certificate per tested source parent; passed false remains uncertified, never an H counterexample',
    'retired':'per-block removed source passed/queue/split indices, node and unique accepted ancestor',
    'checkpoint':'chained registration/certificates/retirement hashes; logical frontier hash and cumulative tested count',
    'prefix_checkpoint':'recoverable partial rows and cursor under frozen block registration',
    'final_state':'all source fields preserved except passed/queue/splits; add task066/status066/fingerprint066/conversion066'},
   'reconstruct':f'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 {sys.executable} -B {Path(__file__).resolve()} --work-dir {work} --reconstruct',
   'resume':f'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 {sys.executable} -B {Path(__file__).resolve()} --work-dir {work} --resume',
   'new_kernel_chain_authorized':False,'pause_file':str(work/'STOP')})
 reg=json.load(open(runreg));assert reg['code_sha256']==codehash and reg['source_sha256']==SOURCE_HASH
 state=json.load(open(SOURCE));assert not state['inflight'] and not state['resume_jobs']
 assert state['epsilon']=='1/64' and state['target']=='1/10000'
 arrays={k:state[k] for k in ('passed','queue','splits')}
 active={k:set(range(len(v))) for k,v in arrays.items()};newleaves=[];tested=set();allaccepted=defaultdict(list)
 cachehash=digest(state['cache']);rootshash=digest(state['roots'])
 partition(state['roots'],arrays['passed'],arrays['queue'],arrays['splits'])
 qcounts=defaultdict(int)
 for r,p in arrays['queue']:
  for depth in range(len(p)):qcounts[(r,p[:depth])]+=1
 ordering=sorted((i for i,s in enumerate(arrays['splits']) if qcounts[key(s['node'])]>0),
  key=lambda i:(len(arrays['splits'][i]['node'][1]),-qcounts[key(arrays['splits'][i]['node'])],*arrays['splits'][i]['node']))
 def frontier_hash():
  return digest({'source_sha256':SOURCE_HASH,'active_source_indices':{k:sorted(v) for k,v in active.items()},
   'accepted_certificates':newleaves,'tested_source_split_indices':sorted(tested)})
 def current_lists():
  return {k:[v[i] for i in sorted(active[k])] for k,v in arrays.items()}
 def apply_rows(rows,expected_retired=None):
  accepted=defaultdict(list)
  for row in rows:
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
 # Replay immutable block deltas, validating the chain and every live frontier.
 while (work/f'batch{block+1:04d}'/'checkpoint.json').exists():
  bdir=work/f'batch{block+1:04d}';cp=json.load(open(bdir/'checkpoint.json'))
  assert cp['stage']=='complete' and cp['parent_checkpoint_sha256']==prevcp
  br=json.load(open(bdir/'registration.json'));assert br['parent_frontier_sha256']==frontier_hash()
  assert br['parent_checkpoint_sha256']==prevcp and cp['registration_sha256']==sha(bdir/'registration.json')
  assert cp['certificates_sha256']==sha(bdir/'certificates.json') and cp['retired_sha256']==sha(bdir/'retired.json')
  apply_rows(json.load(open(bdir/'certificates.json')),json.load(open(bdir/'retired.json')))
  assert cp['frontier_sha256']==frontier_hash()
  prevcp=sha(bdir/'checkpoint.json');block+=1
 def audit():
  lists=current_lists();a=partition(state['roots'],lists['passed']+newleaves,lists['queue'],lists['splits'])
  assert digest(state['cache'])==cachehash and digest(state['roots'])==rootshash
  return a
 if not args.reconstruct:
  while True:
   if (work/'STOP').exists():
    write(work/'stopped.json',{'status':'stopped_before_new_dispatch','pid':os.getpid(),'last_complete_block':block,'checkpoint_sha256':prevcp});return
   selected=[i for i in ordering if i in active['splits'] and i not in tested][:2048]
   if not selected or not active['queue']:break
   bdir=work/f'batch{block+1:04d}';bdir.mkdir(exist_ok=True);brpath=bdir/'registration.json'
   br={'task':'upper_cover066','block':block+1,'source_state':str(SOURCE),'source_state_sha256':SOURCE_HASH,
    'run_registration_sha256':sha(runreg),'code_sha256':codehash,'parent_checkpoint_sha256':prevcp,
    'parent_frontier_sha256':frontier_hash(),'max_tested':2048,'qualification_seconds':120,
    'selected':[{'source_split_index':i,'node':arrays['splits'][i]['node'],'initial_queue_descendants':qcounts[key(arrays['splits'][i]['node'])]} for i in selected]}
   if brpath.exists():assert json.load(open(brpath))==br
   else:write(brpath,br)
   rows=[];cursor=0;accepted=defaultdict(list);pfx=bdir/'prefix_checkpoint.json'
   if pfx.exists():
    prefix=json.load(open(pfx));assert prefix['registration_sha256']==sha(brpath)
    rows=prefix['rows'];cursor=prefix['cursor']
    for row in rows:
     if row['certificate']['passed']:accepted[row['node'][0]].append(row['node'][1])
   blockstart=time.monotonic()
   def saveprefix(stage):write(pfx,{'stage':stage,'registration_sha256':sha(brpath),'pid':os.getpid(),'cursor':cursor,'rows':rows,'seconds':time.monotonic()-blockstart})
   saveprefix('ready')
   try:
    while cursor<len(selected) and len(rows)<2048 and time.monotonic()-blockstart<120:
     if (work/'STOP').exists():saveprefix('stopped');return
     idx=selected[cursor];cursor+=1;n=arrays['splits'][idx]['node']
     if any(n[1].startswith(p) for p in accepted[n[0]]):continue
     row=qualify(arrays['splits'][idx],state,idx);rows.append(row)
     if row['certificate']['passed']:accepted[n[0]].append(n[1])
     if len(rows)%16==0:saveprefix('qualification_prefix')
    qualification_seconds=time.monotonic()-blockstart;saveprefix('qualification_complete')
    retired=apply_rows(rows);write(bdir/'certificates.json',rows);write(bdir/'retired.json',retired)
    a=audit();fh=frontier_hash()
    summary={'task':'upper_cover066','block':block+1,'pid':os.getpid(),'tested':len(rows),
     'accepted':sum(r['certificate']['passed'] for r in rows),'failed_uncertified':sum(not r['certificate']['passed'] for r in rows),
     'points':sum(r['certificate']['checked_points'] for r in rows),'selected':len(selected),'cursor':cursor,
     'retired':{k:len(v) for k,v in retired.items()},'passed':len(active['passed'])+len(newleaves),
     'split':len(active['splits']),'queue':len(active['queue']),'qualification_seconds':qualification_seconds,
     'cumulative_tested':len(tested),'cumulative_accepted':len(newleaves),'cache_entries':len(state['cache']),
     'cache_canonical_sha256':cachehash,'frontier_sha256':fh,'kernel_calls':0,'audit':a,
     'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    write(bdir/'summary.json',summary)
    cp={'stage':'complete','block':block+1,'pid':os.getpid(),'parent_checkpoint_sha256':prevcp,
     'registration_sha256':sha(brpath),'certificates_sha256':sha(bdir/'certificates.json'),
     'retired_sha256':sha(bdir/'retired.json'),'summary_sha256':sha(bdir/'summary.json'),
     'frontier_sha256':fh,'cumulative_tested':len(tested),'cumulative_accepted':len(newleaves)}
    write(bdir/'checkpoint.json',cp);prevcp=sha(bdir/'checkpoint.json');block+=1
    write(work/'CURRENT.json',{'block':block,'checkpoint_path':str(bdir/'checkpoint.json'),'checkpoint_sha256':prevcp,'pid':os.getpid()})
    print(json.dumps({'event':'block_complete','directory':str(bdir),'checkpoint_sha256':prevcp,**summary}),flush=True)
   except BaseException as exc:
    saveprefix('actual_fault_recoverable');write(work/'fault.json',{'pid':os.getpid(),'block':block+1,'error':repr(exc),'checkpoint':str(pfx)});raise
 # A single final materialization, with all unchanged cache and metadata retained.
 audit();lists=current_lists()
 for k in arrays:state[k]=lists[k]+(newleaves if k=='passed' else [])
 state.update({'task066':'upper_cover066','status066':'zero_kernel_conversion_complete_waiting_main_review',
  'fingerprint066':digest({'registration':sha(runreg),'last_checkpoint':prevcp}),
  'conversion066':{'source_state':str(SOURCE),'source_state_sha256':SOURCE_HASH,'source_fingerprint':state['fingerprint'],
   'run_registration_sha256':sha(runreg),'last_block':block,'last_checkpoint_sha256':prevcp,
   'frontier_sha256':frontier_hash(),'tested':len(tested),'accepted':len(newleaves),
   'cache_canonical_sha256':cachehash,'new_kernel_chain_started':False}})
 final=work/'final';final.mkdir(exist_ok=True);write(final/'state.json',state)
 assert sha(SOURCE)==SOURCE_HASH
 for p,h in LOCKS.items():assert sha(p)==h
 summary={'task':'upper_cover066','status':'complete_waiting_main_review','pid':os.getpid(),'blocks':block,
  'tested':len(tested),'accepted':len(newleaves),'failed_uncertified':len(tested)-len(newleaves),
  'passed':len(state['passed']),'split':len(state['splits']),'queue':len(state['queue']),
  'remaining_untested_eligible':sum(i in active['splits'] and i not in tested for i in ordering),
  'state_path':str(final/'state.json'),'state_sha256':sha(final/'state.json'),
  'source_sha256_after':sha(SOURCE),'last_checkpoint_sha256':prevcp,'cache_canonical_sha256':cachehash,
  'kernel_calls':0,'wall_seconds':time.monotonic()-began,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
  'exit_code':0,'global_g_or_C5_claim':False}
 write(final/'summary.json',summary);print(json.dumps({'event':'final_complete',**summary}),flush=True)

if __name__=='__main__':run()

