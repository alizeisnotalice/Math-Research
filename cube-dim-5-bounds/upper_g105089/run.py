#!/usr/bin/env python3
"""Finite g105 first batch after exact zero-H pre-split; no automatic continuation.

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
import gzip
import importlib.util
import multiprocessing as mp
from flint import ctx
import sys
import time

SOURCE=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_g105086/final/state.json')
SOURCE_HASH='8647a4543ef6a869e747ec06d93f5137cf11feb7b3a64a002995d25a681c4f9a'
LOCKS={
 '/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-lower/cube-dim-5-bounds/audit_frame075/report.md':'3308d1e5c660cd11d73d47a013acf794bd1f596afdb5a113573a5bfd7481f225',
 '/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-lower/cube-dim-5-bounds/audit_frame075/manifest.json':'1d5e270eb1664d714ec80f0b371420d40028b24c3943196b9ab41a0b3de1f9d1',
 '/Users/zhengzhihao/Desktop/T/Math-Research-cube5-upper/cube-dim-5-bounds/upper_g105086/manifest.json':'10461718679541d040e06241896c6f01da54165355f57ad6441084d0ddd9cbab',
 '/Users/zhengzhihao/Desktop/T/Math-Research-cube5-lower-audit-cont/cube-dim-5-bounds/audit_kernel081/report.md':'76276ff907887dd68d984aaabc7a91677f6c92f607d9e51c3e4d2d749fbd15ce',
 '/Users/zhengzhihao/Desktop/T/Math-Research-cube5-upper/cube-dim-5-bounds/upper_cover046/arb_cover.py':'df3fb4165aea29e4aa3d1fb7ecacdb9cc5ad4a531bc03e3cc43a0820671fe028',
 '/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-upper/cube-dim-5-bounds/astra_vertex017/boundary_triangulation.json.gz':'8dd2e9895d7465ac599d43feb1b5099506f4e43d1248215eba0a77e11f4d13f9',
 '/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_g105086/pending_order089.json':'68f4ef52ed6dd1360d8487266e9a744a82d4298000677092abd67b026bf0c8c8',
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

API_PATH=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-upper/cube-dim-5-bounds/upper_cover046/arb_cover.py')
TABLE_PATH=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-upper/cube-dim-5-bounds/astra_vertex017/boundary_triangulation.json.gz')
ORDER_PATH=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_g105086/pending_order089.json')
API=None;KERNEL=None
def load_api():
 spec=importlib.util.spec_from_file_location('frozen046api',API_PATH);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def worker_init(ready):
 global API,KERNEL
 assert sha(API_PATH)==LOCKS[str(API_PATH)] and sha(TABLE_PATH)==LOCKS[str(TABLE_PATH)]
 API=load_api();KERNEL=API.FirstJetKernel(json.loads(gzip.decompress(TABLE_PATH.read_bytes())),F(1,64))
 ready.put({'pid':os.getpid(),'initialized':True,'api_sha256':sha(API_PATH),'cones':3504,'vertices':280})
def run_job(task):
 began=time.monotonic();point=task['point'];path=Path(task['point_checkpoint']);attempts=[];pi=0;partial=None;calls=0
 def snapshot(stage,**extra):
  write(path,{'stage':stage,'pid':os.getpid(),'point':point,'precision_index':pi,'partial':partial,
   'attempts':attempts,'kernel_segment_calls':calls,'task_sha256':digest(task),**extra})
 def result(kind,**extra):return {'kind':kind,'pid':os.getpid(),'seconds':time.monotonic()-began,'calls':calls,'attempts':attempts,'point_checkpoint':str(path),**extra}
 try:
  assert KERNEL.eligible(point)['eligible'];snapshot('point_registered_before_kernel')
  while pi<3:
   bits=(256,384,512)[pi]
   # A segment deadline only writes a resumable cone prefix; it never stops
   # this finite job or omits a cone. There is no program wall deadline.
   got=KERNEL.evaluate(point,bits,partial,time.monotonic()+1.0,128);calls+=1
   if not got['complete']:
    partial=got['partial'];snapshot('cone_prefix')
    if Path(task['stop_file']).exists():return result('partial',partial=partial,precision_index=pi)
    continue
   values=got['values'];packed=list(map(API.pack,values));quality=API.numeric_quality(values)
   attempt={k:got[k] for k in ('bits','seconds','cancellation_ratio','minimum_integer_projection_gap','common_denominator','all_cones')}
   attempt['numeric_quality']=quality;attempts.append(attempt)
   counter=values[0].is_finite() and API.endpoint(values[0],True)<F(21,200000)
   separated=values[0].is_finite() and API.endpoint(values[0])>F(21,200000)
   if counter or (quality and separated):
    record={'point':point,'values':packed,'bits':bits,'candidate':task['candidate'],'geometry_tests':task['geometry_tests'],'attempts':attempts}
    snapshot('counterexample' if counter else 'complete',record=record)
    return result('counterexample' if counter else 'complete',record=record,peak_rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
   attempts[-1]['values']=packed;pi+=1;partial=None;snapshot('numeric_precision_retry')
  snapshot('numerical_diagnostic');return result('numerical_diagnostic',precision_index=pi,partial=partial)
 except BaseException as exc:
  snapshot('actual_worker_error',exception=repr(exc));return result('worker_error',exception=repr(exc),precision_index=pi,partial=partial)

def certificate(api,box,record):
 ctx.prec=record['bits'];values=[tuple(map(F,v)) for v in record['values']]
 Hlo,Hhi=values[0];assert F(21,200000)<Hlo<=Hhi
 gradient=[tuple(map(F,api.pack(api.unpack(v)/api.unpack(record['values'][0])))) for v in record['values'][1:]]
 for (a,b),(gl,gu) in zip(values[1:],gradient):
  ratios=[a/Hlo,a/Hhi,b/Hlo,b/Hhi];assert gl<=min(ratios)<=max(ratios)<=gu
 Delta=api.endpoint((api.ab(Hlo)/api.ab(F(21,200000))).log());assert Delta>0
 point=list(map(F,record['point']));assert point[0]==1;lo,hi=box
 assert all(a<=p<=b for a,p,b in zip(lo,point[1:],hi))
 rows=[]
 for bits in product((0,1),repeat=4):
  v=[(hi[j] if bits[j] else lo[j])-point[j+1] for j in range(4)]
  rows.append((Delta+sum(min(gl*x,gu*x) for x,(gl,gu) in zip(v,gradient))-256*sum(x*x for x in v),bits))
 margin,bits=min(rows)
 return {'mechanism':'full_frame_16','quadratic':'Q_full=I4/16','checked_points':16,
  'passed':margin>0,'minimum_log_margin':str(margin),'worst_corner_bits':list(bits),
  'Delta_lower':str(Delta),'loggradient_directed':[[str(a),str(b)] for a,b in gradient],'anchor':api.pointkey(point)}

def apply_event(overlay,event):
 body=event['body'];kind=event['kind'];inflight=overlay['inflight'];c=overlay['counts']
 if kind=='dispatch':
  assert overlay['queue'].pop(0)==body['node'];assert str(body['ordinal']) not in inflight
  inflight[str(body['ordinal'])]=body;c['new_anchor_starts']+=body['point'] is not None and not body['cached']
 elif kind=='received':
  ent=inflight[str(body['ordinal'])];assert 'returned' not in ent;ent['returned']=body['result']
  c['kernel_segment_calls']+=body['result'].get('calls',0)
 elif kind=='merge':
  ent=inflight[str(body['ordinal'])];got=ent['returned'];node=ent['node']
  assert body['ordinal']==overlay['next_merge'];overlay['next_merge']+=1
  if got['kind']=='complete':
   record=got['record'];anchor=','.join(record['point'])
   if ent['cached']:c['cache_hits']+=1
   else:assert anchor not in overlay['new_cache'];overlay['new_cache'][anchor]=record;c['new_complete']+=1
   cert=body['certificate'];overlay['hardness'].pop(str(node[0])+':'+node[1],None)
   if cert['passed']:overlay['new_passed'].append({'node':node,'certificate':cert})
   else:
    children=[[node[0],node[1]+'0'],[node[0],node[1]+'1']]
    overlay['new_splits'].append({'node':node,'axis':len(node[1])%4,'children':children,'failed_certificate':cert})
    overlay['queue'].extend(children)
    for child in children:overlay['hardness'][str(child[0])+':'+child[1]]=[record['values'][0][0],cert['minimum_log_margin']]
   del inflight[str(body['ordinal'])]
  elif got['kind']=='counterexample':
   overlay['counterexample']={'point':ent['point'],'H_interval':got['record']['values'][0],'node':node,
    'target':'21/200000','record':got['record'],'scope':'H_epsilon target only; no fullG or C5 counterexample claim'}
  else:overlay['faults'].append({'node':node,'point':ent['point'],'returned':got})
 elif kind=='finish':
  overlay['status']=body['status'];overlay['pool_evidence']=body['pool_evidence']
  overlay['queue'].sort(key=lambda n:(*map(F,overlay['hardness'][str(n[0])+':'+n[1]]),*n))
 else:raise ValueError(kind)
 overlay['journal_seq']=event['seq'];overlay['journal_hash']=event['hash']

def run():
 began=time.monotonic();ap=argparse.ArgumentParser();ap.add_argument('--work-dir',type=Path,required=True);ap.add_argument('--recover-only',type=Path)
 ap.add_argument('--authorize-next-batch',action='store_true')
 args=ap.parse_args();work=args.work_dir.resolve();work.mkdir(parents=True,exist_ok=True);codehash=sha(__file__)
 assert sha(SOURCE)==SOURCE_HASH
 for p,h in LOCKS.items():assert sha(p)==h
 api=load_api();source=json.load(open(SOURCE));assert [len(source[k]) for k in ('passed','splits','queue','cache')]==[96104,45600,1856,180040]
 assert source['epsilon']==source['cache_model_epsilon']=='1/64' and source['target']=='21/200000'
 assert not source['inflight'] and not source['resume_jobs'];basecachehash=digest(source['cache'])
 source_other_hashes={k:digest(v) for k,v in source.items() if k not in ('passed','splits','queue','cache','roots')}
 runregpath=work/'registration.json'
 if not runregpath.exists():
  write(runregpath,{'task':'upper_g105089','source':str(SOURCE),'source_sha256':SOURCE_HASH,'source_counts':{'passed':96104,'split':45600,'queue':1856,'cache':180040},
   'code_path':str(Path(__file__).resolve()),'code_sha256':codehash,'locks':LOCKS,'epsilon':'1/64','target':'21/200000',
   'model':'alpha1;16 unnormalised body Cauchy directions;inverse Fourier(2pi)^-5','mathematics':'075 Q_full=I4/16,16 closed corners',
   'python':sys.executable,'workers':12,'max_new_anchors_per_finite_batch':4096,'max_selected_nodes':4096,'program_wall_deadline':None,'alarm':False,
   'anchor_rule':'frozen046 midpoint then three eta=minside/4096 generic perturbations beta=1/8,1/16,1/32; exact nonzero/distinct280 integer projections',
   'precision':'256,384,512 only for nonfinite/numeric width/H target straddle; never for geometric margin failure',
   'kernel':'locked046 FirstJetKernel all3504 cones, all5 components; cone-prefix snapshot each cooperative1s segment, stride128; no cones omitted',
   'ordering':'frozen086 failed pending certificates -> exact zero-H pre-split1856; first finite selected3712 closedchildren; subsequent children sorted low cached parent H_lower,low failed fullQ margin,root,path',
   'stops':['queue empty','finite selected prefix exhausted','strict H_upper<target','actual numerical/process/program/state fault','explicit STOP'],
   'schema':{'overlay':'source references plus cumulative new_cache/new_passed/new_splits/current queue+hardness/inflight/counterexample/faults/counters',
    'journal':'hash-chained dispatch/received/merge/finish events; replay from initial_overlay without kernel',
    'point_checkpoint':'point,precision_index,exact partial cone prefix/attempts before every segment/retry; retain on fault',
    'pre_split':'consume all1856 already-qualified new-target failed closed boxes verbatim;split longest side into3712 children;no H/log/gradient or same-anchor requalification',
    'final_state':'old source fields retained; append new cache/passed/splits,replace queue/inflight/resume_jobs; add task089/status089/fingerprint089/conversion089'},
   'first_selected_children':3712,'expected_initial_cache_hits':1124,'expected_initial_new_anchor_starts':2588,
   'later_batch_authorization':'main must review first-batch residual before --authorize-next-batch',
   'final_claim':'closed cover candidate awaiting main homogeneity/global qualification, no C5 or integral claim',
   'registered_utc':datetime.now(timezone.utc).isoformat()})
 runreg=json.load(open(runregpath));assert runreg['code_sha256']==codehash and runreg['source_sha256']==SOURCE_HASH
 currentpath=work/'CURRENT.json';parenthash=None;parentoverlay=None
 if currentpath.exists():
  assert args.recover_only or args.authorize_next_batch,'Main review required before another production batch'
  cur=json.load(open(currentpath));assert sha(cur['state_path'])==cur['state_sha256'];parentoverlay=Path(cur['state_path']);parenthash=cur['state_sha256'];ov=json.load(open(parentoverlay))
  assert ov['status']=='finite_prefix_complete' and not ov['inflight'] and not ov['faults'] and ov['counterexample'] is None
  number=cur['batch']+1
 else:
  records=json.load(open(ORDER_PATH))['records'];assert {key(r['node']) for r in records}=={key(n) for n in source['queue']}
  certificates={key(r['node']):r for r in source['pending_certificates086']};assert len(certificates)==1856
  presplitreg={'task':'upper_g105089_zero_H_presplit','source_sha256':SOURCE_HASH,'code_sha256':codehash,'target':'21/200000','epsilon':'1/64',
   'source_pending_certificates_sha256':digest(source['pending_certificates086']),'pending_order_sha256':sha(ORDER_PATH),
   'failed_parents':1856,'closed_children':3712,'method':'consume exact false fullQ certificates,split longest side;tied first coordinate','new_H_gradient_kernel_log_calls':0}
  if not (work/'pre_split_registration.json').exists():write(work/'pre_split_registration.json',presplitreg)
  else:assert json.load(open(work/'pre_split_registration.json'))==presplitreg
  queue=[];splits=[];hardness={};pre_map=[]
  for record in records:
   node=record['node'];entry=certificates[key(node)];cert=entry['certificate']
   assert cert['target']=='21/200000' and cert['epsilon']=='1/64' and cert['passed'] is False and F(cert['minimum_log_margin'])<=0
   children=[[node[0],node[1]+'0'],[node[0],node[1]+'1']]
   splits.append({'node':node,'axis':len(node[1])%4,'children':children,'failed_certificate':cert,'source086_pending_certificate':True})
   queue.extend(children)
   for child in children:hardness[str(child[0])+':'+child[1]]=[record['parent_cached_H_lower'],record['parent_full_frame_minimum_log_margin']]
   pre_map.append({'source086_pending_node':node,'source070_leaf_index':entry['source070_leaf_index'],'failed_certificate_sha256':digest(cert),'children':children})
  assert len(queue)==3712 and len({key(n) for n in queue})==3712 and len(splits)==1856
  queue.sort(key=lambda n:(*map(F,hardness[str(n[0])+':'+n[1]]),*n))
  write(work/'pre_split_mapper.json',pre_map)
  ov={'source_sha256':SOURCE_HASH,'new_cache':{},'new_passed':[],'new_splits':splits,'queue':queue,'hardness':hardness,
   'inflight':{},'resume_jobs':{},'counterexample':None,'faults':[],'status':'initialized089_after_zero_H_presplit','pre_split_count':1856,
   'counts':{'new_anchor_starts':0,'new_complete':0,'kernel_segment_calls':0,'cache_hits':0},'batch_number':0}
  number=1
 def sort_queue():ov['queue'].sort(key=lambda n:(*map(F,ov['hardness'][str(n[0])+':'+n[1]]),*n))
 def audit():
  pending=ov['queue']+[e['node'] for e in ov['inflight'].values()]
  return partition(source['roots'],source['passed']+ov['new_passed'],pending,source['splits']+ov['new_splits'])
 if number==1 and not args.recover_only:
  tree=audit();assert len(ov['queue'])==3712 and len(ov['new_splits'])==1856
  write(work/'pre_split_validation.json',{'status':'passed','source_sha256':SOURCE_HASH,'tree':tree,'accepted_reused':96104,
   'pre_split_parents':1856,'children':3712,'source_cache_preserved':True,'source_cache_canonical_sha256':basecachehash,
   'pre_split_mapper_sha256':sha(work/'pre_split_mapper.json'),'new_H_physical_gradient_kernel_log_calls':0})
 if args.recover_only:
  bdir=args.recover_only;ov=json.load(open(bdir/'initial_overlay.json'));previous='0'*64;seq=0
  for line in open(bdir/'journal.jsonl'):
   e=json.loads(line);h=e.pop('hash');assert e['previous']==previous and e['seq']==seq+1 and digest(e)==h;e['hash']=h;apply_event(ov,e);previous=h;seq+=1
  audit();write(bdir/'recovered_overlay.json',ov);print(json.dumps({'event':'recovery_only','events':seq,'path':str(bdir/'recovered_overlay.json'),'sha256':sha(bdir/'recovered_overlay.json'),'kernel_calls':0}));return
 assert ov['queue'],'queue already empty; no new batch'
 if number>1:sort_queue()
 audit();selected=ov['queue'][:4096];bdir=work/f'batch{number:04d}';bdir.mkdir(exist_ok=True);assert not (bdir/'registration.json').exists()
 parentcounts={k:len(ov[k]) for k in ('new_cache','new_passed','new_splits')};beforecounters=dict(ov['counts'])
 br={'task':'upper_g105089','batch':number,'master_pid':os.getpid(),'registered_utc':datetime.now(timezone.utc).isoformat(),
  'source_state':str(SOURCE),'source_state_sha256':SOURCE_HASH,'run_registration_sha256':sha(runregpath),'code_sha256':codehash,
  'parent_overlay':str(parentoverlay) if parentoverlay else None,'parent_overlay_sha256':parenthash,'base_cache_canonical_sha256':basecachehash,
  'selected_nodes':selected,'selected_sha256':digest(selected),'max_new_anchor_starts':4096,'workers':12,'wall_deadline':None,'locks':LOCKS,
  'pre_split_count':ov['pre_split_count'],'pre_split_registration_sha256':sha(work/'pre_split_registration.json'),'pre_split_mapper_sha256':sha(work/'pre_split_mapper.json')}
 write(bdir/'registration.json',br);ov.update({'batch_number':number,'journal_seq':0,'journal_hash':'0'*64,'next_merge':0,'status':'running_finite_batch'})
 write(bdir/'initial_overlay.json',ov);initialhash=sha(bdir/'initial_overlay.json');journal=open(bdir/'journal.jsonl','ab');events=0
 def emit(kind,body):
  nonlocal events
  e={'seq':ov['journal_seq']+1,'previous':ov['journal_hash'],'kind':kind,'body':body};e['hash']=digest(e)
  journal.write(json.dumps(e,sort_keys=True,separators=(',',':')).encode()+b'\n');journal.flush();os.fsync(journal.fileno());apply_event(ov,e);events+=1
 def checkpoint():audit();write(bdir/'state.json',ov)
 checkpoint();kernel=api.FirstJetKernel(json.loads(gzip.decompress(TABLE_PATH.read_bytes())),F(1,64))
 context=mp.get_context('spawn');ready=context.Queue();pool=context.Pool(12,worker_init,(ready,));processes=list(pool._pool)
 worker_readies=[ready.get(timeout=60) for _ in range(12)];pids=sorted(r['pid'] for r in worker_readies)
 write(bdir/'pool_registration.json',{'master_pid':os.getpid(),'worker_pids':pids,'workers':12,'registered_before_dispatch':True,'readies':worker_readies})
 print(json.dumps({'event':'pool_ready','master_pid':os.getpid(),'worker_pids':pids,'batch':number,'selected':len(selected),'directory':str(bdir)}),flush=True)
 pending={};dispatch=0;stop=None;terminated=False;seenkeys=set();received=[]
 try:
  while True:
   if stop is None and (work/'STOP').exists():stop='explicit_STOP'
   while stop is None and len(pending)<12 and dispatch<len(selected):
    node=selected[dispatch];point,candidate,tests=api.choose_anchor(api.box_of(node,source['roots']),kernel)
    point=None if point is None else list(map(str,point));anchor=None if point is None else api.pointkey(point)
    cached=anchor in source['cache'] or anchor in ov['new_cache'] if anchor else False
    task={'ordinal':dispatch,'node':node,'point':point,'candidate':candidate,'geometry_tests':tests,'cached':cached,
     'point_checkpoint':str(bdir/f'point{dispatch:04d}.json'),'stop_file':str(work/'STOP')}
    assert anchor is None or anchor not in seenkeys
    if anchor:seenkeys.add(anchor)
    emit('dispatch',task)
    if point is None:emit('received',{'ordinal':dispatch,'result':{'kind':'no_generic_anchor','calls':0}});stop='no_generic_anchor'
    elif cached:
     record=source['cache'][anchor] if anchor in source['cache'] else ov['new_cache'][anchor]
     emit('received',{'ordinal':dispatch,'result':{'kind':'complete','record':record,'calls':0,'cache_reused':True}})
    else:pending[dispatch]=pool.apply_async(run_job,(task,))
    dispatch+=1
   for ordinal,handle in list(pending.items()):
    if not handle.ready():continue
    try:got=handle.get()
    except BaseException as exc:got={'kind':'worker_error','exception':repr(exc),'calls':0}
    emit('received',{'ordinal':ordinal,'result':got});received.append(got);del pending[ordinal]
    if got['kind']!='complete' and stop is None:stop=got['kind']
   if any(p.exitcode is not None for p in processes):
    stop='worker_process_fault';pool.terminate();terminated=True
    for ordinal in list(pending):
     emit('received',{'ordinal':ordinal,'result':{'kind':'worker_process_fault','calls':0,'diagnostic':'pool stopped on actual process death; own point checkpoint retained'}});del pending[ordinal]
   while str(ov['next_merge']) in ov['inflight'] and 'returned' in ov['inflight'][str(ov['next_merge'])]:
    ordinal=ov['next_merge'];entry=ov['inflight'][str(ordinal)];body={'ordinal':ordinal}
    if entry['returned']['kind']=='complete':body['certificate']=certificate(api,api.box_of(entry['node'],source['roots']),entry['returned']['record'])
    emit('merge',body)
    if ov['next_merge']%64==0:
     checkpoint();print(json.dumps({'event':'progress','batch':number,'merged':ov['next_merge'],'new_passed':len(ov['new_passed'])-parentcounts['new_passed'],
      'new_split':len(ov['new_splits'])-parentcounts['new_splits'],'queue':len(ov['queue']),'inflight':len(ov['inflight']),'seconds':time.monotonic()-began}),flush=True)
   if dispatch==len(selected) and stop is None:stop='finite_prefix_complete'
   if not pending and stop is not None:break
   time.sleep(.002)
 finally:
  if not terminated:pool.close()
  pool.join();ready.close();ready.join_thread()
 exits=[{'pid':p.pid,'exitcode':p.exitcode,'is_alive':p.is_alive()} for p in processes]
 assert all(not e['is_alive'] for e in exits)
 if not terminated:assert all(e['exitcode']==0 for e in exits)
 evidence={'master_pid':os.getpid(),'worker_pids':pids,'close_called':not terminated,'join_returned':True,'terminated_on_actual_fault':terminated,'exits':exits}
 emit('finish',{'status':stop,'pool_evidence':evidence});sort_queue();checkpoint();journal.close()
 assert digest(source['cache'])==basecachehash and all(k not in source['cache'] for k in ov['new_cache'])
 assert {k:digest(source[k]) for k in source_other_hashes}==source_other_hashes
 assert sha(SOURCE)==SOURCE_HASH
 for p,h in LOCKS.items():assert sha(p)==h
 validation={'status':'passed','tree':audit(),'base_cache_preserved':True,'base_cache_canonical_sha256':basecachehash,
  'source_other_fields_preserved':True,'new_cache_entries':len(ov['new_cache']),'inflight':len(ov['inflight']),
  'journal_sha256':sha(bdir/'journal.jsonl'),'state_sha256':sha(bdir/'state.json'),'density_calls_in_validation':0}
 write(bdir/'validation.json',validation)
 summary={'task':'upper_g105089','batch':number,'status':stop,'master_pid':os.getpid(),'selected':len(selected),'dispatched':dispatch,
  'new_starts':ov['counts']['new_anchor_starts']-beforecounters['new_anchor_starts'],'new_complete':ov['counts']['new_complete']-beforecounters['new_complete'],
  'new_passed':len(ov['new_passed'])-parentcounts['new_passed'],'new_split':len(ov['new_splits'])-parentcounts['new_splits'],
  'zero_H_pre_split':ov['pre_split_count'],'actual_cache_hits':ov['counts']['cache_hits']-beforecounters['cache_hits'],
  'passed':len(source['passed'])+len(ov['new_passed']),'split':len(source['splits'])+len(ov['new_splits']),'queue':len(ov['queue']),'inflight':len(ov['inflight']),
  'counts':ov['counts'],'counterexample':ov['counterexample'],'faults':len(ov['faults']),'pool':evidence,
  'state_path':str(bdir/'state.json'),'state_sha256':sha(bdir/'state.json'),'journal_sha256':sha(bdir/'journal.jsonl'),
  'wall_seconds':time.monotonic()-began,'main_peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
  'worker_outcomes':dict((kind,sum(r['kind']==kind for r in received)) for kind in set(r['kind'] for r in received)),'global_g_or_C5_claim':False}
 write(bdir/'summary.json',summary)
 write(work/'CURRENT.json',{'batch':number,'state_path':str(bdir/'state.json'),'state_sha256':sha(bdir/'state.json'),'summary_sha256':sha(bdir/'summary.json'),'status':stop})
 if not ov['queue'] and not ov['inflight'] and not ov['counterexample'] and not ov['faults']:
  for k,v in ov['new_cache'].items():source['cache'][k]=v
  source['passed']+=ov['new_passed'];source['splits']+=ov['new_splits'];source['queue']=[]
  source.update({'task089':'upper_g105089','status089':'closed_root_cover_candidate_waiting_main_review',
   'fingerprint089':digest({'run_registration':sha(runregpath),'last_overlay':sha(bdir/'state.json')}),
   'conversion089':{'source_state':str(SOURCE),'source_state_sha256':SOURCE_HASH,'last_overlay':str(bdir/'state.json'),'last_overlay_sha256':sha(bdir/'state.json'),
    'new_counts':ov['counts'],'new_passed':len(ov['new_passed']),'new_splits':len(ov['new_splits']),'math':'075 Q_full=I4/16',
    'closed_root_cover_candidate':True,'global_g_or_C5_claim':False}})
  final=work/'final';final.mkdir(exist_ok=True);write(final/'state.json',source)
  write(final/'summary.json',{**summary,'status':'closed_root_cover_candidate_waiting_main_review','full_state_path':str(final/'state.json'),'full_state_sha256':sha(final/'state.json')})
  print(json.dumps({'event':'closed_cover_candidate','full_state':str(final/'state.json'),'sha256':sha(final/'state.json')}),flush=True)
 print(json.dumps({'event':'batch_complete',**summary}),flush=True)
if __name__=='__main__':run()
