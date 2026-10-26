#!/usr/bin/env python3
"""No-H verification of 070 journals, all new exact margins and materialization.

This is author verification; external independent review remains separate.
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from collections import defaultdict,Counter
import gzip,hashlib,importlib.util,json,os,resource,time
from flint import ctx

WORK=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_cover070')
ENC=json.JSONEncoder(sort_keys=True,separators=(',',':'))
def digest(obj,newline=False):
 h=hashlib.sha256()
 for chunk in ENC.iterencode(obj):h.update(chunk.encode())
 if newline:h.update(b'\n')
 return h.hexdigest()
def sha(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  while data:=f.read(8*1024*1024):h.update(data)
 return h.hexdigest()
def dump(path,obj):path.write_text(json.dumps(obj,sort_keys=True,indent=2)+'\n')
def physical(node,roots):
 r,p=node;lo=[];hi=[]
 for j,k in enumerate(roots[r]):
  bits=p[j::4];d=len(bits);b=int(bits,2) if bits else 0
  a=F(k*(1<<d)+b,32*(1<<d));lo.append(a);hi.append(a+F(1,32*(1<<d)))
 return lo,hi
def verify():
 began=time.monotonic();reg=json.load(open(WORK/'registration.json'));sourcepath=Path(reg['source'])
 assert sha(sourcepath)==reg['source_sha256']
 for p,h in reg['locks'].items():assert sha(p)==h
 spec=importlib.util.spec_from_file_location('verify_locked046',next(p for p in reg['locks'] if p.endswith('/arb_cover.py')));api=importlib.util.module_from_spec(spec);spec.loader.exec_module(api)
 tablepath=next(p for p in reg['locks'] if p.endswith('boundary_triangulation.json.gz'))
 kernel=api.FirstJetKernel(json.loads(gzip.decompress(Path(tablepath).read_bytes())),F(1,64))
 source=json.load(open(sourcepath));current=json.load(open(WORK/'CURRENT.json'));ov=json.load(open(current['state_path']))
 assert sha(current['state_path'])==current['state_sha256'] and not ov['queue'] and not ov['inflight'] and not ov['faults'] and ov['counterexample'] is None
 assert len(ov['new_cache'])==1318 and len(ov['new_passed'])==1001 and len(ov['new_splits'])==317
 journal_events=0;returns={};tasks={};allmerges={};pids=[];batchrows=[];previous_overlay=None
 for b in sorted(WORK.glob('batch[0-9][0-9][0-9][0-9]')):
  br=json.load(open(b/'registration.json'));initial=json.load(open(b/'initial_overlay.json'));final=json.load(open(b/'state.json'));summary=json.load(open(b/'summary.json'))
  assert br['source_state_sha256']==reg['source_sha256'] and br['code_sha256']==reg['code_sha256']
  assert br['selected_nodes']==initial['queue'][:4096] and digest(br['selected_nodes'])==br['selected_sha256']
  if previous_overlay is not None:
   assert sha(br['parent_overlay'])==br['parent_overlay_sha256']
   expected={**previous_overlay,'batch_number':br['batch'],'journal_seq':0,'journal_hash':'0'*64,'next_merge':0,'status':'running_finite_batch'}
   assert initial==expected
  else:
   order=json.load(open(next(p for p in reg['locks'] if p.endswith('/remaining_queue.json'))))['records']
   assert initial['queue']==[r['node'] for r in order]
  previous='0'*64;seq=0;localtasks={};localreturns={};localmerges={}
  for line in open(b/'journal.jsonl'):
   event=json.loads(line);h=event.pop('hash');assert event['previous']==previous and event['seq']==seq+1 and digest(event)==h
   previous=h;seq+=1;kind=event['kind'];body=event['body']
   if kind=='dispatch':
    ordinal=body['ordinal'];assert ordinal==len(localtasks);assert body['node']==br['selected_nodes'][ordinal]
    point,candidate,tests=api.choose_anchor(physical(body['node'],source['roots']),kernel)
    assert list(map(str,point))==body['point'] and candidate==body['candidate'] and tests==body['geometry_tests']
    localtasks[ordinal]=body
   elif kind=='received':localreturns[body['ordinal']]=body['result']
   elif kind=='merge':localmerges[body['ordinal']]=body
   elif kind=='finish':assert body['status']=='finite_prefix_complete'
   else:raise AssertionError(kind)
  assert seq==3*len(localtasks)+1 and previous==final['journal_hash'] and seq==final['journal_seq']
  assert set(localtasks)==set(localreturns)==set(localmerges)
  for i,task in localtasks.items():
   got=localreturns[i];assert got['kind']=='complete' and got['calls']==1 and len(got['attempts'])==1
   pointcp=json.load(open(task['point_checkpoint']));assert pointcp['stage']=='complete' and pointcp['task_sha256']==digest(task)
   assert pointcp['record']==got['record'] and got['record']==ov['new_cache'][','.join(task['point'])]
   assert got['record']['attempts'][0]['all_cones']==3504 and got['record']['bits']==256
   n=tuple(task['node']);assert n not in tasks;tasks[n]=task;returns[n]=got;allmerges[n]=localmerges[i]
  assert summary['pool']['close_called'] and summary['pool']['join_returned'] and not summary['pool']['terminated_on_actual_fault']
  assert len(summary['pool']['exits'])==12 and all(x['exitcode']==0 and not x['is_alive'] for x in summary['pool']['exits'])
  pids.extend([summary['master_pid'],*summary['pool']['worker_pids']]);journal_events+=seq
  assert digest(final,newline=True)==sha(b/'state.json')==sha(b/'recovered_overlay.json')
  dump(b/'recovery_validation.json',{'status':'passed','state_sha256':sha(b/'state.json'),'recovered_sha256':sha(b/'recovered_overlay.json'),'kernel_calls':0})
  batchrows.append({k:summary[k] for k in ('batch','selected','new_starts','new_passed','new_split','queue','wall_seconds')});previous_overlay=final
 allnew=[*ov['new_passed'],*ov['new_splits']];minimum=None;minH=None;precisions=Counter();depths=Counter()
 for row in allnew:
  node=row['node'];cert=row.get('certificate',row.get('failed_certificate'));n=tuple(node);got=returns[n];record=got['record'];point=list(map(F,record['point']))
  assert allmerges[n]['certificate']==cert and cert['anchor']==','.join(record['point'])
  values=[tuple(map(F,pair)) for pair in record['values']];hl,hu=values[0];assert F(1,10000)<hl<=hu
  for lo,hi in values:assert lo<=hi and hi-lo<=max(abs(lo),abs(hi),F(1,10**6))/2**40
  gl=[tuple(map(F,pair)) for pair in cert['loggradient_directed']]
  for (a,b),(l,u) in zip(values[1:],gl):assert l<=min(a/hl,a/hu,b/hl,b/hu)<=max(a/hl,a/hu,b/hl,b/hu)<=u
  ctx.prec=512;assert F(cert['Delta_lower'])<=api.endpoint((api.ab(hl)/api.ab(F(1,10000))).log())
  lo,hi=physical(node,source['roots']);assert all(a<=p<=b for a,p,b in zip(lo,point[1:],hi))
  margins=[]
  for bits in product((0,1),repeat=4):
   v=[(hi[j] if bits[j] else lo[j])-point[j+1] for j in range(4)]
   margins.append((F(cert['Delta_lower'])+sum(min(a*x,b*x) for x,(a,b) in zip(v,gl))-256*sum(x*x for x in v),bits))
  m,bits=min(margins);assert str(m)==cert['minimum_log_margin'] and list(bits)==cert['worst_corner_bits'] and (m>0)==cert['passed']
  assert cert['mechanism']=='full_frame_16' and cert['checked_points']==16
  if cert['passed']:minimum=m if minimum is None else min(minimum,m)
  minH=hl if minH is None else min(minH,hl);precisions[record['bits']]+=1;depths[len(node[1])]+=1
 # Reconstruct exact canonical final bytes from immutable077+small overlay.
 assert all(k not in source['cache'] for k in ov['new_cache'])
 for k,v in ov['new_cache'].items():source['cache'][k]=v
 source['passed']+=ov['new_passed'];source['splits']+=ov['new_splits'];source['queue']=[]
 source.update({'task070':'upper_cover070','status070':'closed_root_cover_candidate_waiting_main_review',
  'fingerprint070':digest({'run_registration':sha(WORK/'registration.json'),'last_overlay':sha(current['state_path'])}),
  'conversion070':{'source_state':str(sourcepath),'source_state_sha256':reg['source_sha256'],'last_overlay':current['state_path'],'last_overlay_sha256':sha(current['state_path']),
   'new_counts':ov['counts'],'new_passed':len(ov['new_passed']),'new_splits':len(ov['new_splits']),'math':'075 Q_full=I4/16','closed_root_cover_candidate':True,'global_g_or_C5_claim':False}})
 expectedhash=digest(source,newline=True);assert expectedhash==sha(WORK/'final/state.json')
 leaves=defaultdict(list);splits=set()
 for row in source['passed']:leaves[row['node'][0]].append(row['node'][1])
 for row in source['splits']:
  r,p=row['node'];assert (r,p) not in splits;splits.add((r,p));assert row['children']==[[r,p+'0'],[r,p+'1']] and row['axis']==len(p)%4
 for r in range(52360):
  ps=sorted(leaves[r]);assert ps and all(not b.startswith(a) for a,b in zip(ps,ps[1:]));d=max(map(len,ps));assert sum(1<<(d-len(p)) for p in ps)==1<<d
  for p in ps:
   for i in range(len(p)):assert (r,p[:i]) in splits
 assert len(source['passed'])==52360+len(source['splits'])==97960
 assert sha(sourcepath)==reg['source_sha256']
 for p,h in reg['locks'].items():assert sha(p)==h
 result={'status':'author_no_H_verification_passed','pid':os.getpid(),'wall_seconds':time.monotonic()-began,
  'journal_events':journal_events,'verified_new_records':len(ov['new_cache']),'verified_new_passes':len(ov['new_passed']),'verified_new_splits':len(ov['new_splits']),
  'new_corners':16*len(allnew),'precisions':dict(precisions),'new_node_depths':dict(depths),'minimum_new_passed_margin':str(minimum),
  'minimum_new_passed_margin_floor_1e12':str(F(minimum.numerator*10**12//minimum.denominator,10**12)),
  'minimum_new_cached_H_lower':str(minH),'minimum_new_cached_H_lower_floor_1e12':str(F(minH.numerator*10**12//minH.denominator,10**12)),
  'full_state_sha256':expectedhash,'full_state_bytes':(WORK/'final/state.json').stat().st_size,'cache':len(source['cache']),
  'passed':len(source['passed']),'split':len(source['splits']),'queue':0,'inflight':0,'roots_exact_Kraft':52360,
  'source_cache_and_old_objects_preserved_by_reconstructed_file_hash':True,'batches':batchrows,'all_production_pids':pids,
  'H_gradient_kernel_calls_in_this_verification':0,'new_log_calls_in_this_verification':len(allnew),
  'main_peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'global_g_or_C5_claim':False}
 dump(WORK/'final/author_validation.json',result);print(json.dumps(result,indent=2))
if __name__=='__main__':verify()
