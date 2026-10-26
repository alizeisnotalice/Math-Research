#!/usr/bin/env python3
"""Finite author gate for one frozen123 batch; no numerical H calls."""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict,Counter
import argparse,hashlib,json,os,time,resource
from flint import arb,fmpq,ctx
H=Path(__file__).resolve().parent;W=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_continue123');T=F(11,100000);K='128d98c94e55632336edafe024aa3a27c1486898670a52be89cb39cc22ed6cc0'
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  while b:=f.read(8*1024*1024):h.update(b)
 return h.hexdigest()
def digest(x):
 h=hashlib.sha256()
 for b in json.JSONEncoder(sort_keys=True,separators=(',',':')).iterencode(x):h.update(b.encode())
 return h.hexdigest()
def write(p,x):
 with open(p,'x') as f:json.dump(x,f,sort_keys=True,separators=(',',':'));f.write('\n');f.flush();os.fsync(f.fileno())
def ab(q):q=F(q);return arb(fmpq(q.numerator,q.denominator))
def box(n,roots):
 lo=[F(x,32) for x in roots[n[0]]];hi=[x+F(1,32) for x in lo]
 for i,bit in enumerate(n[1]):
  j=i%4;m=(lo[j]+hi[j])/2
  if bit=='0':hi[j]=m
  else:assert bit=='1';lo[j]=m
 return lo,hi
def main():
 began=time.monotonic();ap=argparse.ArgumentParser();ap.add_argument('--batch',type=int,choices=(1,2),required=True);args=ap.parse_args();B=W/f'batch{args.batch:04d}';plan=json.load(open(W/'finite_plan.json'));selection=plan['batches'][args.batch-1];summary=json.load(open(B/'summary.json'));assert summary['status']=='finite_original_prefix_complete' and summary['faults']==summary['inflight']==summary['obstructions_H128']==0
 write(W/f'author_batch{args.batch:04d}_registration.json',{'pid':os.getpid(),'code_sha256':sha(__file__),'batch':args.batch,'batch_manifest_sha256':sha(B/'manifest.json'),'selected':selection['count'],'new_H_gradient_kernel_calls':0,'fresh_log_verification_bits':768})
 package=json.load(open(H/'production_prepared_manifest.json'));assert sha(H/'production_prepared_manifest.json')=='1edbc7912626737020fec146c80247f3f5876603732d97b0d6ecdd45aa46165f'
 for n,h in package['repository_files'].items():assert sha(H/n)==h
 for p,h in package['external_files'].items():assert sha(p)==h
 for p,h in package['source_locks'].items():assert sha(p)==h
 seed=json.load(open(W/'seed.json'));assert sha(W/'seed.json')==plan['seed_sha256'];reg=json.load(open(B/'registration.json'));initial=json.load(open(B/'initial_overlay.json'));ov=json.load(open(B/'state.json'));tasks=json.load(open(selection['tasks_path']));assert len(tasks)==selection['count'] and reg['worker_count']==8 and reg['target']=='11/100000'
 assert digest(seed['cache_by_epsilon'])==reg['source_cache_canonical_sha256'] and digest(seed['passed'])==reg['source_passed_canonical_sha256'] and len(seed['passed'])==98007 and len(seed['splits'])==51884
 if args.batch==1:assert not initial['new_passed'] and not initial['new_splits'] and not initial['new_cache_by_epsilon']['1/128']
 else:
  prior=json.load(open(W/'author_batch0001_validation.json'));parent=reg['parent_checkpoint'];assert prior['status']=='author_finite_batch_checks_passed' and prior['state_sha256']==parent['state_sha256']
  assert sha(parent['state'])==parent['state_sha256'] and sha(parent['journal'])==parent['journal_sha256'];parentstate=json.load(open(parent['state']))
  for k in ('new_passed','new_splits','new_cache_by_epsilon','queue','inflight','obstructions_H128','faults','counts'):assert initial[k]==parentstate[k]
 assert initial['queue'][:len(tasks)]==[t['node'] for t in tasks] and not ov['inflight'] and not ov['faults'] and not ov['obstructions_H128'] and all(not t['cached'] for t in tasks)
 manifest=json.load(open(B/'manifest.json'))
 for n,h in manifest['files'].items():assert sha(B/n)==h
 dispatched=[];received={};certs={};merged=[];previous='0'*64;events=0;finish=None
 for line in open(B/'journal.jsonl'):
  e=json.loads(line);h=e.pop('hash');assert e['seq']==events+1 and e['previous']==previous and digest(e)==h;events+=1;previous=h;b=e['body']
  if e['kind']=='dispatch':assert b['task']==tasks[len(dispatched)];dispatched.append(b['task'])
  elif e['kind']=='received':i=b['ordinal'];assert i<len(dispatched) and i not in received;received[i]=b['returned'];assert received[i]['kind']=='complete'
  elif e['kind']=='merge':i=b['ordinal'];assert i==len(merged) and i in received;merged.append(i);certs[i]=b['certificate']
  elif e['kind']=='finish':assert finish is None;finish=b
  else:raise AssertionError('unexpected fault/status journal event')
 count=len(tasks);assert len(dispatched)==len(received)==len(merged)==count and events==3*count+1==ov['journal_seq'] and previous==ov['journal_hash'] and finish['status']==ov['status']=='finite_original_prefix_complete'
 points=json.load(open(B/'point_manifest.json'))['records'];assert len(points)==count
 newcache=dict(initial['new_cache_by_epsilon']['1/128']);newpassed=list(initial['new_passed']);newsplits=list(initial['new_splits']);children=[];calls=0;multisegment=0;bits=Counter();minHL=None;minpass=None;minfail=None;maxfail=None;ctx.prec=768
 for i,(task,point) in enumerate(zip(tasks,points)):
  assert task['batch_ordinal']==point['ordinal']==i and point['global_ordinal']==task['ordinal'] and point['node']==task['node'] and point['point']==task['point'] and point['pair_cache_key']==task['pair_cache_key'] and not point['cached']
  assert point['stage']=='complete' and point['target']=='11/100000' and point['epsilon']=='1/128' and sha(point['point_checkpoint'])==point['point_checkpoint_sha256'];cp=json.load(open(point['point_checkpoint']));record=point['record'];got=received[i]
  assert cp['record']==record==got['record'] and cp['partial'] is None and cp['point']==record['point']==task['point'] and cp['pair_cache_key']==record['pair_cache_key']==got['pair_cache_key']==task['pair_cache_key']
  assert cp['epsilon']==record['epsilon']==got['epsilon']=='1/128' and cp['target']==record['qualification_target']==got['target']=='11/100000' and cp['kernel_sha256']==record['kernel_sha256']==K
  assert cp['calls']==got['calls']>=len(cp['attempts']) and cp['attempts']==record['attempts']==got['attempts'] and cp['precision_index']==len(cp['attempts'])-1 and cp['attempts'][-1]['numeric_quality'] and all(t['all_cones']==3504 for t in cp['attempts']);calls+=cp['calls'];multisegment+=cp['calls']>len(cp['attempts']);bits[record['bits']]+=1
  assert record['geometry_tests']==task['geometry_tests'] and record['geometry_tests'][-1]['eligible'] and record['geometry_tests'][-1]['zero_count']==0 and record['geometry_tests'][-1]['distinct_count']==280
  if task['candidate']!=task['original_H64_candidate']:assert task['same_point_H64_guard'] is None
  guard=task['same_point_H64_guard']
  if guard is not None:assert guard['point']==record['point'] and F(record['values'][0][1])>=F(guard['H_interval'][0])
  p=','.join(map(str,map(F,record['point'])));assert p==task['pair_cache_key'][1] and p not in seed['cache_by_epsilon']['1/128'] and p not in newcache;newcache[p]=record
  values=[tuple(map(F,x)) for x in record['values']];hl,hu=values[0];assert T<hl<=hu;minHL=hl if minHL is None else min(minHL,hl);c=certs[i]
  assert c['target']=='11/100000' and c['leaf_epsilon']=='1/128' and c['penalty_coefficient']=='1024' and c['point_key']==p and c['checked_points']==16;gradient=[tuple(map(F,x)) for x in c['loggradient_directed']]
  for (a,b),(l,u) in zip(values[1:],gradient):assert l<=min(a/hl,a/hu,b/hl,b/hu)<=max(a/hl,a/hu,b/hl,b/hu)<=u
  delta=F(c['Delta_lower']);assert 0<delta<=F(str((ab(hl)/ab(T)).log().lower().fmpq()));lo,hi=box(task['node'],seed['roots']);assert list(map(str,lo))==task['box_lower'] and list(map(str,hi))==task['box_upper'];anchor=list(map(F,record['point']));assert all(a<=x<=b for a,x,b in zip(lo,anchor[1:],hi))
  minimum=delta;worst=[]
  for j,(l,u) in enumerate(gradient):
   offsets=[lo[j]-anchor[j+1],hi[j]-anchor[j+1]];terms=[min(l*x,u*x)-1024*x*x for x in offsets];k=min(range(2),key=lambda k:(terms[k],k));minimum+=terms[k];worst.append(k)
  assert minimum==F(c['minimum_log_margin']) and worst==c['worst_corner_bits'] and c['passed']==(minimum>0)
  if c['passed']:newpassed.append({'node':task['node'],'leaf_epsilon':'1/128','lower_bound_object':'G','target':'11/100000','certificate':c});minpass=minimum if minpass is None else min(minpass,minimum)
  else:
   node=task['node'];child=[[node[0],node[1]+'0'],[node[0],node[1]+'1']];newsplits.append({'node':node,'axis':len(node[1])%4,'children':child,'leaf_epsilon':'1/128','target':'11/100000','failed_certificate':c});children.extend(child);minfail=minimum if minfail is None else min(minfail,minimum);maxfail=minimum if maxfail is None else max(maxfail,minimum)
 assert newcache==ov['new_cache_by_epsilon']['1/128'] and newpassed==ov['new_passed'] and newsplits==ov['new_splits'] and ov['queue']==initial['queue'][count:]+children
 expected={k:v+({'new_anchor_starts':count,'new_complete':count,'kernel_segment_calls':calls,'cache_hits':0}[k]) for k,v in initial['counts'].items()};assert ov['counts']==expected
 by=defaultdict(list);inner=set()
 for leaf in seed['passed']+ov['new_passed']:r,p=leaf['node'];by[r].append(p)
 for r,p in ov['queue']:by[r].append(p)
 for s in seed['splits']+ov['new_splits']:
  r,p=s['node'];assert (r,p) not in inner;inner.add((r,p));assert s['axis']==len(p)%4 and s['children']==[[r,p+'0'],[r,p+'1']]
 for r in range(52360):
  ps=sorted(by[r]);assert ps and all(not b.startswith(a) for a,b in zip(ps,ps[1:]));d=max(map(len,ps));assert sum(1<<(d-len(p)) for p in ps)==1<<d
  for p in ps:
   for j in range(len(p)):assert (r,p[:j]) in inner
 assert 98007+len(newpassed)+len(ov['queue'])==52360+len(inner)
 pool=json.load(open(B/'pool_registration.json'));pids=[pool['master_pid']]+[r['pid'] for r in pool['workers']];assert len(pool['workers'])==8 and all(r['nice']==5 and r['epsilon']=='1/128' and r['target']=='11/100000' for r in pool['workers']) and all(r['exitcode']==0 and not r['alive'] for r in ov['pool_evidence'])
 for pid in pids:
  try:os.kill(pid,0);raise AssertionError(('production PID live',pid))
  except ProcessLookupError:pass
 assert sha(W/'seed.json')==plan['seed_sha256']
 for p,h in package['source_locks'].items():assert sha(p)==h
 result={'status':'author_finite_batch_checks_passed','pid':os.getpid(),'batch':args.batch,'selected':count,'new_H128_records':count,'fresh_Delta_768bit_checks':count,'exact_corner_checks':16*count,'journal_events':events,'new_passed':len(newpassed)-len(initial['new_passed']),'new_splits':len(newsplits)-len(initial['new_splits']),'total_passed':98007+len(newpassed),'queue':len(ov['queue']),'inflight':0,'faults':0,'obstructions_H128':0,'new_counts':{k:ov['counts'][k]-initial['counts'][k] for k in ov['counts']},'all_complete_checkpoint_partial_None':True,'multisegment_complete_points':multisegment,'bits_counts':dict(bits),'minimum_H128_lower':str(minHL),'minimum_accepted_margin':str(minpass) if minpass is not None else None,'minimum_failed_margin':str(minfail) if minfail is not None else None,'maximum_failed_margin':str(maxfail) if maxfail is not None else None,'exact_closed_Kraft_roots':52360,'source_98007_passed_unchanged':True,'source189086_typed_cache_unchanged':True,'original_prefix_only_no_newchild_dispatch':True,'all_production_pids_absent':pids,'state_sha256':sha(B/'state.json'),'point_manifest_sha256':sha(B/'point_manifest.json'),'batch_manifest_sha256':sha(B/'manifest.json'),'new_H_gradient_kernel_calls':0,'wall_seconds':time.monotonic()-began,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'S_credit':0}
 write(W/f'author_batch{args.batch:04d}_validation.json',result);print(json.dumps(result),flush=True)
if __name__=='__main__':main()
