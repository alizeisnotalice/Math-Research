#!/usr/bin/env python3
"""Author no-H verification of finite089 prefixes, including exact pre-split."""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict,Counter
from math import lcm
import gzip,hashlib,json,os,time,resource
from flint import arb,fmpq,ctx
WORK=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_g105089')
ENC=json.JSONEncoder(sort_keys=True,separators=(',',':'))
T=F(21,200000)
def digest(x):
 h=hashlib.sha256()
 for c in ENC.iterencode(x):h.update(c.encode())
 return h.hexdigest()
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  while b:=f.read(8*1024*1024):h.update(b)
 return h.hexdigest()
def write(p,x):p.write_text(json.dumps(x,sort_keys=True,indent=2)+'\n')
def box(node,roots):
 r,p=node;lo=[];hi=[]
 for j,k in enumerate(roots[r]):
  bits=p[j::4];n=len(bits);v=int(bits,2) if bits else 0;a=F(k*(1<<n)+v,32*(1<<n));lo.append(a);hi.append(a+F(1,32*(1<<n)))
 return lo,hi
def ab(x):
 x=F(x);return arb(fmpq(x.numerator,x.denominator))
def verify():
 began=time.monotonic();reg=json.load(open(WORK/'registration.json'));sourcepath=Path(reg['source']);assert sha(sourcepath)==reg['source_sha256']
 for p,h in reg['locks'].items():assert sha(p)==h
 source=json.load(open(sourcepath));assert source['epsilon']==source['cache_model_epsilon']=='1/64' and source['target']=='21/200000'
 current=json.load(open(WORK/'CURRENT.json'));ov=json.load(open(current['state_path']));assert sha(current['state_path'])==current['state_sha256']
 assert not ov['inflight'] and not ov['faults'] and ov['counterexample'] is None and ov['pre_split_count']==1856
 pending={tuple(x['node']):x for x in source['pending_certificates086']};assert len(pending)==1856
 pre=ov['new_splits'][:1856];pre_map=json.load(open(WORK/'pre_split_mapper.json'));assert len(pre_map)==len(pre)==1856
 expectedchildren=set()
 for row,m in zip(pre,pre_map):
  n=tuple(row['node']);c=pending[n]['certificate'];assert row['failed_certificate']==c and c['passed'] is False and c['target']=='21/200000' and F(c['minimum_log_margin'])<=0
  assert row['axis']==len(n[1])%4 and row['children']==[[n[0],n[1]+'0'],[n[0],n[1]+'1']] and m['children']==row['children'] and m['failed_certificate_sha256']==digest(c)
  expectedchildren.update(tuple(x) for x in row['children'])
 assert {tuple(x['node']) for x in pre}=={tuple(x) for x in source['queue']} and len(expectedchildren)==3712
 tablepath=next(p for p in reg['locks'] if p.endswith('boundary_triangulation.json.gz'));vertices=json.loads(gzip.decompress(Path(tablepath).read_bytes()))['vertices'];assert len(vertices)==280
 def geometry(point):
  D=lcm(*(x.denominator for x in point));ints=[int(x*D) for x in point];projections=[sum(a*b for a,b in zip(ints,v)) for v in vertices]
  return {'common_denominator':str(D),'distinct_count':len(set(projections)),'zero_count':projections.count(0),'eligible':0 not in projections and len(set(projections))==280}
 def choose(node):
  lo,hi=box(node,source['roots']);mid=[(a+b)/2 for a,b in zip(lo,hi)];eta=min(b-a for a,b in zip(lo,hi))/4096
  candidates=[mid]+[[x+eta*beta**j for j,x in enumerate(mid)] for beta in (F(1,8),F(1,16),F(1,32))];tests=[]
  for i,tail in enumerate(candidates):
   p=[F(1),*tail];test={'candidate':i,**geometry(p)};tests.append(test)
   if test['eligible']:return list(map(str,p)),i,tests
  raise AssertionError('no generic candidate')
 returns={};tasks={};merges={};events=0;cachehits=0;newrecords=0;previous_overlay=None;allpids=[]
 for b in sorted(WORK.glob('batch[0-9][0-9][0-9][0-9]')):
  if not (b/'summary.json').exists():continue
  br=json.load(open(b/'registration.json'));initial=json.load(open(b/'initial_overlay.json'));final=json.load(open(b/'state.json'));summary=json.load(open(b/'summary.json'))
  assert br['source_state_sha256']==reg['source_sha256'] and br['code_sha256']==reg['code_sha256'] and br['selected_nodes']==initial['queue'][:4096]
  if previous_overlay is None:
   assert {tuple(n) for n in initial['queue']}==expectedchildren and len(initial['queue'])==3712 and initial['new_splits']==pre
  else:
   assert sha(br['parent_overlay'])==br['parent_overlay_sha256'];expected={**previous_overlay,'batch_number':br['batch'],'journal_seq':0,'journal_hash':'0'*64,'next_merge':0,'status':'running_finite_batch'};assert initial==expected
  oldcache={**source['cache'],**initial['new_cache']};localtasks={};localreturns={};localmerges={};previous='0'*64;seq=0
  for line in open(b/'journal.jsonl'):
   e=json.loads(line);h=e.pop('hash');assert e['previous']==previous and e['seq']==seq+1 and digest(e)==h;previous=h;seq+=1;body=e['body']
   if e['kind']=='dispatch':
    i=body['ordinal'];assert i==len(localtasks) and body['node']==br['selected_nodes'][i]
    point,candidate,tests_=choose(body['node']);assert point==body['point'] and candidate==body['candidate'] and tests_==body['geometry_tests']
    assert body['cached']==(','.join(point) in oldcache);localtasks[i]=body
   elif e['kind']=='received':localreturns[body['ordinal']]=body['result']
   elif e['kind']=='merge':localmerges[body['ordinal']]=body['certificate']
   elif e['kind']=='finish':assert body['status']=='finite_prefix_complete'
   else:raise AssertionError(e['kind'])
  assert seq==3*len(localtasks)+1 and seq==final['journal_seq'] and previous==final['journal_hash'] and set(localtasks)==set(localreturns)==set(localmerges)
  for i,task in localtasks.items():
   n=tuple(task['node']);assert n not in tasks;tasks[n]=task;got=localreturns[i];assert got['kind']=='complete';returns[n]=got;merges[n]=localmerges[i];anchor=','.join(task['point'])
   if task['cached']:assert got['record']==oldcache[anchor] and got['calls']==0;cachehits+=1
   else:
    cp=json.load(open(task['point_checkpoint']));assert cp['stage']=='complete' and cp['task_sha256']==digest(task) and cp['record']==got['record']==ov['new_cache'][anchor]
    assert all(a['all_cones']==3504 for a in got['record']['attempts']);newrecords+=1
  assert len(summary['pool']['exits'])==12 and summary['pool']['close_called'] and summary['pool']['join_returned'] and all(x['exitcode']==0 and not x['is_alive'] for x in summary['pool']['exits'])
  assert sha(b/'state.json')==sha(b/'recovered_overlay.json')
  write(b/'recovery_validation.json',{'status':'passed','state_sha256':sha(b/'state.json'),'recovered_sha256':sha(b/'recovered_overlay.json'),'kernel_calls':0})
  allpids.extend([summary['master_pid'],*summary['pool']['worker_pids']]);previous_overlay=final;events+=seq
 ctx.prec=768;minmargin=None;minH=None;bits=Counter();allnew=[*ov['new_passed'],*ov['new_splits'][1856:]]
 for row in allnew:
  n=tuple(row['node']);cert=row.get('certificate',row.get('failed_certificate'));record=returns[n]['record'];assert cert==merges[n] and cert['anchor']==','.join(record['point'])
  vals=[tuple(map(F,p)) for p in record['values']];hl,hu=vals[0];assert T<hl<=hu
  for a,b in vals:assert a<=b and b-a<=max(abs(a),abs(b),F(1,10**6))/2**40
  gradient=[tuple(map(F,p)) for p in cert['loggradient_directed']]
  for (a,b),(gl,gu) in zip(vals[1:],gradient):assert gl<=min(a/hl,a/hu,b/hl,b/hu)<=max(a/hl,a/hu,b/hl,b/hu)<=gu
  delta=F(cert['Delta_lower']);assert delta<=F(str((ab(hl)/ab(T)).log().lower().fmpq()))
  lo,hi=box(row['node'],source['roots']);point=list(map(F,record['point']));m=delta;worst=[]
  for j,(gl,gu) in enumerate(gradient):
   v=[lo[j]-point[j+1],hi[j]-point[j+1]];assert v[0]<=0<=v[1];z=[min(gl*x,gu*x)-256*x*x for x in v];k=min(range(2),key=lambda k:(z[k],k));m+=z[k];worst.append(k)
  assert m==F(cert['minimum_log_margin']) and worst==cert['worst_corner_bits'] and (m>0)==cert['passed'] and cert['checked_points']==16 and cert['quadratic']=='Q_full=I4/16'
  if cert['passed']:minmargin=m if minmargin is None else min(minmargin,m)
  minH=hl if minH is None else min(minH,hl);bits[record['bits']]+=1
 assert len(allnew)==len(tasks) and newrecords==len(ov['new_cache'])==ov['counts']['new_complete'] and cachehits==ov['counts']['cache_hits']
 paths=defaultdict(list);internal=set()
 for row in source['passed']+ov['new_passed']:paths[row['node'][0]].append(row['node'][1])
 for r,p in ov['queue']:paths[r].append(p)
 for row in source['splits']+ov['new_splits']:
  r,p=row['node'];assert (r,p) not in internal;internal.add((r,p));assert row['children']==[[r,p+'0'],[r,p+'1']] and row['axis']==len(p)%4
 for r in range(52360):
  ps=sorted(paths[r]);assert ps and all(not b.startswith(a) for a,b in zip(ps,ps[1:]));d=max(map(len,ps));assert sum(1<<(d-len(p)) for p in ps)==1<<d
  for p in ps:
   for i in range(len(p)):assert (r,p[:i]) in internal
 assert sum(map(len,paths.values()))==52360+len(internal)
 assert sha(sourcepath)==reg['source_sha256']
 for p,h in reg['locks'].items():assert sha(p)==h
 result={'status':'author_no_H_prefix_verification_passed','pid':os.getpid(),'batch':current['batch'],'source086_sha256':reg['source_sha256'],'target':'21/200000','epsilon':'1/64',
  'zero_H_pre_split':1856,'children':3712,'journal_events':events,'verified_child_certificates':len(allnew),'verified_child_corners':16*len(allnew),
  'new_records':newrecords,'cache_hits':cachehits,'record_precisions':dict(bits),'minimum_passed_margin':str(minmargin),'minimum_used_cached_H_lower':str(minH),
  'passed':len(source['passed'])+len(ov['new_passed']),'split':len(internal),'queue':len(ov['queue']),'inflight':0,'cache':len(source['cache'])+len(ov['new_cache']),
  'root_Kraft_exact':52360,'state_sha256':sha(current['state_path']),'old_source_cache_and_leaf_objects_read_only':True,'all_production_pids':allpids,
  'new_H_physical_gradient_kernel_calls':0,'wall_seconds':time.monotonic()-began,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'g105_closed_cover':not ov['queue'],'global_g_or_C5_claim':False}
 write(WORK/f"author_validation_batch{current['batch']:04d}.json",result);print(json.dumps(result,indent=2))
if __name__=='__main__':verify()
