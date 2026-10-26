#!/usr/bin/env python3
"""Author zero-H terminal check, integer log lower bound and direct source compare.

No production kernel/Arb or author transition functions are imported.
Independent numerical re-evaluation remains audit_mixed100's responsibility.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
from itertools import product
import hashlib,json,os,time,resource
HERE=Path(__file__).resolve().parent
WORK=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_mixed098')
B=WORK/'batch0001';TARGET=F(21,200000)
S='19294a9bd493524fdc23ec2abe2bb3e70577f04b762b08d4fc57e78af5377480'
FINAL='73bc873f23f7b68f6f70132e590d28ebbe4cd619012ac9c81980f9fa500314af'
PREP='fe03bfa691b6947a24f13a2f94885962332a77edfce087e2f9e01bc49a3ace15'
BM='0333c44e9d8bc29f2f1c2ac16ab38bbc935af8ae738f334a767812aed4a61928'
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  while b:=f.read(8*1024*1024):h.update(b)
 return h.hexdigest()
def digest(x):return hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def write(p,x):
 t=p.with_suffix(p.suffix+'.tmp')
 with open(t,'w') as f:json.dump(x,f,sort_keys=True,separators=(',',':'));f.write('\n');f.flush();os.fsync(f.fileno())
 os.replace(t,p)
def log_lower(q):
 # For q>1 reduce to y in[1,2). Every fixed-point multiplication rounds
 # downward. Positive atanh series truncation is therefore a strict lower
 # bound, without any floating point or transcendental library call.
 assert q>1;k=0
 while q>=2:q/=2;k+=1
 scale=1<<768
 def series(y):
  z=(y-1)/(y+1);a=z.numerator;b=z.denominator
  power=a*scale//b;total=0
  for j in range(180):
   total+=power//(2*j+1);power=power*a*a//(b*b)
  return 2*total
 return F(series(q)+k*series(F(2)),scale)
def box(node,roots):
 lo=[F(x,32) for x in roots[node[0]]];hi=[x+F(1,32) for x in lo]
 for i,bit in enumerate(node[1]):
  j=i%4;mid=(lo[j]+hi[j])/2
  if bit=='0':hi[j]=mid
  else:assert bit=='1';lo[j]=mid
 return lo,hi
def main():
 began=time.monotonic();out=WORK/'author_terminal_verification.json'
 assert not out.exists(),'Do not overwrite frozen author result'
 write(WORK/'author_terminal_verification_registration.json',{'pid':os.getpid(),'code_sha256':sha(__file__),'source_sha256':S,'final_state_sha256':FINAL,'batch_manifest_sha256':BM,
  'finite_new_points':194,'finite_old_cache':184796,'finite_old_accepted':100804,'new_H_gradient_kernel_Arb_log_calls':0,'integer_log':'768bit lower-rounded180terms atanh after[1,2) range reduction'})
 assert sha(WORK/'initial_mixed_state.json')==S and sha(WORK/'final_state.json')==FINAL and sha(HERE/'manifest.json')==PREP and sha(B/'manifest.json')==BM
 prep=json.load(open(HERE/'manifest.json'))
 for name,h in prep['repository_files'].items():assert sha(HERE/name)==h
 for path,h in prep['external_files'].items():assert sha(path)==h
 bm=json.load(open(B/'manifest.json'))
 for name,h in bm['files'].items():assert sha(B/name)==h
 for name,h in bm['code_sha256'].items():assert sha(HERE/name)==h
 source=json.load(open(WORK/'initial_mixed_state.json'));final=json.load(open(WORK/'final_state.json'));ov=json.load(open(B/'state.json'));initial=json.load(open(B/'initial_overlay.json'))
 assert [len(source[k]) for k in ('passed','splits','queue')]==[100804,48638,194]
 assert [len(final[k]) for k in ('passed','splits','queue')]==[100998,48638,0]
 assert not final['inflight'] and not final['resume_jobs'] and not ov['inflight'] and not ov['queue'] and not ov['faults'] and not ov['obstructions_H128'] and not ov['new_splits']
 assert final['cache_by_epsilon']['1/64']==source['cache_by_epsilon']['1/64'] and len(source['cache_by_epsilon']['1/64'])==184796
 assert final['cache_by_epsilon']['1/128']==ov['new_cache_H128'] and len(ov['new_cache_H128'])==194 and not source['cache_by_epsilon']['1/128']
 assert final['passed'][:100804]==source['passed'] and final['passed'][100804:]==ov['new_passed'] and final['splits']==source['splits']
 assert all(x['leaf_epsilon']=='1/64' for x in source['passed']) and all(x['leaf_epsilon']=='1/128' for x in ov['new_passed'])
 changed={'task','status','passed','splits','queue','inflight','resume_jobs','cache_by_epsilon','accepted_leaf_epsilon','g105_closed_cover','new_H128_records','production_authorized'}
 for k,v in source.items():
  if k not in changed:assert final[k]==v
 assert set(final)-set(source)=={'final098_source_sha256','final098_overlay_sha256','final098_registration_sha256','final098_scope','fingerprint098_final'}
 assert final['target']=='21/200000' and final['lower_bound_object']=='G' and final['g105_closed_cover'] is True
 expectedroots=[[a,b,c,d] for a in range(32) for b in range(a+1) for c in range(b+1) for d in range(c+1)]
 assert source['roots']==final['roots']==expectedroots
 splitkeys=set();byroot=defaultdict(list)
 for row in final['splits']:
  r,p=row['node'];assert (r,p) not in splitkeys;splitkeys.add((r,p));assert row['axis']==len(p)%4 and row['children']==[[r,p+'0'],[r,p+'1']]
 for leaf in final['passed']:
  r,p=leaf['node'];byroot[r].append(p)
  for j in range(len(p)):assert (r,p[:j]) in splitkeys
 assert len(byroot)==52360 and len(final['passed'])==52360+len(splitkeys)
 for r in range(52360):
  ps=sorted(byroot[r]);assert all(not b.startswith(a) for a,b in zip(ps,ps[1:]));d=max(map(len,ps));assert sum(1<<(d-len(p)) for p in ps)==1<<d
 tasks=json.load(open(WORK/'execution_tasks_194.json'));points=json.load(open(B/'point_manifest.json'))['records'];assert len(points)==194
 assert [t['node'] for t in tasks]==source['queue']==[x['node'] for x in ov['new_passed']]
 dispatch=[];received={};merged=[];previous='0'*64;count=0;finish=None
 for line in open(B/'journal.jsonl'):
  e=json.loads(line);h=e.pop('hash');assert e['previous']==previous and e['seq']==count+1 and digest(e)==h;previous=h;count+=1
  b=e['body'];kind=e['kind']
  if kind=='dispatch':assert b['task']==tasks[len(dispatch)];dispatch.append(b['task'])
  elif kind=='received':
   i=b['ordinal'];assert i<len(dispatch) and i not in received;received[i]=b['returned'];assert received[i]['kind']=='complete'
  elif kind=='merge':
   i=b['ordinal'];assert i==len(merged) and i in received and b['certificate']==ov['new_passed'][i]['certificate'];merged.append(i)
  elif kind=='finish':assert finish is None;finish=b
  else:raise AssertionError(kind)
 assert len(dispatch)==len(received)==len(merged)==194 and previous==ov['journal_hash'] and count==ov['journal_seq']==583
 assert finish['status']==ov['status']=='finite194_complete' and finish['pool_evidence']==ov['pool_evidence']
 assert initial['queue']==source['queue'] and not initial['new_passed'] and not initial['new_cache_H128']
 minHL=None;minmargin=None;maxmargin=None;mapper=[];maxwidth=F(0)
 for i,(task,point,leaf) in enumerate(zip(tasks,points,ov['new_passed'])):
  assert point['ordinal']==task['ordinal']==i and point['node']==leaf['node']==task['node'] and point['stage']=='complete'
  assert sha(point['point_checkpoint'])==point['point_checkpoint_sha256'];cp=json.load(open(point['point_checkpoint']));record=point['record'];got=received[i]
  assert cp['record']==record==got['record'] and cp['point']==record['point']==task['point'] and cp['epsilon']==record['epsilon']==got['epsilon']=='1/128'
  assert cp['pair_cache_key']==got['pair_cache_key']==task['pair_cache_key'] and cp['kernel_sha256']==record['kernel_sha256']==bm['code_sha256']['epsilon_kernel.py']
  assert record['bits']==256 and len(record['attempts'])==1 and record['attempts'][0]['all_cones']==3504 and record['attempts'][0]['numeric_quality'] is True
  assert record['geometry_tests'][-1]['eligible'] and record['geometry_tests'][-1]['distinct_count']==280 and record['geometry_tests'][-1]['zero_count']==0
  key=','.join(map(str,map(F,record['point'])));assert final['cache_by_epsilon']['1/128'][key]==record
  vals=[tuple(map(F,p)) for p in record['values']];assert len(vals)==5 and all(a<=b for a,b in vals)
  hl,hu=vals[0];assert TARGET<hl<=hu;minHL=hl if minHL is None else min(minHL,hl);maxwidth=max(maxwidth,hu-hl)
  c=leaf['certificate'];assert c['leaf_epsilon']=='1/128' and c['penalty_coefficient']=='1024' and c['target']=='21/200000' and c['point_key']==key and c['checked_points']==16 and c['passed'] is True
  delta=F(c['Delta_lower']);assert 0<delta<log_lower(hl/TARGET)
  gradients=[tuple(map(F,p)) for p in c['loggradient_directed']]
  for (a,b),(l,u) in zip(vals[1:],gradients):assert l<=min(a/hl,a/hu,b/hl,b/hu)<=max(a/hl,a/hu,b/hl,b/hu)<=u
  lo,hi=box(leaf['node'],expectedroots);anchor=list(map(F,record['point']));assert anchor[0]==1 and all(a<=x<=b for a,x,b in zip(lo,anchor[1:],hi))
  rows=[]
  for bits in product((0,1),repeat=4):
   v=[(hi[j] if bits[j] else lo[j])-anchor[j+1] for j in range(4)]
   rows.append((delta+sum(min(a*x,b*x) for (a,b),x in zip(gradients,v))-1024*sum(x*x for x in v),bits))
  margin,bits=min(rows);assert margin==F(c['minimum_log_margin'])>0 and list(bits)==c['worst_corner_bits']
  minmargin=margin if minmargin is None else min(minmargin,margin);maxmargin=margin if maxmargin is None else max(maxmargin,margin)
  mapper.append({'new_index':i,'source_initial_queue_index':i,'node':leaf['node'],'epsilon':'1/128','record_sha256':digest(record),'certificate_sha256':digest(c),'box_lower':list(map(str,lo)),'box_upper':list(map(str,hi))})
 pools=json.load(open(B/'pool_registration.json'));assert len(pools['workers'])==12 and all(x['exitcode']==0 and not x['alive'] for x in ov['pool_evidence'])
 pids=[pools['master_pid']]+[x['pid'] for x in pools['workers']]
 for pid in pids:
  try:os.kill(pid,0)
  except ProcessLookupError:continue
  raise AssertionError(('production PID still exists',pid))
 write(WORK/'new_H128_leaf_mapper.json',mapper)
 result={'status':'author_zero_H_terminal_checks_passed','pid':os.getpid(),'source_sha256':S,'final_state_sha256':FINAL,'batch_manifest_sha256':BM,
  'direct_old_H64_cache_objects':184796,'direct_old_passed_objects':100804,'new_H128_values':194,'new_physical_gradient_components':776,
  'new_Delta_integer_log_checks':194,'new_exact_corners':3104,'ordered_roots':52360,'exact_Kraft_roots':52360,'final_passed':100998,'final_splits':48638,'queue':0,'inflight':0,
  'journal_events':583,'source_and_final_old_fields_direct_equal':True,'source_tree_mapper_complete':True,'minimum_H128_lower':str(minHL),'minimum_new_margin':str(minmargin),'maximum_new_margin':str(maxmargin),
  'maximum_H128_interval_width':str(maxwidth),'production_pids_absent':pids,'mapper_sha256':sha(WORK/'new_H128_leaf_mapper.json'),'new_H_gradient_kernel_Arb_log_calls':0,
  'wall_seconds':time.monotonic()-began,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'independent_numeric_audit_pending':True,'C5_or_I_claim':False}
 assert sha(WORK/'initial_mixed_state.json')==S and sha(WORK/'final_state.json')==FINAL and sha(B/'manifest.json')==BM
 write(out,result);print(json.dumps(result),flush=True)
if __name__=='__main__':main()
