#!/usr/bin/env python3
"""Verify and replay the fixed-H64 obstruction without evaluating H."""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
from math import lcm
import gzip,hashlib,importlib.util,json,os,subprocess,time
from flint import arb,fmpq,ctx
W=Path('/Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_g105089');T=F(21,200000)
ENC=json.JSONEncoder(sort_keys=True,separators=(',',':'))
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  while b:=f.read(8*1024*1024):h.update(b)
 return h.hexdigest()
def digest(x):return hashlib.sha256(ENC.encode(x).encode()).hexdigest()
def write(p,x):p.write_text(json.dumps(x,sort_keys=True,indent=2)+'\n')
def ab(x):
 x=F(x);return arb(fmpq(x.numerator,x.denominator))
def box(node,roots):
 lo=[F(x,32) for x in roots[node[0]]];hi=[x+F(1,32) for x in lo]
 for bit in node[1]:
  axis=max(range(4),key=lambda j:hi[j]-lo[j]);mid=(lo[axis]+hi[axis])/2
  if bit=='0':hi[axis]=mid
  else:assert bit=='1';lo[axis]=mid
 return lo,hi
def main():
 began=time.monotonic();b=W/'batch0009';reg=json.load(open(W/'registration.json'));br=json.load(open(b/'registration.json'))
 assert sha(Path(reg['source']))==reg['source_sha256'];source=json.load(open(reg['source']));assert source['epsilon']=='1/64' and source['target']=='21/200000'
 for p,h in reg['locks'].items():assert sha(Path(p))==h
 production=Path(reg['code_path']);assert sha(production)==reg['code_sha256']==br['code_sha256']
 spec=importlib.util.spec_from_file_location('frozen089_replay',production);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
 state=json.load(open(b/'initial_overlay.json'));assert br['parent_overlay_sha256']==sha(Path(br['parent_overlay']))
 previous='0'*64;seq=0;tasks={};returns={};merges={}
 for line in open(b/'journal.jsonl'):
  e=json.loads(line);h=e.pop('hash');assert e['previous']==previous and e['seq']==seq+1 and digest(e)==h;e['hash']=h;module.apply_event(state,e);seq+=1;previous=h;body=e['body']
  if e['kind']=='dispatch':tasks[body['ordinal']]=dict(body)
  elif e['kind']=='received':returns[body['ordinal']]=body['result']
  elif e['kind']=='merge':merges[body['ordinal']]=body
 assert seq==37 and state['status']=='counterexample' and len(tasks)==len(returns)==len(merges)==12
 # Event application mutates dispatch dictionaries; reconstruct the original
 # task hash by dropping only the received-result field.
 for task in tasks.values():task.pop('returned',None)
 expected=hashlib.sha256()
 for chunk in ENC.iterencode(state):expected.update(chunk.encode())
 expected.update(b'\n');assert expected.hexdigest()==sha(b/'state.json')
 (b/'recovered_counter_overlay.json').write_text(ENC.encode(state)+'\n');assert sha(b/'recovered_counter_overlay.json')==sha(b/'state.json')
 tablepath=next(p for p in reg['locks'] if p.endswith('boundary_triangulation.json.gz'));vertices=json.loads(gzip.decompress(Path(tablepath).read_bytes()))['vertices']
 obstacles=[];normal=0;ctx.prec=768
 for i in range(12):
  task=tasks[i];got=returns[i];record=got['record'];point=list(map(F,record['point']));assert task['node']==br['selected_nodes'][i] and task['point']==record['point']
  assert point[0]==1 and all(0<=x<=1 for x in point[1:]);lo,hi=box(task['node'],source['roots']);assert all(a<=p<=z for a,p,z in zip(lo,point[1:],hi))
  D=lcm(*(p.denominator for p in point));ints=[int(D*p) for p in point];kk=[sum(a*z for a,z in zip(ints,v)) for v in vertices];assert 0 not in kk and len(set(kk))==280
  cp=json.load(open(task['point_checkpoint']));assert cp['record']==record and cp['task_sha256']==digest(task) and record['bits']==256 and got['calls']==1
  assert len(record['attempts'])==1 and record['attempts'][0]['all_cones']==3504 and record['attempts'][0]['numeric_quality'] is True
  vals=[tuple(map(F,v)) for v in record['values']];hl,hu=vals[0]
  for a,z in vals:assert a<=z and z-a<=max(abs(a),abs(z),F(1,10**6))/2**40
  if got['kind']=='counterexample':
   assert hu<T and str(i) in state['inflight'] and cp['stage']=='counterexample'
   obstacles.append({'ordinal':i,'node':task['node'],'point':record['point'],'H_interval':record['values'][0],'strict_target_minus_H_upper':str(T-hu),'point_checkpoint_sha256':sha(Path(task['point_checkpoint']))})
  else:
   assert got['kind']=='complete' and hl>T and str(i) not in state['inflight'];normal+=1;c=merges[i]['certificate'];assert c['passed'] is False
   gradient=[tuple(map(F,v)) for v in c['loggradient_directed']]
   for (a,z),(l,u) in zip(vals[1:],gradient):assert l<=min(a/hl,a/hu,z/hl,z/hu)<=max(a/hl,a/hu,z/hl,z/hu)<=u
   delta=F(c['Delta_lower']);assert delta<=F(str((ab(hl)/ab(T)).log().lower().fmpq()))
   m=delta;bits=[]
   for j,(l,u) in enumerate(gradient):
    vv=[lo[j]-point[j+1],hi[j]-point[j+1]];zz=[min(l*v,u*v)-256*v*v for v in vv];k=min(range(2),key=lambda k:(zz[k],k));m+=zz[k];bits.append(k)
   assert m==F(c['minimum_log_margin'])<=0 and bits==c['worst_corner_bits']
 assert normal==8 and len(obstacles)==4 and len(state['queue'])==190 and len(state['inflight'])==4
 virtual=state['queue']+[x['node'] for x in state['inflight'].values()];assert len({tuple(n) for n in virtual})==194
 paths=defaultdict(list);internal=set()
 for row in source['passed']+state['new_passed']:paths[row['node'][0]].append(row['node'][1])
 for r,p in virtual:paths[r].append(p)
 for row in source['splits']+state['new_splits']:
  r,p=row['node'];assert (r,p) not in internal;internal.add((r,p));assert row['children']==[[r,p+'0'],[r,p+'1']] and row['axis']==len(p)%4
 for r in range(52360):
  ps=sorted(paths[r]);assert ps and all(not z.startswith(a) for a,z in zip(ps,ps[1:]));d=max(map(len,ps));assert sum(1<<(d-len(p)) for p in ps)==1<<d
  for p in ps:
   for j in range(len(p)):assert (r,p[:j]) in internal
 assert sum(map(len,paths.values()))==52360+len(internal)==100998
 parent=json.load(open(br['parent_overlay']));assert all(state['new_cache'][k]==v for k,v in parent['new_cache'].items())
 assert state['new_passed']==parent['new_passed'] and state['new_splits'][:len(parent['new_splits'])]==parent['new_splits'] and len(state['new_splits'])==len(parent['new_splits'])+8
 pool=state['pool_evidence'];assert pool['close_called'] and pool['join_returned'] and all(x['exitcode']==0 and not x['is_alive'] for x in pool['exits'])
 pids=[pool['master_pid'],*pool['worker_pids']];ps=subprocess.run(['ps','-p',','.join(map(str,pids)),'-o','pid,ppid,command'],capture_output=True,text=True);assert len(ps.stdout.strip().splitlines())==1
 assert sha(Path(reg['source']))==reg['source_sha256']
 result={'status':'author_strict_fixed_H64_obstruction_and_stop_verified','pid':os.getpid(),'source086_sha256':reg['source_sha256'],'target':'21/200000','epsilon':'1/64',
  'selected':186,'dispatched':12,'new_records':12,'positive_complete':8,'strict_H64_obstructions':4,'obstacles':obstacles,
  'journal_events':37,'state_sha256':sha(b/'state.json'),'recovered_counter_state_sha256':sha(b/'recovered_counter_overlay.json'),
  'passed':100804,'split':48638,'queue':190,'retained_counter_inflight':4,'pending_frontier_total':194,'all_root_exact_Kraft':52360,
  'normal_positive_cache_total':180040+len(state['new_cache']),'full_raw_H64_cache_if_counter_included':180040+len(state['new_cache'])+4,
  'old_source_and_previous_positive_prefix_preserved':True,'all13_production_pids_OS_absent':True,'wall_seconds':time.monotonic()-began,
  'new_H_physical_gradient_kernel_calls':0,'fixed_H64_target_impossible':True,'not_G_or_C5_refutation':True,'no_new_dispatch_after_stop':True}
 write(W/'author_stop_validation_batch0009.json',result);print(json.dumps(result,indent=2))
if __name__=='__main__':main()
