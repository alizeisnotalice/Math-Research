#!/usr/bin/env python3
"""Finite saved-function hidden-seam chain, stopping before any unpaid jump."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations
import json,hashlib,importlib.util,sys,time,resource,signal,os
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;OUT=HERE/'results.json'
if OUT.exists():raise FileExistsError('Refuse overwrite '+str(OUT))
resource.setrlimit(resource.RLIMIT_CPU,(20,21))
def wallcap(*_):raise TimeoutError('180 second wall timeout')
signal.signal(signal.SIGALRM,wallcap);signal.alarm(180)
CPU,WALL=time.process_time(),time.monotonic()
LOWER=HERE.parents[1];PREV=LOWER/'round066/post_s2_local_event/results.json'
BASE=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds')
ARITH=BASE/'lower/round001_rectangular_certificate.py';IFACE=BASE/'lower/round008/spacing_chamber/chamber_interface.json'
configs,nodes=0,0
out={'status':'running','pid':os.getpid(),'full_union_calls':0,'new_C6_bound':False,'automatic_reruns':0,'events':[],'best_new_candidate':None}
def budget():
 if time.process_time()-CPU>=19.5:raise TimeoutError('CPU prelimit 19.5 seconds')
 if time.monotonic()-WALL>=180:raise TimeoutError('wall limit')
 if configs>32 or nodes>4000:raise RuntimeError('configuration/local node cap')
try:
 spec=importlib.util.spec_from_file_location('accepted_field',ARITH);A=importlib.util.module_from_spec(spec);spec.loader.exec_module(A)
 A.D=10**80;A.RU,A.RV=A.floor_root(2*A.D**6,6),A.floor_root(3*A.D**6,6)
 A.UL,A.UH=Q(A.RU,A.D),Q(A.RU+1,A.D);A.VL,A.VH=Q(A.RV,A.D),Q(A.RV+1,A.D);A.bounds.cache_clear()
 F,ZERO=A.F,A.F(0)
 def field(z):
  a=[Q(0)]*36
  for ij,x in z.items():i,j=map(int,ij.split(','));a[6*i+j]=Q(x)
  return F(a=a)
 def enc(x):return {f'{i},{j}':str(x.a[6*i+j]) for i in range(6) for j in range(6) if x.a[6*i+j]}
 def view(x):
  lo,hi=A.bounds(x.a);return {'coefficients':enc(x),'exact_interval':[str(lo),str(hi)],'display':float((lo+hi)/2)}
 d=json.loads(PREV.read_text());iface=json.loads(IFACE.read_text())
 defs=d['full_functions_reused_from_064'];funcs=[{'c':field(f['c']),'m':Q(f['m'])} for f in defs]
 fixed=[tuple(field(iface['endpoints'][eid]['c0']) for axis in [0,2,3,4,5] for eid in b[axis]) for b in iface['boxes']]
 moving={z['raw_index']:i for i,f in enumerate(defs) for z in f['labels'] if 'raw_index' in z}
 original=field(d['s3']['coefficients']);g=field(d['slope']['coefficients']);baseline=field(d['candidate']['coefficients'])
 old_upper=Q(d['candidate']['exact_interval'][1]);assert ZERO<g
 out.update({'baseline_s3':view(original),'baseline_complete_066_interval':d['candidate']['exact_interval'],'slope':d['slope'],'complete_functions_reused':defs,'old_passed_pair_comparisons':0,'inherited_E066_manifest_sha256':'89ffcf55fcf3a6c92dc860b498709099d8be783bc22f67bf6d5f4a813fcb8b41'})
 def state(event,side):
  global configs
  configs+=1;budget();keys=[(f['c']+f['m']*event,side*f['m']) for f in funcs]
  sels={};active=set()
  for k,ell in enumerate(iface['box_lengths']):
   lo=max((0,) if ell[1]%2 else (0,moving[k]),key=lambda i:keys[i]);hi=min((1,moving[k]) if ell[1]%2 else (1,),key=lambda i:keys[i])
   if keys[lo]<keys[hi]:sels[str(k)]=[lo,hi];active.update((lo,hi))
  return {'positive_boxes':len(sels),'effective_order':sorted(active,key=lambda i:keys[i]),'selections':sels}
 event=original;order=list(d['complete_s2_right_order']);tie_edges=d['s3_tying_adjacent_ids'];left,right=d['s3_states']['left'],d['s3_states']['right']
 for step in range(8):
  budget();rec={'index':step,'parameter':view(event),'left_state':left,'right_state':right,'tie_edges':tie_edges,'left_complete_order':list(order),'projection_pair_checks':[]};out['events'].append(rec)
  if left['selections']!=right['selections'] or left['positive_boxes']!=right['positive_boxes']:
   rec['classification']='stop_branch_or_positive_width_change';out['status']='stopped_branch';break
  # Find tied blocks using only saved/new exact zero gaps at this seam.
  edge_set={tuple(z) for z in tie_edges};blocks=[];block=[order[0]]
  for i,j in zip(order,order[1:]):
   if (i,j) in edge_set:block.append(j)
   else:blocks.append(block);block=[j]
  blocks.append(block)
  failed=False
  for block in blocks:
   if len(block)<2:continue
   for i,j in combinations(block,2):
    gi=[int(k) for k,p in left['selections'].items() if i in p];gj=[int(k) for k,p in left['selections'].items() if j in p]
    for a in gi:
     for b in gj:
      nodes+=1;budget();B,C=fixed[a],fixed[b]
      I=[(max(B[2*q],C[2*q]),min(B[2*q+1],C[2*q+1])) for q in range(5)]
      sep=[q for q,(lo,hi) in enumerate(I) if hi<lo]
      p={'function_ids':[i,j],'raw_indices':[a,b],'lengths':[iface['box_lengths'][a],iface['box_lengths'][b]]}
      if sep:
       q=sep[0];lo,hi=I[q];p.update({'strictly_disjoint':True,'axis':[1,3,4,5,6][q],'gap':view(lo-hi)})
      else:
       p.update({'strictly_disjoint':False,'intersection':[[enc(lo),enc(hi)] for lo,hi in I],'all_five_intervals_positive':all(lo<hi for lo,hi in I)})
       failed=True
      rec['projection_pair_checks'].append(p)
      if failed:break
     if failed:break
    if failed:break
   if failed:break
  if failed:
   rec['classification']='stop_first_nonempty_or_touching_projection';out['status']='stopped_overlap';break
  rec['classification']='hidden_strict_fixed_projection_disjointness'
  order=[i for block in blocks for i in sorted(block,key=lambda k:funcs[k]['m'])]
  rec['right_complete_order']=list(order)
  if step==7:
   out['status']='finished_event_cap';break
  # NEW 36 adjacent checks only at this newly classified event.
  candidates=[]
  for i,j in zip(order,order[1:]):
   nodes+=1;budget();f,h=funcs[i],funcs[j]
   if f['m']<=h['m']:continue
   t=(h['c']-f['c'])/(f['m']-h['m']);assert event<t;candidates.append((t,i,j))
  nxt=min(candidates,key=lambda z:z[0])[0]
  gaps=[]
  for i,j in zip(order,order[1:]):
   f,h=funcs[i],funcs[j];l=h['c']-f['c']+(h['m']-f['m'])*event;r=h['c']-f['c']+(h['m']-f['m'])*nxt
   assert ZERO<=l and ZERO<=r and (ZERO<l or ZERO<r)
   gaps.append({'ids':[i,j],'left_gap':enc(l),'right_gap':enc(r)})
  nextleft=state(nxt,-1);nextright=state(nxt,1)
  if nextleft!=right:
   rec['classification']='stop_unpaid_interior_state_change';out['status']='stopped_unpaid';break
  candidate=baseline+g*(nxt-original);cl=A.bounds(candidate.a)[0];assert old_upper<cl
  target={'event_parameter':view(nxt),'ratio':view(candidate),'lower_strictly_above_complete_066_upper':True,'exact_margin_over_complete_066_upper':str(cl-old_upper)}
  rec['next_interval']={'right_parameter':view(nxt),'adjacent_checks':36,'adjacent_future_candidates':[{'ids':[i,j],'location':view(t)} for t,i,j in candidates],'complete_order_gaps':gaps,'candidate':target}
  out['best_new_candidate']=target
  event=nxt;tie_edges=[z['ids'] for z in gaps if not z['right_gap']];left,right=nextleft,nextright
 else:out['status']='finished_event_cap'
except Exception as e:out.update({'status':'stopped_exception','error':repr(e)})
finally:
 signal.alarm(0)
 out.update({'event_classifications':len(out['events']),'symbolic_configurations':configs,'local_nodes':nodes,'cpu_seconds':time.process_time()-CPU,'wall_seconds':time.monotonic()-WALL,'input_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [PREV,IFACE,ARITH,HERE/'registration.json',Path(__file__)]}})
 with OUT.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
 print(json.dumps({k:out[k] for k in ['status','pid','event_classifications','symbolic_configurations','local_nodes','cpu_seconds','wall_seconds']}))
 if 'error' in out:print(out['error']);sys.exit(1)
