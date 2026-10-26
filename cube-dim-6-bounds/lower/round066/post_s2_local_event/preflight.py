#!/usr/bin/env python3
"""Saved-064 local seam plus one new adjacency chamber; no union code."""
from pathlib import Path
from fractions import Fraction as Q
import json,hashlib,importlib.util,sys,time,resource,signal,os
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;OUT=HERE/'results.json'
if OUT.exists():raise FileExistsError('Refuse to overwrite '+str(OUT))
resource.setrlimit(resource.RLIMIT_CPU,(20,21))
def wallcap(*_):raise TimeoutError('180 second wall limit')
signal.signal(signal.SIGALRM,wallcap);signal.alarm(180)
CPU,WALL=time.process_time(),time.monotonic()
LOWER=HERE.parents[1]
PREV=LOWER/'round064/a2_one_more_event/results.json'
BASE=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds')
ARITH=BASE/'lower/round001_rectangular_certificate.py'
IFACE=BASE/'lower/round008/spacing_chamber/chamber_interface.json'
configs,nodes=0,0
out={'status':'running','pid':os.getpid(),'full_union_calls':0,'new_C6_bound':False,'automatic_reruns':0}
def budget():
 assert time.process_time()-CPU<19.5 and time.monotonic()-WALL<180
 assert configs<=12 and nodes<=1000
try:
 spec=importlib.util.spec_from_file_location('accepted_field',ARITH)
 A=importlib.util.module_from_spec(spec);spec.loader.exec_module(A)
 A.D=10**80
 A.RU,A.RV=A.floor_root(2*A.D**6,6),A.floor_root(3*A.D**6,6)
 A.UL,A.UH=Q(A.RU,A.D),Q(A.RU+1,A.D)
 A.VL,A.VH=Q(A.RV,A.D),Q(A.RV+1,A.D)
 A.bounds.cache_clear()
 F,ZERO=A.F,A.F(0)
 def field(z):
  a=[Q(0)]*36
  for ij,x in z.items():i,j=map(int,ij.split(','));a[6*i+j]=Q(x)
  return F(a=a)
 def enc(x):return {f'{i},{j}':str(x.a[6*i+j]) for i in range(6) for j in range(6) if x.a[6*i+j]}
 def view(x):
  lo,hi=A.bounds(x.a)
  return {'coefficients':enc(x),'exact_interval':[str(lo),str(hi)],'display':float((lo+hi)/2)}
 d=json.loads(PREV.read_text());iface=json.loads(IFACE.read_text())
 funcs=[{'c':field(f['c']),'m':Q(f['m']),'labels':f['labels']} for f in d['complete_functions']]
 s2=field(d['s2']['coefficients']);g=field(d['slope']['coefficients'])
 baseline=field(d['candidate']['coefficients'])
 # Keep fixed endpoint coordinates from the accepted 127-box interface.
 fixed=[]
 for b in iface['boxes']:
  fixed.append(tuple(field(iface['endpoints'][eid]['c0']) for axis in [0,2,3,4,5] for eid in b[axis]))
 swapped=d['s2_crossing_adjacent_ids'];assert swapped==[[19,14]]
 groups=[[x['raw_index'] for x in funcs[i]['labels'] if 'raw_index' in x] for i in swapped[0]]
 assert len(groups[0])==1 and len(groups[1])==18
 pairs=[]
 for a in groups[0]:
  for b in groups[1]:
   nodes+=1;budget();B,C=fixed[a],fixed[b]
   intersections=[(max(B[2*j],C[2*j]),min(B[2*j+1],C[2*j+1])) for j in range(5)]
   separated=[j for j,(lo,hi) in enumerate(intersections) if hi<=lo]
   if not separated:
    pairs.append({'a':a,'b':b,'nonempty':True,'intersection':[[enc(lo),enc(hi)] for lo,hi in intersections]})
   else:
    j=separated[0];lo,hi=intersections[j]
    assert ZERO<lo-hi
    pairs.append({'a':a,'b':b,'nonempty':False,'separating_original_axis':[1,3,4,5,6][j],'positive_gap':view(lo-hi),'lengths':[iface['box_lengths'][a],iface['box_lengths'][b]]})
 out['s2_swapped_projection_pairs']=pairs
 if any(p['nonempty'] for p in pairs):
  out['status']='stopped_nonempty_projection'
 else:
  out['s2_hidden_by_fixed_projection_disjointness']=True
  out['right_slope_preserved_by_accepted_inclusion_exclusion_lemma']=True
  # Recover s2 full right order from saved left order and ONLY saved seam tie.
  configs+=1;budget();order=list(d['all_function_order_s1_right'])
  for i,j in swapped:
   a,b=order.index(i),order.index(j);assert b==a+1
   assert funcs[i]['m']>funcs[j]['m'];order[a],order[b]=order[b],order[a]
  candidates=[]
  for i,j in zip(order,order[1:]):
   nodes+=1;budget();f,h=funcs[i],funcs[j]
   if f['m']<=h['m']:continue
   t=(h['c']-f['c'])/(f['m']-h['m'])
   assert s2<t;candidates.append((t,i,j))
  s3=min(candidates,key=lambda z:z[0])[0]
  out.update({'s2':view(s2),'s3':view(s3),'delta':view(s3-s2),'slope':view(g),'baseline_complete_064_interval':d['candidate']['exact_interval'],'full_functions_reused_from_064':d['complete_functions'],'complete_s2_right_order':order,'old_inventory_all_pairs_comparisons':0,'next_adjacent_pair_checks':36,'future_adjacent_candidates':[{'ids':[i,j],'location':view(t)} for t,i,j in candidates]})
  gaps=[]
  for i,j in zip(order,order[1:]):
   f,h=funcs[i],funcs[j]
   l=h['c']-f['c']+(h['m']-f['m'])*s2
   r=h['c']-f['c']+(h['m']-f['m'])*s3
   assert ZERO<=l and ZERO<=r and (ZERO<l or ZERO<r)
   gaps.append({'ids':[i,j],'s2_gap':enc(l),'s3_gap':enc(r)})
  out['complete_order_gap_certificate']=gaps
  out['s3_tying_adjacent_ids']=[z['ids'] for z in gaps if not z['s3_gap']]
  # Select min/max endpoints using the complete frozen function definitions.
  moving={z['raw_index']:i for i,f in enumerate(funcs) for z in f['labels'] if 'raw_index' in z}
  def state(event,side):
   keys=[(f['c']+f['m']*event,side*f['m']) for f in funcs]
   sels={};active=set()
   for k,ell in enumerate(iface['box_lengths']):
    lo=max((0,) if ell[1]%2 else (0,moving[k]),key=lambda i:keys[i])
    hi=min((1,moving[k]) if ell[1]%2 else (1,),key=lambda i:keys[i])
    if keys[lo]<keys[hi]:sels[str(k)]=[lo,hi];active.update((lo,hi))
   return {'positive_boxes':len(sels),'effective_order':sorted(active,key=lambda i:keys[i]),'selections':sels}
  out['s2_right_effective_state_reused_from_064']=d['s2_states']['right']
  out['s3_states']={}
  for side,label in [(-1,'left'),(1,'right')]:
   configs+=1;budget();out['s3_states'][label]=state(s3,side)
  assert out['s3_states']['left']==d['s2_states']['right']
  assert ZERO<g
  candidate=baseline+g*(s3-s2)
  cl=A.bounds(candidate.a)[0];oldhi=Q(d['candidate']['exact_interval'][1]);assert oldhi<cl
  out.update({'candidate':view(candidate),'candidate_lower_strictly_above_complete_064_upper':True,'candidate_margin_over_complete_064_upper':str(cl-oldhi),'proof_status':'single exact ideal candidate, same mechanism; finite-L1 and independent review of new candidate not done','status':'finished'})
except Exception as e:out.update({'status':'stopped','error':repr(e)})
finally:
 signal.alarm(0)
 out.update({'symbolic_configurations':configs,'local_nodes':nodes,'cpu_seconds':time.process_time()-CPU,'wall_seconds':time.monotonic()-WALL,'input_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [PREV,IFACE,ARITH,HERE/'registration.json',Path(__file__)]}})
 with OUT.open('x') as fp:json.dump(out,fp,indent=2);fp.write('\n')
 print(json.dumps({k:out[k] for k in ['status','pid','symbolic_configurations','local_nodes','cpu_seconds','wall_seconds']}))
 if out['status']!='finished':print(out.get('error',out['status']));sys.exit(1)
