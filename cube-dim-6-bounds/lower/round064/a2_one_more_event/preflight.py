#!/usr/bin/env python3
"""One nearest-event continuation; every output is exclusive-create."""
from pathlib import Path
from fractions import Fraction as Q
from math import prod
import importlib.util,json,hashlib,sys,resource,time,signal,os
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
OUT=HERE/'results.json'
if OUT.exists(): raise FileExistsError('refuse to overwrite frozen/new output: '+str(OUT))
resource.setrlimit(resource.RLIMIT_CPU,(20,21))
signal.signal(signal.SIGALRM,lambda *_: (_ for _ in ()).throw(TimeoutError('180 second wall cap')))
signal.alarm(180)
CPU,WALL=time.process_time(),time.monotonic()
LOWER=HERE.parents[1]
PREV=LOWER/'round058/a2_next_event_preflight/results.json'
ROOT=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds')
SRC=ROOT/'lower/round008/spacing_chamber/chamber.py'
IFACE=SRC.parent/'chamber_interface.json'
ARITH=ROOT/'lower/round001_rectangular_certificate.py'
configs,nodes=0,0
out={'status':'running','pid':os.getpid(),'full_union_calls':0,'new_C6_bound':False,'automatic_reruns':0}
def budget():
 assert time.process_time()-CPU<19.5 and time.monotonic()-WALL<180
 assert configs<=12 and nodes<=1000
try:
 spec=importlib.util.spec_from_file_location('accepted_chamber',SRC)
 C=importlib.util.module_from_spec(spec); spec.loader.exec_module(C)
 C.budget=budget
 F,ZERO=C.F,C.ZERO
 d=json.loads(PREV.read_text()); iface=json.loads(IFACE.read_text())
 def field(z):return C.field(z)
 def enc(x):return C.coeff(x)
 def view(x):
  lo,hi=C.A.bounds(x.a)
  return {'coefficients':enc(x),'exact_interval':[str(lo),str(hi)],'display':float((lo+hi)/2)}
 s1=field(d['next_raw_location']['coefficients'])
 baseline=field(d['next_affine_ideal_candidate']['coefficients'])
 g=field(d['right_slope']['coefficients'])
 # Saved 058 result has only the two examined crossing definitions, not
 # all 37 functions/event locations. Reconstruct function definitions only
 # from the frozen 127 lengths; no old geometry/source checks or 666 pairs.
 raw=[{'ell':tuple(ell),'radius':C.root_radius(prod(ell))} for ell in iface['box_lengths']]
 funcs=C.make_functions(raw)
 assert len(raw)==127 and len(funcs)==37
 def key(i,event,side):
  f=funcs[i]; return (f['c']+f['m']*event,side*f['m'])
 configs+=1; budget()
 all_order=sorted(range(len(funcs)),key=lambda i:key(i,s1,1))
 candidates=[]
 for i,j in zip(all_order,all_order[1:]):
  nodes+=1;budget()
  f,h=funcs[i],funcs[j]
  if f['m']<=h['m']: continue
  crossing=(h['c']-f['c'])/(f['m']-h['m'])
  assert s1<crossing
  candidates.append((crossing,i,j))
 s2=min(candidates,key=lambda z:z[0])[0]
 out.update({'s1':view(s1),'s2':view(s2),'delta':view(s2-s1),'slope':view(g),'baseline_058_ideal_interval':d['next_affine_ideal_candidate']['exact_interval'],'function_inventory_reconstructed':True,'old_all_pairs_comparisons':0,'nearest_search_adjacent_pairs':36,'future_adjacent_candidates':len(candidates),'complete_functions':[{'id':i,'c':enc(f['c']),'m':str(f['m']),'labels':f['labels']} for i,f in enumerate(funcs)],'all_function_order_s1_right':all_order,'adjacent_future_candidates':[{'ids':[i,j],'location':view(x)} for x,i,j in candidates]})
 # All 36 gaps certify the complete affine-line order on the chamber.
 gaps=[]
 for i,j in zip(all_order,all_order[1:]):
  f,h=funcs[i],funcs[j]
  gap1=h['c']-f['c']+(h['m']-f['m'])*s1
  gap2=h['c']-f['c']+(h['m']-f['m'])*s2
  assert ZERO<=gap1 and ZERO<=gap2 and (ZERO<gap1 or ZERO<gap2)
  gaps.append({'ids':[i,j],'s1_gap':enc(gap1),'s2_gap':enc(gap2)})
 out['complete_order_gap_certificate']=gaps
 saved=d['next_states']['right']
 savedsel={int(k):tuple(pair) for k,pair in saved['selections'].items()}
 out['s1_right_state_reused_from_058']=saved
 states={}
 for side,label in [(-1,'left'),(1,'right')]:
  configs+=1;budget(); sel,order=C.state_at_event(raw,funcs,s2,side)
  states[label]={'positive_boxes':len(sel),'effective_order':list(order),'selections':{str(k):list(pair) for k,pair in sel.items()}}
  if side==-1:
   assert sel==savedsel and list(order)==saved['effective_order']
 out['s2_states']=states
 # Extract ties from the complete order certificate; no all-pair replay.
 out['s2_crossing_adjacent_ids']=[z['ids'] for z in gaps if not any(Q(c) for c in z['s2_gap'].values())]
 out['closed_chamber_effective_gap_certificate']=[]
 for i,j in zip(saved['effective_order'],saved['effective_order'][1:]):
  f,h=funcs[i],funcs[j]
  gap1=h['c']-f['c']+(h['m']-f['m'])*s1
  gap2=h['c']-f['c']+(h['m']-f['m'])*s2
  assert ZERO<=gap1 and ZERO<=gap2 and (ZERO<gap1 or ZERO<gap2)
  out['closed_chamber_effective_gap_certificate'].append({'ids':[i,j],'s1_gap':enc(gap1),'s2_gap':enc(gap2)})
 # Parent accepted E058's hidden s1 exchange, hence right derivative is g.
 # With only coordinate 2 moving and no order/branch event in the open
 # chamber, full union ratio stays affine. Continuity permits s2 endpoint.
 assert ZERO<g
 candidate=baseline+g*(s2-s1)
 lower=C.A.bounds(candidate.a)[0]
 old_upper=Q(d['next_affine_ideal_candidate']['exact_interval'][1])
 assert old_upper<lower
 out.update({'candidate':view(candidate),'candidate_lower_strictly_above_complete_058_upper':True,'candidate_margin_over_complete_058_upper':str(lower-old_upper),'proof_status':'exact ideal continuation candidate conditional on inherited accepted E058 right slope; no new finite-L1 certification','status':'finished'})
except Exception as e:
 out.update({'status':'stopped','error':repr(e)})
finally:
 signal.alarm(0)
 out.update({'symbolic_configurations':configs,'local_nodes':nodes,'cpu_seconds':time.process_time()-CPU,'wall_seconds':time.monotonic()-WALL,'input_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [PREV,SRC,IFACE,ARITH,HERE/'registration.json',Path(__file__)]}})
 with OUT.open('x') as fp:json.dump(out,fp,indent=2);fp.write('\n')
 print(json.dumps({k:out[k] for k in ['status','pid','symbolic_configurations','local_nodes','cpu_seconds','wall_seconds']}))
 if out['status']!='finished': print(out['error']);sys.exit(1)
