#!/usr/bin/env python3
"""Finite exact endpoint/coverage preflight; deliberately has no union routine."""
import os, sys, json, hashlib, importlib.util, time, resource
from pathlib import Path
from fractions import Fraction as Q
from math import prod
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
ROOT=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6')
BASE=ROOT/'cube-dim-6-bounds'
SRC=BASE/'lower/round008/spacing_chamber/chamber.py'
IFACE=SRC.parent/'chamber_interface.json'
OLD=SRC.parent/'chamber_results.json'
CAND=BASE/'audit/lower008/parametric_union/candidate_results.json'
CPU=time.process_time(); WALL=time.monotonic()
resource.setrlimit(resource.RLIMIT_CPU,(20,21))
spec=importlib.util.spec_from_file_location('old_chamber',SRC)
C=importlib.util.module_from_spec(spec); spec.loader.exec_module(C)
F,u,v,ZERO,ONE=C.F,C.u,C.v,C.ZERO,C.ONE
configs=0; coverage_nodes=0
out={'status':'running','pid':os.getpid(),'full_union_calls':0,'new_C6_bound':False,'symbolic_configurations':0}
def budget():
 assert time.process_time()-CPU<19.5 and time.monotonic()-WALL<30
 assert configs<=12 and coverage_nodes<=1000
C.budget=budget

def enc(x): return C.coeff(x)
def view(x):
 lo,hi=C.A.bounds(x.a)
 return {'coefficients':enc(x),'exact_interval':[str(lo),str(hi)],'display':float((lo+hi)/2)}
def fixed(r): return r['base_box'][:2]+r['base_box'][4:]
try:
 d=json.loads(IFACE.read_text()); old=json.loads(OLD.read_text()); candidate=json.loads(CAND.read_text())
 current=u**2*v**2/3+u**2*v**4/6
 start=C.field(old['current_spacing2_coefficients'])
 spacing=(ONE+u**4/2,current,*((ONE+u)/2 for _ in range(4)))
 raw=[]
 for ell0,b in zip(d['box_lengths'],d['boxes']):
  ell=tuple(ell0); r=C.root_radius(prod(ell)); box=[]
  for j,pair in enumerate(b):
   vals=[C.field(d['endpoints'][eid]['c0'])+C.field(d['endpoints'][eid]['c1']) for eid in pair]
   box.extend(vals)
   q=-(ell[j]//2) if ell[j]%2 else 1-ell[j]//2
   for x in vals:
    for z in (q,q+ell[j]-1): assert -r<=x-z*spacing[j]<=r
  assert all(box[2*j]<box[2*j+1] for j in range(6))
  raw.append({'ell':ell,'radius':r,'base_box':tuple(box)})
 assert len(raw)==127
 funcs=C.make_functions(raw)
 events,evpairs=C.event_inventory(funcs,start)
 assert len(events)==217 and events[0]==current
 future=[x for x in events if current<x]; nxt=future[0]
 out.update({'current':view(current),'old_start':view(start),'raw_boxes':len(raw),'functions':len(funcs),'old_future_locations':len(events),'locations_equal_current':sum(x==current for x in events),'locations_strictly_right':len(future),'next_raw_location':view(nxt),'gap':view(nxt-current),'current_candidate_interval':candidate['candidate_ratio_interval']})
 states={}
 for name,event in [('current',current),('next',nxt)]:
  states[name]={}
  for side,label in [(-1,'left'),(1,'right')]:
   configs+=1; budget(); sel,order=C.state_at_event(raw,funcs,event,side)
   states[name][label]=(sel,order)
  out[name+'_states']={label:{'positive_boxes':len(sel),'effective_order':list(order),'selections':{str(k):list(z) for k,z in sel.items()}} for label,(sel,order) in states[name].items()}
  out[name+'_crossings']=[{'ids':[i,j],'functions':[{'c':enc(funcs[k]['c']),'slope':str(funcs[k]['m']),'labels':funcs[k]['labels']} for k in (i,j)]} for i,j in evpairs[event]]
 # Exact order chamber closure check from right germ to next left germ.
 sel,order=states['current']['right']
 assert (sel,order)==states['next']['left']
 order_checks=[]
 for i,j in zip(order,order[1:]):
  f,g=funcs[i],funcs[j]
  a=g['c']-f['c']+(g['m']-f['m'])*current
  b=g['c']-f['c']+(g['m']-f['m'])*nxt
  assert ZERO<=a and ZERO<=b and (ZERO<a or ZERO<b)
  order_checks.append({'ids':[i,j],'left_gap':enc(a),'right_gap':enc(b)})
 out['closed_chamber_order_checks']=order_checks
 # At each crossing, changed cross sections A intersect B can be hidden by
 # witnesses strictly spanning its x2-height. Only single-box containment
 # is used here: a successful test is an exact certificate, never a union.
 for name,event in [('current',current),('next',nxt)]:
  changes=[]
  left,leftorder=states[name]['left']; right,rightorder=states[name]['right']
  for i,j in evpairs[event]:
   if i not in leftorder or j not in leftorder or i not in rightorder or j not in rightorder: continue
   height=funcs[i]['c']+funcs[i]['m']*event
   group_i=[k for k,pair in left.items() if i in pair]
   group_j=[k for k,pair in left.items() if j in pair]
   background=[]
   for k,pair in left.items():
    if i in pair or j in pair: continue
    lo=funcs[pair[0]]['c']+funcs[pair[0]]['m']*event
    hi=funcs[pair[1]]['c']+funcs[pair[1]]['m']*event
    if lo<height<hi: background.append(k)
   pairs=[]
   for a in group_i:
    for b in group_j:
     coverage_nodes+=1; budget()
     A,B=fixed(raw[a]),fixed(raw[b]); I=tuple(z for q in range(5) for z in (max(A[2*q],B[2*q]),min(A[2*q+1],B[2*q+1])))
     if any(I[2*q+1]<=I[2*q] for q in range(5)):
      q=next(q for q in range(5) if I[2*q+1]<=I[2*q])
      pairs.append({'a':a,'b':b,'empty':True,'original_axis':[0,2,3,4,5][q]+1,'separation':view(I[2*q]-I[2*q+1])}); continue
     covers=[k for k in background if all(fixed(raw[k])[2*q]<=I[2*q] and I[2*q+1]<=fixed(raw[k])[2*q+1] for q in range(5))]
     pairs.append({'a':a,'b':b,'empty':False,'single_cover':covers[0] if covers else None})
   changes.append({'ids':[i,j],'background':background,'pair_certificates':pairs,'all_pair_intersections_hidden':all(p.get('empty') or p.get('single_cover') is not None for p in pairs)})
  out[name+'_local_coverage']=changes
 # Derivative source is a frozen exact certificate; no new volume computed.
 slope=C.field(old['spacing_derivative_coefficients'])
 delta=C.field(d['spacing2']['c1'])
 polynomial=candidate['ratio_polynomial_coefficients']
 assert len(polynomial)==2
 assert C.field(polynomial[1])==slope*delta
 assert sum((C.field(x) for x in polynomial),ZERO)==C.field(candidate['candidate_ratio_coefficients'])
 out['inherited_slope_independently_matches_E008_polynomial']=True
 out['inherited_left_slope']=view(slope)
 out['right_slope_equals_left_by_local_coverage']=bool(out['current_local_coverage']) and all(x['all_pair_intersections_hidden'] for x in out['current_local_coverage'])
 if out['right_slope_equals_left_by_local_coverage']:
  out['right_slope']=view(slope)
  assert ZERO<slope
  predicted=C.field(candidate['candidate_ratio_coefficients'])+slope*(nxt-current)
  old_upper=Q(candidate['candidate_ratio_interval'][1])
  predicted_lower=C.A.bounds(predicted.a)[0]
  out['next_affine_ideal_candidate']=view(predicted)
  out['exact_above_complete_current_ideal_upper']=old_upper<predicted_lower
  out['candidate_margin_over_complete_current_upper']=str(predicted_lower-old_upper)
 out['status']='finished'
except Exception as e:
 out['status']='stopped'; out['error']=repr(e)
finally:
 out.update({'symbolic_configurations':configs,'local_coverage_nodes':coverage_nodes,'cpu_seconds':time.process_time()-CPU,'wall_seconds':time.monotonic()-WALL,'input_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [SRC,IFACE,OLD,CAND,HERE/'registration.json',Path(__file__)]}})
 (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({k:out[k] for k in ['status','pid','symbolic_configurations','local_coverage_nodes','cpu_seconds','wall_seconds']}))
