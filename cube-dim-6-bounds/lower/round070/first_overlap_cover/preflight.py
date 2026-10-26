#!/usr/bin/env python3
"""One finite seam coverage proof and one next endpoint; no union."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations
from math import prod
import json,hashlib,importlib.util,sys,time,resource,signal,os
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;OUT=HERE/'results.json'
if OUT.exists():raise FileExistsError('Refuse overwrite '+str(OUT))
resource.setrlimit(resource.RLIMIT_CPU,(20,21))
def wallcap(*_):raise TimeoutError('180 second wall cap')
signal.signal(signal.SIGALRM,wallcap);signal.alarm(180)
CPU,WALL=time.process_time(),time.monotonic()
LOWER=HERE.parents[1];PREV=LOWER/'round068/adjacent_hidden_chain/results.json'
BASE=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds')
ARITH=BASE/'lower/round001_rectangular_certificate.py';IFACE=BASE/'lower/round008/spacing_chamber/chamber_interface.json'
configs,nodes=0,0
out={'status':'running','pid':os.getpid(),'full_union_calls':0,'new_C6_bound':False,'automatic_reruns':0}
def budget():
 if time.process_time()-CPU>=19.5:raise TimeoutError('CPU prelimit')
 if time.monotonic()-WALL>=180:raise TimeoutError('wall cap')
 if configs>12 or nodes>1500:raise RuntimeError('configuration/local node cap')
try:
 spec=importlib.util.spec_from_file_location('accepted_field',ARITH);A=importlib.util.module_from_spec(spec);spec.loader.exec_module(A)
 A.D=10**80;A.RU,A.RV=A.floor_root(2*A.D**6,6),A.floor_root(3*A.D**6,6)
 A.UL,A.UH=Q(A.RU,A.D),Q(A.RU+1,A.D);A.VL,A.VH=Q(A.RV,A.D),Q(A.RV+1,A.D);A.bounds.cache_clear()
 F,ZERO,u=A.F,A.F(0),A.u
 def field(z):
  a=[Q(0)]*36
  for ij,x in z.items():i,j=map(int,ij.split(','));a[6*i+j]=Q(x)
  return F(a=a)
 def enc(x):return {f'{i},{j}':str(x.a[6*i+j]) for i in range(6) for j in range(6) if x.a[6*i+j]}
 def view(x):
  lo,hi=A.bounds(x.a);return {'coefficients':enc(x),'exact_interval':[str(lo),str(hi)],'display':float((lo+hi)/2)}
 d=json.loads(PREV.read_text());iface=json.loads(IFACE.read_text());seam=d['events'][-1]
 defs=d['complete_functions_reused'];funcs=[{'c':field(f['c']),'m':Q(f['m'])} for f in defs]
 event=field(d['best_new_candidate']['event_parameter']['coefficients']);baseline=field(d['best_new_candidate']['ratio']['coefficients']);g=field(d['slope']['coefficients'])
 oldupper=Q(d['best_new_candidate']['ratio']['exact_interval'][1]);left,right=seam['left_state'],seam['right_state']
 assert left['selections']==right['selections'] and left['positive_boxes']==right['positive_boxes']==127
 fixed=[tuple(field(iface['endpoints'][eid]['c0']) for j in [0,2,3,4,5] for eid in b[j]) for b in iface['boxes']]
 lengths=[tuple(z) for z in iface['box_lengths']];index={ell:i for i,ell in enumerate(lengths)}
 half_faces=[(F(1)+u**4/2)/2]+[(F(1)+u)/4]*4
 def val(i,t):return funcs[i]['c']+funcs[i]['m']*t
 configs+=1;budget()
 out.update({'current_s6':view(event),'baseline_complete_E068_interval':d['best_new_candidate']['ratio']['exact_interval'],'current_slope':d['slope'],'current_left_state_reused':left,'current_right_state_reused':right,'complete_functions_reused':defs,'tie_diagnostics':[]})
 all_paid=True
 for i,j in seam['tie_edges']:
  gi=[int(k) for k,p in left['selections'].items() if i in p];gj=[int(k) for k,p in left['selections'].items() if j in p]
  rec={'ids':[i,j],'effective_groups':[gi,gj]};out['tie_diagnostics'].append(rec)
  if not gi and not gj:
   checks=[]
   for fid in [i,j]:
    for label in defs[fid]['labels']:
     if 'raw_index' not in label:continue
     nodes+=1;budget();k=label['raw_index'];selected=left['selections'][str(k)]
     assert selected[1]==1 and fid not in selected
     margin=val(fid,event)-event/2;assert ZERO<margin
     checks.append({'function':fid,'raw_index':k,'selected_endpoints':selected,'unused_upper_minus_halfcell':view(margin)})
   rec.update({'classification':'both_unclipped_functions_inactive_strictly_above_halfcell','checks':checks});continue
  if not gi or not gj:
   rec['classification']='unpaid_one_group_active';all_paid=False;break
  covers=[];rec['pair_covers']=covers
  for a in gi:
   for b in gj:
    nodes+=1;budget();B,C=fixed[a],fixed[b]
    I=[(max(B[2*q],C[2*q]),min(B[2*q+1],C[2*q+1])) for q in range(5)]
    # This seam's suggested candidate is generated from actual A tail lengths.
    ell=(1,2,*[2 if z==4 else 1 for z in lengths[a][2:]])
    p={'raw_indices':[a,b],'lengths':[lengths[a],lengths[b]],'intersection':[[enc(lo),enc(hi)] for lo,hi in I],'cover_lengths':ell};covers.append(p)
    if ell not in index:p['covered']=False;all_paid=False;break
    c=index[ell];D=fixed[c];p['cover_raw_index']=c
    assert prod(ell)==8 and (2*(u**3/2))**6==F(8)
    face_checks=[];paid=True
    for q,(lo,hi) in enumerate(I):
     assert lo<hi
     lowgap=lo-D[2*q];highgap=D[2*q+1]-hi
     if lowgap<ZERO or highgap<ZERO:paid=False;break
     lowrelative=lowgap==ZERO and lo==ZERO and D[2*q]==ZERO
     highrelative=highgap==ZERO and hi==half_faces[q] and D[2*q+1]==half_faces[q]
     if not (ZERO<lowgap or lowrelative) or not (ZERO<highgap or highrelative):paid=False;break
     face_checks.append({'original_axis':[1,3,4,5,6][q],'low_gap':view(lowgap),'high_gap':view(highgap),'low_equality_is_cell_face':lowrelative,'high_equality_is_cell_face':highrelative})
    selected=left['selections'][str(c)];height=val(i,event);assert height==val(j,event)
    normal_lo=val(selected[0],event);normal_hi=val(selected[1],event)
    collarlo=height-normal_lo;collarhi=normal_hi-height
    p.update({'cover_radius':enc(u**3/2),'cover_selected_normal_endpoints':selected,'fixed_face_checks':face_checks,'normal_height':view(height),'normal_collar_below':view(collarlo),'normal_collar_above':view(collarhi)})
    if not (ZERO<collarlo and ZERO<collarhi):paid=False
    p['covered_with_strict_relative_cell_collar']=paid
    if not paid:all_paid=False;break
   if not all_paid:break
  rec['classification']='covered_by_existing_admissible_box_with_collar' if all_paid else 'unpaid_cover';
  if not all_paid:break
 out['all_simultaneous_effective_strips_paid']=all_paid
 if not all_paid:out['status']='stopped_uncovered'
 else:
  # The strict normal collar plus fixed relative-cell cover masks every
  # swapped intersection for a two-sided neighborhood. Same affine slope.
  out['seam_hidden_and_slope_preserved']=True
  order=list(seam['left_complete_order']);edges={tuple(z) for z in seam['tie_edges']};blocks=[];block=[order[0]]
  for i,j in zip(order,order[1:]):
   if (i,j) in edges:block.append(j)
   else:blocks.append(block);block=[j]
  blocks.append(block);order=[i for b in blocks for i in sorted(b,key=lambda k:funcs[k]['m'])]
  candidates=[]
  for i,j in zip(order,order[1:]):
   nodes+=1;budget();f,h=funcs[i],funcs[j]
   if f['m']<=h['m']:continue
   t=(h['c']-f['c'])/(f['m']-h['m']);assert event<t;candidates.append((t,i,j))
  nxt=min(candidates,key=lambda z:z[0])[0];gaps=[]
  for i,j in zip(order,order[1:]):
   f,h=funcs[i],funcs[j];lo=h['c']-f['c']+(h['m']-f['m'])*event;hi=h['c']-f['c']+(h['m']-f['m'])*nxt
   assert ZERO<=lo and ZERO<=hi and (ZERO<lo or ZERO<hi);gaps.append({'ids':[i,j],'left_gap':enc(lo),'right_gap':enc(hi)})
  moving={z['raw_index']:i for i,f in enumerate(defs) for z in f['labels'] if 'raw_index' in z}
  def state(t,side):
   global configs
   configs+=1;budget();keys=[(val(i,t),side*funcs[i]['m']) for i in range(37)];sels={};active=set()
   for k,ell in enumerate(lengths):
    lo=max((0,) if ell[1]%2 else (0,moving[k]),key=lambda i:keys[i]);hi=min((1,moving[k]) if ell[1]%2 else (1,),key=lambda i:keys[i])
    if keys[lo]<keys[hi]:sels[str(k)]=[lo,hi];active.update((lo,hi))
   return {'positive_boxes':len(sels),'effective_order':sorted(active,key=lambda i:keys[i]),'selections':sels}
  nl,nr=state(nxt,-1),state(nxt,1);assert nl==right
  assert ZERO<g;candidate=baseline+g*(nxt-event);cl=A.bounds(candidate.a)[0];assert oldupper<cl
  out.update({'next_s7':view(nxt),'delta':view(nxt-event),'right_complete_order':order,'new_adjacent_checks':36,'new_future_adjacent_candidates':[{'ids':[i,j],'location':view(t)} for t,i,j in candidates],'complete_order_gap_certificate':gaps,'next_tie_edges':[z['ids'] for z in gaps if not z['right_gap']],'next_states':{'left':nl,'right':nr},'candidate':view(candidate),'candidate_lower_strictly_above_complete_E068_upper':True,'candidate_margin_over_complete_E068_upper':str(cl-oldupper),'status':'finished'})
except Exception as e:out.update({'status':'stopped_exception','error':repr(e)})
finally:
 signal.alarm(0)
 out.update({'symbolic_configurations':configs,'local_nodes':nodes,'cpu_seconds':time.process_time()-CPU,'wall_seconds':time.monotonic()-WALL,'input_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [PREV,IFACE,ARITH,HERE/'registration.json',Path(__file__)]}})
 with OUT.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
 print(json.dumps({k:out[k] for k in ['status','pid','symbolic_configurations','local_nodes','cpu_seconds','wall_seconds']}))
 if 'error' in out:print(out['error']);sys.exit(1)
