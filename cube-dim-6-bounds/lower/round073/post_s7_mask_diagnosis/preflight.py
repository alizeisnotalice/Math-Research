#!/usr/bin/env python3
"""Frozen-s7 local mask diagnosis and one adjacent endpoint. No union."""
from pathlib import Path
from fractions import Fraction as Q
import json,hashlib,importlib.util,sys,time,resource,signal,os
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;OUT=HERE/'results.json'
if OUT.exists():raise FileExistsError('Refuse overwrite '+str(OUT))
resource.setrlimit(resource.RLIMIT_CPU,(20,21))
def wallcap(*_):raise TimeoutError('180-second total program wall cap including cleanup')
signal.signal(signal.SIGALRM,wallcap);signal.alarm(180)
CPU,WALL=time.process_time(),time.monotonic()
LOWER=HERE.parents[1];PREV=LOWER/'round070/first_overlap_cover/results.json'
BASE=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds')
ARITH=BASE/'lower/round001_rectangular_certificate.py';IFACE=BASE/'lower/round008/spacing_chamber/chamber_interface.json'
configs,nodes=0,0
out={'status':'running','pid':os.getpid(),'full_union_calls':0,'new_C6_bound':False,'automatic_reruns':0,'existing_cover_candidates_tested':0}
def budget():
 if time.process_time()-CPU>=19.5:raise TimeoutError('CPU prelimit')
 if time.monotonic()-WALL>=179:raise TimeoutError('wall prelimit with cleanup reserve')
 if configs>12 or nodes>1000:raise RuntimeError('configuration/spatial cap')
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
 d=json.loads(PREV.read_text());iface=json.loads(IFACE.read_text());defs=d['complete_functions_reused']
 funcs=[{'c':field(f['c']),'m':Q(f['m'])} for f in defs]
 event=field(d['next_s7']['coefficients']);baseline=field(d['candidate']['coefficients']);g=field(d['current_slope']['coefficients'])
 oldupper=Q(d['candidate']['exact_interval'][1]);left,right=d['next_states']['left'],d['next_states']['right']
 def val(i,t):return funcs[i]['c']+funcs[i]['m']*t
 fixed=[tuple(field(iface['endpoints'][eid]['c0']) for q in [0,2,3,4,5] for eid in b[q]) for b in iface['boxes']]
 out.update({'current_s7':view(event),'baseline_complete_E070_interval':d['candidate']['exact_interval'],'current_slope':d['current_slope'],'complete_functions_reused':defs,'current_left_state_reused':left,'current_right_state_reused':right,'tie_edges':d['next_tie_edges'],'tie_diagnostics':[]})
 configs+=1;budget();paid=left['selections']==right['selections'] and left['positive_boxes']==right['positive_boxes']
 if not paid:out['stop_reason']='current endpoint selection or positive-width branch changed'
 for i,j in d['next_tie_edges'] if paid else []:
  rec={'ids':[i,j],'all_labels':[],'pair_diagnostics':[]};out['tie_diagnostics'].append(rec)
  for fid in [i,j]:
   for label in defs[fid]['labels']:
    if 'raw_index' not in label:continue
    nodes+=1;budget();k=label['raw_index'];sel=left['selections'].get(str(k));rec['all_labels'].append({'function_id':fid,'label':label,'selected_endpoints':sel,'function_is_effectively_selected':sel is not None and fid in sel})
  gi=[int(k) for k,p in left['selections'].items() if i in p];gj=[int(k) for k,p in left['selections'].items() if j in p];rec['effective_groups']=[gi,gj]
  for a in gi:
   for b in gj:
    nodes+=1;budget();B,C=fixed[a],fixed[b];I=[(max(B[2*q],C[2*q]),min(B[2*q+1],C[2*q+1])) for q in range(5)]
    p={'raw_indices':[a,b],'lengths':[iface['box_lengths'][a],iface['box_lengths'][b]],'intersection':[[enc(lo),enc(hi)] for lo,hi in I]};rec['pair_diagnostics'].append(p)
    sep=[q for q,(lo,hi) in enumerate(I) if hi<lo]
    if sep:
     q=sep[0];lo,hi=I[q];p.update({'classification':'strictly_disjoint','separating_original_axis':[1,3,4,5,6][q],'positive_gap':view(lo-hi)});continue
    # Only for a nonempty/touching projection, test existing single boxes.
    cover=None
    for ck,pair in left['selections'].items():
     nodes+=1;budget();out['existing_cover_candidates_tested']+=1;c=int(ck);D=fixed[c]
     if not all(D[2*q]<=lo and hi<=D[2*q+1] for q,(lo,hi) in enumerate(I)):continue
     h=val(i,event);assert h==val(j,event)
     lower=h-val(pair[0],event);upper=val(pair[1],event)-h
     if not (ZERO<lower and ZERO<upper):continue
     llo=A.bounds(lower.a)[0];ulo=A.bounds(upper.a)[0]
     # Rational event-neighborhood radius from Lipschitz slopes. Actual
     # min/max endpoints have slopes at most 2 in magnitude, and movers <=2.
     delta=min(llo,ulo)/16
     cover={'raw_index':c,'lengths':iface['box_lengths'][c],'normal_collar_lower':view(lower),'normal_collar_upper':view(upper),'explicit_rational_parameter_neighborhood':str(delta),'projection_nonnegative_face_gaps':[{'low':enc(lo-D[2*q]),'high':enc(D[2*q+1]-hi)} for q,(lo,hi) in enumerate(I)]};break
    if cover:p.update({'classification':'single_existing_box_cover','cover':cover})
    else:
     p['classification']='uncovered_by_all_tested_single_existing_boxes';paid=False;out['stop_reason']='needs a new finite 5D multi-box cover or mask integral; not a true-kink proof';break
   if not paid:break
  if not paid:break
 out['all_effective_strips_paid']=paid
 if not paid:out['status']='stopped_unpaid'
 else:
  out['seam_hidden_and_g_preserved']=True
  order=list(d['right_complete_order']);edges={tuple(z) for z in d['next_tie_edges']};blocks=[];block=[order[0]]
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
   f,h=funcs[i],funcs[j];l=h['c']-f['c']+(h['m']-f['m'])*event;r=h['c']-f['c']+(h['m']-f['m'])*nxt
   assert ZERO<=l and ZERO<=r and (ZERO<l or ZERO<r);gaps.append({'ids':[i,j],'left_gap':enc(l),'right_gap':enc(r)})
  moving={z['raw_index']:i for i,f in enumerate(defs) for z in f['labels'] if 'raw_index' in z}
  def state(t,side):
   global configs
   configs+=1;budget();keys=[(val(i,t),side*funcs[i]['m']) for i in range(37)];sels={};active=set()
   for k,ell in enumerate(iface['box_lengths']):
    lo=max((0,) if ell[1]%2 else (0,moving[k]),key=lambda i:keys[i]);hi=min((1,moving[k]) if ell[1]%2 else (1,),key=lambda i:keys[i])
    if keys[lo]<keys[hi]:sels[str(k)]=[lo,hi];active.update((lo,hi))
   return {'positive_boxes':len(sels),'effective_order':sorted(active,key=lambda i:keys[i]),'selections':sels}
  nl,nr=state(nxt,-1),state(nxt,1);assert nl==right
  assert ZERO<g;candidate=baseline+g*(nxt-event);cl=A.bounds(candidate.a)[0];assert oldupper<cl
  out.update({'next_s8':view(nxt),'delta':view(nxt-event),'right_complete_order':order,'new_adjacent_checks':36,'new_future_adjacent_candidates':[{'ids':[i,j],'location':view(t)} for t,i,j in candidates],'complete_order_gap_certificate':gaps,'next_tie_edges':[z['ids'] for z in gaps if not z['right_gap']],'next_states':{'left':nl,'right':nr},'candidate':view(candidate),'candidate_lower_strictly_above_complete_E070_upper':True,'candidate_margin_over_complete_E070_upper':str(cl-oldupper),'status':'finished'})
except Exception as e:out.update({'status':'stopped_exception','error':repr(e)})
finally:
 out.update({'symbolic_configurations':configs,'spatial_nodes':nodes,'cpu_seconds':time.process_time()-CPU,'wall_seconds_before_exclusive_save':time.monotonic()-WALL,'input_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [PREV,IFACE,ARITH,HERE/'registration.json',Path(__file__)]}})
 with OUT.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
 print(json.dumps({k:out[k] for k in ['status','pid','symbolic_configurations','spatial_nodes','cpu_seconds','wall_seconds_before_exclusive_save']}))
 signal.alarm(0)
 if 'error' in out:print(out['error']);sys.exit(1)
