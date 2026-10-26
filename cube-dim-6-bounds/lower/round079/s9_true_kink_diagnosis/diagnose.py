#!/usr/bin/env python3
"""New five-dimensional background mask at s9; no six-dimensional union."""
from pathlib import Path
from fractions import Fraction as Q
from functools import lru_cache
import json,hashlib,importlib.util,sys,time,resource,signal,os
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;OUT=HERE/'results.json'
if OUT.exists():raise FileExistsError('Refuse overwrite '+str(OUT))
resource.setrlimit(resource.RLIMIT_CPU,(20,21))
def wallcap(*_):raise TimeoutError('180 second total wall cap including cleanup')
signal.signal(signal.SIGALRM,wallcap);signal.alarm(180)
CPU,WALL=time.process_time(),time.monotonic()
LOWER=HERE.parents[1];PREV=LOWER/'round075/post_s8_overlap/results.json'
BASE=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds')
ARITH=BASE/'lower/round001_rectangular_certificate.py';IFACE=BASE/'lower/round008/spacing_chamber/chamber_interface.json'
configs,nodes,unionnodes=0,0,0
out={'status':'running','pid':os.getpid(),'new_C6_bound':False,'full_6D_union_calls':0,'new_5D_union_calls':0,'automatic_reruns':0}
def budget():
 if time.process_time()-CPU>=19.5:raise TimeoutError('CPU prelimit')
 if time.monotonic()-WALL>=179:raise TimeoutError('wall prelimit')
 if configs>8 or nodes>1000:raise RuntimeError('configuration/combined node cap')
try:
 spec=importlib.util.spec_from_file_location('accepted_field',ARITH);A=importlib.util.module_from_spec(spec);spec.loader.exec_module(A)
 A.D=10**80;A.RU,A.RV=A.floor_root(2*A.D**6,6),A.floor_root(3*A.D**6,6)
 A.UL,A.UH=Q(A.RU,A.D),Q(A.RU+1,A.D);A.VL,A.VH=Q(A.RV,A.D),Q(A.RV+1,A.D);A.bounds.cache_clear()
 F,ZERO,HALF=A.F,A.F(0),A.F(Q(1,2))
 def field(z):
  a=[Q(0)]*36
  for ij,x in z.items():i,j=map(int,ij.split(','));a[6*i+j]=Q(x)
  return F(a=a)
 def enc(x):return {f'{i},{j}':str(x.a[6*i+j]) for i in range(6) for j in range(6) if x.a[6*i+j]}
 def view(x):
  lo,hi=A.bounds(x.a);return {'coefficients':enc(x),'exact_interval':[str(lo),str(hi)],'display':float((lo+hi)/2)}
 def boxenc(B):return [[enc(B[2*q]),enc(B[2*q+1])] for q in range(len(B)//2)]
 def vol(B):
  z=F(1)
  for q in range(len(B)//2):z=z*(B[2*q+1]-B[2*q])
  return z
 d=json.loads(PREV.read_text());iface=json.loads(IFACE.read_text());defs=d['complete_functions_reused']
 funcs=[{'c':field(f['c']),'m':Q(f['m'])} for f in defs]
 event=field(d['next_s9']['coefficients']);gleft=field(d['current_slope']['coefficients'])
 left,right=d['next_states']['left'],d['next_states']['right'];assert d['next_tie_edges']==[[7,2]]
 assert left['selections']==right['selections'] and left['positive_boxes']==right['positive_boxes']==127
 assert left['selections']['16']==[7,1] and left['selections']['0']==[0,2]
 assert funcs[7]['m']-funcs[2]['m']==1
 def val(i,t):return funcs[i]['c']+funcs[i]['m']*t
 height=val(2,event);assert height==HALF and val(7,event)==height
 fixed=[tuple(field(iface['endpoints'][eid]['c0']) for q in [0,2,3,4,5] for eid in b[q]) for b in iface['boxes']]
 P=tuple(z for q in range(5) for z in (max(fixed[0][2*q],fixed[16][2*q]),min(fixed[0][2*q+1],fixed[16][2*q+1])))
 assert all(P[2*q]==ZERO and P[2*q+1]==HALF for q in range(5))
 configs+=1;budget()
 out.update({'s9':view(event),'baseline_complete_E075_ideal_interval':d['candidate']['exact_interval'],'left_slope':view(gleft),'swapped_labels':{'lower':defs[7],'upper':defs[2]},'tie_edges':d['next_tie_edges'],'current_left_state_reused':left,'current_right_state_reused':right,'projection_P':boxenc(P),'projection_P_volume':view(vol(P)),'jump_coefficient_verified':str(funcs[7]['m']-funcs[2]['m']),'all_background_normal_classifications':[]})
 boxes=[];positive_indices=[];normalmargins=[]
 for sk,pair in left['selections'].items():
  k=int(sk)
  if k in (0,16):continue
  nodes+=1;budget();lo,hi=val(pair[0],event),val(pair[1],event)
  rec={'raw_index':k,'lengths':iface['box_lengths'][k],'normal_endpoint_ids':pair}
  if lo<height<hi:
   mlo,mhi=height-lo,hi-height;normalmargins.extend((A.bounds(mlo.a)[0],A.bounds(mhi.a)[0]))
   B=tuple(z for q in range(5) for z in (max(ZERO,fixed[k][2*q]),min(HALF,fixed[k][2*q+1])))
   positive=all(B[2*q]<B[2*q+1] for q in range(5))
   rec.update({'classification':'strictly_active_at_height','strict_normal_collar_below':view(mlo),'strict_normal_collar_above':view(mhi),'clipped_cross_section_in_P':boxenc(B),'positive_5D_cross_section':positive})
   if positive:boxes.append(B);positive_indices.append(k)
  elif hi<height:
   gap=height-hi;normalmargins.append(A.bounds(gap.a)[0]);rec.update({'classification':'strictly_below_height','strict_normal_distance':view(gap)})
  elif height<lo:
   gap=lo-height;normalmargins.append(A.bounds(gap.a)[0]);rec.update({'classification':'strictly_above_height','strict_normal_distance':view(gap)})
  else:raise RuntimeError('unpaid extra selected endpoint equals seam height')
  out['all_background_normal_classifications'].append(rec)
 assert len(out['all_background_normal_classifications'])==125
 assert normalmargins and min(normalmargins)>0
 out['explicit_rational_common_normal_neighborhood']=str(min(normalmargins)/16)
 unique=sorted(set(boxes));kept=[B for B in unique if not any(C!=B and all(C[2*q]<=B[2*q] and B[2*q+1]<=C[2*q+1] for q in range(5)) for C in unique)]
 out.update({'positive_background_indices':positive_indices,'positive_background_box_count':len(boxes),'unique_positive_background_boxes':len(unique),'retained_noncontained_background_boxes':len(kept),'retained_cross_sections':[boxenc(B) for B in kept]})
 # Explicit origin cube excluded by every background: each has some
 # strictly positive lower coordinate. This is evidence, not the full jump.
 thresholds=[max(B[::2]) for B in kept]
 assert thresholds and all(ZERO<x for x in thresholds)
 threshold=min(thresholds);epsilon=A.bounds(threshold.a)[0]/2
 assert 0<epsilon<Q(1,2)
 witnesses=[]
 for B in kept:
  q=max(range(5),key=lambda j:B[2*j]);assert F(epsilon)<B[2*q]
  witnesses.append({'box':boxenc(B),'separating_original_axis':[1,3,4,5,6][q],'box_lower':enc(B[2*q])})
 out['open_uncovered_origin_cube']={'interval_each_coordinate':['0',str(epsilon)],'axes':[1,3,4,5,6],'verified_separating_box_witnesses':witnesses,'not_used_as_full_uncovered_volume':True}
 # Safe fallback bounds are saved before any new exact union computation.
 maxsingle=max((vol(B) for B in kept),default=ZERO);sumvol=sum((vol(B) for B in kept),ZERO)
 out['pre_union_volume_bounds']={'background_lower':view(maxsingle),'background_upper':view(min(vol(P),sumvol))}
 @lru_cache(None)
 def new_section_union(Bs):
  global nodes,unionnodes
  nodes+=1;unionnodes+=1;budget()
  if not Bs:return ZERO
  if len(Bs[0])==2:
   ordered=sorted(Bs);lo,hi=ordered[0];total=ZERO
   for a,b in ordered[1:]:
    if hi<a:total=total+hi-lo;lo,hi=a,b
    else:hi=max(hi,b)
   return total+hi-lo
  cuts=sorted({x for B in Bs for x in B[:2]});total=ZERO
  for a,b in zip(cuts,cuts[1:]):
   child=tuple(sorted({B[2:] for B in Bs if B[0]<=a and b<=B[1]}))
   if child:total=total+(b-a)*new_section_union(child)
  return total
 out['new_5D_union_calls']=1
 background=new_section_union(tuple(sorted(kept)));uncovered=vol(P)-background
 assert ZERO<uncovered
 gright=gleft-64*uncovered
 out.update({'exact_background_5D_union':view(background),'exact_uncovered_5D_volume':view(uncovered),'exact_ratio_slope_jump':view(-64*uncovered),'right_slope':view(gright),'new_5D_union_cache_states':new_section_union.cache_info().currsize})
 # Only one NEW right chamber, no full6D union, using saved full left order.
 order=list(d['right_complete_order']);pos=order.index(7);assert order[pos+1]==2;order[pos],order[pos+1]=order[pos+1],order[pos]
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
 out.update({'right_chamber_next_endpoint':view(nxt),'right_complete_order':order,'right_chamber_full_order_gaps':gaps,'next_tie_edges':[z['ids'] for z in gaps if not z['right_gap']],'next_states':{'left':nl,'right':nr},'new_adjacent_checks':36})
 if gright<ZERO:
  assert ZERO<gleft
  out.update({'status':'finished_strict_local_maximum','strict_local_maximum_family':'same127 unit-mass common-radius standard central witness union, fixed a1 and all4 tail spacings, absolute a2 varied only','strict_local_maximum_parameter':enc(event),'left_comparison_chamber_start':d['current_s8']['coefficients'],'right_comparison_chamber_end':enc(nxt),'new_improving_candidate':None})
 elif ZERO<gright:
  value=field(d['candidate']['coefficients'])+gright*(nxt-event);lower=A.bounds(value.a)[0];oldhi=Q(d['candidate']['exact_interval'][1])
  out.update({'status':'finished_positive_right_slope','new_endpoint_ideal_candidate':view(value),'candidate_lower_above_complete_E075_upper':oldhi<lower,'candidate_margin':str(lower-oldhi)})
 else:out['status']='finished_zero_right_slope'
except Exception as e:out.update({'status':'stopped_exception','error':repr(e)})
finally:
 out.update({'symbolic_configurations':configs,'combined_nodes':nodes,'new_5D_union_nodes':unionnodes,'cpu_seconds':time.process_time()-CPU,'wall_seconds_before_exclusive_save':time.monotonic()-WALL,'input_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [PREV,IFACE,ARITH,HERE/'registration.json',Path(__file__)]}})
 with OUT.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
 print(json.dumps({k:out[k] for k in ['status','pid','symbolic_configurations','combined_nodes','new_5D_union_nodes','cpu_seconds','wall_seconds_before_exclusive_save']}))
 signal.alarm(0)
 if 'error' in out:print(out['error']);sys.exit(1)
