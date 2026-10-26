#!/usr/bin/env python3
"""Continue once from frozen NEW section boxes; no old-node replay."""
from pathlib import Path
from fractions import Fraction as Q
from functools import lru_cache
import json,hashlib,importlib.util,sys,time,resource,signal,os
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;OUT=HERE/'section_results.json'
if OUT.exists():raise FileExistsError('Refuse overwrite '+str(OUT))
resource.setrlimit(resource.RLIMIT_CPU,(20,21))
def wallcap(*_):raise TimeoutError('180 second total program wall cap')
signal.signal(signal.SIGALRM,wallcap);signal.alarm(180)
CPU,WALL=time.process_time(),time.monotonic()
SAVED=HERE/'results.json';LOWER=HERE.parents[1];PREV=LOWER/'round075/post_s8_overlap/results.json'
BASE=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds')
ARITH=BASE/'lower/round001_rectangular_certificate.py';IFACE=BASE/'lower/round008/spacing_chamber/chamber_interface.json'
s=json.loads(SAVED.read_text());configs=s['symbolic_configurations'];nodes=s['combined_nodes'];unionnodes=0
out={'status':'running','pid':os.getpid(),'new_C6_bound':False,'full_6D_union_calls':0,'new_5D_union_calls':0,'old_phase_reruns':0,'old_background_nodes_replayed':0}
def budget():
 if time.process_time()-CPU+s['cpu_seconds']>=19.5:raise TimeoutError('combined CPU prelimit')
 if time.monotonic()-WALL>=179:raise TimeoutError('wall prelimit')
 if configs>8 or nodes>1000:raise RuntimeError('combined node/configuration cap')
try:
 assert s['status']=='stopped_exception' and s['error']=='AssertionError()' and s['new_5D_union_nodes']==0
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
 def boxdecode(B):return tuple(field(z) for pair in B for z in pair)
 Bs=tuple(sorted(boxdecode(B) for B in s['retained_cross_sections']))
 def vol(B):
  z=F(1)
  for q in range(len(B)//2):z=z*(B[2*q+1]-B[2*q])
  return z
 Pvolume=F(Q(1,32))
 # Fallback saved-only geometry bounds before union; no old-node evaluation.
 maxsingle=max((vol(B) for B in Bs),default=ZERO);sumvol=sum((vol(B) for B in Bs),ZERO)
 out['pre_union_volume_bounds']={'background_lower':view(maxsingle),'background_upper':view(min(Pvolume,sumvol))}
 @lru_cache(None)
 def new_union(B):
  global nodes,unionnodes
  nodes+=1;unionnodes+=1;budget()
  dim=len(B[0])//2 if B else 1
  if not B:return ZERO,(ZERO,HALF)
  if dim==1:
   ordered=sorted(B);lo,hi=ordered[0];total=ZERO;sample=(ZERO,lo) if ZERO<lo else None
   for a,b in ordered[1:]:
    if hi<a:
     if sample is None:sample=(hi,a)
     total=total+hi-lo;lo,hi=a,b
    else:hi=max(hi,b)
   if sample is None and hi<HALF:sample=(hi,HALF)
   return total+hi-lo,sample
  cuts=sorted({ZERO,HALF,*[x for b in B for x in b[:2]]});total=ZERO;sample=None
  for a,b in zip(cuts,cuts[1:]):
   child=tuple(sorted({z[2:] for z in B if z[0]<=a and b<=z[1]}))
   if child:
    value,childgap=new_union(child);total=total+(b-a)*value
    if sample is None and childgap is not None:sample=(a,b)+childgap
   elif sample is None:sample=(a,b)+tuple(z for _ in range(dim-1) for z in (ZERO,HALF))
  return total,sample
 out['new_5D_union_calls']=1
 bg,sample=new_union(Bs);uncovered=Pvolume-bg;assert ZERO<uncovered and sample is not None
 gleft=field(s['left_slope']['coefficients']);gright=gleft-64*uncovered
 rational=[]
 for q in range(5):
  lo,hi=sample[2*q:2*q+2];assert lo<hi
  aa=A.bounds(lo.a)[1];bb=A.bounds(hi.a)[0];assert aa<bb;rational.append([str(aa),str(bb)])
 witnesses=[]
 for k,B in enumerate(Bs):
  disjoint=[q for q in range(5) if sample[2*q+1]<=B[2*q] or B[2*q+1]<=sample[2*q]]
  assert disjoint
  witnesses.append({'retained_box_position':k,'separating_original_axis':[1,3,4,5,6][disjoint[0]]})
 out.update({'exact_background_5D_union':view(bg),'exact_uncovered_5D_volume':view(uncovered),'exact_ratio_slope_jump':view(-64*uncovered),'right_slope':view(gright),'new_5D_union_cache_states':new_union.cache_info().currsize,'open_uncovered_region':{'exact_endpoints':boxenc(sample),'rational_inward_open_intervals':rational,'original_axes':[1,3,4,5,6],'all_21_retained_boxes_disjoint_witnesses':witnesses,'not_used_as_full_volume':True}})
 # Only new right chamber using already frozen37 definitions/current states.
 d=json.loads(PREV.read_text());iface=json.loads(IFACE.read_text());defs=d['complete_functions_reused'];funcs=[{'c':field(f['c']),'m':Q(f['m'])} for f in defs]
 event=field(s['s9']['coefficients']);order=list(d['right_complete_order']);pos=order.index(7);assert order[pos+1]==2;order[pos],order[pos+1]=order[pos+1],order[pos]
 def val(i,t):return funcs[i]['c']+funcs[i]['m']*t
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
 nl,nr=state(nxt,-1),state(nxt,1);assert nl==s['current_right_state_reused']
 out.update({'right_chamber_next_endpoint':view(nxt),'right_complete_order':order,'right_chamber_full_order_gaps':gaps,'next_tie_edges':[z['ids'] for z in gaps if not z['right_gap']],'next_states':{'left':nl,'right':nr},'new_adjacent_checks':36})
 if gright<ZERO:
  assert ZERO<gleft
  out.update({'status':'finished_strict_local_maximum','strict_local_maximum_family':'same127 unit-mass common-radius standard central witness union, fixed a1 and fourtail spacings, absolute a2 varied only','strict_local_maximum_parameter':s['s9']['coefficients'],'left_comparison_chamber_start':d['current_s8']['coefficients'],'right_comparison_chamber_end':enc(nxt),'new_improving_candidate':None})
 elif ZERO<gright:
  candidate=field(d['candidate']['coefficients'])+gright*(nxt-event);cl=A.bounds(candidate.a)[0];oldhi=Q(d['candidate']['exact_interval'][1])
  out.update({'status':'finished_positive_right_slope','new_endpoint_ideal_candidate':view(candidate),'candidate_lower_above_complete_E075_upper':oldhi<cl,'candidate_margin':str(cl-oldhi)})
 else:out['status']='finished_zero_right_slope'
except Exception as e:out.update({'status':'stopped_exception','error':repr(e)})
finally:
 out.update({'combined_symbolic_configurations':configs,'combined_nodes':nodes,'new_5D_union_nodes':unionnodes,'continuation_cpu_seconds':time.process_time()-CPU,'combined_cpu_seconds':s['cpu_seconds']+time.process_time()-CPU,'continuation_wall_seconds_before_exclusive_save':time.monotonic()-WALL,'input_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [SAVED,PREV,IFACE,ARITH,HERE/'continuation_registration.json',Path(__file__)]}})
 with OUT.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
 print(json.dumps({k:out[k] for k in ['status','pid','combined_symbolic_configurations','combined_nodes','new_5D_union_nodes','continuation_cpu_seconds','combined_cpu_seconds','continuation_wall_seconds_before_exclusive_save']}))
 signal.alarm(0)
 if 'error' in out:print(out['error']);sys.exit(1)
