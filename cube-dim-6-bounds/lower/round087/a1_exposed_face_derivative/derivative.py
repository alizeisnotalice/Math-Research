#!/usr/bin/env python3
"""New a1 Fubini derivative via complete NEW five-tail cross sections."""
from pathlib import Path
from fractions import Fraction as Q
from functools import lru_cache
from itertools import permutations
from math import prod
import json,hashlib,importlib.util,sys,time,resource,signal,os
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;OUT=HERE/'results.json'
if OUT.exists():raise FileExistsError('Refuse overwrite '+str(OUT))
resource.setrlimit(resource.RLIMIT_CPU,(20,21))
def wallcap(*_):raise TimeoutError('180 second total wall cap')
signal.signal(signal.SIGALRM,wallcap);signal.alarm(180)
CPU,WALL=time.process_time(),time.monotonic()
LOWER=HERE.parents[1];PREV=LOWER/'round075/post_s8_overlap/results.json'
BASE=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds')
ARITH=BASE/'lower/round001_rectangular_certificate.py';IFACE=BASE/'lower/round008/spacing_chamber/chamber_interface.json'
nodes,configs,canon_calls=0,0,0
out={'status':'running','pid':os.getpid(),'static_6D_union_calls':0,'old_program_replays':0,'new_C6_bound':False,'automatic_reruns':0,'sides':{}}
def budget():
 if time.process_time()-CPU>=19.5 or time.monotonic()-WALL>=179:raise TimeoutError('CPU/wall prelimit')
 if nodes>1000 or configs>8:raise RuntimeError('DP/configuration cap')
try:
 spec=importlib.util.spec_from_file_location('accepted_field',ARITH);A=importlib.util.module_from_spec(spec);spec.loader.exec_module(A)
 A.D=10**80;A.RU,A.RV=A.floor_root(2*A.D**6,6),A.floor_root(3*A.D**6,6)
 A.UL,A.UH=Q(A.RU,A.D),Q(A.RU+1,A.D);A.VL,A.VH=Q(A.RV,A.D),Q(A.RV+1,A.D);A.bounds.cache_clear()
 F,ZERO,ONE,u=A.F,A.F(0),A.F(1),A.u
 def field(z):
  a=[Q(0)]*36
  for ij,x in z.items():i,j=map(int,ij.split(','));a[6*i+j]=Q(x)
  return F(a=a)
 def enc(x):return {f'{i},{j}':str(x.a[6*i+j]) for i in range(6) for j in range(6) if x.a[6*i+j]}
 def view(x):
  lo,hi=A.bounds(x.a);return {'coefficients':enc(x),'exact_interval':[str(lo),str(hi)],'display':float((lo+hi)/2)}
 d=json.loads(PREV.read_text());iface=json.loads(IFACE.read_text());ells=[tuple(z) for z in iface['box_lengths']]
 current=ONE+u**4/2;tail=(ONE+u)/2;tailvals=[]
 for k,b in enumerate(iface['boxes']):
  B=[]
  for axis in range(1,6):
   if axis==1:
    pair=d['next_states']['left']['selections'][str(k)]
    B.extend(field(d['complete_functions_reused'][i]['c'])+Q(d['complete_functions_reused'][i]['m'])*tail for i in pair)
   else:B.extend(field(iface['endpoints'][eid]['c0']) for eid in b[axis])
  tailvals.append(tuple(B))
 cuts=sorted({x for B in tailvals for x in B});rank={x:i for i,x in enumerate(cuts)}
 projections=[tuple(rank[x] for x in B) for B in tailvals]
 funcs=[];idx={};records=[]
 def add(c,m,label):
  key=(c,Q(m))
  if key not in idx:idx[key]=len(funcs);funcs.append({'c':c,'m':Q(m),'labels':[]})
  i=idx[key];funcs[i]['labels'].append(label);return i
 zero=add(ZERO,0,{'face':'zero'});cap=add(ZERO,Q(1,2),{'face':'moving_halfcell'})
 for k,ell in enumerate(ells):
  radius=A.root(prod(ell))/2;l=ell[0]
  if l%2:
   moving=add(radius,-Q(l-1,2),{'raw_index':k,'ell':ell,'mass':prod(ell),'face':'odd_upper'})
   opts=((zero,),(cap,moving))
  else:
   moving=add(-radius,Q(l,2),{'raw_index':k,'ell':ell,'mass':prod(ell),'face':'even_lower'})
   opts=((zero,moving),(cap,))
  records.append({'ell':ell,'options':opts,'projection':projections[k]})
 def val(i,t):return funcs[i]['c']+funcs[i]['m']*t
 def state(t,side):
  global configs
  configs+=1;budget();keys=[(val(i,t),side*f['m']) for i,f in enumerate(funcs)];sel={};used=set()
  for k,r in enumerate(records):
   lo=max(r['options'][0],key=lambda i:keys[i]);hi=min(r['options'][1],key=lambda i:keys[i])
   if keys[lo]<keys[hi]:sel[k]=(lo,hi);used.update((lo,hi))
  return sel,sorted(used,key=lambda i:keys[i]),keys
 out.update({'a1_current':view(current),'fixed_tail_spacing':view(tail),'old_full_E075_ideal_interval':d['candidate']['exact_interval'],'normal_functions':[{'id':i,'c':enc(f['c']),'m':str(f['m']),'labels':f['labels']} for i,f in enumerate(funcs)],'tail_rank_coordinates':[enc(x) for x in cuts],'raw_projection_ranks':projections,'radius_rule':'prod(ell)^(1/6)/2 constant along a1; same radius all6coordinates'})
 # Canonicalize ENTIRE box family under all d! global coordinate permutations.
 # No symmetry of individual boxes is assumed; volume is permutation invariant.
 perms={j:list(permutations(range(j))) for j in range(1,6)}
 @lru_cache(None)
 def canonical(Bs):
  global canon_calls
  canon_calls+=1;budget()
  if not Bs:return ()
  dim=len(Bs[0])//2
  return min(tuple(sorted(tuple(z for j in p for z in B[2*j:2*j+2]) for B in Bs)) for p in perms[dim])
 @lru_cache(None)
 def area(Bs):
  global nodes
  nodes+=1;budget()
  if not Bs:return ZERO
  if len(Bs[0])==2:
   ordered=sorted(Bs);lo,hi=ordered[0];v=ZERO
   for a,b in ordered[1:]:
    if hi<a:v=v+cuts[hi]-cuts[lo];lo,hi=a,b
    else:hi=max(hi,b)
   return v+cuts[hi]-cuts[lo]
  cutr=sorted({x for B in Bs for x in B[:2]});v=ZERO
  for a,b in zip(cutr,cutr[1:]):
   child=canonical(tuple(sorted({B[2:] for B in Bs if B[0]<=a and b<=B[1]})))
   if child:v=v+(cuts[b]-cuts[a])*area(child)
  return v
 savedstates={}
 for side,label in [(-1,'left'),(1,'right')]:
  sel,order,keys=state(current,side);assert len(sel)==127
  pos={k:i for i,k in enumerate(order)};slabs=[];g=ZERO
  for lo,hi in zip(order,order[1:]):
   slope=funcs[hi]['m']-funcs[lo]['m']
   if not slope:continue
   active=[k for k,(a,b) in sel.items() if pos[a]<=pos[lo] and pos[hi]<=pos[b]]
   B=canonical(tuple(sorted(set(projections[k] for k in active))))
   ar=area(B);term=64*slope*ar;g=g+term
   slabs.append({'low_id':lo,'high_id':hi,'width_slope':str(slope),'active_raw_indices':active,'canonical_projected_boxes':B,'exact_5D_area':view(ar),'ratio_derivative_contribution':view(term),'zero_baseline_width':val(lo,current)==val(hi,current)})
  ties=[]
  allorder=sorted(range(len(funcs)),key=lambda i:keys[i]);group=[allorder[0]]
  for i,j in zip(allorder,allorder[1:]):
   if val(i,current)==val(j,current):group.append(j)
   else:
    if len(group)>1:ties.append(group)
    group=[j]
  if len(group)>1:ties.append(group)
  out['sides'][label]={'positive_boxes':len(sel),'selections':{str(k):list(v) for k,v in sel.items()},'effective_order':order,'all_function_order':allorder,'all_tie_groups':ties,'slabs_with_nonzero_derivative_weight':slabs,'exact_ratio_derivative':view(g)}
  savedstates[label]=(sel,order,allorder,g)
 # All full normal functions affine => one adjacent event in either direction.
 # No further candidate unless exact sided sign gives an improving direction.
 gl,gr=savedstates['left'][3],savedstates['right'][3]
 if gr>ZERO:direction=1;slope=gr
 elif gl<ZERO:direction=-1;slope=gl
 else:direction=0
 out['new_5D_cache_states']=area.cache_info().currsize;out['canonicalization_states']=canonical.cache_info().currsize
 if direction:
  label='right' if direction>0 else 'left';order=savedstates[label][2];events=[]
  for i,j in zip(order,order[1:]):
   f,h=funcs[i],funcs[j]
   if f['m']==h['m']:continue
   t=(h['c']-f['c'])/(f['m']-h['m'])
   if (current<t if direction>0 else t<current):events.append((t,i,j))
  nxt=(min if direction>0 else max)(events,key=lambda z:z[0])[0]
  gaps=[]
  for i,j in zip(order,order[1:]):
   a=val(j,current)-val(i,current);b=val(j,nxt)-val(i,nxt)
   assert ZERO<=a and ZERO<=b and (ZERO<a or ZERO<b);gaps.append({'ids':[i,j],'current_gap':enc(a),'next_gap':enc(b)})
  ns=state(nxt,-direction);assert ns[0]==savedstates[label][0]
  base=field(d['candidate']['coefficients']);candidate=base+slope*(nxt-current);lower=A.bounds(candidate.a)[0];oldupper=Q(d['candidate']['exact_interval'][1]);assert oldupper<lower
  out.update({'status':'finished_improving_candidate','direction':direction,'next_endpoint':view(nxt),'new_chamber_full_order_gaps':gaps,'candidate':view(candidate),'candidate_lower_above_complete_E075_upper':True,'candidate_margin_over_complete_E075_upper':str(lower-oldupper)})
 elif gl>ZERO and gr<ZERO:out.update({'status':'finished_strict_local_maximum','scope':'fixed127 witness union, fixedfive tail spacings, a1 only; not all-family global optimality'})
 else:out.update({'status':'finished_no_strict_improving_sided_derivative','scope':'Only exact sided derivative signs; zero case unresolved higher order'})
except Exception as e:out.update({'status':'stopped_exception','error':repr(e)})
finally:
 out.update({'new_configurations':configs,'new_DP_nodes':nodes,'canonical_calls':canon_calls,'CPU_seconds':time.process_time()-CPU,'wall_seconds_before_exclusive_save':time.monotonic()-WALL,'input_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [PREV,IFACE,ARITH,HERE/'registration.json',Path(__file__)]}})
 with OUT.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
 print(json.dumps({k:out[k] for k in ['status','pid','new_configurations','new_DP_nodes','canonical_calls','CPU_seconds','wall_seconds_before_exclusive_save']}))
 signal.alarm(0)
 if 'error' in out:print(out['error']);sys.exit(1)
