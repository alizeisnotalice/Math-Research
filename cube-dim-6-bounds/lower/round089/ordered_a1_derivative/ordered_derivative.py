#!/usr/bin/env python3
"""One bounded new batch: exact a1 derivatives by S5 order-statistic complements."""
from pathlib import Path
from fractions import Fraction as Q
from functools import cmp_to_key, lru_cache
from math import factorial,prod
from itertools import permutations
import json,hashlib,time,resource,signal,os,sys
sys.dont_write_bytecode=True
H=Path(__file__).resolve().parent; OUT=H/'results.json'
if OUT.exists():raise FileExistsError('Exclusive result already exists')
CPU,WALL=time.process_time(),time.monotonic()
resource.setrlimit(resource.RLIMIT_CPU,(20,20))
def alarm(*_):raise TimeoutError('170 second soft wall cap, cleanup reserve remains')
signal.signal(signal.SIGALRM,alarm);signal.alarm(170)
import field as A
A.bounds=lru_cache(None)(A.bounds)
Z=A.Z; C=A.constant; add=A.add; sub=A.sub; scale=A.scale; mul=A.mul
parse=A.parse; encode=A.encode; sign=A.sign; view=A.view
admitted=0;refused=0;requested=0;configs=0;cache={};cell_receipts=[]
out={'status':'running','pid':os.getpid(),'sides':{},'section_families':[], 'static_6D_union_calls':0,'old_DP_replays':0,'old_main_run_executions':0,'automatic_reruns':0,'new_C6_claim':False}
def budget():
 if time.process_time()-CPU>=18.5:raise TimeoutError('18.5 second algorithm prelimit; no retry')
 if time.monotonic()-WALL>=169:raise TimeoutError('wall algorithm prelimit; no retry')
def admit():
 global admitted,requested,refused
 budget();requested+=1
 if admitted>=1000:
  refused+=1;raise RuntimeError('1001st integration-state admission refused')
 admitted+=1

def run():
 global configs
 reg=json.loads((H/'registration.json').read_text())
 for path,sha in reg['source_pins'].items():assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==sha
 paths=list(reg['source_pins']); iface=json.loads(Path(paths[2]).read_text());prev=json.loads(Path(paths[3]).read_text());d=json.loads(Path(paths[5]).read_text());st=json.loads(Path(paths[6]).read_text());e=json.loads(Path(paths[7]).read_text())
 assert prev['status']=='finished' and e['status']=='pass_new_s8_disjoint_seam_chamber_and_s9_finite_L1'
 assert parse(prev['candidate']['coefficients'])==parse(e['new_ideal_ratio']['coefficients'])
 ell=[tuple(x) for x in iface['box_lengths']];assert len(ell)==len(set(ell))==127
 t=add(C(1),A.mono(4,0,Q(1,2)));a=scale(add(C(1),A.mono(1,0)),Q(1,2));cap=scale(a,Q(1,2))
 assert t==parse(d['a1_current']['coefficients']) and a==parse(d['fixed_tail_spacing']['coefficients'])
 cuts=[parse(x) for x in d['tail_rank_coordinates']];assert cuts[0]==Z and cuts[-1]==cap and len(cuts)==21
 for x,y in zip(cuts,cuts[1:]):assert sign(sub(y,x))>0
 ranks={x:i for i,x in enumerate(cuts)};projections=[]
 for k,ls in enumerate(ell):
  radius=A.rad(ls);p=[]
  for l in ls[1:]:
   q=-(l//2) if l%2 else 1-l//2
   lo=A.maximum(Z,sub(scale(a,q+l-1),radius));hi=A.minimum(cap,add(scale(a,q),radius))
   assert sign(sub(hi,lo))>0 and (lo==Z or hi==cap)
   p.extend((ranks[lo],ranks[hi]))
  assert p==d['raw_projection_ranks'][k];projections.append(tuple(p))
 funcs=[(parse(z['c']),Q(z['m'])) for z in d['normal_functions']];assert len(funcs)==24
 moving={z['raw_index']:i for i,f in enumerate(d['normal_functions']) for z in f['labels'] if 'raw_index' in z}
 for k,ls in enumerate(ell):
  r=A.rad(ls);l=ls[0];want=(r,-Q(l-1,2)) if l%2 else (scale(r,-1),Q(l,2))
  assert funcs[moving[k]]==want
 assert funcs[0]==(Z,Q(0)) and funcs[1]==(Z,Q(1,2))
 vals=[add(c,scale(t,m)) for c,m in funcs];states={};all_families={};full5=cap
 for _ in range(4):full5=mul(full5,cap)
 out.update({'baseline':view(parse(prev['candidate']['coefficients'])),'a1':view(t),'tail_spacing':view(a),'tail_cap':view(cap),'tail_cuts':[encode(x) for x in cuts],'geometry_all_127_projections_regenerated':True,'root_enclosures':{'denominator':A.DEN,'u':[str(A.UL),str(A.UH)],'v':[str(A.VL),str(A.VH)]}})
 for side,label in [(-1,'left'),(1,'right')]:
  configs+=1;assert configs<=8;budget()
  def cmp(i,j):
   s=sign(sub(vals[i],vals[j]));return s if s else ((side*funcs[i][1]>side*funcs[j][1])-(side*funcs[i][1]<side*funcs[j][1]))
  key=cmp_to_key(cmp);sel={};used=set()
  for k,ls in enumerate(ell):
   m=moving[k];lo=0 if ls[0]%2 else max((0,m),key=key);hi=min((1,m),key=key) if ls[0]%2 else 1
   assert cmp(lo,hi)<0;sel[str(k)]=[lo,hi];used.update((lo,hi))
  order=sorted(used,key=key);allorder=sorted(range(24),key=key)
  assert sel==st['sides'][label]['selections'] and order==st['sides'][label]['effective_order'] and allorder==st['sides'][label]['all_function_order']
  pos={v:j for j,v in enumerate(order)};slabs=[]
  for lo,hi in zip(order,order[1:]):
   weight=funcs[hi][1]-funcs[lo][1]
   if not weight:continue
   raw=[k for k in range(127) if pos[sel[str(k)][0]]<=pos[lo] and pos[hi]<=pos[sel[str(k)][1]]]
   labels={ell[k] for k in raw};closures=0
   for ls in labels:
    for j in range(1,5):
     swap=list(ls);swap[j],swap[j+1]=swap[j+1],swap[j];assert tuple(swap) in labels;closures+=1
   orbits={}
   for k in raw:orbits.setdefault((ell[k][0],tuple(sorted(ell[k][1:]))),[]).append(k)
   ordered={}
   for orbit,ids in sorted(orbits.items()):
    k=ids[0];pairs=[tuple(projections[k][2*j:2*j+2]) for j in range(5)]
    pref=sorted((x,y) for x,y in pairs if x==0 and y<20);suf=sorted((x,y) for x,y in pairs if x>0);full=[(0,20)]*(5-len(pref)-len(suf));box=tuple(pref+full+suf)
    for h in ids:assert sorted(tuple(projections[h][2*j:2*j+2]) for j in range(5))==sorted(pairs)
    expected=len(set(permutations(orbit[1])));assert len(ids)==expected
    ordered.setdefault(box,[]).append({'normal_length':orbit[0],'sorted_tail_lengths':list(orbit[1]),'raw_indices':ids,'orbit_size':expected})
   boxes=tuple(sorted(ordered));keep=tuple(B for B in boxes if not any(B!=D and all(c<=x and y<=z for (x,y),(c,z) in zip(B,D)) for D in boxes))
   removed=[{'box':B,'covering_retained_box':next(D for D in keep if all(c<=x and y<=z for (x,y),(c,z) in zip(B,D)))} for B in boxes if B not in keep]
   events=[];L=[0]*5;U=[20]*5
   for B in keep:
    pref=[(i,y) for i,(x,y) in enumerate(B) if y<20];suf=[(i,x) for i,(x,y) in enumerate(B) if x>0]
    assert bool(pref)!=bool(suf),'every retained event must be a single prefix or suffix threshold'
    if pref:
     assert [i for i,y in pref]==list(range(len(pref))) and len({y for i,y in pref})==1
     i,y=pref[-1];events.append({'kind':'prefix','order_index':i+1,'cut_rank':y});L[i]=max(L[i],y)
    else:
     assert [i for i,x in suf]==list(range(5-len(suf),5)) and len({x for i,x in suf})==1
     i,x=suf[0];events.append({'kind':'suffix','order_index':i+1,'cut_rank':x});U[i]=min(U[i],x)
   for i in range(1,5):L[i]=max(L[i],L[i-1])
   for i in range(3,-1,-1):U[i]=min(U[i],U[i+1])
   ck=(tuple(L),tuple(U));fid=all_families.setdefault(ck,len(all_families))
   slabs.append({'ids':[lo,hi],'width_slope':str(weight),'active_raw_indices':raw,'S5_adjacent_swap_checks':closures,'label_orbits':list(ordered.values()),'all_ordered_boxes':boxes,'retained_ordered_boxes':keep,'containment_deletions':removed,'threshold_events':events,'complement_lower_ranks':L,'complement_upper_ranks':U,'family_id':fid,'zero_baseline_width':vals[lo]==vals[hi]})
  assert [{'ids':x['ids'],'width_slope':x['width_slope'],'active_raw_indices':x['active_raw_indices'],'zero_baseline_width':x['zero_baseline_width'],'5D_union_area':'NOT_COMPUTED'} for x in slabs]==st['sides'][label]['required_weighted_slabs']
  states[label]={'order':order,'allorder':allorder,'selections':sel};out['sides'][label]={'slabs':slabs,'effective_order':order,'all_function_order':allorder,'selections':sel}
 planned=sum(5*(len(set((0,20)+L+U))-1) for L,U in all_families if all(x<y for x,y in zip(L,U)))
 out['planned_unique_families']=len(all_families);out['planned_integration_states']=planned;assert planned<=1000
 for (L,U),fid in all_families.items():
  budget();knots=sorted(set((0,20)+L+U));initial=admitted;rows=[]
  if any(x>=y for x,y in zip(L,U)):vol=Z
  else:
   prev_polys=[[C(1)] for _ in range(len(knots)-1)]
   for i in range(5):
    carry=Z;next_polys=[];row=[]
    for j,(l,h) in enumerate(zip(knots,knots[1:])):
     admit();active=L[i]<=l and h<=U[i];p=[carry]+[scale(c,Q(1,k+1)) for k,c in enumerate(prev_polys[j])] if active else [carry]
     w=sub(cuts[h],cuts[l]);v=Z
     for c in reversed(p):v=add(mul(v,w),c)
     row.append({'interval_ranks':[l,h],'active':active,'local_polynomial':[encode(c) for c in p],'right_value':encode(v)})
     carry=v;next_polys.append(p)
    rows.append(row);prev_polys=next_polys
   vol=carry
  area=sub(full5,scale(vol,120));assert sign(area)>=0 and sign(sub(full5,area))>=0
  cache[fid]=area
  out['section_families'].append({'family_id':fid,'complement_lower_ranks':L,'complement_upper_ranks':U,'knots':knots,'integration_states':admitted-initial,'recurrence_rows':rows,'ordered_complement_volume':view(vol),'full_section_union_area':view(area)})
 assert admitted==planned
 for label in ('left','right'):
  deriv=Z
  for slab in out['sides'][label]['slabs']:
   term=scale(cache[slab['family_id']],64*Q(slab['width_slope']));deriv=add(deriv,term);slab['weighted_derivative_term']=view(term)
  out['sides'][label]['exact_derivative']=view(deriv)
  out['sides'][label]['derivative_sign']=sign(deriv)
 # Adjacent affine-order event gives the complete a1 chamber, without union rerun.
 for direction,label in [(-1,'left'),(1,'right')]:
  order=states[label]['allorder'];events=[]
  for i,j in zip(order,order[1:]):
   c,m=funcs[i];h,n=funcs[j]
   if m==n:continue
   x=scale(sub(h,c),1/(m-n))
   if sign(sub(x,t))==direction:events.append((x,i,j))
  assert events
  endpoint=(min if direction>0 else max)(events,key=lambda item:cmp_to_key(A.cmp)(item[0]))[0]
  gaps=[]
  for i,j in zip(order,order[1:]):
   c,m=funcs[i];h,n=funcs[j];g0=sub(vals[j],vals[i]);g1=add(sub(h,c),scale(endpoint,n-m))
   assert sign(g0)>=0 and sign(g1)>=0 and (sign(g0)>0 or sign(g1)>0)
   gaps.append({'ids':[i,j],'at_current':encode(g0),'at_adjacent_endpoint':encode(g1)})
  out['sides'][label]['adjacent_endpoint']=view(endpoint);out['sides'][label]['complete_order_gap_certificate']=gaps
 gl=out['sides']['left']['derivative_sign'];gr=out['sides']['right']['derivative_sign']
 if gl>0 and gr<0:out['status']='finished_strict_a1_local_maximum'
 elif gl<0 or gr>0:
  direction=1 if gr>0 else -1;label='right' if direction>0 else 'left';side=out['sides'][label]
  candidate=add(parse(prev['candidate']['coefficients']),mul(parse(side['exact_derivative']['coefficients']),sub(parse(side['adjacent_endpoint']['coefficients']),t)))
  assert A.bounds(candidate)[0]>Q(prev['candidate']['exact_interval'][1]);out.update({'status':'finished_improving_a1_candidate','candidate_direction':direction,'candidate':view(candidate),'candidate_endpoint':side['adjacent_endpoint']})
 else:out['status']='finished_no_strict_improving_derivative'
 out['scope']='fixed 127 witnesses, fixed five equal tail spacings, a1 adjacent chambers only; no all-family global optimum or new C6 registration'
 return reg

try:run()
except Exception as exc:out.update({'status':'stopped_exception','error':repr(exc)})
finally:
 out['resource']={'configurations':configs,'integration_state_requests':requested,'admitted_new_states':admitted,'refused_requests':refused,'completed_family_count':len(cache),'algorithm_CPU_seconds':time.process_time()-CPU,'wall_seconds_before_save':time.monotonic()-WALL,'MATLAB_calls':0}
 out['input_sha256']={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [H/'registration.json',H/'field.py',Path(__file__)]}
 with OUT.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
 receipt={'status':out['status'],'pid':out['pid'],'CPU_seconds_including_save':time.process_time()-CPU,'wall_seconds_including_save':time.monotonic()-WALL,'resource':out['resource'],'result_sha256':hashlib.sha256(OUT.read_bytes()).hexdigest(),'exit_code':1 if 'error' in out else 0}
 with (H/'receipt.json').open('x') as f:json.dump(receipt,f,indent=2);f.write('\n')
 print(json.dumps(receipt));signal.alarm(0)
 if 'error' in out:print(out['error']);sys.exit(1)
