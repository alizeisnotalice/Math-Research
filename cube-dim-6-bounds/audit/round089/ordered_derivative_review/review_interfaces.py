#!/usr/bin/env python3
"""Two pure exact interface audits; never generate a union/integral state."""
from pathlib import Path
from fractions import Fraction as Q
from functools import cmp_to_key,lru_cache
from math import prod,factorial
from collections import Counter
from decimal import Decimal,localcontext
from types import SimpleNamespace
import ast,hashlib,json,os,resource,signal,sys,time,traceback
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;OUT=HERE/'interface_results.json'
if OUT.exists():raise FileExistsError('Refuse overwrite '+str(OUT))
CPU=time.process_time();WALL=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(15,15))
def timeout(*_):raise TimeoutError('180-second total wall including cleanup')
signal.signal(signal.SIGALRM,timeout);signal.alarm(180)
AUTHOR=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-lower-orderstat-astra/cube-dim-6-bounds/lower/round089/ordered_a1_derivative')
MAIN=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds')
DROOT=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-lower-clean/cube-dim-6-bounds/lower')
ARITH=MAIN/'audit/round058/a2_next_event_review/review.py';IFACE=MAIN/'lower/round008/spacing_chamber/chamber_interface.json'
EPREV=HERE.parents[1]/'round075/post_s8_overlap_review/review_results.json'
configs=0;checks=Counter()
out={'status':'running','pid':os.getpid(),'new_integration_states':0,'old087_DP_replays':0,'static6D_union_calls':0}
def check(v,kind):
 checks[kind]+=1;assert v,kind
 assert configs<=2 and time.process_time()-CPU<14 and time.monotonic()-WALL<179

def run():
 global configs
 gate=json.loads((HERE/'source_gate.json').read_text());hashes=dict(gate['source_sha256'])
 for p,h in hashes.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h
 for p in [HERE/'registration.json',HERE/'source_gate.json',Path(__file__)]:hashes[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
 # Load only explicitly whitelisted accepted pure AST definitions, never old module code.
 names={'constant','mono','parse','encode','add','scale','sub','mul','rootbox','bounds','sign','cmp','minimum','maximum','dec','view','rad'}
 body=[]
 for n in ast.parse(ARITH.read_text()).body:
  if isinstance(n,ast.FunctionDef) and n.name in names:n.decorator_list=[];body.append(n)
 assert {n.name for n in body}==names
 ns={'Q':Q,'prod':prod,'Decimal':Decimal,'localcontext':localcontext,'DEN':10**100,'Z':(Q(0),)*36}
 exec(compile(ast.Module(body=body,type_ignores=[]),str(ARITH)+'::pure-whitelist','exec'),ns)
 ns['UL'],ns['UH']=ns['rootbox'](2);ns['VL'],ns['VH']=ns['rootbox'](3)
 ns['POWER_BOUNDS']=[(ns['UL']**(k//6)*ns['VL']**(k%6),ns['UH']**(k//6)*ns['VH']**(k%6)) for k in range(36)]
 ns['bounds']=lru_cache(None)(ns['bounds']);A=SimpleNamespace(**ns)
 r=json.loads((AUTHOR/'results.json').read_text());iface=json.loads(IFACE.read_text());old=json.loads((DROOT/'round075/post_s8_overlap/results.json').read_text());ep=json.loads(EPREV.read_text());d=json.loads((DROOT/'round087/a1_exposed_face_derivative/results.json').read_text());st=json.loads((DROOT/'round087/a1_exposed_face_derivative/structure_results.json').read_text())
 assert r['status']=='finished_strict_a1_local_maximum' and ep['new_endpoint_all_127_boxes_valid']
 t=A.add(A.constant(1),A.mono(4,0,Q(1,2)));a=A.scale(A.add(A.constant(1),A.mono(1,0)),Q(1,2));cap=A.scale(a,Q(1,2))
 assert t==A.parse(r['a1']['coefficients'])==A.parse(d['a1_current']['coefficients'])
 assert a==A.parse(ep['s9']['coefficients'])==A.parse(r['tail_spacing']['coefficients'])==A.parse(d['fixed_tail_spacing']['coefficients'])
 assert A.parse(r['baseline']['coefficients'])==A.parse(ep['new_ideal_ratio']['coefficients'])==A.parse(old['candidate']['coefficients'])
 ell=[tuple(x) for x in iface['box_lengths']];assert len(ell)==len(set(ell))==127
 configs+=1
 cuts=[A.parse(x) for x in r['tail_cuts']];check(cuts==[A.parse(x) for x in d['tail_rank_coordinates']] and cuts[0]==A.Z and cuts[-1]==cap and len(cuts)==21,'complete_tail_cut_binding')
 for x,y in zip(cuts,cuts[1:]):check(A.sign(A.sub(y,x))>0,'strict_cut_order')
 ranks={x:i for i,x in enumerate(cuts)};proj=[];oldvals=[A.add(A.parse(f['c']),A.scale(a,Q(f['m']))) for f in old['complete_functions_reused']]
 for k in range(127):
  B=[]
  for j in range(1,6):
   if j==1:ends=[oldvals[i] for i in old['next_states']['left']['selections'][str(k)]]
   else:
    ends=[]
    for eid in iface['boxes'][k][j]:
     z=iface['endpoints'][eid];check(A.parse(z['c1'])==A.Z,'fixed_tail_source');ends.append(A.parse(z['c0']))
   lo,hi=[ranks[x] for x in ends];check(0<=lo<hi<=20 and (lo==0 or hi==20),'tail_prefix_suffix_or_full');B.append((lo,hi))
  check([q for pair in B for q in pair]==d['raw_projection_ranks'][k],'tail_projection_matches_accepted_E075');proj.append(tuple(B))
 def canonical(B):
  pref=sorted((lo,hi) for lo,hi in B if lo==0 and hi<20)
  suff=sorted((lo,hi) for lo,hi in B if lo>0)
  return tuple(pref+[(0,20)]*(5-len(pref)-len(suff))+suff)
 def tobox(B):return tuple(tuple(p) for p in B)
 def contains(D,B):return all(c<=x and y<=z for (x,y),(c,z) in zip(B,D))
 families=r['section_families'];check([f['family_id'] for f in families]==list(range(11)),'family_ids')
 check(len({(tuple(f['complement_lower_ranks']),tuple(f['complement_upper_ranks'])) for f in families})==11,'distinct_complement_families')
 structural=[];emptyfamilies=[];all_orbits=0
 for label in ('left','right'):
  for s in r['sides'][label]['slabs']:
   raw=s['active_raw_indices'];check(raw==sorted(set(raw)) and all(0<=k<127 for k in raw),'active_set_format')
   ls={ell[k] for k in raw};swaps=0
   for x in ls:
    for j in range(1,5):
     y=list(x);y[j],y[j+1]=y[j+1],y[j];check(tuple(y) in ls,'S5_adjacent_generator_closure');swaps+=1
   check(swaps==s['S5_adjacent_swap_checks'],'S5_count')
   expected={}
   for k in raw:expected.setdefault((ell[k][0],tuple(sorted(ell[k][1:]))),[]).append(k)
   seen={};canons=set()
   for group in s['label_orbits']:
    groupboxes=set()
    for orb in group:
     key=(orb['normal_length'],tuple(orb['sorted_tail_lengths']));ids=orb['raw_indices'];check(key not in seen and expected[key]==ids,'complete_orbit_partition');seen[key]=ids
     n=factorial(5)//prod(factorial(q) for q in Counter(key[1]).values());check(n==orb['orbit_size']==len(ids),'orbit_multinomial_cardinality')
     B=canonical(proj[ids[0]])
     for k in ids:check(sorted(proj[k])==sorted(proj[ids[0]]) and canonical(proj[k])==B,'orbit_interval_multiset_and_canonical_box')
     groupboxes.add(B);all_orbits+=1
    check(len(groupboxes)==1,'canonical_dedup_group');canons.update(groupboxes)
   check(seen==expected,'all_active_orbits_accounted')
   allboxes=[tobox(B) for B in s['all_ordered_boxes']];kept=[tobox(B) for B in s['retained_ordered_boxes']]
   check(set(allboxes)==canons and len(allboxes)==len(canons),'canonical_whole_union')
   check(len(kept)==len(set(kept)) and set(kept)<=canons and 6<=len(kept)<=8,'retained_boxes_from_true_union')
   deleted=[]
   for z in s['containment_deletions']:
    B,D=tobox(z['box']),tobox(z['covering_retained_box']);check(D in kept and contains(D,B),'saved_final_containment_witness');deleted.append(B)
   check(set(deleted).isdisjoint(kept) and len(deleted)==len(set(deleted)) and set(deleted)|set(kept)==canons,'containment_deletion_complete')
   events=[];L=[0]*5;U=[20]*5
   for B in kept:
    prefix=[(i,hi) for i,(lo,hi) in enumerate(B) if hi<20];suffix=[(i,lo) for i,(lo,hi) in enumerate(B) if lo>0]
    check(bool(prefix)!=bool(suffix),'one_order_statistic_event')
    if prefix:
     check([i for i,c in prefix]==list(range(len(prefix))) and len({c for i,c in prefix})==1,'equal_prefix_thresholds')
     i,c=prefix[-1];events.append({'kind':'prefix','order_index':i+1,'cut_rank':c});L[i]=max(L[i],c)
    else:
     check([i for i,c in suffix]==list(range(5-len(suffix),5)) and len({c for i,c in suffix})==1,'equal_suffix_thresholds')
     i,c=suffix[0];events.append({'kind':'suffix','order_index':i+1,'cut_rank':c});U[i]=min(U[i],c)
   for i in range(1,5):L[i]=max(L[i],L[i-1])
   for i in range(3,-1,-1):U[i]=min(U[i],U[i+1])
   f=families[s['family_id']]
   check(events==s['threshold_events'] and L==s['complement_lower_ranks']==f['complement_lower_ranks'] and U==s['complement_upper_ranks']==f['complement_upper_ranks'],'DeMorgan_ordered_complement_complete')
   structural.append({'side':label,'slab_ids':s['ids'],'active_raw_count':len(raw),'complete_orbits':len(expected),'S5_generator_checks':swaps,'canonical_boxes':len(canons),'retained_boxes':len(kept),'deleted_boxes':len(deleted),'family_id':s['family_id']})
 full=cap
 for _ in range(4):full=A.mul(full,cap)
 for f in families:
  if any(x>=y for x,y in zip(f['complement_lower_ranks'],f['complement_upper_ranks'])):
   check(f['recurrence_rows']==[] and f['integration_states']==0 and A.parse(f['ordered_complement_volume']['coefficients'])==A.Z and A.parse(f['full_section_union_area']['coefficients'])==full,'empty_complement_zero_certificate');emptyfamilies.append(f['family_id'])
 check(emptyfamilies==[0,1],'both_empty_families_preserved')
 # Phase 2: independently establish all normal slabs and the complete adjacent chambers.
 configs+=1
 defs=d['normal_functions'];funcs=[(A.parse(x['c']),Q(x['m'])) for x in defs]
 check(len(funcs)==len(set(funcs))==24 and [x['id'] for x in defs]==list(range(24)),'normal_function_distinct_complete')
 moving={}
 for i,f in enumerate(defs):
  for x in f['labels']:
   if 'raw_index' in x:
    k=x['raw_index'];check(k not in moving and tuple(x['ell'])==ell[k] and x['mass']==prod(ell[k]),'normal_source_label');moving[k]=i
 check(set(moving)==set(range(127)) and funcs[0]==(A.Z,Q(0)) and funcs[1]==(A.Z,Q(1,2)),'normal_faces_and_labels_complete')
 want={(A.Z,Q(0)),(A.Z,Q(1,2))}
 for k,ls in enumerate(ell):
  radius=A.rad(ls);l=ls[0];v=(radius,-Q(l-1,2)) if l%2 else (A.scale(radius,-1),Q(l,2))
  check(funcs[moving[k]]==v,'common_radius_normal_source_formula');want.add(v)
 check(want==set(funcs),'no_missing_normal_function')
 values=[A.add(c,A.scale(t,m)) for c,m in funcs];side_results={}
 for side,label in [(-1,'left'),(1,'right')]:
  s=r['sides'][label]
  def comp(i,j):
   z=A.cmp(values[i],values[j]);return z if z else ((side*funcs[i][1]>side*funcs[j][1])-(side*funcs[i][1]<side*funcs[j][1]))
  key=cmp_to_key(comp);allorder=sorted(range(24),key=key);sel={};used=set()
  for k,ls in enumerate(ell):
   fid=moving[k];lo=0 if ls[0]%2 else max((0,fid),key=key);hi=min((1,fid),key=key) if ls[0]%2 else 1
   check(comp(lo,hi)<0 and A.sign(A.sub(values[hi],values[lo]))>0,'normal_germ_positive_source_box');sel[str(k)]=[lo,hi];used.update((lo,hi))
  order=sorted(used,key=key);check(sel==s['selections']==st['sides'][label]['selections'] and order==s['effective_order']==st['sides'][label]['effective_order'] and allorder==s['all_function_order']==st['sides'][label]['all_function_order'],'normal_germ_exact_binding')
  check(order[0]==0 and order[-1]==1,'halfcell_normal_partition_endpoints')
  currentties=[[i,j] for i,j in zip(allorder,allorder[1:]) if values[i]==values[j]];check(currentties==([[13,6]] if side<0 else [[6,13]]),'unique_current_tie')
  positions={i:j for j,i in enumerate(order)};weighted=[];zero_weights=0;total=A.Z;zero_width=0
  for lo,hi in zip(order,order[1:]):
   speed=funcs[hi][1]-funcs[lo][1]
   if speed==0:zero_weights+=1;continue
   raw=[k for k in range(127) if positions[sel[str(k)][0]]<=positions[lo] and positions[hi]<=positions[sel[str(k)][1]]]
   zero=values[lo]==values[hi];zero_width+=int(zero)
   weighted.append({'ids':[lo,hi],'width_slope':str(speed),'active_raw_indices':raw,'zero_baseline_width':zero})
  check(weighted==[{k:x[k] for k in ('ids','width_slope','active_raw_indices','zero_baseline_width')} for x in s['slabs']],'all_mutually_disjoint_weighted_slabs_complete')
  check(len(weighted)==(13 if side<0 else 15),'required_slab_count')
  for slab in s['slabs']:
   area=A.parse(families[slab['family_id']]['full_section_union_area']['coefficients']);term=A.scale(area,64*Q(slab['width_slope']))
   check(term==A.parse(slab['weighted_derivative_term']['coefficients']),'saved_area_weight_exact');total=A.add(total,term)
  check(total==A.parse(s['exact_derivative']['coefficients']) and A.sign(total)==-side and list(map(str,A.bounds(total)))==s['exact_derivative']['exact_interval'],'complete_derivative_coefficients_and_sign')
  candidates=[]
  for i,j in zip(allorder,allorder[1:]):
   c,m=funcs[i];h,n=funcs[j]
   if m!=n:
    x=A.scale(A.sub(h,c),1/(m-n))
    if A.sign(A.sub(x,t))==side:candidates.append(x)
  endpoint=(max if side<0 else min)(candidates,key=cmp_to_key(A.cmp))
  explicit=A.scale(A.add(A.mono(1,0),A.mono(3,3)),Q(1,2)) if side<0 else A.scale(A.add(A.mono(1,4),A.mono(2,0)),Q(1,2))
  check(endpoint==explicit==A.parse(s['adjacent_endpoint']['coefficients']),'nearest_full_affine_event')
  gaps=[]
  for i,j in zip(allorder,allorder[1:]):
   c,m=funcs[i];h,n=funcs[j];g0=A.sub(values[j],values[i]);g1=A.add(A.sub(h,c),A.scale(endpoint,n-m))
   check(A.sign(g0)>=0 and A.sign(g1)>=0 and (A.sign(g0)>0 or A.sign(g1)>0),'whole_chamber_gap_not_identically_zero')
   gaps.append({'ids':[i,j],'at_current':A.encode(g0),'at_adjacent_endpoint':A.encode(g1)})
  check(gaps==s['complete_order_gap_certificate'],'all23_gap_coefficients')
  loss=A.scale(A.mul(total,A.sub(endpoint,t)),-1);check(A.sign(loss)>0,'strict_loss_to_adjacent_endpoint')
  side_results[label]={'complete_derivative':A.view(total),'adjacent_endpoint':A.view(endpoint),'all_function_count':24,'adjacent_gaps':gaps,'current_ties':currentties,'boundary_ties':[x['ids'] for x in gaps if not x['at_adjacent_endpoint']],'weighted_slabs':len(weighted),'zero_speed_slabs_omitted':zero_weights,'zero_baseline_width_slabs_retained':zero_width,'endpoint_ratio_loss':A.view(loss)}
 check(r['sides']['left']['selections']==r['sides']['right']['selections'],'no_minmax_change_at_current')
 for p,h in hashes.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h
 return {'status':'pass_two_exact_interfaces_pending_saved_polynomial_check','input_sha256':hashes,'a1':r['a1'],'tail_spacing':r['tail_spacing'],'tail_cap':r['tail_cap'],'baseline_complete_E075':ep['new_ideal_ratio'],'tail_geometry_binding':'127 projections tied to independently accepted E075 endpoints; no old source-main or union execution','tail_projection_count':127,'tail_coordinate_intervals_verified':635,'whole_S5_slab_reduction_certificates':structural,'orbit_count_across_slabs':all_orbits,'empty_complement_family_ids':emptyfamilies,'nonempty_complement_families':9,'total_distinct_complement_families':11,'piecewise_polynomial_rows_generated':0,'saved_polynomial_cells_regenerated':0,'normal_source_formula_verified_for_all127':True,'sides':side_results,'uncrossing_proof_scope':'S5-complete interval orbit, all tails prefix/suffix/full, ordered simplex; containment deletes only witnessed subsets','all_family_overlap_handling':'De Morgan complement, not event-area addition','resource':{'new_exact_interface_configurations':configs,'checks_by_type':dict(checks),'checks_total':sum(checks.values()),'new_integration_states':0,'new_section_integrals':0,'static6D_union_calls':0,'old087_DP_replays':0,'author_main_executions':0,'MATLAB_calls':0},'proof_boundary':'Fixed127 witnesses and five equal fixed tails only; neither new C6 nor joint/global optimality','new_C6_bound':False}

try:out.update(run())
except Exception as e:out.update({'status':'stopped','error':repr(e),'traceback':traceback.format_exc(),'resource':{'new_exact_interface_configurations':configs,'checks_by_type':dict(checks),'checks_total':sum(checks.values())}})
finally:
 out.update({'CPU_seconds_before_save':time.process_time()-CPU,'wall_seconds_before_save':time.monotonic()-WALL,'wall_policy':'180-second alarm covers exclusive save and stdout'})
 with OUT.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
 print(json.dumps({'status':out['status'],'pid':out['pid'],'configurations':configs,'checks':sum(checks.values()),'CPU_seconds_including_save':time.process_time()-CPU,'wall_seconds_including_save':time.monotonic()-WALL,'error':out.get('error')},indent=2))
 signal.alarm(0)
 if out['status']!='pass_two_exact_interfaces_pending_saved_polynomial_check':sys.exit(1)
