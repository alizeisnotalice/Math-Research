#!/usr/bin/env python3
"""E064: new adjacent event, geometry and finite L1 only; no old union replay."""
from pathlib import Path
from fractions import Fraction as Q
from functools import cmp_to_key
from math import prod
import importlib.util,hashlib,json,os,resource,signal,time,sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
OUT=HERE/'review_results.json'
if OUT.exists():raise FileExistsError('Frozen output exists: '+str(OUT))
CPU=time.process_time();WALL=time.monotonic()
resource.setrlimit(resource.RLIMIT_CPU,(30,31))
def alarm(*_):raise TimeoutError('180 second wall cap')
signal.signal(signal.SIGALRM,alarm);signal.alarm(180)
MAIN=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds')
AUTHOR=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-lower-clean/cube-dim-6-bounds/lower/round064/a2_one_more_event')
PREV=MAIN/'audit/round058/a2_next_event_review'
DPREV=MAIN/'lower/round058/a2_next_event_preflight'
IFACE=MAIN/'lower/round008/spacing_chamber/chamber_interface.json'
configs=0;nodes={};source_checks=0
out={'status':'running','pid':os.getpid(),'automatic_reruns':0}
def node(kind):
 nodes[kind]=nodes.get(kind,0)+1
 assert sum(nodes.values())<=1800 and configs<=12
 assert time.process_time()-CPU<29 and time.monotonic()-WALL<180

def run():
 global configs,source_checks
 paths=[AUTHOR/'manifest.json',AUTHOR/'results.json',AUTHOR/'preflight.py',AUTHOR/'report.md',AUTHOR/'registration.json',
 PREV/'review.py',PREV/'review_results.json',PREV/'report.md',PREV/'manifest.json',
 DPREV/'results.json',IFACE,HERE/'frozen_lower_ledger.json',HERE/'registration.json',Path(__file__)]
 hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
 manifest=json.loads((AUTHOR/'manifest.json').read_text())
 for item in manifest['artifacts']:
  assert hashlib.sha256((AUTHOR/item['path']).read_bytes()).hexdigest()==item['sha256']
 for item in manifest['dependencies']:
  assert item['match'] and hashlib.sha256(Path(item['path']).read_bytes()).hexdigest()==item['frozen_sha256']==item['current_sha256']
 assert manifest['execution']['status']=='finished' and manifest['saved_only_checker']['pass']
 assert hashes[str(AUTHOR/'manifest.json')]=='7f1c3cf6187ae91a2102e4a429c944432eeed65a0497ed21ca3e461c7ebe89d9'
 assert hashes[str(PREV/'review.py')]=='3e32128be4c14330a7e8f96c7c9317c3576a143cb3d5e57f94d0773f5133e40d'
 assert hashes[str(PREV/'review_results.json')]=='adf704b2927cfaf7172ce73830371e17c8e20d17ec7b7c028e41e7a6ebaf400a'
 # Reuse previously accepted independent arithmetic only. Its main/run,
 # old inventory, geometry, crossings and union baseline are never invoked.
 spec=importlib.util.spec_from_file_location('accepted_E058_arithmetic',PREV/'review.py')
 A=importlib.util.module_from_spec(spec);spec.loader.exec_module(A)
 Z=A.Z;parse=A.parse;add=A.add;sub=A.sub;scale=A.scale;mono=A.mono;constant=A.constant
 mul=A.mul;sign=A.sign;cmp=A.cmp;bounds=A.bounds;view=A.view;encode=A.encode
 d=json.loads((AUTHOR/'results.json').read_text())
 prev=json.loads((PREV/'review_results.json').read_text())
 dp=json.loads((DPREV/'results.json').read_text())
 ledger=json.loads((HERE/'frozen_lower_ledger.json').read_text())['lower']
 iface=json.loads(IFACE.read_text())
 assert d['status']=='finished' and prev['status']=='pass_independent_local_and_finite_L1_certificate'
 assert ledger['id']=='R6_rectangular_a2_next_event_round058'
 assert ledger['certificate_sha256']==hashes[str(DPREV/'results.json')]
 assert ledger['exact_interval']==dp['next_affine_ideal_candidate']['exact_interval']==d['baseline_058_ideal_interval']
 assert Q(ledger['strict_simple_lower'])==Q(2067759428416,10**12)
 s1=parse(prev['s1']['coefficients']);g=parse(prev['inherited_slope']['coefficients'])
 baseline=parse(prev['new_ideal_ratio']['coefficients'])
 assert s1==parse(d['s1']['coefficients'])==parse(dp['next_raw_location']['coefficients'])
 assert g==parse(d['slope']['coefficients'])==parse(dp['right_slope']['coefficients'])
 assert baseline==parse(dp['next_affine_ideal_candidate']['coefficients'])
 s2=add(constant(Q(1,2)),mono(2,3,Q(1,4)))
 assert s2==parse(d['s2']['coefficients']) and sign(sub(s2,s1))>0
 # Bind every saved function to the exact source-block definition, not labels.
 items=d['complete_functions'];assert [x['id'] for x in items]==list(range(37))
 funcs=[(parse(x['c']),Q(x['m'])) for x in items]
 assert len(set(funcs))==37 and funcs[:2]==[(Z,Q(0)),(Z,Q(1,2))]
 fmap={x:i for i,x in enumerate(funcs)};expected=set(funcs[:2]);raw=[]
 assert len(iface['box_lengths'])==127 and len({tuple(x) for x in iface['box_lengths']})==127
 for ell in iface['box_lengths']:
  node('source_block_function_definition')
  assert len(ell)==6 and all(1<=l<=4 for l in ell)
  r=A.rad(ell);l=ell[1]
  f=(r,-Q(l-1,2)) if l%2 else (scale(r,-1),Q(l,2))
  expected.add(f);assert f in fmap
  raw.append({'ell':ell,'r':r,'moving':fmap[f]})
 assert expected==set(funcs)
 saved=dp['next_states']['right']
 assert d['s1_right_state_reused_from_058']==saved
 savedsel=[tuple(saved['selections'][str(k)]) for k in range(127)]
 # One new complete-order configuration; 36 adjacent pairs are sufficient.
 configs+=1
 order=d['all_function_order_s1_right'];assert sorted(order)==list(range(37))
 candidates=[];gaps=[]
 for i,j in zip(order,order[1:]):
  node('new_adjacent_pair')
  c,m=funcs[i];h,n=funcs[j]
  gap1=add(sub(h,c),scale(s1,n-m));gap2=add(sub(h,c),scale(s2,n-m))
  assert sign(gap1)>=0 and sign(gap2)>=0 and (sign(gap1)>0 or sign(gap2)>0)
  if sign(gap1)==0:assert n>m
  if m>n:
   crossing=scale(sub(h,c),1/(m-n));assert cmp(crossing,s1)>0
   candidates.append((crossing,i,j))
  gaps.append({'ids':[i,j],'s1_gap':encode(gap1),'s2_gap':encode(gap2)})
 assert len(candidates)==15 and min((x[0] for x in candidates),key=cmp_to_key(cmp))==s2
 assert gaps==d['complete_order_gap_certificate']
 assert [{'ids':[i,j],'location':encode(x)} for x,i,j in candidates]==[{'ids':x['ids'],'location':x['location']['coefficients']} for x in d['adjacent_future_candidates']]
 ties=[z['ids'] for z in gaps if not z['s2_gap']]
 assert ties==d['s2_crossing_adjacent_ids']==[[19,14]]
 assert funcs[19]==(constant(-1),Q(1)) and funcs[14]==(mono(2,3,Q(1,2)),Q(-1))
 assert [sum(pair[0]==19 or pair[1]==19 for pair in savedsel),sum(pair[0]==14 or pair[1]==14 for pair in savedsel)]==[1,18]
 effective=[]
 for i,j in zip(saved['effective_order'],saved['effective_order'][1:]):
  node('new_effective_order_gap')
  c,m=funcs[i];h,n=funcs[j]
  gap1=add(sub(h,c),scale(s1,n-m));gap2=add(sub(h,c),scale(s2,n-m))
  assert sign(gap1)>=0 and sign(gap2)>=0 and (sign(gap1)>0 or sign(gap2)>0)
  effective.append({'ids':[i,j],'s1_gap':encode(gap1),'s2_gap':encode(gap2)})
 assert len(effective)==21 and effective==d['closed_chamber_effective_gap_certificate']
 states={}
 for side,label in [(-1,'left'),(1,'right')]:
  configs+=1
  values=[add(c,scale(s2,m)) for c,m in funcs]
  def compare(i,j):
   z=cmp(values[i],values[j])
   return z if z else ((side*funcs[i][1]>side*funcs[j][1])-(side*funcs[i][1]<side*funcs[j][1]))
  key=cmp_to_key(compare);selected=[];active=set()
  for record in raw:
   node('new_s2_germ_box')
   odd=record['ell'][1]%2;k=record['moving']
   lo=0 if odd else max((0,k),key=key);hi=min((1,k),key=key) if odd else 1
   assert compare(lo,hi)<0
   selected.append((lo,hi));active.update((lo,hi))
  eo=sorted(active,key=key)
  assert selected==savedsel and len(eo)==22
  if side==-1:assert eo==saved['effective_order']
  claimed=d['s2_states'][label]
  assert claimed['positive_boxes']==127 and claimed['effective_order']==eo
  assert claimed['selections']=={str(k):list(p) for k,p in enumerate(selected)}
  states[label]={'positive_boxes':127,'effective_order':eo,'selections_unchanged_from_accepted_s1_right':True}
 # New endpoint source geometry, with one common radius for all six axes.
 configs+=1
 sp=[add(constant(1),mono(4,0,Q(1,2))),s2]+[add(constant(Q(1,2)),mono(1,0,Q(1,2)))]*4
 minwidth=None;maxindex=0;boxes=[]
 for k,record in enumerate(raw):
  node('new_s2_source_box')
  ell,r=record['ell'],record['r'];b=[]
  for j,(l,a) in enumerate(zip(ell,sp)):
   q=-(l//2) if l%2 else 1-l//2;maxindex=max(maxindex,abs(q),abs(q+l-1))
   lo=A.maximum(Z,sub(scale(a,q+l-1),r));hi=A.minimum(scale(a,Q(1,2)),add(scale(a,q),r))
   width=sub(hi,lo);assert sign(width)>0 and sign(lo)>=0 and sign(sub(scale(a,Q(1,2)),hi))>=0
   assert sign(sub(constant(1),width))>0
   minwidth=width if minwidth is None else A.minimum(minwidth,width)
   for x in (lo,hi):
    for n in (q,q+l-1):
     v=sub(x,scale(a,n));assert sign(add(r,v))>=0 and sign(sub(r,v))>=0;source_checks+=1
   if j==1:
    assert [lo,hi]==[add(funcs[z][0],scale(s2,funcs[z][1])) for z in savedsel[k]]
   else:
    frozen=[]
    for eid in iface['boxes'][k][j]:
     ep=iface['endpoints'][eid];assert parse(ep['c1'])==Z;frozen.append(parse(ep['c0']))
    assert [lo,hi]==frozen
   b.extend((lo,hi))
  boxes.append(b)
 # Full new value from accepted right slope and proven new open chamber.
 result=add(baseline,mul(g,sub(s2,s1)))
 assert result==parse(d['candidate']['coefficients'])
 low,high=bounds(result);author_lo,author_hi=map(Q,d['candidate']['exact_interval'])
 assert author_lo<=low<=high<=author_hi
 old_upper=Q(ledger['exact_interval'][1])
 assert low>old_upper
 # Erosion and finite L1 arithmetic are newly evaluated at s2.
 configs+=1
 eps=Q(1,10**22);eta=eps/2;theta=1-Q(1,10**16);K=10**16
 target=Q(2068316441022,10**12)
 assert sign(sub(minwidth,constant(2*eps)))>0 and maxindex==2
 assert all(sign(sub(constant(1),scale(a,Q(1,2))))>0 for a in sp)
 loss=64*127*12*eps;inward=low-loss;factor=Q(2*K+1,2*K+5)**6;finite=theta*inward*factor
 assert 0<theta<1 and finite>target>old_upper
 assert configs<=12
 return {'status':'pass_new_interval_and_finite_L1','input_sha256':hashes,
 'scope':'only new s1-to-s2 interval; accepted E058 exact baseline/right-slope/hidden-s1 proof inherited without replay',
 'old_baseline_coefficients_recomputed':False,'old_666_inventory_replayed':False,'old_22_crossing_checks_replayed':False,
 'accepted_s1':prev['s1'],'s2':view(s2),'delta':view(sub(s2,s1)),
 'all_function_order_s1_right':order,'complete_order_gaps':gaps,'nearest_adjacent_candidates':[{'ids':[i,j],'location':view(x)} for x,i,j in candidates],
 'effective_order_gaps':effective,'s2_crossing_pair':ties[0],'s2_crossing_group_sizes':[1,18],'s2_states':states,
 'new_endpoint_all_127_boxes_valid':True,'new_ideal_ratio':view(result),'complete_current058_upper':str(old_upper),
 'new_ideal_margin_over_complete058_upper':str(low-old_upper),
 'finite_L1':{'epsilon':str(eps),'bump_halfwidth':str(eta),'threshold':str(theta),'receiver_index_radius':K,'source_index_radius':K+2,'common_radius':'(product ell)^(1/6)/2','max_relative_source_index':maxindex,'minimum_new_box_width':view(minwidth),'all_half_cell_widths_below_one':True,'union_loss_upper':str(loss),'inward_ratio_lower':str(inward),'boundary_factor':str(factor),'weak_ratio_lower':str(finite),'weak_ratio_display_lower':A.dec(finite),'strict_rational_target':str(target),'strict_rational_target_unreduced':'2068316441022/1000000000000','margin_over_target':str(finite-target),'margin_over_target_display':A.dec(finite-target),'margin_over_complete058_upper':str(finite-old_upper)},
 'root_denominator':A.DEN,'resource':{'configurations':configs,'local_nodes':sum(nodes.values()),'nodes_by_type':nodes,'new_scalar_source_endpoint_checks':source_checks,'full_union_calls':0,'old_passed_baseline_recomputations':0,'MATLAB_calls':0},
 'C6_registration_performed':False,'post_s2_slope':'not checked; no extension beyond s2 is asserted',
 'proof_status':'complete new endpoint ideal union and explicit finite-L1 implication, using accepted E058 dependency; parent owns acceptance/ledger'}

try:
 out.update(run())
except Exception as e:
 out.update({'status':'stopped','error':repr(e),'resource':{'configurations':configs,'local_nodes':sum(nodes.values()),'nodes_by_type':nodes,'new_scalar_source_endpoint_checks':source_checks}})
finally:
 signal.alarm(0)
 out.update({'cpu_seconds':time.process_time()-CPU,'wall_seconds':time.monotonic()-WALL})
 with OUT.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
 print(json.dumps({'status':out['status'],'pid':out['pid'],'cpu_seconds':out['cpu_seconds'],'wall_seconds':out['wall_seconds'],'resource':out['resource'],'new_ideal_display':out.get('new_ideal_ratio',{}).get('display_lower'),'finite_display':out.get('finite_L1',{}).get('weak_ratio_display_lower'),'target':out.get('finite_L1',{}).get('strict_rational_target_unreduced'),'margin':out.get('finite_L1',{}).get('margin_over_target_display'),'error':out.get('error')},indent=2))
 if out['status']!='pass_new_interval_and_finite_L1':sys.exit(1)
