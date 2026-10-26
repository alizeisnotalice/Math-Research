#!/usr/bin/env python3
"""E068: three NEW hidden seams/chambers and terminal finite L1; zero union."""
from pathlib import Path
from fractions import Fraction as Q
from functools import cmp_to_key
from itertools import combinations
import hashlib,importlib.util,json,os,resource,signal,sys,time
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;OUT=HERE/'review_results.json'
if OUT.exists():raise FileExistsError('Refuse overwrite '+str(OUT))
CPU=time.process_time();WALL=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(30,31))
def alarm(*_):raise TimeoutError('180 second wall cap')
signal.signal(signal.SIGALRM,alarm);signal.alarm(180)
MAIN=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds')
AUTHOR=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-lower-clean/cube-dim-6-bounds/lower/round068/adjacent_hidden_chain')
DPREV=AUTHOR.parents[1]/'round066/post_s2_local_event';EPREV=HERE.parents[1]/'round066/post_s2_local_event_review'
ARITH=MAIN/'audit/round058/a2_next_event_review/review.py';IFACE=MAIN/'lower/round008/spacing_chamber/chamber_interface.json'
configs=0;nodes={};source_checks=0
out={'status':'running','pid':os.getpid(),'automatic_reruns':0}
def node(kind):
 nodes[kind]=nodes.get(kind,0)+1
 assert sum(nodes.values())<=2500 and configs<=16
 assert time.process_time()-CPU<29 and time.monotonic()-WALL<180

def run():
 global configs,source_checks
 paths=[AUTHOR/'manifest.json',AUTHOR/'results.json',AUTHOR/'preflight.py',AUTHOR/'report.md',AUTHOR/'registration.json',
 DPREV/'results.json',EPREV/'review_results.json',EPREV/'report.md',EPREV/'manifest.json',ARITH,IFACE,
 MAIN/'lower/round066/primary/results.json',MAIN/'lower/round066/primary/run_receipt.json',
 HERE/'frozen_066_baseline.json',HERE/'registration.json',Path(__file__)]
 hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
 assert hashes[str(AUTHOR/'manifest.json')]=='810cee8830932c3005454de1ab96b07ca19be3e9790f62d896868f481ac73c6f'
 assert hashes[str(AUTHOR/'results.json')]=='06876806d7fb31fe12692450c74212aab661a96cb23335c84ec15f25765a84c7'
 m=json.loads((AUTHOR/'manifest.json').read_text())
 for x in m['artifacts']:assert hashlib.sha256((AUTHOR/x['path']).read_bytes()).hexdigest()==x['sha256']
 for x in m['dependencies']:assert x['matches'] and hashlib.sha256(Path(x['path']).read_bytes()).hexdigest()==x['frozen_sha256']==x['current_sha256']
 assert m['stop_is_successful_protocol_exit'] and m['saved_only_checker']['pass']
 assert hashes[str(EPREV/'review_results.json')]=='3afb4d73216932d5d1f7720df9c80881ad0df8c8f7f02a515ad7eb9ded54086a'
 assert hashes[str(DPREV/'results.json')]=='2a7d51b6b102e4515be4fbbde3e4a59831b7874f48e2f98713de085d290eb327'
 assert hashes[str(ARITH)]=='3e32128be4c14330a7e8f96c7c9317c3576a143cb3d5e57f94d0773f5133e40d'
 spec=importlib.util.spec_from_file_location('accepted_independent_arithmetic',ARITH)
 A=importlib.util.module_from_spec(spec);spec.loader.exec_module(A)
 Z=A.Z;parse=A.parse;add=A.add;sub=A.sub;scale=A.scale;mono=A.mono;constant=A.constant
 cmp=A.cmp;sign=A.sign;mul=A.mul;bounds=A.bounds;view=A.view;encode=A.encode
 d=json.loads((AUTHOR/'results.json').read_text());dp=json.loads((DPREV/'results.json').read_text());ep=json.loads((EPREV/'review_results.json').read_text())
 iface=json.loads(IFACE.read_text());freeze=json.loads((HERE/'frozen_066_baseline.json').read_text())
 assert d['status']=='stopped_overlap' and len(d['events'])==4
 assert ep['status']=='pass_new_s2_seam_new_interval_and_finite_L1' and freeze['accepted_constructed_round']==66
 for p,sha in freeze['input_sha256'].items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==sha
 original=parse(ep['s3']['coefficients']);g=parse(dp['slope']['coefficients']);baseline=parse(ep['new_ideal_ratio']['coefficients'])
 assert original==parse(dp['s3']['coefficients'])==parse(d['baseline_s3']['coefficients'])
 assert baseline==parse(dp['candidate']['coefficients'])==parse(freeze['complete_ideal_coefficients'])
 assert g==parse(d['slope']['coefficients'])
 oldupper=Q(freeze['comparison_upper'])
 assert oldupper==max(Q(v[1]) for v in freeze['complete_ideal_intervals'].values())
 assert oldupper==Q(dp['candidate']['exact_interval'][1])==Q(d['baseline_complete_066_interval'][1])
 defs=dp['full_functions_reused_from_064'];assert d['complete_functions_reused']==defs and len(defs)==37
 funcs=[(parse(x['c']),Q(x['m'])) for x in defs]
 moving={label['raw_index']:i for i,f in enumerate(defs) for label in f['labels'] if 'raw_index' in label}
 assert set(moving)==set(range(127))
 fixed={}
 def projection(k):
  if k not in fixed:
   node('new_projection_box_read');coords=[]
   for j in (0,2,3,4,5):
    for eid in iface['boxes'][k][j]:
     end=iface['endpoints'][eid];assert parse(end['c1'])==Z;coords.append(parse(end['c0']))
   fixed[k]=coords
  return fixed[k]
 def state(event,side):
  global configs
  configs+=1;vals=[add(c,scale(event,m)) for c,m in funcs]
  def compare(i,j):
   z=cmp(vals[i],vals[j]);return z if z else ((side*funcs[i][1]>side*funcs[j][1])-(side*funcs[i][1]<side*funcs[j][1]))
  key=cmp_to_key(compare);sels={};active=set()
  for k,ell in enumerate(iface['box_lengths']):
   node('new_endpoint_germ_box');f=moving[k]
   lo=0 if ell[1]%2 else max((0,f),key=key);hi=min((1,f),key=key) if ell[1]%2 else 1
   assert compare(lo,hi)<0;sels[str(k)]=[lo,hi];active.update((lo,hi))
  return {'positive_boxes':127,'effective_order':sorted(active,key=key),'selections':sels}
 # Exactly the accepted s3 boundary data; no old gap or old state replay.
 event=original;order=list(dp['complete_s2_right_order']);ties=dp['s3_tying_adjacent_ids']
 left,right=dp['s3_states']['left'],dp['s3_states']['right']
 assert ep['s3_tying_adjacent_ids']==ties==[[36,16]]
 assert ep['s3_states']['left']['selections_unchanged_from_accepted_s2_right'] and ep['s3_states']['right']['selections_unchanged_from_accepted_s2_right']
 locations=[original,add(mono(2,4,Q(1,4)),mono(4,0,Q(1,4))),add(mono(1,5,Q(1,4)),mono(3,0,Q(1,4))),sub(mono(1,2),mono(1,0,Q(1,2)))]
 expectedties=[[[36,16]],[[10,27]],[[9,28]],[[31,7],[3,29]]]
 expectedseparations=[sub(sub(scale(add(constant(1),mono(1,0)),Q(3,2)),mono(4,1)),mono(0,5,Q(1,2))),sub(sub(add(constant(1),mono(1,0)),mono(2,4,Q(1,2))),mono(4,0,Q(1,2))),sub(sub(add(constant(1),mono(1,0)),mono(1,5,Q(1,2))),mono(3,0,Q(1,2)))]
 chain=[];candidate=baseline
 for step in range(3):
  rec=d['events'][step];assert rec['index']==step
  assert event==locations[step]==parse(rec['parameter']['coefficients'])
  assert rec['left_complete_order']==order and rec['tie_edges']==ties==expectedties[step]
  assert rec['left_state']==left and rec['right_state']==right
  assert left['positive_boxes']==right['positive_boxes']==127 and left['selections']==right['selections']
  assert rec['classification']=='hidden_strict_fixed_projection_disjointness'
  configs+=1
  edge_set={tuple(x) for x in ties};blocks=[];block=[order[0]];used=set()
  for i,j in zip(order,order[1:]):
   if (i,j) in edge_set:block.append(j);used.add((i,j))
   else:blocks.append(block);block=[j]
  blocks.append(block);assert used==edge_set and sorted(order)==list(range(37))
  # Complete coverage: every pair within every tied block, every selected box
  # carrying either endpoint. Identical min/max choices exclude branch changes.
  paircerts=[];coverage=[]
  for block in blocks:
   if len(block)<2:continue
   for i,j in combinations(block,2):
    assert add(funcs[i][0],scale(event,funcs[i][1]))==add(funcs[j][0],scale(event,funcs[j][1]))
    gi=[int(k) for k,p in left['selections'].items() if i in p]
    gj=[int(k) for k,p in left['selections'].items() if j in p]
    coverage.append({'function_pair':[i,j],'left_boxes':gi,'right_boxes':gj,'required_pair_count':len(gi)*len(gj)})
    for a in gi:
     for b in gj:
      node('new_hidden_projection_pair');B=projection(a);C=projection(b);found=None
      for q in range(5):
       gap=sub(A.maximum(B[2*q],C[2*q]),A.minimum(B[2*q+1],C[2*q+1]))
       if sign(gap)>0:found=(q,gap);break
      assert found is not None,'unpaid crossing projection'
      q,gap=found;assert gap==expectedseparations[step]
      paircerts.append({'function_ids':[i,j],'raw_indices':[a,b],'lengths':[iface['box_lengths'][a],iface['box_lengths'][b]],'axis':[1,3,4,5,6][q],'gap':view(gap)})
  assert len(paircerts)==[1,16,6][step]==sum(x['required_pair_count'] for x in coverage)
  assert len(rec['projection_pair_checks'])==len(paircerts)
  for a,b in zip(paircerts,rec['projection_pair_checks']):
   assert b['strictly_disjoint'] and all(a[k]==b[k] for k in ['function_ids','raw_indices','lengths','axis'])
   assert a['gap']['coefficients']==b['gap']['coefficients']
  order=[i for block in blocks for i in sorted(block,key=lambda k:funcs[k][1])]
  assert order==rec['right_complete_order']
  configs+=1;nxt=locations[step+1];interval=rec['next_interval']
  assert nxt==parse(interval['right_parameter']['coefficients']) and sign(sub(nxt,event))>0
  adjacent=[];gaps=[]
  for i,j in zip(order,order[1:]):
   node('new_adjacent_pair_and_gap');c,m=funcs[i];h,n=funcs[j]
   gl=add(sub(h,c),scale(event,n-m));gr=add(sub(h,c),scale(nxt,n-m))
   assert sign(gl)>=0 and sign(gr)>=0 and (sign(gl)>0 or sign(gr)>0)
   if sign(gl)==0:assert n>m
   if m>n:
    x=scale(sub(h,c),1/(m-n));assert cmp(x,event)>0;adjacent.append((x,i,j))
   gaps.append({'ids':[i,j],'left_gap':encode(gl),'right_gap':encode(gr)})
  assert min((x[0] for x in adjacent),key=cmp_to_key(cmp))==nxt
  assert gaps==interval['complete_order_gaps']
  assert [{'ids':[i,j],'location':encode(x)} for x,i,j in adjacent]==[{'ids':x['ids'],'location':x['location']['coefficients']} for x in interval['adjacent_future_candidates']]
  effective=[]
  for i,j in zip(right['effective_order'],right['effective_order'][1:]):
   node('new_effective_gap');c,m=funcs[i];h,n=funcs[j]
   gl=add(sub(h,c),scale(event,n-m));gr=add(sub(h,c),scale(nxt,n-m))
   assert sign(gl)>=0 and sign(gr)>=0 and (sign(gl)>0 or sign(gr)>0)
   effective.append({'ids':[i,j],'left_gap':encode(gl),'right_gap':encode(gr)})
  assert len(effective)==21
  nextleft=state(nxt,-1);assert nextleft==right==d['events'][step+1]['left_state']
  nextright=state(nxt,1) if step<2 else None
  if step<2:assert nextright==d['events'][step+1]['right_state']
  candidate=add(candidate,mul(g,sub(nxt,event)))
  assert candidate==add(baseline,mul(g,sub(nxt,original)))==parse(interval['candidate']['ratio']['coefficients'])
  lo,hi=bounds(candidate);authorlo,authorhi=map(Q,interval['candidate']['ratio']['exact_interval'])
  assert authorlo<=lo<=hi<=authorhi and lo>oldupper
  nextties=[z['ids'] for z in gaps if not z['right_gap']]
  assert nextties==expectedties[step+1]==d['events'][step+1]['tie_edges']
  chain.append({'index':step,'left_parameter':view(event),'hidden_tie_edges':ties,'all_tied_function_pairs_and_box_coverage':coverage,'projection_pair_certificates':paircerts,'all_minmax_choices_unchanged':True,'right_parameter':view(nxt),'complete_right_order':order,'new_complete_order_gaps':gaps,'new_effective_gaps':effective,'next_boundary_ties':nextties,'next_left_state_equals_current_right':True,'new_endpoint_ideal_ratio':view(candidate)})
  event=nxt;ties=nextties;left=nextleft;right=nextright
 assert event==parse(d['best_new_candidate']['event_parameter']['coefficients'])==locations[3]
 assert candidate==parse(d['best_new_candidate']['ratio']['coefficients'])
 assert ties==[[31,7],[3,29]]
 # Last endpoint: no right state, projection or slope classification is needed.
 configs+=1
 sp=[add(constant(1),mono(4,0,Q(1,2))),event]+[add(constant(Q(1,2)),mono(1,0,Q(1,2)))]*4
 minwidth=None;maxindex=0
 for k,ell in enumerate(iface['box_lengths']):
  node('new_final_source_box');r=A.rad(ell)
  for j,(l,a) in enumerate(zip(ell,sp)):
   q=-(l//2) if l%2 else 1-l//2;maxindex=max(maxindex,abs(q),abs(q+l-1))
   lo=A.maximum(Z,sub(scale(a,q+l-1),r));hi=A.minimum(scale(a,Q(1,2)),add(scale(a,q),r));w=sub(hi,lo)
   assert sign(w)>0 and sign(lo)>=0 and sign(sub(scale(a,Q(1,2)),hi))>=0 and sign(sub(constant(1),w))>0
   minwidth=w if minwidth is None else A.minimum(minwidth,w)
   for x in (lo,hi):
    for n in (q,q+l-1):
     v=sub(x,scale(a,n));assert sign(add(r,v))>=0 and sign(sub(r,v))>=0;source_checks+=1
   if j==1:assert [lo,hi]==[add(funcs[i][0],scale(event,funcs[i][1])) for i in left['selections'][str(k)]]
   else:
    actual=[]
    for eid in iface['boxes'][k][j]:
     end=iface['endpoints'][eid];assert parse(end['c1'])==Z;actual.append(parse(end['c0']))
    assert [lo,hi]==actual
 configs+=1;low,high=bounds(candidate)
 eps=Q(1,10**22);eta=eps/2;theta=1-Q(1,10**16);K=10**16;target=Q(2071938494109,10**12)
 assert sign(sub(minwidth,constant(2*eps)))>0 and maxindex==2
 assert all(sign(sub(constant(1),scale(a,Q(1,2))))>0 for a in sp)
 loss=64*127*12*eps;inside=low-loss;factor=Q(2*K+1,2*K+5)**6;finite=theta*inside*factor
 assert 0<theta<1 and finite>target>oldupper and configs<=16
 return {'status':'pass_three_new_hidden_seams_chambers_and_terminal_finite_L1','input_sha256':hashes,'accepted_baseline_round':66,
 'baseline_s3':ep['s3'],'complete_current066_upper':str(oldupper),'new_hidden_chain':chain,
 'new_hidden_pair_count':sum(len(z['projection_pair_certificates']) for z in chain),'complete_tie_and_minmax_coverage_at_every_crossed_seam':True,
 'final_parameter':view(event),'final_boundary_ties':ties,'final_right_state_or_slope_classification_performed':False,'final_projection_checks':0,
 'new_endpoint_all_127_boxes_valid':True,'new_ideal_ratio':view(candidate),'new_ideal_margin_over_complete066_upper':str(low-oldupper),
 'finite_L1':{'epsilon':str(eps),'bump_halfwidth':str(eta),'threshold':str(theta),'receiver_index_radius':K,'source_index_radius':K+2,'common_radius':'(product ell)^(1/6)/2','max_relative_source_index':maxindex,'minimum_new_box_width':view(minwidth),'all_half_cell_widths_below_one':True,'union_loss_upper':str(loss),'inward_ratio_lower':str(inside),'boundary_factor':str(factor),'weak_ratio_lower':str(finite),'weak_ratio_display_lower':A.dec(finite),'strict_rational_target':str(target),'strict_rational_target_unreduced':'2071938494109/1000000000000','margin_over_target':str(finite-target),'margin_over_target_display':A.dec(finite-target),'margin_over_complete066_upper':str(finite-oldupper)},
 'root_denominator':A.DEN,'resource':{'configurations':configs,'local_nodes':sum(nodes.values()),'nodes_by_type':nodes,'new_scalar_source_endpoint_checks':source_checks,'full_union_calls':0,'old_gap_inventory_seam_or_baseline_replays':0,'MATLAB_calls':0},'C6_registration_performed':False,
 'proof_status':'complete new best endpoint union and finite-L1 implication from accepted E066 and three fully covered new hidden seams; terminal right slope explicitly unknown'}

try:out.update(run())
except Exception as e:out.update({'status':'stopped','error':repr(e),'resource':{'configurations':configs,'local_nodes':sum(nodes.values()),'nodes_by_type':nodes,'new_scalar_source_endpoint_checks':source_checks}})
finally:
 signal.alarm(0);out.update({'cpu_seconds':time.process_time()-CPU,'wall_seconds':time.monotonic()-WALL})
 with OUT.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
 print(json.dumps({'status':out['status'],'pid':out['pid'],'cpu_seconds':out['cpu_seconds'],'wall_seconds':out['wall_seconds'],'resource':out['resource'],'ideal_display':out.get('new_ideal_ratio',{}).get('display_lower'),'finite_display':out.get('finite_L1',{}).get('weak_ratio_display_lower'),'target':out.get('finite_L1',{}).get('strict_rational_target_unreduced'),'margin':out.get('finite_L1',{}).get('margin_over_target_display'),'error':out.get('error')},indent=2))
 if out['status']!='pass_three_new_hidden_seams_chambers_and_terminal_finite_L1':sys.exit(1)
