#!/usr/bin/env python3
"""E075 new seam, chamber and finite L1 only; alarm includes cleanup."""
from pathlib import Path
from fractions import Fraction as Q
from functools import cmp_to_key
import hashlib,importlib.util,json,os,resource,signal,sys,time
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;OUT=HERE/'review_results.json'
if OUT.exists():raise FileExistsError('Refuse overwrite '+str(OUT))
CPU=time.process_time();WALL=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(30,31))
def alarm(*_):raise TimeoutError('180 second total wall cap including cleanup')
signal.signal(signal.SIGALRM,alarm);signal.alarm(180)
MAIN=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds')
AUTHOR=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-lower-clean/cube-dim-6-bounds/lower/round075/post_s8_overlap')
DPREV=AUTHOR.parents[1]/'round073/post_s7_mask_diagnosis';EPREV=HERE.parents[1]/'round073/post_s7_mask_review'
ARITH=MAIN/'audit/round058/a2_next_event_review/review.py';IFACE=MAIN/'lower/round008/spacing_chamber/chamber_interface.json'
configs=0;nodes={};source_checks=0
out={'status':'running','pid':os.getpid(),'automatic_reruns':0}
def node(kind):
 nodes[kind]=nodes.get(kind,0)+1
 assert sum(nodes.values())<=1500 and configs<=10
 assert time.process_time()-CPU<29 and time.monotonic()-WALL<179

def run():
 global configs,source_checks
 paths=[AUTHOR/'manifest.json',AUTHOR/'results.json',AUTHOR/'preflight.py',AUTHOR/'report.md',AUTHOR/'registration.json',DPREV/'results.json',EPREV/'review_results.json',EPREV/'report.md',EPREV/'manifest.json',ARITH,IFACE,HERE/'frozen_073_baseline.json',HERE/'registration.json',Path(__file__)]
 hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
 assert hashes[str(AUTHOR/'manifest.json')]=='abc96c935bda6b50528ad0d12138a4476defbef8143da06c907d3c9ecf4e45da'
 assert hashes[str(AUTHOR/'results.json')]=='dec62274abcf8204cf698cd7665dd6c754db3fee87a55beb5bfa1ece676e2817'
 m=json.loads((AUTHOR/'manifest.json').read_text())
 for x in m['artifacts']:assert hashlib.sha256((AUTHOR/x['path']).read_bytes()).hexdigest()==x['sha256']
 for x in m['frozen_dependencies']:assert x['matches'] and hashlib.sha256(Path(x['path']).read_bytes()).hexdigest()==x['frozen_sha256']==x['current_sha256']
 assert m['execution']['status']=='finished' and m['checker_receipt']['pass']
 assert hashes[str(EPREV/'review_results.json')]=='760e0d0d3930acd4a2ffcede7fca464889ba27a814462759933901e994c04c1d'
 assert hashes[str(ARITH)]=='3e32128be4c14330a7e8f96c7c9317c3576a143cb3d5e57f94d0773f5133e40d'
 spec=importlib.util.spec_from_file_location('accepted_independent_arithmetic',ARITH);A=importlib.util.module_from_spec(spec);spec.loader.exec_module(A)
 Z=A.Z;parse=A.parse;add=A.add;sub=A.sub;scale=A.scale;mono=A.mono;constant=A.constant
 cmp=A.cmp;sign=A.sign;mul=A.mul;bounds=A.bounds;view=A.view;encode=A.encode
 d=json.loads((AUTHOR/'results.json').read_text());dp=json.loads((DPREV/'results.json').read_text());ep=json.loads((EPREV/'review_results.json').read_text());iface=json.loads(IFACE.read_text());frozen=json.loads((HERE/'frozen_073_baseline.json').read_text())
 assert d['status']=='finished' and ep['status']=='pass_new_s7_disjoint_seam_chamber_and_s8_finite_L1' and frozen['constructed_round']==73
 for p,sha in frozen['input_sha256'].items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==sha
 event=parse(ep['s8']['coefficients']);baseline=parse(ep['new_ideal_ratio']['coefficients']);g=parse(dp['current_slope']['coefficients'])
 assert event==parse(d['current_s8']['coefficients'])==parse(dp['next_s8']['coefficients'])
 assert baseline==parse(dp['candidate']['coefficients'])==parse(frozen['complete_ideal_coefficients']) and g==parse(d['current_slope']['coefficients'])
 oldupper=Q(frozen['comparison_upper']);assert oldupper==max(Q(v[1]) for v in frozen['complete_ideal_intervals'].values())==Q(d['baseline_complete_E073_interval'][1])
 defs=dp['complete_functions_reused'];assert defs==d['complete_functions_reused'] and len(defs)==37
 funcs=[(parse(x['c']),Q(x['m'])) for x in defs]
 moving={x['raw_index']:i for i,f in enumerate(defs) for x in f['labels'] if 'raw_index' in x};assert set(moving)==set(range(127))
 def val(i,t):return add(funcs[i][0],scale(t,funcs[i][1]))
 def state(t,side):
  global configs
  configs+=1;values=[val(i,t) for i in range(37)]
  def compare(i,j):
   z=cmp(values[i],values[j]);return z if z else ((side*funcs[i][1]>side*funcs[j][1])-(side*funcs[i][1]<side*funcs[j][1]))
  key=cmp_to_key(compare);sels={};active=set()
  for k,ell in enumerate(iface['box_lengths']):
   node('new_germ_box');f=moving[k];lo=0 if ell[1]%2 else max((0,f),key=key);hi=min((1,f),key=key) if ell[1]%2 else 1
   assert compare(lo,hi)<0;sels[str(k)]=[lo,hi];active.update((lo,hi))
  return {'positive_boxes':127,'effective_order':sorted(active,key=key),'selections':sels}
 left=dp['next_states']['left'];assert ep['s8_left_state_equals_s7_right'] and left==d['current_left_state_reused']
 assert ep['next_ties']==dp['next_tie_edges']==d['tie_edges']==[[31,2]]
 # This right germ was deliberately not audited in E073: new calculation.
 right=state(event,1);assert right==dp['next_states']['right']==d['current_right_state_reused']
 assert left['positive_boxes']==right['positive_boxes']==127 and left['selections']==right['selections']
 configs+=1
 groups=[[int(k) for k,p in left['selections'].items() if fid in p] for fid in (31,2)]
 assert groups==[[110,111,112,114,115,117],[0]]
 assert len(d['tie_diagnostics'])==1 and d['tie_diagnostics'][0]['ids']==[31,2] and d['tie_diagnostics'][0]['effective_groups']==groups
 labels=[]
 for fid,group in zip((31,2),groups):
  actual=[]
  for x in defs[fid]['labels']:
   if 'raw_index' not in x:continue
   node('new_effective_label');k=x['raw_index'];assert fid in left['selections'][str(k)] and fid in right['selections'][str(k)]
   assert list(x['lengths'])==iface['box_lengths'][k];actual.append(k)
   labels.append({'function':fid,'raw_index':k,'lengths':iface['box_lengths'][k],'selected_endpoints':left['selections'][str(k)]})
  assert actual==group
 assert len(labels)==7 and all(x['function_is_effectively_selected'] for x in d['tie_diagnostics'][0]['all_labels'])
 assert val(31,event)==val(2,event)
 fixed={}
 for k in groups[0]+groups[1]:
  node('new_fixed_projection_read');b=[]
  for j in (0,2,3,4,5):
   for eid in iface['boxes'][k][j]:
    e=iface['endpoints'][eid];assert parse(e['c1'])==Z;b.append(parse(e['c0']))
  fixed[k]=b
 expectedgap=sub(add(constant(Q(1,2)),mono(1,0)),mono(1,2))
 assert sign(expectedgap)>0
 pairs=[]
 for a in groups[0]:
  for b in groups[1]:
   node('new_disjoint_projection_pair');B=fixed[a];C=fixed[b];I=[(A.maximum(B[2*q],C[2*q]),A.minimum(B[2*q+1],C[2*q+1])) for q in range(5)];sep=None
   for q,(lo,hi) in enumerate(I):
    gap=sub(lo,hi)
    if sign(gap)>0:sep=(q,gap);break
   assert sep is not None;q,gap=sep;assert gap==expectedgap
   pairs.append({'raw_indices':[a,b],'lengths':[iface['box_lengths'][a],iface['box_lengths'][b]],'intersection':[[encode(lo),encode(hi)] for lo,hi in I],'axis':[1,3,4,5,6][q],'gap':view(gap)})
 assert len(pairs)==6==len(d['tie_diagnostics'][0]['pair_diagnostics'])
 for ours,theirs in zip(pairs,d['tie_diagnostics'][0]['pair_diagnostics']):
  assert theirs['classification']=='strictly_disjoint' and ours['raw_indices']==theirs['raw_indices'] and ours['intersection']==theirs['intersection']
  assert ours['axis']==theirs['separating_original_axis'] and ours['gap']['coefficients']==theirs['positive_gap']['coefficients']
 configs+=1
 order=list(dp['right_complete_order']);i,j=order.index(31),order.index(2);assert j==i+1 and funcs[31][1]>funcs[2][1];order[i],order[j]=order[j],order[i]
 assert order==d['right_complete_order'] and sorted(order)==list(range(37))
 nxt=scale(add(constant(1),mono(1,0)),Q(1,2));assert nxt==parse(d['next_s9']['coefficients']) and sign(sub(nxt,event))>0
 gaps=[];candidates=[]
 for i,j in zip(order,order[1:]):
  node('new_adjacent_pair_and_gap');c,m=funcs[i];h,n=funcs[j];gl=add(sub(h,c),scale(event,n-m));gr=add(sub(h,c),scale(nxt,n-m))
  assert sign(gl)>=0 and sign(gr)>=0 and (sign(gl)>0 or sign(gr)>0)
  if sign(gl)==0:assert n>m
  if m>n:
   x=scale(sub(h,c),1/(m-n));assert cmp(x,event)>0;candidates.append((x,i,j))
  gaps.append({'ids':[i,j],'left_gap':encode(gl),'right_gap':encode(gr)})
 assert min((x[0] for x in candidates),key=cmp_to_key(cmp))==nxt and gaps==d['complete_order_gap_certificate']
 assert [{'ids':[i,j],'location':encode(x)} for x,i,j in candidates]==[{'ids':x['ids'],'location':x['location']['coefficients']} for x in d['new_future_adjacent_candidates']]
 ties=[x['ids'] for x in gaps if not x['right_gap']];assert ties==d['next_tie_edges']==[[7,2]]
 effective=[]
 for i,j in zip(right['effective_order'],right['effective_order'][1:]):
  node('new_effective_gap');c,m=funcs[i];h,n=funcs[j];gl=add(sub(h,c),scale(event,n-m));gr=add(sub(h,c),scale(nxt,n-m))
  assert sign(gl)>=0 and sign(gr)>=0 and (sign(gl)>0 or sign(gr)>0);effective.append({'ids':[i,j],'left_gap':encode(gl),'right_gap':encode(gr)})
 nl=state(nxt,-1);assert nl==right==d['next_states']['left']
 configs+=1
 sp=[add(constant(1),mono(4,0,Q(1,2))),nxt]+[scale(add(constant(1),mono(1,0)),Q(1,2))]*4
 minwidth=None;maxindex=0
 for k,ell in enumerate(iface['box_lengths']):
  node('new_s9_source_box');r=A.rad(ell)
  for j,(l,a) in enumerate(zip(ell,sp)):
   q=-(l//2) if l%2 else 1-l//2;maxindex=max(maxindex,abs(q),abs(q+l-1))
   lo=A.maximum(Z,sub(scale(a,q+l-1),r));hi=A.minimum(scale(a,Q(1,2)),add(scale(a,q),r));w=sub(hi,lo)
   assert sign(w)>0 and sign(lo)>=0 and sign(sub(scale(a,Q(1,2)),hi))>=0 and sign(sub(constant(1),w))>0
   minwidth=w if minwidth is None else A.minimum(minwidth,w)
   for x in (lo,hi):
    for z in (q,q+l-1):
     v=sub(x,scale(a,z));assert sign(add(r,v))>=0 and sign(sub(r,v))>=0;source_checks+=1
   if j==1:assert [lo,hi]==[val(i,nxt) for i in nl['selections'][str(k)]]
   else:
    frozen=[]
    for eid in iface['boxes'][k][j]:
     e=iface['endpoints'][eid];assert parse(e['c1'])==Z;frozen.append(parse(e['c0']))
    assert [lo,hi]==frozen
 result=add(baseline,mul(g,sub(nxt,event)));assert result==parse(d['candidate']['coefficients'])
 lo,hi=bounds(result);alo,ahi=map(Q,d['candidate']['exact_interval']);assert alo<=lo<=hi<=ahi and lo>oldupper
 configs+=1;eps=Q(1,10**22);eta=eps/2;theta=1-Q(1,10**16);K=10**16;target=Q(2073015624830,10**12)
 assert sign(sub(minwidth,constant(2*eps)))>0 and maxindex==2
 loss=64*127*12*eps;inside=lo-loss;factor=Q(2*K+1,2*K+5)**6;finite=theta*inside*factor
 assert 0<theta<1 and finite>target>oldupper and configs<=10
 return {'status':'pass_new_s8_disjoint_seam_chamber_and_s9_finite_L1','input_sha256':hashes,'baseline_round':73,'s8':ep['s8'],'all_simultaneous_ties':[[31,2]],'new_s8_right_germ_verified':True,'all_127_minmax_choices_unchanged':True,'all_seven_tie_labels_effective':labels,'complete_effective_groups':groups,'six_disjoint_projection_pairs':pairs,'uniform_separation':view(expectedgap),'existing_cover_candidates_tested':0,
 's9':view(nxt),'new_complete_order_gaps':gaps,'new_effective_gaps':effective,'next_ties':ties,'s9_left_state_equals_s8_right':True,'new_endpoint_all_127_boxes_valid':True,'new_ideal_ratio':view(result),'complete_current073_upper':str(oldupper),'new_ideal_margin_over_complete073_upper':str(lo-oldupper),
 'finite_L1':{'epsilon':str(eps),'bump_halfwidth':str(eta),'threshold':str(theta),'receiver_index_radius':K,'source_index_radius':K+2,'common_radius':'(product ell)^(1/6)/2','max_relative_source_index':maxindex,'minimum_new_box_width':view(minwidth),'all_box_widths_below_one':True,'union_loss_upper':str(loss),'inward_ratio_lower':str(inside),'boundary_factor':str(factor),'weak_ratio_lower':str(finite),'weak_ratio_display_lower':A.dec(finite),'strict_rational_target':str(target),'strict_rational_target_unreduced':'2073015624830/1000000000000','margin_over_target':str(finite-target),'margin_over_target_display':A.dec(finite-target),'margin_over_complete073_upper':str(finite-oldupper)},
 'root_denominator':A.DEN,'resource':{'configurations':configs,'local_nodes':sum(nodes.values()),'nodes_by_type':nodes,'new_scalar_source_endpoint_checks':source_checks,'full_union_calls':0,'old073_070_068_projection_gap_inventory_or_baseline_replays':0,'MATLAB_calls':0},'C6_registration_performed':False,'post_s9_slope':'not checked; no extrapolation beyond s9','proof_status':'complete new closed-endpoint union and finite-L1 from E073 dependency; main owns final acceptance and ledger'}

try:out.update(run())
except Exception as e:out.update({'status':'stopped','error':repr(e),'resource':{'configurations':configs,'local_nodes':sum(nodes.values()),'nodes_by_type':nodes,'new_scalar_source_endpoint_checks':source_checks}})
finally:
 out.update({'cpu_seconds_before_save':time.process_time()-CPU,'wall_seconds_before_save':time.monotonic()-WALL,'wall_policy':'180 second signal remains active through exclusive save and stdout; cancelled only afterwards'})
 with OUT.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
 print(json.dumps({'status':out['status'],'pid':out['pid'],'cpu_seconds_before_save':out['cpu_seconds_before_save'],'wall_seconds_before_save':out['wall_seconds_before_save'],'resource':out['resource'],'ideal_display':out.get('new_ideal_ratio',{}).get('display_lower'),'finite_display':out.get('finite_L1',{}).get('weak_ratio_display_lower'),'target':out.get('finite_L1',{}).get('strict_rational_target_unreduced'),'margin':out.get('finite_L1',{}).get('margin_over_target_display'),'error':out.get('error'),'cpu_seconds_including_save':time.process_time()-CPU,'wall_seconds_including_save':time.monotonic()-WALL},indent=2))
 signal.alarm(0)
 if out['status']!='pass_new_s8_disjoint_seam_chamber_and_s9_finite_L1':sys.exit(1)
