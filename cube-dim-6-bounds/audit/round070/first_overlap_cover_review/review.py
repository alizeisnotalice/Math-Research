#!/usr/bin/env python3
"""E070: exact two-sided strip mask, new chamber and finite L1; no union."""
from pathlib import Path
from fractions import Fraction as Q
from functools import cmp_to_key
from math import prod
import hashlib,importlib.util,json,os,resource,signal,sys,time
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;OUT=HERE/'review_results.json'
if OUT.exists():raise FileExistsError('Refuse overwrite '+str(OUT))
CPU=time.process_time();WALL=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(30,31))
def alarm(*_):raise TimeoutError('180 second wall cap')
signal.signal(signal.SIGALRM,alarm);signal.alarm(180)
MAIN=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds')
AUTHOR=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-lower-clean/cube-dim-6-bounds/lower/round070/first_overlap_cover')
DPREV=AUTHOR.parents[1]/'round068/adjacent_hidden_chain';EPREV=HERE.parents[1]/'round068/adjacent_hidden_chain_review'
ARITH=MAIN/'audit/round058/a2_next_event_review/review.py';IFACE=MAIN/'lower/round008/spacing_chamber/chamber_interface.json'
configs=0;nodes={};source_checks=0
out={'status':'running','pid':os.getpid(),'automatic_reruns':0}
def node(kind):
 nodes[kind]=nodes.get(kind,0)+1
 assert sum(nodes.values())<=1500 and configs<=12
 assert time.process_time()-CPU<29 and time.monotonic()-WALL<180

def run():
 global configs,source_checks
 paths=[AUTHOR/'manifest.json',AUTHOR/'results.json',AUTHOR/'preflight.py',AUTHOR/'report.md',AUTHOR/'registration.json',DPREV/'results.json',EPREV/'review_results.json',EPREV/'report.md',EPREV/'manifest.json',ARITH,IFACE,HERE/'frozen_068_baseline.json',HERE/'registration.json',Path(__file__)]
 hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
 assert hashes[str(AUTHOR/'manifest.json')]=='e6e836f74d8a6284e22b82657baa86eb6e750f50736dbfd78929cf6de0bb2e67'
 assert hashes[str(AUTHOR/'results.json')]=='19fd191c40f522bfb8d41647a31886935b285fa080135e8754c6531440dea8cb'
 m=json.loads((AUTHOR/'manifest.json').read_text())
 for x in m['artifacts']:assert hashlib.sha256((AUTHOR/x['path']).read_bytes()).hexdigest()==x['sha256']
 for x in m['frozen_dependencies']:assert x['matches'] and hashlib.sha256(Path(x['path']).read_bytes()).hexdigest()==x['frozen_sha256']==x['current_sha256']
 assert m['execution']['status']=='finished' and m['saved_only_check']['pass']
 assert hashes[str(EPREV/'review_results.json')]=='56fc910446b041f64f0b244d0f6d6e2ab340b895a0e39677945b1fff96e227e3'
 assert hashes[str(ARITH)]=='3e32128be4c14330a7e8f96c7c9317c3576a143cb3d5e57f94d0773f5133e40d'
 spec=importlib.util.spec_from_file_location('accepted_independent_arithmetic',ARITH);A=importlib.util.module_from_spec(spec);spec.loader.exec_module(A)
 Z=A.Z;parse=A.parse;add=A.add;sub=A.sub;scale=A.scale;mono=A.mono;constant=A.constant
 cmp=A.cmp;sign=A.sign;mul=A.mul;bounds=A.bounds;view=A.view;encode=A.encode
 d=json.loads((AUTHOR/'results.json').read_text());dp=json.loads((DPREV/'results.json').read_text());ep=json.loads((EPREV/'review_results.json').read_text());iface=json.loads(IFACE.read_text());frozen=json.loads((HERE/'frozen_068_baseline.json').read_text())
 assert d['status']=='finished' and ep['status']=='pass_three_new_hidden_seams_chambers_and_terminal_finite_L1' and frozen['constructed_round']==68
 for p,sha in frozen['input_sha256'].items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==sha
 event=parse(ep['final_parameter']['coefficients']);baseline=parse(ep['new_ideal_ratio']['coefficients']);g=parse(dp['slope']['coefficients'])
 assert event==parse(d['current_s6']['coefficients'])==parse(dp['best_new_candidate']['event_parameter']['coefficients'])
 assert baseline==parse(frozen['complete_ideal_coefficients'])==parse(dp['best_new_candidate']['ratio']['coefficients'])
 assert g==parse(d['current_slope']['coefficients'])
 oldupper=Q(frozen['comparison_upper']);assert oldupper==max(Q(v[1]) for v in frozen['complete_ideal_intervals'].values())==Q(d['baseline_complete_E068_interval'][1])
 defs=dp['complete_functions_reused'];assert defs==d['complete_functions_reused'] and len(defs)==37
 funcs=[(parse(x['c']),Q(x['m'])) for x in defs]
 moving={z['raw_index']:i for i,f in enumerate(defs) for z in f['labels'] if 'raw_index' in z};assert set(moving)==set(range(127))
 lengths=[tuple(x) for x in iface['box_lengths']];index={x:k for k,x in enumerate(lengths)}
 seam=dp['events'][-1];left=seam['left_state']
 assert ep['final_boundary_ties']==seam['tie_edges']==[[31,7],[3,29]]
 assert left==d['current_left_state_reused']
 def val(i,t):return add(funcs[i][0],scale(t,funcs[i][1]))
 def state(t,side):
  global configs
  configs+=1;values=[val(i,t) for i in range(37)]
  def compare(i,j):
   s=cmp(values[i],values[j]);return s if s else ((side*funcs[i][1]>side*funcs[j][1])-(side*funcs[i][1]<side*funcs[j][1]))
  key=cmp_to_key(compare);sels={};active=set()
  for k,ell in enumerate(lengths):
   node('new_germ_box');f=moving[k];lo=0 if ell[1]%2 else max((0,f),key=key);hi=min((1,f),key=key) if ell[1]%2 else 1
   assert compare(lo,hi)<0;sels[str(k)]=[lo,hi];active.update((lo,hi))
  return {'positive_boxes':127,'effective_order':sorted(active,key=key),'selections':sels}
 # E068 did not certify this right germ; it is new work here.
 right=state(event,1);assert right==seam['right_state']==d['current_right_state_reused']
 assert left['positive_boxes']==right['positive_boxes']==127 and left['selections']==right['selections']
 fixed={}
 def projection(k):
  if k not in fixed:
   node('new_fixed_projection_read');b=[]
   for j in (0,2,3,4,5):
    for eid in iface['boxes'][k][j]:
     e=iface['endpoints'][eid];assert parse(e['c1'])==Z;b.append(parse(e['c0']))
   fixed[k]=b
  return fixed[k]
 diagnostics={tuple(x['ids']):x for x in d['tie_diagnostics']};assert set(diagnostics)=={(31,7),(3,29)}
 groups={}
 for i,j in seam['tie_edges']:
  assert val(i,event)==val(j,event)
  groups[(i,j)]=[[int(k) for k,p in left['selections'].items() if z in p] for z in (i,j)]
  assert groups[(i,j)]==diagnostics[(i,j)]['effective_groups']
 assert [len(x) for x in groups[(31,7)]]==[6,1] and groups[(31,7)][1]==[16]
 assert groups[(3,29)]==[[],[]]
 configs+=1
 covers=[];half=[scale(add(constant(1),mono(4,0,Q(1,2))),Q(1,2))]+[scale(add(constant(1),mono(1,0)),Q(1,4))]*4
 for a in groups[(31,7)][0]:
  for b in groups[(31,7)][1]:
   node('new_cover_pair');B=projection(a);D=projection(b)
   I=[(A.maximum(B[2*q],D[2*q]),A.minimum(B[2*q+1],D[2*q+1])) for q in range(5)]
   ell=(1,2,*[2 if x==4 else 1 for x in lengths[a][2:]])
   assert ell in index and prod(ell)==8;c=index[ell];C=projection(c);r=A.rad(ell);assert r==mono(3,0,Q(1,2))
   margins=[]
   for q,(lo,hi) in enumerate(I):
    assert sign(sub(hi,lo))>0
    lowgap=sub(lo,C[2*q]);highgap=sub(C[2*q+1],hi)
    assert sign(lowgap)>=0 and sign(highgap)>=0
    # Non-strict fixed inclusion is sufficient; equal cell faces are recorded.
    margins.append({'axis':[1,3,4,5,6][q],'low_gap':view(lowgap),'high_gap':view(highgap),'low_equal_cell_face':lowgap==Z and lo==Z and C[2*q]==Z,'high_equal_cell_face':highgap==Z and hi==half[q] and C[2*q+1]==half[q]})
   normal=left['selections'][str(c)];h=val(31,event);below=sub(h,val(normal[0],event));above=sub(val(normal[1],event),h)
   assert normal==[9,1]
   assert below==scale(sub(mono(3,0),mono(1,0)),Q(1,2)) and above==scale(sub(mono(1,0),event),Q(1,2))
   assert sign(below)>0 and sign(above)>0
   covers.append({'raw_indices':[a,b],'cover_raw_index':c,'cover_lengths':list(ell),'cover_mass':8,'cover_common_radius':encode(r),'intersection':[[encode(lo),encode(hi)] for lo,hi in I],'fixed_inclusion_margins':margins,'normal_endpoint_ids':normal,'normal_collar_below':view(below),'normal_collar_above':view(above)})
 assert [x['cover_raw_index'] for x in covers]==[19,21,22,25,26,28]
 assert len(covers)==len(diagnostics[(31,7)]['pair_covers'])==6
 for ours,theirs in zip(covers,diagnostics[(31,7)]['pair_covers']):
  assert ours['raw_indices']==theirs['raw_indices'] and ours['cover_raw_index']==theirs['cover_raw_index'] and ours['intersection']==theirs['intersection']
  assert ours['normal_collar_below']['coefficients']==theirs['normal_collar_below']['coefficients'] and ours['normal_collar_above']['coefficients']==theirs['normal_collar_above']['coefficients']
 configs+=1;inactive=[]
 for fid in (3,29):
  for label in defs[fid]['labels']:
   if 'raw_index' not in label:continue
   node('new_inactive_label');k=label['raw_index'];selected=left['selections'][str(k)]
   assert selected[1]==1 and fid not in selected and right['selections'][str(k)]==selected
   margin=sub(val(fid,event),scale(event,Q(1,2)));assert margin==scale(sub(mono(1,0),event),Q(1,2)) and sign(margin)>0
   inactive.append({'function':fid,'raw_index':k,'selected_endpoints':selected,'margin':view(margin)})
 assert len(inactive)==8
 # Complete new right order, paying both known simultaneous ties.
 configs+=1;order=list(seam['left_complete_order']);edges={tuple(x) for x in seam['tie_edges']};blocks=[];block=[order[0]];used=set()
 for i,j in zip(order,order[1:]):
  if (i,j) in edges:block.append(j);used.add((i,j))
  else:blocks.append(block);block=[j]
 blocks.append(block);assert used==edges
 order=[i for b in blocks for i in sorted(b,key=lambda z:funcs[z][1])]
 assert order==d['right_complete_order'] and sorted(order)==list(range(37))
 nxt=scale(add(mono(3,3),mono(5,0)),Q(1,4));assert nxt==parse(d['next_s7']['coefficients'])
 gaps=[];candidates=[]
 for i,j in zip(order,order[1:]):
  node('new_adjacent_pair_and_gap');c,m=funcs[i];h,n=funcs[j]
  gl=add(sub(h,c),scale(event,n-m));gr=add(sub(h,c),scale(nxt,n-m))
  assert sign(gl)>=0 and sign(gr)>=0 and (sign(gl)>0 or sign(gr)>0)
  if sign(gl)==0:assert n>m
  if m>n:
   x=scale(sub(h,c),1/(m-n));assert cmp(x,event)>0;candidates.append((x,i,j))
  gaps.append({'ids':[i,j],'left_gap':encode(gl),'right_gap':encode(gr)})
 assert min((x[0] for x in candidates),key=cmp_to_key(cmp))==nxt
 assert gaps==d['complete_order_gap_certificate']
 assert [{'ids':[i,j],'location':encode(x)} for x,i,j in candidates]==[{'ids':x['ids'],'location':x['location']['coefficients']} for x in d['new_future_adjacent_candidates']]
 ties=[x['ids'] for x in gaps if not x['right_gap']];assert ties==d['next_tie_edges']==[[11,26]]
 effective=[]
 for i,j in zip(right['effective_order'],right['effective_order'][1:]):
  node('new_effective_gap');c,m=funcs[i];h,n=funcs[j];gl=add(sub(h,c),scale(event,n-m));gr=add(sub(h,c),scale(nxt,n-m))
  assert sign(gl)>=0 and sign(gr)>=0 and (sign(gl)>0 or sign(gr)>0)
  effective.append({'ids':[i,j],'left_gap':encode(gl),'right_gap':encode(gr)})
 nl=state(nxt,-1);assert nl==right==d['next_states']['left']
 # Explicit rational TWO-SIDED collar: all changed endpoint strips lie in C.
 delta=Q(1,10**6);prevleft=parse(ep['new_hidden_chain'][-1]['left_parameter']['coefficients'])
 assert sign(sub(sub(event,prevleft),constant(delta)))>0 and sign(sub(sub(nxt,event),constant(delta)))>0
 collar=[];r8=mono(3,0,Q(1,2))
 def geometry(k,t,kind,expected=None):
  global source_checks
  node(kind);ell=lengths[k];r=A.rad(ell);sp=[add(constant(1),mono(4,0,Q(1,2))),t]+[scale(add(constant(1),mono(1,0)),Q(1,2))]*4;b=[]
  for j,(l,a) in enumerate(zip(ell,sp)):
   q=-(l//2) if l%2 else 1-l//2;lo=A.maximum(Z,sub(scale(a,q+l-1),r));hi=A.minimum(scale(a,Q(1,2)),add(scale(a,q),r))
   assert sign(sub(hi,lo))>0 and sign(lo)>=0 and sign(sub(scale(a,Q(1,2)),hi))>=0
   for x in (lo,hi):
    for z in (q,q+l-1):
     v=sub(x,scale(a,z));assert sign(add(r,v))>=0 and sign(sub(r,v))>=0;source_checks+=1
   if expected is not None and j==1:assert [lo,hi]==[val(i,t) for i in expected]
   if j!=1:
    actual=[]
    for eid in iface['boxes'][k][j]:
     e=iface['endpoints'][eid];assert parse(e['c1'])==Z;actual.append(parse(e['c0']))
    assert [lo,hi]==actual
   b.extend((lo,hi))
  return b
 for side in (-1,1):
  configs+=1;t=add(event,constant(side*delta));endpointchecks=[]
  for record in covers:
   c=record['cover_raw_index'];box=geometry(c,t,'new_cover_source_box_at_collar',[9,1]);assert A.rad(lengths[c])==r8
   margins=[]
   for fid in (31,7):
    lo=sub(val(fid,t),box[2]);hi=sub(box[3],val(fid,t));assert sign(lo)>0 and sign(hi)>0
    margins.append({'function':fid,'from_lower':view(lo),'from_upper':view(hi)})
   endpointchecks.append({'cover_raw_index':c,'both_swapped_endpoints_strictly_inside':margins})
  unused=[]
  for x in inactive:
   node('new_inactive_collar_check');gap=sub(val(x['function'],t),scale(t,Q(1,2)));assert sign(gap)>0;unused.append({'function':x['function'],'raw_index':x['raw_index'],'gap':view(gap)})
  collar.append({'side':side,'parameter':view(t),'cover_checks':endpointchecks,'inactive_checks':unused})
 # Affine margins at both collar endpoints imply the ENTIRE intervening strip
 # is covered for every s in this neighborhood. Fixed I inclusion is exact.
 configs+=1;minwidth=None;maxindex=0
 for k,ell in enumerate(lengths):
  box=geometry(k,nxt,'new_s7_source_box',nl['selections'][str(k)])
  for j,l in enumerate(ell):
   w=sub(box[2*j+1],box[2*j]);assert sign(sub(constant(1),w))>0;minwidth=w if minwidth is None else A.minimum(minwidth,w)
   q=-(l//2) if l%2 else 1-l//2;maxindex=max(maxindex,abs(q),abs(q+l-1))
 result=add(baseline,mul(g,sub(nxt,event)));assert result==parse(d['candidate']['coefficients'])
 lo,hi=bounds(result);alo,ahi=map(Q,d['candidate']['exact_interval']);assert alo<=lo<=hi<=ahi and lo>oldupper
 configs+=1;eps=Q(1,10**22);eta=eps/2;theta=1-Q(1,10**16);K=10**16;target=Q(2071993199356,10**12)
 assert sign(sub(minwidth,constant(2*eps)))>0 and maxindex==2
 loss=64*127*12*eps;inside=lo-loss;factor=Q(2*K+1,2*K+5)**6;finite=theta*inside*factor
 assert finite>target>oldupper and 0<theta<1 and configs<=12
 return {'status':'pass_complete_cover_mask_new_interval_and_finite_L1','input_sha256':hashes,'baseline_round':68,'s6':ep['final_parameter'],'simultaneous_ties':seam['tie_edges'],'effective_groups':{'31,7':groups[(31,7)],'3,29':groups[(3,29)]},'new_s6_right_germ_verified':True,'all_127_minmax_choices_unchanged':True,
 'six_complete_cover_pairs':covers,'eight_inactive_labels':inactive,'two_sided_strip_collar':{'parameter_radius':str(delta),'inside_preceding_and_new_certified_chambers':True,'endpoint_checks':collar,'proof':'Affine strict normal margins at both endpoints cover both swapping endpoints and every point between them throughout the collar; fixed projection inclusion allows equality.'},
 's7':view(nxt),'new_complete_order_gaps':gaps,'new_effective_gaps':effective,'next_ties':ties,'next_left_state_equals_s6_right':True,'new_endpoint_all_127_boxes_valid':True,
 'new_ideal_ratio':view(result),'complete_current068_upper':str(oldupper),'new_ideal_margin_over_complete068_upper':str(lo-oldupper),
 'finite_L1':{'epsilon':str(eps),'bump_halfwidth':str(eta),'threshold':str(theta),'receiver_index_radius':K,'source_index_radius':K+2,'common_radius':'(product ell)^(1/6)/2','max_relative_source_index':maxindex,'minimum_new_box_width':view(minwidth),'all_box_widths_below_one':True,'union_loss_upper':str(loss),'inward_ratio_lower':str(inside),'boundary_factor':str(factor),'weak_ratio_lower':str(finite),'weak_ratio_display_lower':A.dec(finite),'strict_rational_target':str(target),'strict_rational_target_unreduced':'2071993199356/1000000000000','margin_over_target':str(finite-target),'margin_over_target_display':A.dec(finite-target),'margin_over_complete068_upper':str(finite-oldupper)},
 'root_denominator':A.DEN,'resource':{'configurations':configs,'local_nodes':sum(nodes.values()),'nodes_by_type':nodes,'new_scalar_source_endpoint_checks':source_checks,'full_union_calls':0,'old23_pair_or_gap_or_baseline_replays':0,'MATLAB_calls':0},'C6_registration_performed':False,'post_s7_slope':'not checked; no extrapolation beyond s7','proof_status':'complete Fubini strip mask and new endpoint finite L1, with E068 inherited; main owns final acceptance and ledger'}

try:out.update(run())
except Exception as e:out.update({'status':'stopped','error':repr(e),'resource':{'configurations':configs,'local_nodes':sum(nodes.values()),'nodes_by_type':nodes,'new_scalar_source_endpoint_checks':source_checks}})
finally:
 signal.alarm(0);out.update({'cpu_seconds':time.process_time()-CPU,'wall_seconds':time.monotonic()-WALL})
 with OUT.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
 print(json.dumps({'status':out['status'],'pid':out['pid'],'cpu_seconds':out['cpu_seconds'],'wall_seconds':out['wall_seconds'],'resource':out['resource'],'ideal_display':out.get('new_ideal_ratio',{}).get('display_lower'),'finite_display':out.get('finite_L1',{}).get('weak_ratio_display_lower'),'target':out.get('finite_L1',{}).get('strict_rational_target_unreduced'),'margin':out.get('finite_L1',{}).get('margin_over_target_display'),'error':out.get('error')},indent=2))
 if out['status']!='pass_complete_cover_mask_new_interval_and_finite_L1':sys.exit(1)
