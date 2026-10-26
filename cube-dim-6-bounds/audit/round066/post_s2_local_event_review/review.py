#!/usr/bin/env python3
"""E066: new s2 seam, one new chamber, s3 geometry and finite L1; no union."""
from pathlib import Path
from fractions import Fraction as Q
from functools import cmp_to_key
import hashlib,importlib.util,json,os,resource,signal,sys,time
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;OUT=HERE/'review_results.json'
if OUT.exists():raise FileExistsError('Refuse to overwrite frozen output '+str(OUT))
CPU=time.process_time();WALL=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(30,31))
def alarm(*_):raise TimeoutError('180 second wall cap')
signal.signal(signal.SIGALRM,alarm);signal.alarm(180)
MAIN=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds')
AUTHOR=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-lower-clean/cube-dim-6-bounds/lower/round066/post_s2_local_event')
DPREV=AUTHOR.parents[1]/'round064/a2_one_more_event'
EPREV=HERE.parents[1]/'round064/a2_one_more_event_review'
ARITH=MAIN/'audit/round058/a2_next_event_review/review.py'
IFACE=MAIN/'lower/round008/spacing_chamber/chamber_interface.json'
configs=0;nodes={};source_checks=0
out={'status':'running','pid':os.getpid(),'automatic_reruns':0}
def node(kind):
 nodes[kind]=nodes.get(kind,0)+1
 assert sum(nodes.values())<=1500 and configs<=12
 assert time.process_time()-CPU<29 and time.monotonic()-WALL<180

def run():
 global configs,source_checks
 paths=[AUTHOR/'manifest.json',AUTHOR/'results.json',AUTHOR/'preflight.py',AUTHOR/'report.md',AUTHOR/'registration.json',
 DPREV/'results.json',EPREV/'review_results.json',EPREV/'report.md',EPREV/'manifest.json',ARITH,IFACE,
 MAIN/'lower/round064/primary/results.json',MAIN/'lower/round064/primary/run_receipt.json',
 HERE/'frozen_064_baseline.json',HERE/'registration.json',Path(__file__)]
 hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
 manifest=json.loads((AUTHOR/'manifest.json').read_text())
 for item in manifest['artifacts']:assert hashlib.sha256((AUTHOR/item['path']).read_bytes()).hexdigest()==item['sha256']
 for item in manifest['dependencies']:assert item['matches'] and hashlib.sha256(Path(item['path']).read_bytes()).hexdigest()==item['frozen_sha256']==item['current_sha256']
 assert manifest['execution']['status']=='finished' and manifest['saved_only_check']['pass']
 assert hashes[str(AUTHOR/'manifest.json')]=='15f5d02857bd74d3bac268edee7c6144db22883ce903db5457ff65473a0a8954'
 assert hashes[str(EPREV/'review_results.json')]=='0d4ae7b7dc39b80d51e71ccf6044d481caf80523b1d155020cb7a532e20faf9a'
 assert hashes[str(DPREV/'results.json')]=='ab1e3a37722ce9779e8dd81b5004df3d3e47c8905fc36c1e3bf2666cd699bd0b'
 assert hashes[str(ARITH)]=='3e32128be4c14330a7e8f96c7c9317c3576a143cb3d5e57f94d0773f5133e40d'
 # Arithmetic-only import: no E058 run/main, old inventory, geometry or union.
 spec=importlib.util.spec_from_file_location('accepted_E058_arithmetic',ARITH)
 A=importlib.util.module_from_spec(spec);spec.loader.exec_module(A)
 Z=A.Z;parse=A.parse;add=A.add;sub=A.sub;scale=A.scale;mono=A.mono;constant=A.constant
 cmp=A.cmp;sign=A.sign;mul=A.mul;bounds=A.bounds;view=A.view;encode=A.encode
 d=json.loads((AUTHOR/'results.json').read_text());dp=json.loads((DPREV/'results.json').read_text())
 ep=json.loads((EPREV/'review_results.json').read_text());iface=json.loads(IFACE.read_text())
 frozen=json.loads((HERE/'frozen_064_baseline.json').read_text())
 assert d['status']=='finished' and ep['status']=='pass_new_interval_and_finite_L1'
 assert frozen['accepted_constructed_round']==64
 for p,sha in frozen['input_sha256'].items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==sha
 s2=parse(ep['s2']['coefficients']);g=parse(dp['slope']['coefficients']);baseline=parse(ep['new_ideal_ratio']['coefficients'])
 assert s2==parse(dp['s2']['coefficients'])==parse(d['s2']['coefficients'])
 assert g==parse(d['slope']['coefficients']) and baseline==parse(dp['candidate']['coefficients'])==parse(frozen['complete_ideal_coefficients'])
 funcs=[(parse(x['c']),Q(x['m'])) for x in dp['complete_functions']]
 assert len(funcs)==37 and d['full_functions_reused_from_064']==dp['complete_functions']
 # Definitions already bound and accepted in E064: only reuse them here.
 moving={label['raw_index']:i for i,f in enumerate(dp['complete_functions']) for label in f['labels'] if 'raw_index' in label}
 assert set(moving)==set(range(127)) and len(iface['boxes'])==len(iface['box_lengths'])==127
 saved=dp['s2_states']['right'];savedsel=[tuple(saved['selections'][str(k)]) for k in range(127)]
 assert d['s2_right_effective_state_reused_from_064']==saved
 assert dp['s2_crossing_adjacent_ids']==[ep['s2_crossing_pair']]==[[19,14]]
 assert dp['s2_states']['left']['selections']==saved['selections']
 assert ep['s2_states']['left']['selections_unchanged_from_accepted_s1_right'] and ep['s2_states']['right']['selections_unchanged_from_accepted_s1_right']
 # The only new seam: all 1x18 actual selected crossing-box projections.
 configs+=1
 groups=[[k for k,p in enumerate(savedsel) if i in p] for i in (19,14)]
 assert [len(x) for x in groups]==[1,18] and not set(groups[0])&set(groups[1])
 fixed={}
 for k in sorted(set(groups[0]+groups[1])):
  node('new_projection_box_read');coords=[]
  for axis in (0,2,3,4,5):
   for eid in iface['boxes'][k][axis]:
    end=iface['endpoints'][eid];assert parse(end['c1'])==Z;coords.append(parse(end['c0']))
  fixed[k]=coords
 expectedgap=sub(mono(1,0),mono(2,3,Q(1,2)));assert sign(expectedgap)>0
 certs=[]
 for a in groups[0]:
  for b in groups[1]:
   node('new_s2_projection_pair');B=fixed[a];C=fixed[b];separation=None
   for j in range(5):
    gap=sub(A.maximum(B[2*j],C[2*j]),A.minimum(B[2*j+1],C[2*j+1]))
    if sign(gap)>0:separation=(j,gap);break
   assert separation is not None
   j,gap=separation;axis=[1,3,4,5,6][j]
   assert gap==expectedgap
   certs.append({'a':a,'b':b,'axis':axis,'gap':view(gap),'lengths':[iface['box_lengths'][a],iface['box_lengths'][b]]})
 assert len(certs)==18
 for mine,theirs in zip(certs,d['s2_swapped_projection_pairs']):
  assert not theirs['nonempty']
  assert (mine['a'],mine['b'],mine['axis'])==(theirs['a'],theirs['b'],theirs['separating_original_axis'])
  assert mine['gap']['coefficients']==theirs['positive_gap']['coefficients']
 # Recover the new right order by the single already accepted boundary tie.
 configs+=1
 order=list(dp['all_function_order_s1_right']);ia,ib=order.index(19),order.index(14)
 assert ib==ia+1 and funcs[19][1]>funcs[14][1]
 order[ia],order[ib]=order[ib],order[ia]
 assert order==d['complete_s2_right_order'] and sorted(order)==list(range(37))
 s3=add(mono(0,5,Q(1,6)),mono(4,1,Q(1,3)))
 assert s3==parse(d['s3']['coefficients']) and sign(sub(s3,s2))>0
 candidates=[];gaps=[]
 for i,j in zip(order,order[1:]):
  node('new_adjacent_pair_and_gap');c,m=funcs[i];h,n=funcs[j]
  gap2=add(sub(h,c),scale(s2,n-m));gap3=add(sub(h,c),scale(s3,n-m))
  assert sign(gap2)>=0 and sign(gap3)>=0 and (sign(gap2)>0 or sign(gap3)>0)
  if sign(gap2)==0:assert n>m
  if m>n:
   x=scale(sub(h,c),1/(m-n));assert cmp(x,s2)>0;candidates.append((x,i,j))
  gaps.append({'ids':[i,j],'s2_gap':encode(gap2),'s3_gap':encode(gap3)})
 assert min((x[0] for x in candidates),key=cmp_to_key(cmp))==s3
 assert gaps==d['complete_order_gap_certificate']
 assert [{'ids':[i,j],'location':encode(x)} for x,i,j in candidates]==[{'ids':x['ids'],'location':x['location']['coefficients']} for x in d['future_adjacent_candidates']]
 ties=[x['ids'] for x in gaps if not x['s3_gap']];assert ties==d['s3_tying_adjacent_ids']==[[36,16]]
 assert funcs[36]==(scale(mono(4,1),-1),Q(2)) and funcs[16]==(mono(0,5,Q(1,2)),Q(-1))
 effective=[]
 for i,j in zip(saved['effective_order'],saved['effective_order'][1:]):
  node('new_effective_gap');c,m=funcs[i];h,n=funcs[j]
  gap2=add(sub(h,c),scale(s2,n-m));gap3=add(sub(h,c),scale(s3,n-m))
  assert sign(gap2)>=0 and sign(gap3)>=0 and (sign(gap2)>0 or sign(gap3)>0)
  effective.append({'ids':[i,j],'s2_gap':encode(gap2),'s3_gap':encode(gap3)})
 assert len(effective)==21
 states={}
 for side,label in [(-1,'left'),(1,'right')]:
  configs+=1;values=[add(c,scale(s3,m)) for c,m in funcs]
  def compare(i,j):
   z=cmp(values[i],values[j]);return z if z else ((side*funcs[i][1]>side*funcs[j][1])-(side*funcs[i][1]<side*funcs[j][1]))
  key=cmp_to_key(compare);selected=[];active=set()
  for k,ell in enumerate(iface['box_lengths']):
   node('new_s3_germ_box');f=moving[k]
   lo=0 if ell[1]%2 else max((0,f),key=key);hi=min((1,f),key=key) if ell[1]%2 else 1
   assert compare(lo,hi)<0;selected.append((lo,hi));active.update((lo,hi))
  eo=sorted(active,key=key);assert selected==savedsel and len(eo)==22
  if side==-1:assert eo==saved['effective_order']
  actual={'positive_boxes':127,'effective_order':eo,'selections':{str(k):list(p) for k,p in enumerate(selected)}}
  assert actual==d['s3_states'][label]
  states[label]={'positive_boxes':127,'effective_order':eo,'selections_unchanged_from_accepted_s2_right':True}
 # New s3 source geometry; no source geometry is rerun at s2 or before.
 configs+=1
 sp=[add(constant(1),mono(4,0,Q(1,2))),s3]+[add(constant(Q(1,2)),mono(1,0,Q(1,2)))]*4
 minwidth=None;maxindex=0
 for k,ell in enumerate(iface['box_lengths']):
  node('new_s3_source_box');r=A.rad(ell)
  for j,(l,a) in enumerate(zip(ell,sp)):
   q=-(l//2) if l%2 else 1-l//2;maxindex=max(maxindex,abs(q),abs(q+l-1))
   lo=A.maximum(Z,sub(scale(a,q+l-1),r));hi=A.minimum(scale(a,Q(1,2)),add(scale(a,q),r));w=sub(hi,lo)
   assert sign(w)>0 and sign(lo)>=0 and sign(sub(scale(a,Q(1,2)),hi))>=0 and sign(sub(constant(1),w))>0
   minwidth=w if minwidth is None else A.minimum(minwidth,w)
   for x in (lo,hi):
    for n in (q,q+l-1):
     delta=sub(x,scale(a,n));assert sign(add(r,delta))>=0 and sign(sub(r,delta))>=0;source_checks+=1
   if j==1:assert [lo,hi]==[add(funcs[i][0],scale(s3,funcs[i][1])) for i in savedsel[k]]
   else:
    freeze=[]
    for eid in iface['boxes'][k][j]:
     end=iface['endpoints'][eid];assert parse(end['c1'])==Z;freeze.append(parse(end['c0']))
    assert [lo,hi]==freeze
 # Inherited g extends across the newly proved hidden s2 seam.
 result=add(baseline,mul(g,sub(s3,s2)));assert result==parse(d['candidate']['coefficients'])
 lo,hi=bounds(result);alo,ahi=map(Q,d['candidate']['exact_interval']);assert alo<=lo<=hi<=ahi
 oldupper=Q(frozen['comparison_upper'])
 assert oldupper==max(Q(v[1]) for v in frozen['complete_ideal_intervals'].values())
 assert oldupper==Q(dp['candidate']['exact_interval'][1])==Q(d['baseline_complete_064_interval'][1])
 assert lo>oldupper
 configs+=1
 eps=Q(1,10**22);eta=eps/2;theta=1-Q(1,10**16);K=10**16;target=Q(2070186602851,10**12)
 assert sign(sub(minwidth,constant(2*eps)))>0 and maxindex==2
 assert all(sign(sub(constant(1),scale(a,Q(1,2))))>0 for a in sp)
 loss=64*127*12*eps;inside=lo-loss;factor=Q(2*K+1,2*K+5)**6;finite=theta*inside*factor
 assert finite>target>oldupper and 0<theta<1 and configs<=12
 return {'status':'pass_new_s2_seam_new_interval_and_finite_L1','input_sha256':hashes,
 'scope':'only new 18 s2 projection separations, one new chamber, new s3 geometry and finite-L1; accepted E064 baseline inherited',
 's2':ep['s2'],'s2_projection_pair_certificate':certs,'s2_projection_group_sizes':[1,18],
 's2_uniform_separation':view(expectedgap),'s2_hidden_by_complete_crossing_pair_disjointness':True,'s2_right_slope_equals_accepted_left_slope':True,
 's3':view(s3),'delta':view(sub(s3,s2)),'complete_s2_right_order':order,'new_complete_order_gaps':gaps,
 'new_adjacent_candidates':[{'ids':[i,j],'location':view(x)} for x,i,j in candidates],
 'new_effective_order_gaps':effective,'s3_tying_adjacent_ids':ties,'s3_states':states,'new_endpoint_all_127_boxes_valid':True,
 'new_ideal_ratio':view(result),'complete_current064_upper':str(oldupper),'new_ideal_margin_over_complete064_upper':str(lo-oldupper),
 'finite_L1':{'epsilon':str(eps),'bump_halfwidth':str(eta),'threshold':str(theta),'receiver_index_radius':K,'source_index_radius':K+2,'common_radius':'(product ell)^(1/6)/2','max_relative_source_index':maxindex,'minimum_new_box_width':view(minwidth),'all_half_cell_widths_below_one':True,'union_loss_upper':str(loss),'inward_ratio_lower':str(inside),'boundary_factor':str(factor),'weak_ratio_lower':str(finite),'weak_ratio_display_lower':A.dec(finite),'strict_rational_target':str(target),'strict_rational_target_unreduced':'2070186602851/1000000000000','margin_over_target':str(finite-target),'margin_over_target_display':A.dec(finite-target),'margin_over_complete064_upper':str(finite-oldupper)},
 'root_denominator':A.DEN,'resource':{'configurations':configs,'local_nodes':sum(nodes.values()),'nodes_by_type':nodes,'new_scalar_source_endpoint_checks':source_checks,'full_union_calls':0,'old36_gap_chamber_replays':0,'old666_inventory_replays':0,'old22_crossing_replays':0,'old_passed_baseline_recomputations':0,'MATLAB_calls':0},
 'C6_registration_performed':False,'post_s3_slope':'not checked; no extrapolation beyond s3',
 'proof_status':'complete new endpoint ideal union and explicit finite-L1 implication, using accepted E064 dependency; parent owns acceptance/ledger'}

try:out.update(run())
except Exception as e:out.update({'status':'stopped','error':repr(e),'resource':{'configurations':configs,'local_nodes':sum(nodes.values()),'nodes_by_type':nodes,'new_scalar_source_endpoint_checks':source_checks}})
finally:
 signal.alarm(0);out.update({'cpu_seconds':time.process_time()-CPU,'wall_seconds':time.monotonic()-WALL})
 with OUT.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
 print(json.dumps({'status':out['status'],'pid':out['pid'],'cpu_seconds':out['cpu_seconds'],'wall_seconds':out['wall_seconds'],'resource':out['resource'],'ideal_display':out.get('new_ideal_ratio',{}).get('display_lower'),'finite_display':out.get('finite_L1',{}).get('weak_ratio_display_lower'),'target':out.get('finite_L1',{}).get('strict_rational_target_unreduced'),'margin':out.get('finite_L1',{}).get('margin_over_target_display'),'error':out.get('error')},indent=2))
 if out['status']!='pass_new_s2_seam_new_interval_and_finite_L1':sys.exit(1)
