#!/usr/bin/env python3
"""E079: independent new-section audit; no author code, old union or old gaps."""
from pathlib import Path
from fractions import Fraction as Q
from functools import cmp_to_key
import hashlib, importlib.util, json, os, resource, signal, sys, time, traceback
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent; OUT=HERE/'review_results.json'
if OUT.exists():raise FileExistsError('Refuse overwrite '+str(OUT))
CPU=time.process_time();WALL=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(30,31))
def alarm(*_):raise TimeoutError('180 seconds including cleanup')
signal.signal(signal.SIGALRM,alarm);signal.alarm(180)
MAIN=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds')
DROOT=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-lower-clean/cube-dim-6-bounds/lower')
AUTHOR=DROOT/'round079/s9_true_kink_diagnosis';DPREV=DROOT/'round075/post_s8_overlap'
EPREV=HERE.parents[1]/'round075/post_s8_overlap_review'
ARITH=MAIN/'audit/round058/a2_next_event_review/review.py';IFACE=MAIN/'lower/round008/spacing_chamber/chamber_interface.json'
configs=0;nodes={};containment_tests=0
out={'status':'running','pid':os.getpid(),'automatic_reruns':0}
def node(kind):
 nodes[kind]=nodes.get(kind,0)+1
 assert sum(nodes.values())<=1500 and configs<=10
 assert time.process_time()-CPU<29 and time.monotonic()-WALL<179

def run():
 global configs,containment_tests
 gate=json.loads((HERE/'source_gate.json').read_text()); hashes=dict(gate['input_sha256'])
 for p,h in hashes.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h
 for p in [ARITH,IFACE,EPREV/'review_results.json',EPREV/'manifest.json',EPREV/'report.md',HERE/'registration.json',HERE/'source_gate.json',Path(__file__)]:hashes[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
 assert hashes[str(ARITH)]=='3e32128be4c14330a7e8f96c7c9317c3576a143cb3d5e57f94d0773f5133e40d'
 assert hashes[str(EPREV/'review_results.json')]=='03ee9356878432d22cf88c42889bddbc19a7bd021bb99da654249b088d0f0590'
 assert hashes[str(EPREV/'manifest.json')]=='537ef0370bf9b850f5627d54edfd5c704861e59525da5fc48e6ba989d430b410'
 spec=importlib.util.spec_from_file_location('accepted_independent_exact_arithmetic',ARITH);A=importlib.util.module_from_spec(spec);spec.loader.exec_module(A)
 Z=A.Z;parse=A.parse;encode=A.encode;add=A.add;sub=A.sub;scale=A.scale;mul=A.mul;mono=A.mono;constant=A.constant;cmp=A.cmp;sign=A.sign;bounds=A.bounds;view=A.view
 d=json.loads((AUTHOR/'results.json').read_text());s=json.loads((AUTHOR/'section_results.json').read_text());dp=json.loads((DPREV/'results.json').read_text());ep=json.loads((EPREV/'review_results.json').read_text());iface=json.loads(IFACE.read_text());manifest=json.loads((AUTHOR/'manifest.json').read_text())
 assert d['status']=='stopped_exception' and d['error']=='AssertionError()' and d['new_5D_union_nodes']==0
 assert s['status']=='finished_strict_local_maximum' and s['new_5D_union_cache_states']==s['new_5D_union_nodes']==63
 assert s['old_phase_reruns']==s['old_background_nodes_replayed']==0
 assert ep['status']=='pass_new_s8_disjoint_seam_chamber_and_s9_finite_L1'
 assert Q(iface['normalization_factor'])==64
 canonical=hashlib.sha256(json.dumps(d['retained_cross_sections'],sort_keys=True,separators=(',',':')).encode()).hexdigest()
 assert canonical==gate['saved21_canonical_sha256']=='988b3a0d05ea18d3ce43e22d76be0ff8d283cef6d4c6e986869d221d93bf1ea4'
 event=parse(ep['s9']['coefficients']);start=parse(ep['s8']['coefficients']);g=parse(dp['current_slope']['coefficients'])
 assert event==parse(dp['next_s9']['coefficients'])==parse(d['s9']['coefficients'])==scale(add(constant(1),mono(1,0)),Q(1,2))
 assert parse(ep['new_ideal_ratio']['coefficients'])==parse(dp['candidate']['coefficients'])
 assert g==parse(d['left_slope']['coefficients']) and sign(g)>0
 defs=dp['complete_functions_reused'];assert len(defs)==37 and [x['id'] for x in defs]==list(range(37))
 funcs=[(parse(x['c']),Q(x['m'])) for x in defs]
 moving={x['raw_index']:i for i,f in enumerate(defs) for x in f['labels'] if 'raw_index' in x};assert set(moving)==set(range(127))
 assert all(abs(m)<=2 for c,m in funcs)
 def val(i,t):return add(funcs[i][0],scale(t,funcs[i][1]))
 def state(t,side):
  global configs
  configs+=1;values=[val(i,t) for i in range(37)]
  def compare(i,j):
   r=cmp(values[i],values[j]);return r if r else ((side*funcs[i][1]>side*funcs[j][1])-(side*funcs[i][1]<side*funcs[j][1]))
  key=cmp_to_key(compare);sels={};active=set()
  for k,ell in enumerate(iface['box_lengths']):
   node('new_germ_box');f=moving[k];lo=0 if ell[1]%2 else max((0,f),key=key);hi=min((1,f),key=key) if ell[1]%2 else 1
   assert compare(lo,hi)<0;sels[str(k)]=[lo,hi];active.update((lo,hi))
  return {'positive_boxes':127,'effective_order':sorted(active,key=key),'selections':sels}
 left=dp['next_states']['left'];assert ep['s9_left_state_equals_s8_right'] and ep['next_ties']==dp['next_tie_edges']==d['tie_edges']==[[7,2]]
 assert left==d['current_left_state_reused']
 right=state(event,1);assert right==dp['next_states']['right']==d['current_right_state_reused'] and right['selections']==left['selections']
 configs+=1
 assert [x for x in defs[7]['labels'] if 'raw_index' in x]==[{'raw_index':16,'lengths':[1,2,1,1,1,1],'mass':2,'face':'even_lower'}]
 assert [x for x in defs[2]['labels'] if 'raw_index' in x]==[{'raw_index':0,'lengths':[1,1,1,1,1,1],'mass':1,'face':'odd_upper'}]
 assert left['selections']['16']==[7,1] and left['selections']['0']==[0,2]
 assert funcs[7]==(mono(1,0,Q(-1,2)),Q(1)) and funcs[2]==(constant(Q(1,2)),Q(0))
 assert [k for k,p in left['selections'].items() if 7 in p]==['16'] and [k for k,p in left['selections'].items() if 2 in p]==['0']
 half=constant(Q(1,2));assert val(7,event)==val(2,event)==half
 def fixed(k):
  B=[]
  for j in (0,2,3,4,5):
   ends=[]
   for eid in iface['boxes'][k][j]:
    e=iface['endpoints'][eid];assert parse(e['c1'])==Z;ends.append(parse(e['c0']))
   B.append(tuple(ends))
  return tuple(B)
 P=tuple((A.maximum(a[0],b[0]),A.minimum(a[1],b[1])) for a,b in zip(fixed(0),fixed(16)))
 assert P==((Z,half),)*5 and [[encode(x),encode(y)] for x,y in P]==d['projection_P']
 records=[];positive=[];margins=[];counts={'active':0,'below':0,'above':0};minmargin=None
 saved={x['raw_index']:x for x in d['all_background_normal_classifications']};assert set(saved)==set(range(127))-{0,16}
 def exactmatch(a,v):
  assert a==parse(v['coefficients']);lo,hi=bounds(a);dl,dh=map(Q,v['exact_interval']);assert dl<=lo<=hi<=dh
 for k in range(127):
  if k in (0,16):continue
  node('new_background_classification');ids=left['selections'][str(k)];lo,hi=(val(i,event) for i in ids);x=saved[k]
  assert x['lengths']==iface['box_lengths'][k] and x['normal_endpoint_ids']==ids
  r={'raw_index':k,'normal_endpoint_ids':ids}
  if cmp(lo,half)<0 and cmp(half,hi)<0:
   cls='active';assert x['classification']=='strictly_active_at_height';ms=[sub(half,lo),sub(hi,half)]
   exactmatch(ms[0],x['strict_normal_collar_below']);exactmatch(ms[1],x['strict_normal_collar_above'])
   B=tuple((A.maximum(Z,a),A.minimum(half,b)) for a,b in fixed(k));positive_flag=all(cmp(a,b)<0 for a,b in B)
   assert [[encode(a),encode(b)] for a,b in B]==x['clipped_cross_section_in_P'] and positive_flag==x['positive_5D_cross_section']
   r.update({'cross_section':[[encode(a),encode(b)] for a,b in B],'positive_5D':positive_flag})
   if positive_flag:positive.append((k,B))
  elif cmp(hi,half)<0:
   cls='below';assert x['classification']=='strictly_below_height';ms=[sub(half,hi)];exactmatch(ms[0],x['strict_normal_distance'])
  elif cmp(half,lo)<0:
   cls='above';assert x['classification']=='strictly_above_height';ms=[sub(lo,half)];exactmatch(ms[0],x['strict_normal_distance'])
  else:raise AssertionError('additional selected endpoint at crossing height')
  counts[cls]+=1;r.update({'classification':cls,'normal_margins':[view(a) for a in ms]});records.append(r);margins.extend(bounds(a)[0] for a in ms)
 assert counts=={'active':76,'below':43,'above':6} and len(positive)==44 and len({b for k,b in positive})==44
 assert [k for k,b in positive]==d['positive_background_indices']
 # Verify exact two-sided union inclusion, without depending on the author's pruning logic.
 retained=[tuple((parse(a),parse(b)) for a,b in B) for B in d['retained_cross_sections']]
 assert len(retained)==len(set(retained))==21 and all(B in {b for k,b in positive} for B in retained)
 coverage=[]
 for k,B in positive:
  node('new_section_reduction_candidate');cover=None
  for j,C in enumerate(retained):
   containment_tests+=1
   if all(cmp(c,a)<=0 and cmp(b,d)<=0 for (a,b),(c,d) in zip(B,C)):cover=j;break
  assert cover is not None;coverage.append({'raw_index':k,'retained_index':cover})
 origin=[j for j,B in enumerate(retained) if all(a==Z and sign(b)>0 for a,b in B)];assert origin
 # Reverse coordinates, then build a DAG with an explicit stack. No recursive author slicer.
 configs+=1
 coords=sorted({x for B in retained for p in B for x in p},key=cmp_to_key(cmp));cid={x:i for i,x in enumerate(coords)}
 root=tuple(sorted(tuple((cid[a],cid[b]) for a,b in reversed(B)) for B in retained))
 states=[];index={};pending=[root]
 while pending:
  key=pending.pop()
  if key in index:continue
  node('new_5D_union_state');idx=len(states);index[key]=idx;dim=len(key[0]);edges=[]
  if dim:
   cuts=sorted({q for B in key for q in B[0]})
   for l,h in zip(cuts,cuts[1:]):
    child=tuple(sorted({B[1:] for B in key if B[0][0]<=l and h<=B[0][1]}))
    if child:edges.append((l,h,child));pending.append(child)
  states.append({'key':key,'dimension':dim,'edges':edges})
 values={};certificate=[]
 for idx in sorted(range(len(states)),key=lambda i:states[i]['dimension']):
  st=states[idx];volume=constant(1) if st['dimension']==0 else Z
  for l,h,child in st['edges']:volume=add(volume,mul(sub(coords[h],coords[l]),values[index[child]]))
  values[idx]=volume
 for i,st in enumerate(states):certificate.append({'id':i,'dimension':st['dimension'],'boxes':st['key'],'edges':[{'lo':l,'hi':h,'child':index[ch]} for l,h,ch in st['edges']],'volume':encode(values[i])})
 covered=values[0];uncovered=sub(constant(Q(1,32)),covered);jump=scale(uncovered,-64);gr=add(g,jump)
 for a,v in [(covered,s['exact_background_5D_union']),(uncovered,s['exact_uncovered_5D_volume']),(jump,s['exact_ratio_slope_jump']),(gr,s['right_slope'])]:exactmatch(a,v)
 assert sign(covered)>0 and sign(uncovered)>0 and sign(gr)<0
 # An independent check of the failed origin claim and the successful exposed open box.
 sample=tuple((parse(a),parse(b)) for a,b in s['open_uncovered_region']['exact_endpoints']);separations=[]
 for (a,b),(aa,bb) in zip(sample,s['open_uncovered_region']['rational_inward_open_intervals']):
  assert cmp(a,b)<0 and Q(aa)<Q(bb) and cmp(a,constant(Q(aa)))<=0 and cmp(constant(Q(bb)),b)<=0
 for j,B in enumerate(retained):
  node('new_exposed_box_separation');axes=[q for q,((a,b),(c,d)) in enumerate(zip(sample,B)) if cmp(b,c)<=0 or cmp(d,a)<=0];assert axes;separations.append({'retained_index':j,'original_axis':[1,3,4,5,6][axes[0]]})
 # Only the 36 new right-chamber adjacent pairs are evaluated.
 configs+=1;order=list(dp['right_complete_order']);i,j=order.index(7),order.index(2);assert j==i+1;order[i],order[j]=order[j],order[i]
 assert order==s['right_complete_order'] and sorted(order)==list(range(37))
 nxt=add(mono(1,5,Q(1,6)),mono(5,0,Q(1,3)));assert nxt==parse(s['right_chamber_next_endpoint']['coefficients'])
 gaps=[];future=[]
 for i,j in zip(order,order[1:]):
  node('new_adjacent_pair');c,m=funcs[i];d0,n=funcs[j];gl=add(sub(d0,c),scale(event,n-m));gh=add(sub(d0,c),scale(nxt,n-m))
  assert sign(gl)>=0 and sign(gh)>=0 and (sign(gl)>0 or sign(gh)>0)
  if not sign(gl):assert n>m and [i,j]==[2,7]
  if m>n:
   x=scale(sub(d0,c),1/(m-n));assert cmp(x,event)>0;future.append(x)
  gaps.append({'ids':[i,j],'left_gap':encode(gl),'right_gap':encode(gh)})
 assert min(future,key=cmp_to_key(cmp))==nxt and gaps==s['right_chamber_full_order_gaps']
 ties=[x['ids'] for x in gaps if not x['right_gap']];assert ties==s['next_tie_edges']==[[33,28]]
 assert state(nxt,-1)==right==s['next_states']['left']
 # Explicit real normal collar, tested at both parameter endpoints for every background.
 self_margins=[sub(half,val(0,event)),sub(val(1,event),half)]
 rho=min([Q(d['explicit_rational_common_normal_neighborhood'])]+[q/16 for q in margins]+[bounds(x)[0]/16 for x in self_margins]+[bounds(sub(event,start))[0]/4,bounds(sub(nxt,event))[0]/4])
 assert rho>0;collars=[]
 for side in (-1,1):
  configs+=1;t=add(event,constant(side*rho));minimum=None
  for r in records:
   node('new_collar_background_box');k=r['raw_index'];lo,hi=(val(i,t) for i in right['selections'][str(k)]);ends=[val(7,t),val(2,t)]
   diffs=([sub(x,lo) for x in ends]+[sub(hi,x) for x in ends]) if r['classification']=='active' else ([sub(x,hi) for x in ends] if r['classification']=='below' else [sub(lo,x) for x in ends])
   for v in diffs:
    assert sign(v)>0;minimum=v if minimum is None else A.minimum(minimum,v)
  for k in (0,16):
   lo,hi=(val(i,t) for i in right['selections'][str(k)]);assert cmp(lo,hi)<0
  # Each exchanging interval extends beyond the thin strip on its own side.
  assert cmp(val(0,t),A.minimum(val(7,t),val(2,t)))<0 and cmp(A.maximum(val(7,t),val(2,t)),val(1,t))<0
  collars.append({'side':side,'parameter':view(t),'minimum_strict_background_strip_margin':view(minimum),'all125_pass':True})
 assert configs==7
 # Verify no source changed during the computation.
 for p,h in hashes.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h
 return {'status':'pass_complete_s9_kink_and_local_maximum','input_sha256':hashes,'inherited_E075_complete_ideal':ep['new_ideal_ratio'],'s8':view(start),'s9':view(event),'s10':view(nxt),'unique_simultaneous_tie':[7,2],'new_s9_right_germ':right,'new_s10_left_germ_equals_s9_right':True,'projection_P':[[encode(a),encode(b)] for a,b in P],'projection_P_volume':view(constant(Q(1,32))),'background_counts':counts,'all125_background_records':records,'positive_section_count':44,'retained_section_count':21,'retained_cross_sections':d['retained_cross_sections'],'saved21_canonical_sha256':canonical,'all44_contained_in_retained21':coverage,'origin_covered_by_retained_indices':origin,'false_origin_trial':'Preserved as failed hypothesis, not excluded construction or slope evidence','author_execution_history':manifest['execution_history'],'independent_5D_algorithm':'Explicit-stack DAG, reversed original axes [6,5,4,3,1], all actual endpoint slabs; exact field arithmetic','union_certificate':{'original_axis_order':[6,5,4,3,1],'coordinate_catalog':[view(x) for x in coords],'root_state':0,'states':certificate},'exact_background_5D_union':view(covered),'exact_uncovered_5D_volume':view(uncovered),'left_slope':view(g),'exact_ratio_slope_jump':view(jump),'right_slope':view(gr),'jump_formula':'right-left=-64*(alpha_prime-beta_prime)*full_uncovered_volume; alpha_prime-beta_prime=1','open_uncovered_box':{'exact_endpoints':s['open_uncovered_region']['exact_endpoints'],'rational_inward_open_intervals':s['open_uncovered_region']['rational_inward_open_intervals'],'separations':separations,'not_used_as_full_jump':True},'new_right_complete_order':order,'new_right_chamber_gap_certificate':gaps,'next_ties':ties,'normal_collar_radius':str(rho),'normal_collar_endpoint_certificates':collars,'strict_local_maximum':{'parameter':view(event),'closed_comparison_interval':[view(start),view(nxt)],'unique':True,'family':'Same frozen 127 common-radius source blocks, a1 and four tail spacings fixed, only a2 varies','post_s10_slope':'not classified'},'resource':{'configurations':configs,'local_nodes':sum(nodes.values()),'nodes_by_type':nodes,'scalar_containment_box_comparisons':containment_tests,'new_5D_union_calls':1,'new_5D_union_states':len(states),'full_6D_union_calls':0,'old_passed_gap_projection_union_inventory_replays':0,'MATLAB_calls':0},'new_C6_bound':False,'global_optimum_claim':False,'proof_status':'Independent complete local kink certificate; main retains final review and ledger ownership'}

try:out.update(run())
except Exception as e:out.update({'status':'stopped','error':repr(e),'traceback':traceback.format_exc(),'resource':{'configurations':configs,'local_nodes':sum(nodes.values()),'nodes_by_type':nodes,'scalar_containment_box_comparisons':containment_tests}})
finally:
 out.update({'cpu_seconds_before_save':time.process_time()-CPU,'wall_seconds_before_save':time.monotonic()-WALL,'wall_policy':'180-second alarm remains active through exclusive output save and stdout'})
 with OUT.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
 print(json.dumps({'status':out['status'],'pid':out['pid'],'resource':out['resource'],'CPU_seconds_including_save':time.process_time()-CPU,'wall_seconds_including_save':time.monotonic()-WALL,'error':out.get('error'),'right_slope':out.get('right_slope',{}).get('display_lower')},indent=2))
 signal.alarm(0)
 if out['status']!='pass_complete_s9_kink_and_local_maximum':sys.exit(1)
