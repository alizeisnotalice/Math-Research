#!/usr/bin/env python3
"""E083 saved root/aggregate-extrema/box certificates; no union or author code."""
from pathlib import Path
from fractions import Fraction as Q
from math import prod
import hashlib,importlib.util,json,os,resource,signal,sys,time,traceback
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;OUT=HERE/'review_results.json'
if OUT.exists():raise FileExistsError('Refuse overwrite '+str(OUT))
CPU=time.process_time();WALL=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(30,31))
def alarm(*_):raise TimeoutError('180-second total wall cap including cleanup')
signal.signal(signal.SIGALRM,alarm);signal.alarm(180)
MAIN=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds')
DROOT=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-lower-clean/cube-dim-6-bounds/lower')
AUTHOR=DROOT/'round083/body_center_satellite_preflight';OLD=DROOT/'round075/post_s8_overlap/results.json'
EPREV=HERE.parents[1]/'round075/post_s8_overlap_review/review_results.json'
ARITH=MAIN/'audit/round058/a2_next_event_review/review.py';IFACE=MAIN/'lower/round008/spacing_chamber/chamber_interface.json'
configs=0;axes=0
out={'status':'running','pid':os.getpid(),'automatic_reruns':0}
def budget():
 assert configs<=6 and axes<=18 and time.process_time()-CPU<29 and time.monotonic()-WALL<179
def box(b):return [tuple(map(Q,p)) for p in b]
def volume(b):return prod(hi-lo for lo,hi in b) if all(lo<hi for lo,hi in b) else Q(0)
def ebox(b):return [[str(lo),str(hi)] for lo,hi in b]

def run():
 global configs,axes
 gate=json.loads((HERE/'source_gate.json').read_text());hashes=dict(gate['input_sha256'])
 for p,h in hashes.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h
 for p in [EPREV,ARITH,IFACE,HERE/'registration.json',HERE/'source_gate.json',Path(__file__)]:hashes[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
 assert hashes[str(EPREV)]=='03ee9356878432d22cf88c42889bddbc19a7bd021bb99da654249b088d0f0590'
 assert hashes[str(ARITH)]=='3e32128be4c14330a7e8f96c7c9317c3576a143cb3d5e57f94d0773f5133e40d'
 spec=importlib.util.spec_from_file_location('accepted_independent_arithmetic',ARITH);A=importlib.util.module_from_spec(spec);spec.loader.exec_module(A)
 d=json.loads((AUTHOR/'results.json').read_text());old=json.loads(OLD.read_text());ep=json.loads(EPREV.read_text());iface=json.loads(IFACE.read_text())
 assert d['status']=='finished' and d['union_calls']==d['old_program_replays']==0 and not d['new_C6_bound']
 assert d['new_source_extremum_nodes']==144 and ep['status']=='pass_new_s8_disjoint_seam_chamber_and_s9_finite_L1'
 configs+=1;budget()
 a=A.scale(A.add(A.constant(1),A.mono(1,0)),Q(1,2));spacing=[A.add(A.constant(1),A.mono(4,0,Q(1,2)))]+[a]*5
 assert [A.encode(x) for x in spacing]==d['spacing_coefficients'] and a==A.parse(ep['s9']['coefficients'])==A.parse(old['next_s9']['coefficients'])
 ab=[A.bounds(x) for x in spacing];cell=A.constant(1)
 for x in spacing:cell=A.mul(cell,x)
 cL,cU=A.bounds(cell);dcL,dcU=map(Q,d['period_cell_volume_interval']);assert dcL<=cL<=cU<=dcU
 ideal=A.parse(ep['new_ideal_ratio']['coefficients']);assert ideal==A.parse(old['candidate']['coefficients'])
 eL,eU=map(Q,ep['new_ideal_ratio']['exact_interval']);oL,oU=map(Q,old['candidate']['exact_interval']);assert oL<=eL<=eU<=oU
 assert d['old_full_075_ideal_interval']==old['candidate']['exact_interval']
 oldL=min(oL,eL);oldU=max(oU,eU);assert oldL==oL and oldU==oU
 assert list(map(Q,d['equal_weight_w1']['level1_receiver_ratio_ceiling_interval']))==[dcL/2,dcU/2] and dcU/2<oldL
 headroom=dcL/oldU-1;w=Q(d['derived_weight']);assert headroom==Q(d['certified_positive_headroom']) and w==Q(1,16) and 0<w<=headroom/2<2*w
 safe_sufficient=((1+w)*oldU-oldL)/64;assert safe_sufficient==Q(d['normalization_gain_requirement_for_any_improvement'])
 # This is a saved discrete label/count comparison, not old source geometry or inventory generation.
 labels=[tuple(x) for x in iface['box_lengths']];eligible=[ell for ell in labels if min(ell)>=2]
 best=max(eligible,key=lambda ell:(Q(prod(l-1 for l in ell),prod(ell)),ell))
 assert best==(3,4,4,4,4,4) and labels[126]==best and labels[54]==(2,)*6 and labels[0]==(1,)*6
 assert len(d['candidate_witnesses'])==3
 expected=[('minimal_interior',54,64,1,Q(1025,16)),('max_relative_satellite_count_inherited',126,3072,486,Q(24819,8)),('unit0_plus_single_bodycenter',0,1,1,Q(17,16))]
 results=[];total=Q(0);index_collar=0
 for z,(label,refidx,mwant,kwant,Mwant) in zip(d['candidate_witnesses'],expected):
  configs+=1;budget();assert z['label']==label and z['reference_index']==refidx
  ul,uh,sl,sh=[z[k] for k in ('unit_lo','unit_hi','sat_lo','sat_hi')]
  assert all(type(x)==int for vals in (ul,uh,sl,sh) for x in vals) and all(len(vals)==6 for vals in (ul,uh,sl,sh))
  m=prod(hi-lo+1 for lo,hi in zip(ul,uh));k=prod(hi-lo+1 for lo,hi in zip(sl,sh));M=Q(m)+k*w
  assert m==mwant==z['unit_count'] and k==kwant==z['satellite_count'] and M==Mwant==Q(z['selected_mass'])
  if refidx in (54,126):
   ell=labels[refidx];q=[-(l//2) if l%2 else 1-l//2 for l in ell]
   assert ul==sl==q and uh==[x+l-1 for x,l in zip(q,ell)] and sh==[x+l-2 for x,l in zip(q,ell)]
  else:assert ul==uh==sl==sh==[0]*6
  rL,rU=Q(z['radius_lower']),Q(z['radius_upper']);eta=Q(z['strict_source_margin'])
  assert 0<rL<rU and rL**6<=M/64<rU**6 and rU-rL==Q(1,10**60) and eta==Q(1,10**35)
  radius_test=rL-eta/2;assert radius_test>0 and (2*radius_test)**6<M
  inn,outer,ref=map(box,[z['rational_inward_box'],z['rational_outer_box'],z['reference_old_inner_box']]);cert=[]
  assert len(inn)==len(outer)==len(ref)==6
  for j in range(6):
   axes+=1;budget()
   # Cartesian extrema are attained. The combined min/max suffices for all atoms.
   entries=[Q(ul[j]),Q(uh[j]),Q(sl[j])+Q(1,2),Q(sh[j])+Q(1,2)];mn,mx=min(entries),max(entries)
   rec=z['true_source_extrema'][j];assert rec['axis']==j+1 and rec['unit_index_extrema']==[ul[j],uh[j]] and rec['satellite_integer_index_extrema']==[sl[j],sh[j]] and rec['actual_combined_normalized_coordinate_extrema']==[str(mn),str(mx)]
   mnB=A.bounds(A.scale(spacing[j],mn));mxB=A.bounds(A.scale(spacing[j],mx))
   lo,hi=inn[j];ol,oh=outer[j];il,ih=ref[j]
   assert 0<=lo<hi<=ab[j][0]/2 and ol<=lo<hi<=oh
   # Two aggregate inequalities replace the author's 8 endpoint/extremum pairs.
   source_lower_slack=lo-(mxB[1]-rL+eta);source_upper_slack=(mnB[0]+rL-eta)-hi
   assert source_lower_slack>=0 and source_upper_slack>=0
   assert ol<=max(Q(0),mxB[0]-rU) and oh>=min(ab[j][1]/2,mnB[1]+rU)
   if j==1:
    ids=old['next_states']['left']['selections'][str(refidx)];ends=[A.add(A.parse(old['complete_functions_reused'][i]['c']),A.scale(a,Q(old['complete_functions_reused'][i]['m']))) for i in ids]
   else:
    ends=[]
    for eid in iface['boxes'][refidx][j]:
     e=iface['endpoints'][eid];assert A.parse(e['c1'])==A.Z;ends.append(A.parse(e['c0']))
   oldlo,oldhi=[A.bounds(x) for x in ends]
   assert oldlo[1]+eta<=il<ih<=oldhi[0]-eta
   unit_reflected=[-uh[j],-ul[j]];sat_reflected=[-sh[j]-1,-sl[j]-1]
   assert unit_reflected[1]-unit_reflected[0]==uh[j]-ul[j] and sat_reflected[1]-sat_reflected[0]==sh[j]-sl[j]
   assert Q(sat_reflected[0])+Q(1,2)==-Q(sh[j])-Q(1,2) and Q(sat_reflected[1])+Q(1,2)==-Q(sl[j])-Q(1,2)
   index_collar=max(index_collar,*map(abs,[ul[j],uh[j],sl[j],sh[j],*unit_reflected,*sat_reflected]))
   cert.append({'axis':j+1,'normalized_source_extrema':[str(mn),str(mx)],'physical_source_min_interval':list(map(str,mnB)),'physical_source_max_interval':list(map(str,mxB)),'inward_box_interval':list(map(str,inn[j])),'outer_box_interval':list(map(str,outer[j])),'aggregate_source_slacks':[str(source_lower_slack),str(source_upper_slack)],'old_reference_exact_endpoints':[A.encode(x) for x in ends],'old_reference_endpoint_intervals':[list(map(str,x)) for x in (oldlo,oldhi)],'old_reference_inner_interval':list(map(str,ref[j])),'unit_reflected_range':unit_reflected,'satellite_reflected_integer_range':sat_reflected})
  intersection=[(max(lo,a),min(hi,b)) for (lo,hi),(a,b) in zip(outer,ref)];du=volume(outer)-volume(intersection)
  assert du>=0 and du==Q(z['outer_minus_reference_inner_volume_upper']) and volume(inn)==Q(z['inward_box_volume'])>0
  total+=du
  results.append({'label':label,'reference_index':refidx,'unit_index_box':[ul,uh],'satellite_integer_index_box':[sl,sh],'unit_count':m,'satellite_count':k,'mass':str(M),'critical_radius_enclosure':[str(rL),str(rU)],'root_sixth_power_residuals':[str(M/64-rL**6),str(rU**6-M/64)],'strict_source_margin':str(eta),'one_common_strict_test_radius':str(radius_test),'strict_level1_mass_residual':str(M-(2*radius_test)**6),'axis_certificates':cert,'rational_inward_box':ebox(inn),'rational_outer_box':ebox(outer),'reference_old_inner_box':ebox(ref),'outer_intersection_reference_inner':ebox(intersection),'outer_volume':str(volume(outer)),'intersection_volume':str(volume(intersection)),'increment_volume_upper':str(du),'increment_display':A.dec(du),'actual_increment_not_computed':True})
 configs+=1;budget()
 assert total==Q(d['sum_new_volume_increment_upper']) and index_collar==d['common_integer_index_collar_for_selected_reflections']==2
 family_upper=(oldU+64*total)/(1+w)
 assert family_upper==Q(d['full_three_mixed_box_family_ratio_upper']) and family_upper-oldL==Q(d['three_box_subfamily_upper_minus_old_lower']) and family_upper<oldL
 assert configs==5 and axes==18
 for p,h in hashes.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h
 return {'status':'pass_explicit_three_mixed_subfamily_upper','input_sha256':hashes,'source_model':{'dimension':6,'unit_phase':[0]*6,'satellite_phase':['1/2']*6,'period_cell_mass':'1+w','periods':'a_j in coordinate j; need not be a minimal lattice basis','halfcell_copy_count':64,'ratio_normalization':'64*volume_of_designated_halfcell_witness_union/(1+w)','reflection_unit_index':'n -> -n','reflection_satellite_index':'n -> -n-1','family_scope':'Translate and reflect the H-clipped boxes; not an assertion about every possible center captured by all source subsets'},'dedup_scope':gate['dedup_scope'],'spacing':[A.view(x) for x in spacing],'period_cell_volume':A.view(cell),'author_period_cell_interval':d['period_cell_volume_interval'],'complete_E075_ideal_coefficients':A.encode(ideal),'complete_E075_intervals':{'author':old['candidate']['exact_interval'],'independent':ep['new_ideal_ratio']['exact_interval']},'comparison_lower':str(oldL),'comparison_upper':str(oldU),'equal_weight_level1_ceiling_upper':str(dcU/2),'equal_weight_level1_ceiling_display':A.dec(dcU/2),'equal_weight_exclusion_scope':'w=1 and lambda=1 only; no assertion about another threshold or complete weak ratio','weight':str(w),'cell_mass':str(1+w),'positive_headroom':str(headroom),'derived_weight_rule_pass':True,'frozen_label_selection':{'labels_read':len(labels),'eligible_labels':len(eligible),'maximum_relative_interior_count_length':list(best),'old_geometry_replays':0},'witnesses':results,'common_reflected_integer_index_collar':index_collar,'sum_increment_upper':str(total),'sum_increment_display':A.dec(total),'family_ratio_upper':str(family_upper),'family_ratio_upper_display':A.dec(family_upper),'margin_below_complete_E075_lower':str(oldL-family_upper),'margin_display':A.dec(oldL-family_upper),'actual_union_volume_computed':False,'normalization_threshold_note':{'exact_actual_increment_condition':'Delta > w*R075/64','necessary_threshold_lower':str(w*oldL/64),'threshold_upper':str(w*oldU/64),'author_saved_safe_sufficient_threshold':str(safe_sufficient),'interpretation':'The author key containing any_improvement is a safe sufficient threshold, not an exact or necessary threshold; the family exclusion uses its direct complete upper instead.'},'resource':{'new_small_configurations':configs,'saved_aggregate_axis_certificates':axes,'old144_source_endpoint_nodes_replayed':0,'author_program_imports':0,'union_calls':0,'new_radius_searches':0,'old_scientific_program_replays':0,'MATLAB_calls':0},'new_C6_bound':False,'excluded_scope':'Exactly w=1/16,lambda=1,old127 U plus the three specified H-clipped mixed boxes and designated reflected/translational copies','not_excluded':['other mixed boxes','other masses','other phases','other thresholds','whole two-coset model'],'proof_status':'Independent bounded saved-certificate verification complete; main retains final acceptance'}

try:out.update(run())
except Exception as e:out.update({'status':'stopped','error':repr(e),'traceback':traceback.format_exc(),'resource':{'new_small_configurations':configs,'saved_aggregate_axis_certificates':axes}})
finally:
 out.update({'CPU_seconds_before_save':time.process_time()-CPU,'wall_seconds_before_save':time.monotonic()-WALL,'wall_policy':'180-second alarm covers exclusive save and stdout; cancelled only afterwards'})
 with OUT.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
 print(json.dumps({'status':out['status'],'pid':out['pid'],'resource':out['resource'],'CPU_seconds_including_save':time.process_time()-CPU,'wall_seconds_including_save':time.monotonic()-WALL,'family_upper':out.get('family_ratio_upper_display'),'margin':out.get('margin_display'),'error':out.get('error')},indent=2))
 signal.alarm(0)
 if out['status']!='pass_explicit_three_mixed_subfamily_upper':sys.exit(1)
