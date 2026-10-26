#!/usr/bin/env python3
"""New two-coset source interface; one derived mass; no union or old replay."""
from pathlib import Path
from fractions import Fraction as Q
from math import prod
import json,hashlib,importlib.util,sys,time,resource,signal,os
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;OUT=HERE/'results.json'
if OUT.exists():raise FileExistsError('Refuse overwrite '+str(OUT))
resource.setrlimit(resource.RLIMIT_CPU,(20,21))
def wallcap(*_):raise TimeoutError('180 second total wall cap')
signal.signal(signal.SIGALRM,wallcap);signal.alarm(180)
CPU,WALL=time.process_time(),time.monotonic()
LOWER=HERE.parents[1];OLD=LOWER/'round075/post_s8_overlap/results.json'
BASE=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds')
ARITH=BASE/'lower/round001_rectangular_certificate.py';IFACE=BASE/'lower/round008/spacing_chamber/chamber_interface.json'
nodes,configs=0,0
out={'status':'running','pid':os.getpid(),'union_calls':0,'old_program_replays':0,'new_C6_bound':False,'automatic_reruns':0}
def budget():
 if time.process_time()-CPU>=19.5 or time.monotonic()-WALL>=179:raise TimeoutError('CPU/wall prelimit')
 if nodes>1000 or configs>8:raise RuntimeError('node/config cap')
try:
 spec=importlib.util.spec_from_file_location('accepted_field',ARITH);A=importlib.util.module_from_spec(spec);spec.loader.exec_module(A)
 A.D=10**80;A.RU,A.RV=A.floor_root(2*A.D**6,6),A.floor_root(3*A.D**6,6)
 A.UL,A.UH=Q(A.RU,A.D),Q(A.RU+1,A.D);A.VL,A.VH=Q(A.RV,A.D),Q(A.RV+1,A.D);A.bounds.cache_clear()
 F,u=A.F,A.u
 def field(z):
  a=[Q(0)]*36
  for ij,x in z.items():i,j=map(int,ij.split(','));a[6*i+j]=Q(x)
  return F(a=a)
 def enc(x):return {f'{i},{j}':str(x.a[6*i+j]) for i in range(6) for j in range(6) if x.a[6*i+j]}
 def interval(x):return A.bounds(x.a)
 def vbox(B):
  if any(b<=a for a,b in B):return Q(0)
  return prod(b-a for a,b in B)
 def intersect(B,C):return [(max(a,c),min(b,d)) for (a,b),(c,d) in zip(B,C)]
 def encodebox(B):return [[str(a),str(b)] for a,b in B]
 old=json.loads(OLD.read_text());iface=json.loads(IFACE.read_text());ells=[tuple(z) for z in iface['box_lengths']]
 spacing=(F(1)+u**4/2,*[(F(1)+u)/2]*5);ab=[interval(z) for z in spacing]
 assert field(old['next_s9']['coefficients'])==spacing[1]
 oldL,oldU=map(Q,old['candidate']['exact_interval'])
 cell=F(1)
 for z in spacing:cell=cell*z
 cellL,cellU=interval(cell)
 configs+=1;budget()
 headroom=cellL/oldU-1;assert 0<headroom
 w=Q(1)
 while headroom/2<w:w/=2
 out.update({'source_model':{'dimension':6,'unit_orbit':'diag(a)*Z^6 mass1','satellite_orbit':'diag(a)*(Z^6+(1/2,...,1/2)) massw','period_cell_mass':'1+w','halfcell_normalization':'64*union_half_volume/(1+w)','distinctness':'integer lattice and half-integer coset disjoint in every coordinate','reflection_satellite_index_rule':'n -> -n-1'},'spacing_coefficients':[enc(z) for z in spacing],'old_full_075_ideal_interval':old['candidate']['exact_interval'],'period_cell_volume_interval':[str(cellL),str(cellU)],'equal_weight_w1':{'mass':'1','level1_receiver_ratio_ceiling_interval':[str(cellL/2),str(cellU/2)],'strictly_below_full_old_lower':cellU/2<oldL},'derived_weight':str(w),'certified_positive_headroom':str(headroom),'weight_derivation':'largest power-of-two inverse <=(cell_volume_lower/old_ratio_upper-1)/2','normalization_gain_requirement_for_any_improvement':str(((1+w)*oldU-oldL)/64)})
 assert cellU/2<oldL
 # Pick minimal interior and maximal satellite/unit count among existing127.
 eligible=[ell for ell in ells if min(ell)>=2]
 best=max(eligible,key=lambda ell:(Q(prod(z-1 for z in ell),prod(ell)),ell))
 assert (2,)*6 in ells
 defs=[]
 for label,ell in [('minimal_interior',(2,)*6),('max_relative_satellite_count_inherited',best)]:
  q=[-(l//2) if l%2 else 1-l//2 for l in ell]
  defs.append({'label':label,'reference_index':ells.index(ell),'unit_lo':q,'unit_hi':[a+l-1 for a,l in zip(q,ell)],'sat_lo':q,'sat_hi':[a+l-2 for a,l in zip(q,ell)]})
 defs.append({'label':'unit0_plus_single_bodycenter','reference_index':0,'unit_lo':[0]*6,'unit_hi':[0]*6,'sat_lo':[0]*6,'sat_hi':[0]*6})
 D=10**60;eta=Q(1,10**35)
 def radius(M):
  target=M*D**6/64;R=A.floor_root(target.numerator//target.denominator,6)
  lo,hi=Q(R,D),Q(R+1,D)
  assert lo**6<=M/64<hi**6
  return lo,hi
 def coord_bounds(index,axis):
  lo,hi=ab[axis];return (index*lo,index*hi) if index>=0 else (index*hi,index*lo)
 def reference_inner(k):
  ret=[]
  for axis in range(6):
   if axis==1:
    pair=old['next_states']['left']['selections'][str(k)];vals=[field(old['complete_functions_reused'][i]['c'])+Q(old['complete_functions_reused'][i]['m'])*spacing[1] for i in pair]
   else:vals=[field(iface['endpoints'][eid]['c0']) for eid in iface['boxes'][k][axis]]
   ret.append((interval(vals[0])[1]+eta,interval(vals[1])[0]-eta))
  assert 0<vbox(ret);return ret
 results=[];delta_upper_sum=Q(0);max_index=0
 for z in defs:
  configs+=1;budget();m=prod(b-a+1 for a,b in zip(z['unit_lo'],z['unit_hi']));k=prod(b-a+1 for a,b in zip(z['sat_lo'],z['sat_hi']));M=Q(m)+k*w;rL,rU=radius(M)
  inn=[];outer=[];source=[]
  for axis in range(6):
   normalized=[Q(z['unit_lo'][axis]),Q(z['unit_hi'][axis]),Q(z['sat_lo'][axis])+Q(1,2),Q(z['sat_hi'][axis])+Q(1,2)]
   mn,mx=min(normalized),max(normalized)
   mnL,mnU=coord_bounds(mn,axis);mxL,mxU=coord_bounds(mx,axis)
   il=max(Q(0),mxU-rL+eta);ih=min(ab[axis][0]/2,mnL+rL-eta)
   ol=max(Q(0),mxL-rU);oh=min(ab[axis][1]/2,mnU+rU)
   assert 0<=il<ih<=ab[axis][0]/2
   for t in normalized:
    tl,th=coord_bounds(t,axis)
    for x in (il,ih):
     nodes+=1;budget();assert -rL+eta<=x-th and x-tl<=rL-eta
   # Actual reflected half-integer index transformation, preserving counts.
   assert -Q(z['sat_hi'][axis])-Q(1,2)==Q(-z['sat_hi'][axis]-1)+Q(1,2)
   assert -Q(z['sat_lo'][axis])-Q(1,2)==Q(-z['sat_lo'][axis]-1)+Q(1,2)
   max_index=max(max_index,*map(abs,[z['unit_lo'][axis],z['unit_hi'][axis],z['sat_lo'][axis],z['sat_hi'][axis],-z['sat_lo'][axis]-1,-z['sat_hi'][axis]-1]))
   inn.append((il,ih));outer.append((ol,oh));source.append({'axis':axis+1,'unit_index_extrema':[z['unit_lo'][axis],z['unit_hi'][axis]],'satellite_integer_index_extrema':[z['sat_lo'][axis],z['sat_hi'][axis]],'actual_combined_normalized_coordinate_extrema':[str(mn),str(mx)]})
  ref=reference_inner(z['reference_index']);inter=intersect(outer,ref);du=vbox(outer)-vbox(inter);assert du>=0;delta_upper_sum+=du
  results.append({**z,'unit_count':m,'satellite_count':k,'selected_mass':str(M),'radius_lower':str(rL),'radius_upper':str(rU),'strict_source_margin':str(eta),'true_source_extrema':source,'rational_inward_box':encodebox(inn),'rational_outer_box':encodebox(outer),'reference_old_inner_box':encodebox(ref),'outer_minus_reference_inner_volume_upper':str(du),'inward_box_volume':str(vbox(inn))})
 family_upper=(oldU+64*delta_upper_sum)/(1+w)
 out.update({'candidate_witnesses':results,'weights_examined':2,'mixed_witness_count':len(results),'old127_preserved':True,'common_integer_index_collar_for_selected_reflections':max_index,'full_three_mixed_box_family_ratio_upper':str(family_upper),'sum_new_volume_increment_upper':str(delta_upper_sum),'upper_bound_logic':'union(oldU,new3) increment <= sum_i |outer_i minus old_reference_inner_i|; overlaps with all oldU/each other can only reduce this upper bound','three_box_subfamily_strictly_below_old_complete_lower':family_upper<oldL,'three_box_subfamily_upper_minus_old_lower':str(family_upper-oldL),'no_complete_new_union_volume_computed':True,'status':'finished'})
except Exception as e:out.update({'status':'stopped_exception','error':repr(e)})
finally:
 out.update({'new_configurations':configs,'new_source_extremum_nodes':nodes,'CPU_seconds':time.process_time()-CPU,'wall_seconds_before_exclusive_save':time.monotonic()-WALL,'input_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [OLD,IFACE,ARITH,HERE/'registration.json',Path(__file__)]}})
 with OUT.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
 print(json.dumps({k:out[k] for k in ['status','pid','new_configurations','new_source_extremum_nodes','CPU_seconds','wall_seconds_before_exclusive_save']}))
 signal.alarm(0)
 if 'error' in out:print(out['error']);sys.exit(1)
