#!/usr/bin/env python3
"""One finite full-spacing structural certificate; no union/integration routine."""
from pathlib import Path
from fractions import Fraction as Q
from functools import lru_cache,cmp_to_key
from itertools import permutations
from math import prod,comb
import json,hashlib,sys,resource,time,signal,os
sys.dont_write_bytecode=True
H=Path(__file__).resolve().parent;OUT=H/'results.json'
if OUT.exists():raise FileExistsError('Exclusive result exists')
CPU,WALL=time.process_time(),time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(15,15))
def alarm(*_):raise TimeoutError('170 second soft wall budget with cleanup reserve')
signal.signal(signal.SIGALRM,alarm);signal.alarm(170)
import field as A
A.bounds=lru_cache(None)(A.bounds)
Z=A.Z;C=A.constant;add=A.add;sub=A.sub;scale=A.scale;parse=A.parse;encode=A.encode;sign=A.sign;view=A.view
configs=0;checks=0;endpoint_checks=0;order_gap_checks=0
out={'status':'running','pid':os.getpid(),'old089_area_or_integral_state_replays':0,'new_integration_states':0,'old_DP_calls':0,'static6D_union_calls':0,'automatic_reruns':0,'MATLAB_calls':0,'new_Cn_claim':False}
def ck(cond):
 global checks
 checks+=1;assert cond
 if time.process_time()-CPU>13.5 or time.monotonic()-WALL>169:raise TimeoutError('Algorithm prelimit')
def run():
 global configs,endpoint_checks,order_gap_checks
 reg=json.loads((H/'registration.json').read_text())
 for p,sha in reg['source_pins'].items():ck(hashlib.sha256(Path(p).read_bytes()).hexdigest()==sha)
 paths=list(reg['source_pins']);iface=json.loads(Path(paths[0]).read_text());e79=json.loads(Path(paths[1]).read_text());r89=json.loads(Path(paths[3]).read_text())
 ck(e79['status']=='pass_complete_s9_kink_and_local_maximum');ck(r89['status']=='finished_strict_a1_local_maximum')
 a=scale(add(C(1),A.mono(1,0)),Q(1,2));t=add(C(1),A.mono(4,0,Q(1,2)));base=[t]+[a]*5
 ck(parse(e79['s9']['coefficients'])==a);ck(parse(r89['a1']['coefficients'])==t);ck(parse(r89['tail_spacing']['coefficients'])==a)
 ck(parse(e79['inherited_E075_complete_ideal']['coefficients'])==parse(r89['baseline']['coefficients']))
 labels=[tuple(x) for x in iface['box_lengths']];labelset=set(labels);ck(len(labels)==len(labelset)==127)
 swaps=[]
 for k,ls in enumerate(labels):
  ck(all(1<=x<=4 for x in ls))
  for j in range(1,5):
   y=list(ls);y[j],y[j+1]=y[j+1],y[j];y=tuple(y);ck(y in labelset and prod(y)==prod(ls));swaps.append([k,j,labels.index(y)])
 orbits={}
 for k,ls in enumerate(labels):orbits.setdefault((ls[0],tuple(sorted(ls[1:]))),[]).append(k)
 orbit_receipts=[]
 for (first,tail),ids in sorted(orbits.items()):
  count=len(set(permutations(tail)));ck(len(ids)==count);orbit_receipts.append({'first_length':first,'sorted_tail_lengths':tail,'raw_indices':ids,'orbit_size':count})
 delta=Q(1,10**17);mu=Q(1,25);slopes={}
 for name,obj,sgn in [('a1_left',r89['sides']['left']['exact_derivative'],1),('a1_right',r89['sides']['right']['exact_derivative'],-1),('tail_left',e79['left_slope'],1),('tail_right',e79['right_slope'],-1)]:
  x=parse(obj['coefficients']);lo,hi=A.bounds(x);oldlo,oldhi=map(Q,obj['exact_interval']);ck(oldlo<=lo<=hi<=oldhi)
  margin=lo if sgn==1 else -hi;ck(margin>mu);slopes[name]={'coefficients':encode(x),'exact_interval':[str(lo),str(hi)],'display':A.dec(lo),'certified_magnitude_lower':str(margin)}
 funcs=[];maps=[]
 for axis in range(6):
  fs=[(Z,Q(0)),(Z,Q(1,2))];idx={f:i for i,f in enumerate(fs)};mov=[]
  for ls in labels:
   radius=A.rad(ls);l=ls[axis];f=(radius,-Q(l-1,2)) if l%2 else (scale(radius,-1),Q(l,2))
   if f not in idx:idx[f]=len(fs);fs.append(f)
   mov.append(idx[f])
  ck(all(-1<=m<=2 for c,m in fs));funcs.append(fs);maps.append(mov)
 ck([len(fs) for fs in funcs]==[24,37,37,37,37,37])
 ck(all(set(fs)==set(funcs[1]) for fs in funcs[2:]))
 values={};geometry=[];center_min=None
 for direction,name in [(0,'center'),(-1,'all_minus_delta'),(1,'all_plus_delta')]:
  configs+=1;ck(configs<=3);sp=[add(x,C(direction*delta)) for x in base];per_axis=[]
  for axis in range(6):
   vals=[add(c,scale(sp[axis],m)) for c,m in funcs[axis]];values[(direction,axis)]=vals
   ck(sign(sp[axis])>0 and sign(sub(C(1),scale(sp[axis],Q(1,2))))>0)
   widths=[];minw=None;minraw=None
   for k,ls in enumerate(labels):
    m=maps[axis][k];lo=Z if ls[axis]%2 else A.maximum(Z,vals[m]);hi=A.minimum(vals[1],vals[m]) if ls[axis]%2 else vals[1]
    w=sub(hi,lo);ck(sign(w)>0);ck(sign(lo)>=0 and sign(sub(vals[1],hi))>=0);endpoint_checks+=1
    if minw is None or A.cmp(w,minw)<0:minw,minraw=w,k
    if direction==0:widths.append(encode(w))
   if direction==0:
    ck(sign(sub(minw,C(3*delta)))>0)
    if center_min is None or A.cmp(minw,center_min)<0:center_min=minw
   record={'axis':axis+1,'spacing':view(sp[axis]),'minimum_width':view(minw),'minimum_width_raw':minraw}
   if direction==0:record['all127_width_coefficients']=widths
   per_axis.append(record)
  geometry.append({'configuration':name,'direction':direction,'per_axis':per_axis})
 chamber=[]
 for axis in range(6):
  fs=funcs[axis];v=values[(0,axis)];record={'axis':axis+1,'functions':[{'c':encode(c),'m':str(m)} for c,m in fs],'single_coordinate_dependence':'c+m*a_axis; all radii fixed','sides':{}}
  for side,name in [(-1,'minus'),(1,'plus')]:
   def cmp(i,j):
    q=sign(sub(v[i],v[j]));return q if q else ((side*fs[i][1]>side*fs[j][1])-(side*fs[i][1]<side*fs[j][1]))
   order=sorted(range(len(fs)),key=cmp_to_key(cmp));end=values[(side,axis)];gaps=[]
   for i,j in zip(order,order[1:]):
    g0=sub(v[j],v[i]);g1=sub(end[j],end[i]);ck(sign(g0)>=0 and sign(g1)>0);order_gap_checks+=1
    gaps.append({'ids':[i,j],'at_center':encode(g0),'at_delta_endpoint':encode(g1)})
   record['sides'][name]={'complete_order':order,'adjacent_gap_certificate':gaps}
  chamber.append(record)
 N=prod(len(fs)-1 for fs in funcs);degree_sum=sum(comb(6,k)*3**k for k in range(2,7));B=64*N*degree_sum
 ck(degree_sum==4077);ck(B*delta<mu/2);ck(delta<=1)
 guaranteed=mu-B*delta;ck(guaranteed>Q(1,50))
 out.update({'status':'proved_full_spacing_strict_local_maximum_from_pinned_derivatives','full127_S5_closure':{'adjacent_swap_count':len(swaps),'adjacent_swap_targets':swaps,'orbits':orbit_receipts,'full_family_not_only_weighted_sections':True},'base_spacing_coefficients':[encode(x) for x in base],'inherited_ideal':r89['baseline'],'inherited_derivatives':slopes,'configurations':geometry,'endpoint_chambers':chamber,'explicit_neighborhood':{'norm':'infinity','delta':str(delta),'positive_base_minimum_width':view(center_min),'width_lipschitz_constant':3,'all127_remain_positive_for_entire_neighborhood':True,'function_counts':[len(fs) for fs in funcs],'possible_orthants':64,'six_dimensional_elementary_cell_count_upper':N,'cell_width_velocity_abs_upper':3,'base_included_cell_width_upper':1,'sum_binomial_velocity_powers_degrees2to6':degree_sum,'uniform_remainder_B':B,'remainder_statement':'abs(R(a0+h)-R(a0)-sum_i g_i^{sign(h_i)}*h_i) <= B*delta*norm1(h), for normInf(h)<=delta','uniform_first_order_decrease_mu':str(mu),'B_times_delta':str(B*delta),'certified_decrease_coefficient':str(guaranteed),'simple_decrease_coefficient':'1/50','conclusion':'0<normInf(h)<=delta implies R(a0+h)<=R(a0)-(1/50)*norm1(h)<R(a0)'},'directional_derivative':'g1^{sign(v1)} v1 + sum_{i=2}^6 gt^{sign(vi)} vi <= -(1/25)*norm1(v) for every nonzero v; zero components contribute zero','acceptance_dependency':'R079 already independently audited; R089 source certificate sealed and awaiting its independent E acceptance by parent','scope':'Only the same 127 fixed labels/common radii/unit-mass lattice family, varying six positive spacings. No new masses, shifted witnesses, radii, labels, or new Cn bound excluded/claimed.'})
 for p,sha in reg['source_pins'].items():ck(hashlib.sha256(Path(p).read_bytes()).hexdigest()==sha)
try:run()
except Exception as exc:out.update({'status':'stopped_exception','error':repr(exc)})
finally:
 out['resource']={'exact_configurations':configs,'scalar_assertions':checks,'box_coordinate_positivity_checks':endpoint_checks,'adjacent_affine_gap_checks':order_gap_checks,'new_integration_states':0,'old089_area_replays':0,'old375_state_replays':0,'old_DP_calls':0,'static6D_union_calls':0,'CPU_seconds_before_save':time.process_time()-CPU,'wall_seconds_before_save':time.monotonic()-WALL}
 out['input_sha256']={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [H/'registration.json',H/'field.py',Path(__file__)]}
 with OUT.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
 receipt={'status':out['status'],'pid':out['pid'],'exit_code':int('error' in out),'resource':out['resource'],'CPU_seconds_including_save':time.process_time()-CPU,'wall_seconds_including_save':time.monotonic()-WALL,'results_sha256':hashlib.sha256(OUT.read_bytes()).hexdigest()}
 with (H/'receipt.json').open('x') as f:json.dump(receipt,f,indent=2);f.write('\n')
 print(json.dumps(receipt));signal.alarm(0)
 if 'error' in out:print(out['error']);sys.exit(1)
