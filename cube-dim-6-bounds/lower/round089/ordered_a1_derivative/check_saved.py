#!/usr/bin/env python3
"""Verify saved algebraic certificate equations only; no union/integration search."""
from pathlib import Path
from fractions import Fraction as Q
import json,hashlib,time,resource,signal,sys,ast
sys.dont_write_bytecode=True
H=Path(__file__).resolve().parent
CPU,WALL=time.process_time(),time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(20,20));signal.alarm(170)
import field as A
checks=0
out={'status':'running','new_integration_states':0,'new_parameter_configurations':0,'replayed_union_algorithms':0}
def check(v):
 global checks
 checks+=1;assert v
 assert time.process_time()-CPU<18.5 and time.monotonic()-WALL<169

def coarse_interval(interval,places=18):
 lo,hi=map(Q,interval);den=10**places
 l=lo.numerator*den//lo.denominator;h=-((-hi.numerator*den)//hi.denominator)
 check(Q(l,den)<=lo<=hi<=Q(h,den))
 return [str(Q(l,den)),str(Q(h,den))]
try:
 r=json.loads((H/'results.json').read_text());receipt=json.loads((H/'receipt.json').read_text());reg=json.loads((H/'registration.json').read_text())
 check(receipt['result_sha256']==hashlib.sha256((H/'results.json').read_bytes()).hexdigest());check(r['status']=='finished_strict_a1_local_maximum')
 for path,sha in {**reg['source_pins'],**r['input_sha256']}.items():check(hashlib.sha256(Path(path).read_bytes()).hexdigest()==sha)
 for m,l,h in [(2,A.UL,A.UH),(3,A.VL,A.VH)]:check(l**6<=m<=h**6 and h-l==Q(1,A.DEN))
 # Confirm pure reused arithmetic definitions exactly match the pinned E058 AST.
 source=next(p for p in reg['source_pins'] if p.endswith('round058/a2_next_event_review/review.py'))
 orig={n.name:ast.dump(n,include_attributes=False) for n in ast.parse(Path(source).read_text()).body if isinstance(n,ast.FunctionDef)}
 for n in ast.parse((H/'field.py').read_text()).body:
  if isinstance(n,ast.FunctionDef):
   n.decorator_list=[];old=next(x for x in ast.parse(Path(source).read_text()).body if isinstance(x,ast.FunctionDef) and x.name==n.name);old.decorator_list=[];check(ast.dump(n,include_attributes=False)==ast.dump(old,include_attributes=False))
 cuts=[A.parse(x) for x in r['tail_cuts']];cap=A.parse(r['tail_cap']['coefficients']);full=cap
 for _ in range(4):full=A.mul(full,cap)
 families={};cells=0
 for f in r['section_families']:
  L,U=f['complement_lower_ranks'],f['complement_upper_ranks'];knots=f['knots'];check(knots==sorted(set([0,20]+L+U)))
  prev=[[A.constant(1)] for _ in range(len(knots)-1)];count=0
  for i,row in enumerate(f['recurrence_rows']):
   check(len(row)==len(knots)-1);carry=A.Z;nxt=[]
   for j,c in enumerate(row):
    l,h=knots[j:j+2];check(c['interval_ranks']==[l,h]);check(c['active']==(L[i]<=l and h<=U[i]))
    p=[A.parse(x) for x in c['local_polynomial']];want=[carry]+[A.scale(x,Q(1,k+1)) for k,x in enumerate(prev[j])] if c['active'] else [carry];check(p==want)
    v=A.Z;w=A.sub(cuts[h],cuts[l])
    for x in reversed(p):v=A.add(A.mul(v,w),x)
    check(v==A.parse(c['right_value']));carry=v;nxt.append(p);count+=1
   prev=nxt
  check(len(f['recurrence_rows'])==5);check(count==f['integration_states']);cells+=count
  comp=A.parse(f['ordered_complement_volume']['coefficients']);check(carry==comp)
  area=A.parse(f['full_section_union_area']['coefficients']);check(A.sub(full,A.scale(comp,120))==area);check(A.sign(area)>0 and A.sign(A.sub(full,area))>0);families[f['family_id']]=area
  for obj in ('ordered_complement_volume','full_section_union_area'):check(list(map(Q,f[obj]['exact_interval']))==list(A.bounds(A.parse(f[obj]['coefficients']))))
 check(cells==375==r['resource']['admitted_new_states']==r['planned_integration_states']);check(len(families)==11)
 summary={}
 for label in ('left','right'):
  s=r['sides'][label];total=A.Z
  for slab in s['slabs']:
   retained=[tuple(tuple(p) for p in B) for B in slab['retained_ordered_boxes']]
   for item in slab['containment_deletions']:
    B,D=item['box'],tuple(tuple(p) for p in item['covering_retained_box']);check(D in retained and all(c<=x and y<=z for (x,y),(c,z) in zip(B,D)))
   f=r['section_families'][slab['family_id']];check(slab['complement_lower_ranks']==f['complement_lower_ranks'] and slab['complement_upper_ranks']==f['complement_upper_ranks'])
   term=A.scale(families[slab['family_id']],64*Q(slab['width_slope']));check(term==A.parse(slab['weighted_derivative_term']['coefficients']));total=A.add(total,term)
  check(total==A.parse(s['exact_derivative']['coefficients']));check(A.sign(total)==(1 if label=='left' else -1));check(list(A.bounds(total))==list(map(Q,s['exact_derivative']['exact_interval'])))
  for gap in s['complete_order_gap_certificate']:check(A.sign(A.parse(gap['at_current']))>=0 and A.sign(A.parse(gap['at_adjacent_endpoint']))>=0)
  dt=A.sub(A.parse(s['adjacent_endpoint']['coefficients']),A.parse(r['a1']['coefficients']));loss=A.scale(A.mul(total,dt),-1);check(A.sign(loss)>0)
  summary[label]={'derivative_display':s['exact_derivative']['display_lower'],'derivative_outward_rational_interval':coarse_interval(s['exact_derivative']['exact_interval']),'adjacent_endpoint_display':s['adjacent_endpoint']['display_lower'],'adjacent_endpoint_coefficients':s['adjacent_endpoint']['coefficients'],'adjacent_endpoint_outward_rational_interval':coarse_interval(s['adjacent_endpoint']['exact_interval']),'ratio_loss_to_adjacent_endpoint':A.view(loss)}
 out.update({'status':'pass_saved_certificate','checked_equations':checks,'verified_saved_integral_cells':cells,'summary':summary,'a1_display':r['a1']['display_lower'],'baseline_display':r['baseline']['display_lower']})
except Exception as exc:out.update({'status':'failed','error':repr(exc),'checked_equations':checks})
finally:
 out.update({'CPU_seconds_before_save':time.process_time()-CPU,'wall_seconds_before_save':time.monotonic()-WALL,'input_sha256':{str(H/n):hashlib.sha256((H/n).read_bytes()).hexdigest() for n in ['results.json','receipt.json','registration.json','field.py','check_saved.py']}})
 with (H/'checker_receipt.json').open('x') as f:json.dump(out,f,indent=2);f.write('\n')
 print(json.dumps({k:out.get(k) for k in ['status','checked_equations','verified_saved_integral_cells','CPU_seconds_before_save','wall_seconds_before_save','error']}));signal.alarm(0)
 if out['status']!='pass_saved_certificate':sys.exit(1)
