#!/usr/bin/env python3
"""Read-only saved scalar/closure/gap certificate verification, no geometric configurations."""
from pathlib import Path
from fractions import Fraction as Q
from math import prod,comb
import json,time,resource,signal,hashlib,sys
sys.dont_write_bytecode=True
H=Path(__file__).resolve().parent;CPU,WALL=time.process_time(),time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(15,15));signal.alarm(170)
import field as A
n=0
out={'status':'running','new_exact_configurations':0,'new_integration_states':0,'old089_area_replays':0,'old375_state_replays':0,'old_DP_calls':0,'static6D_union_calls':0}
def ck(b):
 global n
 n+=1;assert b;assert time.process_time()-CPU<13.5
try:
 reg=json.loads((H/'registration.json').read_text());r=json.loads((H/'results.json').read_text());rec=json.loads((H/'receipt.json').read_text())
 ck(rec['results_sha256']==hashlib.sha256((H/'results.json').read_bytes()).hexdigest())
 for p,sha in {**reg['source_pins'],**r['input_sha256']}.items():ck(hashlib.sha256(Path(p).read_bytes()).hexdigest()==sha)
 labels=json.loads(Path(next(iter(reg['source_pins']))).read_text())['box_lengths']
 swaps=r['full127_S5_closure']['adjacent_swap_targets'];ck(len(swaps)==127*4)
 for k,j,target in swaps:
  x=labels[k].copy();x[j],x[j+1]=x[j+1],x[j];ck(x==labels[target])
 cover=[k for o in r['full127_S5_closure']['orbits'] for k in o['raw_indices']];ck(sorted(cover)==list(range(127)))
 delta=Q(r['explicit_neighborhood']['delta']);base=[A.parse(x) for x in r['base_spacing_coefficients']]
 for axis,c in enumerate(r['endpoint_chambers']):
  fs=[(A.parse(x['c']),Q(x['m'])) for x in c['functions']]
  for side,name in [(-1,'minus'),(1,'plus')]:
   s=c['sides'][name];ck(sorted(s['complete_order'])==list(range(len(fs))))
   ck([x['ids'] for x in s['adjacent_gap_certificate']]==[list(x) for x in zip(s['complete_order'],s['complete_order'][1:])])
   for gap in s['adjacent_gap_certificate']:
    i,j=gap['ids'];g=A.add(A.sub(fs[j][0],fs[i][0]),A.scale(base[axis],fs[j][1]-fs[i][1]));end=A.add(g,A.constant(side*delta*(fs[j][1]-fs[i][1])))
    ck(g==A.parse(gap['at_center']) and end==A.parse(gap['at_delta_endpoint']));ck(A.sign(g)>=0 and A.sign(end)>0)
 for s in r['inherited_derivatives'].values():
  a=A.parse(s['coefficients']);lo,hi=A.bounds(a);ck([lo,hi]==list(map(Q,s['exact_interval'])));ck(lo>Q(1,25) or hi<-Q(1,25))
 for ax in r['configurations'][0]['per_axis']:
  mn=A.parse(ax['minimum_width']['coefficients']);ws=[A.parse(x) for x in ax['all127_width_coefficients']];ck(mn in ws);ck(all(A.sign(A.sub(w,mn))>=0 for w in ws));ck(A.sign(A.sub(mn,A.constant(3*delta)))>0)
 e=r['explicit_neighborhood'];N=prod(x-1 for x in e['function_counts']);B=64*N*sum(comb(6,k)*3**k for k in range(2,7));ck(N==e['six_dimensional_elementary_cell_count_upper'] and B==e['uniform_remainder_B'])
 ck(B*delta<Q(1,50));ck(Q(1,25)-B*delta==Q(e['certified_decrease_coefficient'])>Q(1,50));ck(r['resource']['exact_configurations']==3)
 out.update({'status':'pass_saved_scalar_and_chamber_certificate','checked_assertions':n,'delta':str(delta),'B':B,'simple_decrease_coefficient':'1/50'})
except Exception as exc:out.update({'status':'failed','error':repr(exc),'checked_assertions':n})
finally:
 out.update({'CPU_seconds_before_save':time.process_time()-CPU,'wall_seconds_before_save':time.monotonic()-WALL,'input_sha256':{str(H/p):hashlib.sha256((H/p).read_bytes()).hexdigest() for p in ['results.json','registration.json','receipt.json','field.py','check_saved.py']}})
 with (H/'checker_receipt.json').open('x') as f:json.dump(out,f,indent=2);f.write('\n')
 print(json.dumps({k:out.get(k) for k in ['status','checked_assertions','CPU_seconds_before_save','wall_seconds_before_save','error']}));signal.alarm(0)
 if out['status']!='pass_saved_scalar_and_chamber_certificate':sys.exit(1)
