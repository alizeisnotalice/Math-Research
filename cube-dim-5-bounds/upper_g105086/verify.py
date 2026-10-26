#!/usr/bin/env python3
"""Author zero-H verification of the new target only, by separable endpoints."""
from fractions import Fraction as F
from pathlib import Path
from collections import Counter,defaultdict
import argparse,json,os,time,resource
from flint import arb,fmpq,ctx
from scan import SOURCE,SOURCE_SHA,TARGET,LOCKS,sha,digest,write

def ab(x):
 x=F(x);return arb(fmpq(x.numerator,x.denominator))
def box(node,roots):
 lo=[F(x,32) for x in roots[node[0]]];hi=[x+F(1,32) for x in lo]
 for bit in node[1]:
  widths=[b-a for a,b in zip(lo,hi)];axis=max(range(4),key=lambda j:widths[j]);mid=(lo[axis]+hi[axis])/2
  if bit=='0':hi[axis]=mid
  else:assert bit=='1';lo[axis]=mid
 return lo,hi
def run():
 began=time.monotonic();ap=argparse.ArgumentParser();ap.add_argument('--work-dir',type=Path,required=True);args=ap.parse_args();out=args.work_dir.resolve()
 assert sha(SOURCE)==SOURCE_SHA
 for p,h in LOCKS.items():assert sha(p)==h
 write(out/'verification_registration.json',{'pid':os.getpid(),'source':str(SOURCE),'source_sha256':SOURCE_SHA,
  'code':str(Path(__file__).resolve()),'code_sha256':sha(__file__),'target':'21/200000','finite_cases':97960,'block_size':2048,
  'method':'exact physical-gradient ratio containment;768bit log lower comparison;independent separable endpoint minimum;closed-root/source mapper/complete residual identity',
  'new_H_physical_gradient_kernel_calls':0,'old_g100_replay':False})
 state=json.load(open(SOURCE));expectedroots=[[i,j,k,l] for i in range(32) for j in range(i+1) for k in range(j+1) for l in range(k+1)]
 assert state['roots']==expectedroots and not state['queue'] and not state['inflight']
 splitkeys={tuple(s['node']) for s in state['splits']};assert len(splitkeys)==45600
 for s in state['splits']:
  r,p=s['node'];assert s['axis']==len(p)%4 and s['children']==[[r,p+'0'],[r,p+'1']]
 rootpaths=defaultdict(list);mapper_count=0
 for row in map(json.loads,open(out/'active_leaf_mapper.jsonl')):
  i=row['source_leaf_index'];assert i==mapper_count;mapper_count+=1;leaf=state['passed'][i];c=leaf['certificate'];record=state['cache'][c['anchor']]
  lo,hi=box(leaf['node'],expectedroots)
  assert row['node']==leaf['node'] and row['source_certificate_sha256']==digest(c) and row['source_leaf_type']==c.get('mechanism',c.get('quadratic'))
  assert row['H_interval']==record['values'][0] and row['anchor']==c['anchor'] and row['point']==record['point']
  assert row['box_lower']==list(map(str,lo)) and row['box_upper']==list(map(str,hi))
  r,p=leaf['node'];rootpaths[r].append(p)
  for n in range(len(p)):assert (r,p[:n]) in splitkeys
 assert mapper_count==97960
 for r in range(52360):
  ps=sorted(rootpaths[r]);assert ps and all(not b.startswith(a) for a,b in zip(ps,ps[1:]));d=max(map(len,ps));assert sum(1<<(d-len(p)) for p in ps)==1<<d
 order=json.load(open(out/'risk_order.json'))['source_leaf_indices'];assert set(order)==set(range(97960))
 counts=Counter();types=Counter();depths=Counter();roots=Counter();seen=set();prev=None;rank=0;minimum_passed=None;minimum_failed=None;maximum_failed=None;failedHL=[];failedvolume=F(0)
 failure_iter=iter(open(out/'uncertified_boxes.jsonl'));status_iter=iter(open(out/'all_leaf_status.jsonl'));failedrows=[];ctx.prec=768
 for n,b in enumerate(sorted(out.glob('batch[0-9][0-9][0-9][0-9]')),1):
  reg=json.load(open(b/'registration.json'));cp=json.load(open(b/'checkpoint.json'));ss=json.load(open(b/'summary.json'))
  assert cp['parent_checkpoint_sha256']==prev and cp['registration_sha256']==sha(b/'registration.json') and cp['certificates_sha256']==sha(b/'certificates.jsonl') and cp['summary_sha256']==sha(b/'summary.json')
  selected=order[rank:rank+2048];assert reg['selected_source_leaf_indices']==selected and reg['start_rank']==rank
  local=Counter();index=0
  for row in map(json.loads,open(b/'certificates.jsonl')):
   i=row['source_leaf_index'];assert i==selected[index] and i not in seen;seen.add(i);index+=1;leaf=state['passed'][i];old=leaf['certificate'];record=state['cache'][old['anchor']];cert=row['certificate']
   assert row['node']==leaf['node'] and row['source_leaf_type']==old.get('mechanism',old.get('quadratic')) and row['source_certificate_sha256']==digest(old)
   assert cert['target']=='21/200000' and cert['epsilon']=='1/64' and cert['anchor']==old['anchor'] and cert['quadratic']=='Q_full=I4/16' and cert['checked_points']==16
   vals=[tuple(map(F,p)) for p in record['values']];hl,hu=vals[0];assert TARGET<hl<=hu and row['H_interval']==record['values'][0]
   gl=[tuple(map(F,p)) for p in cert['loggradient_directed']]
   for (a,z),(l,u) in zip(vals[1:],gl):assert l<=min(a/hl,a/hu,z/hl,z/hu)<=max(a/hl,a/hu,z/hl,z/hu)<=u
   delta=F(cert['Delta_lower']);assert delta<=F(str((ab(hl)/ab(TARGET)).log().lower().fmpq()))
   lo,hi=box(leaf['node'],expectedroots);assert row['box_lower']==list(map(str,lo)) and row['box_upper']==list(map(str,hi));point=list(map(F,record['point']))
   minimum=delta;bits=[]
   for j,(a,z) in enumerate(gl):
    vv=[lo[j]-point[j+1],hi[j]-point[j+1]];terms=[min(a*v,z*v)-256*v*v for v in vv]
    best=min(range(2),key=lambda k:(terms[k],k));minimum+=terms[best];bits.append(best)
   assert minimum==F(cert['minimum_log_margin']) and bits==cert['worst_corner_bits'] and (minimum>0)==cert['passed']
   status=json.loads(next(status_iter));assert status=={'source_leaf_index':i,'node':row['node'],'source_leaf_type':row['source_leaf_type'],'passed':cert['passed'],'minimum_log_margin':cert['minimum_log_margin'],'certificate_block':n}
   if cert['passed']:
    counts['passed']+=1;local['passed']+=1;minimum_passed=minimum if minimum_passed is None else min(minimum_passed,minimum)
   else:
    counts['uncertified']+=1;local['uncertified']+=1;assert json.loads(next(failure_iter))==row;failedrows.append(row)
    types[row['source_leaf_type']]+=1;depths[len(row['node'][1])]+=1;roots[row['node'][0]]+=1;failedHL.append(hl);failedvolume+=F(1,32**4*2**len(row['node'][1]))
    minimum_failed=minimum if minimum_failed is None else min(minimum_failed,minimum);maximum_failed=minimum if maximum_failed is None else max(maximum_failed,minimum)
  assert index==len(selected)==ss['tested'] and local['passed']==ss['passed'] and local['uncertified']==ss['uncertified']
  rank+=index;assert cp['next_rank']==rank;prev=sha(b/'checkpoint.json')
  write(out/'verification_checkpoint.json',{'pid':os.getpid(),'verified_blocks':n,'verified_leaves':rank,'passed':counts['passed'],'uncertified':counts['uncertified'],'last_source_checkpoint_sha256':prev})
 assert rank==97960 and len(seen)==97960 and next(failure_iter,None) is None and next(status_iter,None) is None
 summary=json.load(open(out/'summary.json'));assert counts['passed']==summary['passed']==96104 and counts['uncertified']==summary['uncertified']==1856
 residual={'uncertified':counts['uncertified'],'distinct_roots':len(roots),'source_leaf_types':dict(types),'depths':dict(depths),
  'minimum_failed_margin':str(minimum_failed),'maximum_failed_margin':str(maximum_failed),'minimum_uncertified_H_lower':str(min(failedHL)),
  'maximum_uncertified_H_lower':str(max(failedHL)),'exact_physical_box_union_volume':str(failedvolume),
  'all_residual_anchor_H_lower_above_target':True,'not_H_counterexamples':True,'source_sha256':SOURCE_SHA,'target':'21/200000'}
 write(out/'uncertified_summary.json',residual)
 result={'status':'author_new_target_zero_H_verification_passed','pid':os.getpid(),'target':'21/200000','verified_leaves':rank,'verified_768bit_log_lower_comparisons':rank,
  'passed':counts['passed'],'uncertified':counts['uncertified'],'new_g105_closed_root_cover':False,'source_mapper_complete':True,'ordered_root_list_exact':True,
  'closed_binary_root_Kraft':52360,'new_certificate_types_independent_of_old_types':True,'minimum_passed_margin':str(minimum_passed),
  'last_checkpoint_sha256':prev,'uncertified_boxes_sha256':sha(out/'uncertified_boxes.jsonl'),'residual':residual,
  'wall_seconds':time.monotonic()-began,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'new_H_physical_gradient_kernel_calls':0,'old_g100_replay':False,'global_g_or_C5_claim':False}
 assert sha(SOURCE)==SOURCE_SHA
 for p,h in LOCKS.items():assert sha(p)==h
 write(out/'author_validation.json',result);print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':run()
