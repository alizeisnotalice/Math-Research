#!/usr/bin/env python3
"""Independent exact local audit. No imports of author code; no union routine."""
from fractions import Fraction as Q
from functools import lru_cache, cmp_to_key
from math import prod
from pathlib import Path
from decimal import Decimal, localcontext
import hashlib, json, os, resource, time
START_CPU=time.process_time(); START_WALL=time.monotonic()
resource.setrlimit(resource.RLIMIT_CPU,(30,31))
HERE=Path(__file__).resolve().parent
MAIN=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds')
AUTHOR=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-lower-clean/cube-dim-6-bounds/lower/round058/a2_next_event_preflight')
DROOT=MAIN/'lower/round008/spacing_chamber'; EROOT=MAIN/'audit/lower008/parametric_union'
DEN=10**100; Z=(Q(0),)*36
nodes={}; configs=0; scalar_source_checks=0

def node(kind):
 nodes[kind]=nodes.get(kind,0)+1
 assert sum(nodes.values())<=2000
 assert time.process_time()-START_CPU<29

def constant(x):return (Q(x),)+Z[1:]
def mono(i,j,c=Q(1)):
 a=list(Z);a[6*(i%6)+j%6]=Q(c)*2**(i//6)*3**(j//6);return tuple(a)
def parse(d):
 a=list(Z)
 for ij,c in d.items():
  i,j=map(int,ij.split(','));assert 0<=i<6 and 0<=j<6 and type(c) in (str,int)
  a[6*i+j]=Q(c)
 return tuple(a)
def encode(a):return {f'{i//6},{i%6}':str(c) for i,c in enumerate(a) if c}
def add(a,b):return tuple(x+y for x,y in zip(a,b))
def scale(a,c):return tuple(x*c for x in a)
def sub(a,b):return add(a,scale(b,-1))
def mul(a,b):
 r=list(Z)
 for i,x in enumerate(a):
  if x:
   for j,y in enumerate(b):
    if y:
     p,q=i//6+j//6,i%6+j%6
     r[6*(p%6)+q%6]+=x*y*2**(p//6)*3**(q//6)
 return tuple(r)
def rootbox(m):
 lo,hi=DEN,2*DEN;target=m*DEN**6
 assert lo**6<=target<hi**6
 while hi-lo>1:
  mid=(lo+hi)//2
  if mid**6<=target:lo=mid
  else:hi=mid
 assert lo**6<=target<(lo+1)**6
 return Q(lo,DEN),Q(lo+1,DEN)
UL,UH=rootbox(2);VL,VH=rootbox(3)
POWER_BOUNDS=[(UL**(k//6)*VL**(k%6),UH**(k//6)*VH**(k%6)) for k in range(36)]
@lru_cache(None)
def bounds(a):
 lo=hi=Q(0)
 for c,(l,h) in zip(a,POWER_BOUNDS):
  lo+=c*(l if c>=0 else h);hi+=c*(h if c>=0 else l)
 return lo,hi
def sign(a):
 if a==Z:return 0
 lo,hi=bounds(a)
 assert lo>0 or hi<0,'Unresolved exact sign'
 return 1 if lo>0 else -1
def cmp(a,b):return sign(sub(a,b))
key=cmp_to_key(cmp)
def minimum(a,b):return a if cmp(a,b)<=0 else b
def maximum(a,b):return a if cmp(a,b)>=0 else b
def dec(q):
 with localcontext() as c:
  c.prec=65;return str(Decimal(q.numerator)/Decimal(q.denominator))
def view(a):
 lo,hi=bounds(a)
 return {'coefficients':encode(a),'exact_interval':[str(lo),str(hi)],'display_lower':dec(lo)}
def rad(ell):
 n=prod(ell);p=q=0
 while n%2==0:p+=1;n//=2
 while n%3==0:q+=1;n//=3
 assert n==1
 return mono(p,q,Q(1,2))

def run():
 global configs,scalar_source_checks
 paths=[DROOT/'chamber_interface.json',DROOT/'chamber_results.json',DROOT/'chamber.py',
 EROOT/'candidate_results.json',EROOT/'parametric_union.py',EROOT/'interface_proof.md',
 MAIN/'results/bounds_ledger.json',MAIN/'audit/lower008/primary/inward_smallstep.py',
 MAIN/'audit/lower008/primary/inward_smallstep_results.json',
 MAIN/'audit/round008/lower_low_resource_review/report.md',AUTHOR/'results.json',
 AUTHOR/'preflight.py',AUTHOR/'manifest.json',HERE/'registration.json',Path(__file__)]
 hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
 d=json.loads(paths[0].read_text());old=json.loads(paths[1].read_text());e=json.loads(paths[3].read_text())
 ledger=json.loads((MAIN/'results/bounds_ledger.json').read_text());author=json.loads((AUTHOR/'results.json').read_text())
 # Binding to the complete previously accepted union certificate, no recomputation.
 assert hashes[str(paths[0])]=='99e762ad8fb63303f90b56c5d7f90fb855a3335903cae21650a23bfb07fea933'
 assert hashes[str(paths[3])]=='6849d3774ab3db1524969715970bbb3aecd07de4467faee5d77b69775cf5ed3d'
 assert ledger['lower']['id']=='R6-008-spacing2-first-event'
 assert ledger['lower']['certificate_sha256']==hashes[str(paths[3])]
 assert ledger['lower']['exact_interval']==e['candidate_ratio_interval']
 assert e['status']=='certified_interval_union_polynomial' and e['degree']==1 and e['static_check_exact_match']
 assert Q(d['normalization_factor'])==Q(e['normalization_factor'])==64
 s0=add(mono(2,2,Q(1,3)),mono(2,4,Q(1,6)))
 s1=add(mono(3,3,Q(1,6)),mono(4,1,Q(1,3)))
 sstart=parse(d['spacing2']['c0']);delta=parse(d['spacing2']['c1'])
 assert add(sstart,delta)==s0 and sstart==parse(old['current_spacing2_coefficients'])
 spac=[add(constant(1),mono(4,0,Q(1,2))),s0]+[add(constant(Q(1,2)),mono(1,0,Q(1,2)))]*4
 raw=[];funcs=[(Z,Q(0)),(Z,Q(1,2))];fmap={f:i for i,f in enumerate(funcs)}
 ep={n:(parse(v['c0']),parse(v['c1'])) for n,v in d['endpoints'].items()}
 assert len(d['boxes'])==len(d['box_lengths'])==127
 assert len({tuple(l) for l in d['box_lengths']})==127
 for k,(ell,b) in enumerate(zip(d['box_lengths'],d['boxes'])):
  assert all(1<=l<=4 for l in ell)
  r=rad(ell);l=ell[1]
  f=(r,-Q(l-1,2)) if l%2 else (scale(r,-1),Q(l,2))
  if f not in fmap:fmap[f]=len(funcs);funcs.append(f)
  raw.append({'ell':ell,'r':r,'moving':fmap[f],'boxids':b})
 assert len(funcs)==37
 geom={};minwidth=None;max_source_index=0
 for label,s,t in [('old_start',sstart,0),('current',s0,1),('next',s1,None)]:
  configs+=1;sp=list(spac);sp[1]=s;boxes=[]
  for record in raw:
   node('geometric_box_at_parameter')
   ell,r=record['ell'],record['r'];box=[]
   for j,(l,a) in enumerate(zip(ell,sp)):
    q=-(l//2) if l%2 else 1-l//2
    max_source_index=max(max_source_index,abs(q),abs(q+l-1))
    lo=maximum(Z,sub(scale(a,q+l-1),r));hi=minimum(scale(a,Q(1,2)),add(scale(a,q),r))
    assert sign(sub(hi,lo))>0
    assert sign(lo)>=0 and sign(sub(scale(a,Q(1,2)),hi))>=0
    for x in (lo,hi):
     for n in (q,q+l-1):
      v=sub(x,scale(a,n));assert sign(add(r,v))>=0 and sign(sub(r,v))>=0
      scalar_source_checks+=1
    if t is not None:
     want=[add(ep[n][0],scale(ep[n][1],t)) for n in record['boxids'][j]]
     assert want==[lo,hi]
    if label=='next':
     w=sub(hi,lo);assert sign(sub(constant(1),w))>0
     minwidth=w if minwidth is None else minimum(minwidth,w)
    box.extend((lo,hi))
   boxes.append(tuple(box))
  geom[label]=boxes
 # Exhaustive endpoint-pair inventory for this frozen 127-box family.
 events={}
 for i in range(len(funcs)):
  for j in range(i):
   node('unclipped_function_pair')
   a,ma=funcs[j];b,mb=funcs[i]
   if ma==mb:continue
   x=scale(sub(a,b),1/(mb-ma))
   if cmp(x,sstart)>0:events.setdefault(x,[]).append((j,i))
 ordered=sorted(events,key=key)
 assert len(ordered)==217 and ordered[0]==s0 and ordered[1]==s1
 assert sum(cmp(x,s0)>0 for x in ordered)==216
 def state(event,side):
  global configs
  configs+=1
  values=[add(c,scale(event,m)) for c,m in funcs]
  def compare(i,j):
   s=cmp(values[i],values[j])
   return s if s else ((side*funcs[i][1]>side*funcs[j][1])-(side*funcs[i][1]<side*funcs[j][1]))
  fk=cmp_to_key(compare);selected=[];active=set()
  for record in raw:
   node('germ_box_selection');m=record['moving'];odd=record['ell'][1]%2
   lo=0 if odd else max((0,m),key=fk)
   hi=min((1,m),key=fk) if odd else 1
   assert compare(lo,hi)<0
   selected.append((lo,hi));active.update((lo,hi))
  return selected,sorted(active,key=fk)
 states={f'{name}_{side}':state(s,side) for name,s in [('current',s0),('next',s1)] for side in [-1,1]}
 assert states['current_1']==states['next_-1']
 assert states['current_-1'][0]==states['current_1'][0]==states['next_-1'][0]==states['next_1'][0]
 sel,order=states['current_1'];assert len(order)==22
 gaps=[]
 for i,j in zip(order,order[1:]):
  node('adjacent_effective_order_gap');a,ma=funcs[i];b,mb=funcs[j]
  g0=sub(add(b,scale(s0,mb)),add(a,scale(s0,ma)))
  g1=sub(add(b,scale(s1,mb)),add(a,scale(s1,ma)))
  assert sign(g0)>=0 and sign(g1)>=0 and (sign(g0)>0 or sign(g1)>0)
  gaps.append({'ids':[i,j],'at_s0':encode(g0),'at_s1':encode(g1)})
 crossing=[]
 for name,s in [('current',s0),('next',s1)]:
  assert len(events[s])==1
  for i,j in events[s]:
   assert i in states[name+'_-1'][1] and j in states[name+'_-1'][1]
   assert i in states[name+'_1'][1] and j in states[name+'_1'][1]
   A=[k for k,p in enumerate(sel) if i in p];B=[k for k,p in enumerate(sel) if j in p]
   assert set(A).isdisjoint(B)
   certs=[]
   for a in A:
    for b in B:
     node('crossing_box_pair')
     xa,xb=geom[name][a],geom[name][b]
     gap=sub(maximum(xa[0],xb[0]),minimum(xa[1],xb[1]))
     assert sign(gap)>0
     certs.append({'boxes':[a,b],'lengths':[raw[a]['ell'],raw[b]['ell']],'axis':1,'strict_separation':view(gap)})
   crossing.append({'event':name,'function_ids':[i,j],'functions':[{'constant':encode(funcs[k][0]),'slope':str(funcs[k][1])} for k in (i,j)],'group_sizes':[len(A),len(B)],'pair_certificates':certs})
 assert [c['group_sizes'] for c in crossing]==[[4,4],[6,1]]
 # Absolute-spacing slope is bound to the complete E008 polynomial.
 p=[parse(x) for x in e['ratio_polynomial_coefficients']];assert len(p)==2
 g=parse(old['spacing_derivative_coefficients'])
 assert mul(g,delta)==p[1] and sign(g)>0
 baseline=add(p[0],p[1]);assert baseline==parse(e['candidate_ratio_coefficients'])
 result=add(baseline,mul(g,sub(s1,s0)))
 assert result==parse(author['next_affine_ideal_candidate']['coefficients'])
 lo,hi=bounds(result);old_upper=Q(e['candidate_ratio_interval'][1])
 assert lo>old_upper
 alo,ahi=map(Q,author['next_affine_ideal_candidate']['exact_interval'])
 assert alo<=lo<=hi<=ahi
 # Direct finite L1 interface. Uniformly eroded algebraic boxes suffice.
 configs+=1
 eps=Q(1,10**22);eta=eps/2;theta=1-Q(1,10**16);K=10**16
 target=Q(2067759428416,10**12)
 assert sign(sub(minwidth,constant(2*eps)))>0 and max_source_index==2
 assert all(sign(sub(constant(1),scale(a,Q(1,2))))>0 for a in [spac[0],s1,*spac[2:]])
 # For 127 boxes, erasing a 2eps layer per coordinate loses <=127*12eps.
 loss=64*127*12*eps
 rin=lo-loss
 factor=Q(2*K+1,2*K+5)**6
 finite=theta*rin*factor
 assert 0<theta<1 and finite>target>old_upper
 assert configs<=12
 return {'status':'pass_independent_local_and_finite_L1_certificate',
 'input_sha256':hashes,'root_denominator':DEN,'root_enclosures':{'u':[str(UL),str(UH)],'v':[str(VL),str(VH)]},
 'scope':'same frozen 127 central block witnesses; exact inherited E008 union, no new union',
 's0':view(s0),'s1':view(s1),'old_start':view(sstart),'functions':37,
 'event_inventory':{'all_function_pairs':666,'strictly_after_old_start':217,'first_is_s0':True,'next_is_s1':True,'strictly_after_s0':216},
 'raw_boxes':127,'positive_boxes_all_three_parameters':True,'four_germ_selections_identical':True,
 'effective_order':order,'closed_interval_gaps':gaps,'crossings':crossing,
 'inherited_slope':view(g),'E008_t_slope_equals_absolute_slope_times_delta':True,
 'baseline_sum_equals_complete_E008_candidate':True,'new_ideal_ratio':view(result),
 'strict_margin_over_complete_old_upper':str(lo-old_upper),
 'finite_L1':{'epsilon':str(eps),'bump_halfwidth':str(eta),'threshold':str(theta),'receiver_index_radius':K,'source_index_radius':K+2,'common_radius':'(product ell)^(1/6)/2','max_relative_source_index':max_source_index,'minimum_uneroded_width':view(minwidth),'all_half_cell_widths_below_one':True,'union_loss_upper':str(loss),'inward_union_ratio_lower':str(rin),'boundary_factor':str(factor),'weak_ratio_lower':str(finite),'weak_ratio_display_lower':dec(finite),'strict_rational_target':str(target),'margin_over_target':str(finite-target),'margin_over_target_display':dec(finite-target),'margin_over_complete_old_upper':str(finite-old_upper)},
 'resource':{'configurations':configs,'local_nodes':sum(nodes.values()),'nodes_by_type':nodes,'scalar_source_endpoint_checks':scalar_source_checks,'full_union_calls':0,'old_baseline_recomputations':0,'MATLAB_calls':0},
 'C6_registration_performed':False,'proof_status':'complete conditional only on the accepted frozen E008 union certificate; new finite-L1 implication established; parent owns ledger acceptance'}

if __name__=='__main__':
 out=run();out.update({'pid':os.getpid(),'cpu_seconds':time.process_time()-START_CPU,'wall_seconds':time.monotonic()-START_WALL})
 (HERE/'review_results.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({'status':out['status'],'pid':out['pid'],'cpu_seconds':out['cpu_seconds'],'wall_seconds':out['wall_seconds'],'resource':out['resource'],'new_ideal_display':out['new_ideal_ratio']['display_lower'],'finite_lower_display':out['finite_L1']['weak_ratio_display_lower'],'target':out['finite_L1']['strict_rational_target'],'margin':out['finite_L1']['margin_over_target_display']},indent=2))
