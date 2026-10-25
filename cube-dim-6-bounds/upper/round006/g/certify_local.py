#!/usr/bin/env python3
"""Finite Arb body32 H_epsilon batch with rigorous local box extension.

The complete signed-permutation group and exact confluent divided differences
are the previously audited point formula. This file groups exact identical
images, and caches only identical exact node lists. All transcendental and
final sign computations use Arb, not binary64.
"""
from collections import Counter
from fractions import Fraction as Q
from functools import cmp_to_key
from pathlib import Path
from time import perf_counter
import hashlib, itertools, json, math, os
from flint import arb, acb, ctx

HERE=Path(__file__).resolve().parent
MAIN=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6')
PROJECT=MAIN/'cube-dim-6-bounds'
FAN=PROJECT/'upper/round003/direct_kernel/fan_integer_input.txt'
ctx.prec=384
POINTS=[
 ('minimum_tie',[1000,1000,893,807,467,97]),
 ('minimum_interior',[1000,990,893,807,467,97]),
 ('triple_tie_boundary',[1000,1000,1000,750,250,0]),
 ('S2_one_four',[1000,750,750,750,750,750]),
 ('S1_vertex',[1000,1000,1000,1000,1000,1000]),
 ('S2_axis',[1000,0,0,0,0,0]),
 ('S2_two_equal',[1000,1000,0,0,0,0]),
 ('generic_control',[1000,900,800,600,400,200]),
]
EPS=[Q(1,200),Q(1,100),Q(1,50)]
INPUTS=[MAIN/'AGENTS.md',
 PROJECT/'coordination/phase_20261025/protocol.md',
 PROJECT/'coordination/phase_20261025/frozen_thresholds.json',
 PROJECT/'coordination/phase_20261025/roles.json',PROJECT/'STATE.md',
 PROJECT/'results/bounds_ledger.json',
 PROJECT/'upper/round004/spatial/005_spatial_multiplicative_certificate.md',
 PROJECT/'upper/round004/spatial/006_second_order_and_vertex_certificates.md',
 PROJECT/'upper/round004/local_box/arb384_box.py',
 PROJECT/'upper/round005/primary/new_tail_budget.md',
 PROJECT/'upper/round005/primary/coarse_mesh_obstruction.md',FAN]

def A(q): return arb(q.numerator)/q.denominator
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def cmp(a,b): return (a[0]*b[1]>b[0]*a[1])-(a[0]*b[1]<b[0]*a[1])
def eq(a,b): return a[0]*b[1]==b[0]*a[1]
def normnode(n,d):
 g=math.gcd(n,d);return n//g,d//g
def jet(x,e,k):
 if not k:return (x/e).atan()
 i=acb(0,1)
 z=(i**k/acb(e,x)**k-(-i)**k/acb(e,-x)**k)/(2*k*i)
 return (1 if k&1 else -1)*z.real
def dd5(z,e):
 v=[arb(n)/d for n,d in z]
  # The explicit triangular table avoids any ambiguity at repeated nodes.
 t=[[arb(0) for _ in range(6)]for _ in range(6)]
 for j in range(6):t[j][0]=jet(v[j],e,0)
 for k in range(1,6):
  for j in range(k,6):t[j][k]=jet(v[j],e,k) if eq(z[j],z[j-k]) else (t[j][k-1]-t[j-1][k-1])/(v[j]-v[j-k])
 return t[5][5]
def fan():
 it=iter(FAN.read_text().split());cones=[]
 for _ in range(int(next(it))):
  det=int(next(it));rs=[];alpha=Q(det)
  for _ in range(6):
   q=int(next(it));r=tuple(int(next(it))for _ in range(6));rs.append((q,r));alpha/=q
  cones.append((alpha,rs))
 return cones
def orbit(x):
 count=Counter()
 for p in itertools.permutations(range(6)):
  for m in range(64):count[tuple((-1 if m>>j&1 else 1)*x[p[j]]for j in range(6))]+=1
 assert sum(count.values())==46080
 return count
def compute(name,x,cones,started):
 tic=perf_counter();group=orbit(x);weighted={};equal_events=0
 for y,mult in group.items():
  for alpha,rs in cones:
   z=tuple(sorted((normnode(sum(y[j]*r[j]for j in range(6)),1000*q)for q,r in rs),key=cmp_to_key(cmp)))
   equal_events+=mult*sum(eq(z[j-1],z[j])for j in range(1,6))
   # atan is odd and the divided-difference order is odd, so reflection leaves DD5 unchanged.
   neg=tuple((-n,d)for n,d in reversed(z));z=min(z,neg)
   weighted[z]=weighted.get(z,Q(0))+mult*alpha
 totals=[arb(0)for _ in EPS]
 for z,w in weighted.items():
  if perf_counter()-started>1200:raise TimeoutError('registered 20 minute limit')
  aa=A(w)
  for j,e in enumerate(EPS):totals[j]+=aa*dd5(z,A(e))
 vals=[v/(2*arb.pi())**6 for v in totals]
 return {'name':name,'point_numerators':x,'denominator':1000,'unique_group_images':len(group),
  'raw_group_images':46080,'complete_cone_terms':967680,'unique_exact_DD_lists':len(weighted),
  'equal_adjacent_events_with_multiplicity':equal_events,
  'values':[{'epsilon':str(e),'H_arb384':v.str(100,more=True)}for e,v in zip(EPS,vals)],
  'elapsed_seconds':perf_counter()-tic}
def boxes(rows):
 out=[];target=arb(1)/60000
 for row in rows:
  options=[]
  for d in row['values']:
   e=Q(d['epsilon']);v=arb(d['H_arb384']);r=0
   # Exact rational radii at 1e-5 increments, capped at 0.1; every acceptance is an Arb comparison.
   if v>target:
    lo,hi=0,10001
    while hi-lo>1:
     mid=(lo+hi)//2
     if v*(-A(Q(mid,100000))/A(e)).exp()>target:lo=mid
     else:hi=mid
    r=lo
   options.append((r,e,v))
  r,e,v=max(options,key=lambda o:o[0]);radius=Q(r,100000)
  x=[Q(a,1000)for a in row['point_numerators']]
  bounds=[[str(max(Q(0),a-radius)),str(min(Q(1),a+radius))]for a in x[1:]]
  lower=v*(-A(radius)/A(e)).exp();assert r==0 or lower>target
  out.append({'name':row['name'],'face':'x1=1','bounds_x2_to_x6':bounds,'radius':str(radius),
   'epsilon':str(e),'G_lower_arb384':lower.str(95,more=True),'strict_above_1_over_60000':bool(lower>target),
   'extension':'G>=H_epsilon>=exp(-radius/epsilon)*H_epsilon(center)'})
 return out
def main():
 started=perf_counter();cones=fan();rows=[]
 result={'status':'RUNNING_FINITE_REGISTERED_BATCH','pid':os.getpid(),'arithmetic':'python-flint Arb 384 bits',
  'limits':{'max_boxes':8,'max_point_epsilon_pairs':24,'max_seconds':1200},
  'input_sha256':{str(p):sha(p)for p in INPUTS},'source_sha256':sha(Path(__file__)),
  'outputs':rows,'uncovered':'all of [0,1]^5 outside the union of the listed clipped boxes; no global g claim'}
 (HERE/'input_manifest.json').write_text(json.dumps({k:v for k,v in result.items()if k!='outputs'},indent=2)+'\n')
 for name,x in POINTS:
  try:row=compute(name,x,cones,started)
  except TimeoutError as ex:
   result['stop_reason']=str(ex);break
  rows.append(row);result['boxes']=boxes(rows);result['elapsed_seconds']=perf_counter()-started
  (HERE/'local_results.json').write_text(json.dumps(result,indent=2)+'\n')
  print(name,row['elapsed_seconds'],[(d['epsilon'],d['H_arb384'][:28])for d in row['values']],flush=True)
 result['status']='COMPLETE_LOCAL_ONLY' if len(rows)==8 else 'STOPPED_LOCAL_ONLY'
 result['elapsed_seconds']=perf_counter()-started;result['boxes']=boxes(rows)
 result['new_C6_bound']=False;result['global_g_target_certified']=False
 c=Q(4225,4224)
 result['new_ODE_coefficient_exact']=str(c);result['new_ODE_coefficient_arb384']=A(c).str(95,more=True)
 result['old_over_new_radius_factor']=((arb(129)/128)/A(c)).sqrt().str(95,more=True)
 (HERE/'local_results.json').write_text(json.dumps(result,indent=2)+'\n')
 print(result['status'],'nodes',len(rows)*3,'boxes',len(rows),'seconds',result['elapsed_seconds'],flush=True)
if __name__=='__main__':main()
