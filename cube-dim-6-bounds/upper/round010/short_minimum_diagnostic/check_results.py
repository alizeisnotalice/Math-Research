#!/usr/bin/env python3
"""No kernel replay: exact geometry, original-source hashes and strict target arithmetic."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
from fractions import Fraction as F
from time import perf_counter
import hashlib,itertools,json,os
import flint
from flint import arb,fmpq,fmpz,ctx
ctx.prec=384
HERE=Path(__file__).resolve().parent;PROJECT=HERE.parents[2]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def A(x):return arb(x.numerator)/x.denominator
def reconstruct_g(g):return arb(fmpq(fmpz(g['G_times_pi5_numerator']),fmpz(g['G_times_pi5_denominator'])))/arb.pi()**5

def geometry(ns,den):
 x=list(map(lambda n:F(n,den),ns));assert x[0]==1 and all(x[i]>x[i+1]>0 for i in range(5))
 V=[(1,)+s for s in itertools.product((-1,1),repeat=5)];ss=sum(t*t for t in x)
 best=None;which=None
 for i,v in enumerate(V):
  b1=sum(t*a for t,a in zip(x,v))
  for j in range(i):
   w=V[j];b2=sum(t*a for t,a in zip(x,w));g=sum(a*b for a,b in zip(v,w))
   d2=ss-(6*b1*b1-2*g*b1*b2+6*b2*b2)/(36-g*g)
   if best is None or d2<best:best=d2;which=[j,i]
 assert best>F(1,80000)**2
 return {'strict_ordered_positive_face':True,'S2_distance_squared_exact':str(best),'closest_pair_indices':which,'S2_distance_arb':A(best).sqrt().str(95,more=True),'outside_A009_barrier_neighborhood':True}
start=perf_counter()
grad=json.loads((HERE/'gradient_results.json').read_text());p1=json.loads((HERE/'point1_results.json').read_text());oth=json.loads((HERE/'other_time_results.json').read_text());b2=json.loads((HERE/'batch2_results.json').read_text());pt2=json.loads((HERE/'point2.json').read_text())
records=[grad,p1,oth,b2,pt2,json.loads((HERE/'selection_results.json').read_text())]
checked_hashes=set()
for r in records:
 for inp in r.get('inputs',[]):
  key=(inp['path'],inp['sha256'])
  if key not in checked_hashes:assert sha(inp['path'])==inp['sha256'],inp['path'];checked_hashes.add(key)
assert sha(PROJECT/'upper/round003/direct_kernel/fan_integer_input.txt')=='7a57caad9b7713d9ce0896e51ab1dd402566feb1cba50829104c593a6d12a812'
old=json.loads((PROJECT/'upper/round007/g_time_moments/results.json').read_text())
for row in grad['gradient_passes']:
 prev=next(x for x in old['time_nodes']if x['tau']==row['tau'])
 assert (arb(row['H_arb384'])-arb(prev['H_arb384'])).contains(0)
 assert row['complete_group_terms']==967680
assert grad['new_full_passes']==2 and p1['new_full_passes_this_point']==2 and oth['new_H_passes']==1 and b2['new_full_passes']==2
assert all(r['exit_code']==0 for r in [grad,p1,oth,b2])
G1=reconstruct_g(p1['G']);G2=reconstruct_g(b2['G'])
for g,G in [(p1['G'],G1),(b2['G'],G2)]:
 assert g['complete_group_terms']==967680 and g['summed_divergent_Laurent_coefficients']==['0','0','0']
 assert (G-arb(g['G_arb384']['ball'])).contains(0)
K64=G1-arb(p1['H_rows'][0]['H_arb384']);K32first=G1-arb(oth['H_arb384']);K32second=G2-arb(b2['H_arb384'])
assert K64<arb(31)/2000000<arb(1)/64000
assert K32first>arb(781516)/10**11>arb(1)/128000
assert K32first<arb(781517)/10**11<arb(1)/125000
assert K32second<arb(1)/128000
xx=[F(n,p1['denominator'])for n in p1['point_numerators']];cc=[F(n,grad['denominator'])for n in grad['center_numerators']];yy=[F(n,b2['denominator'])for n in b2['point_numerators']]
assert yy[0]==1 and yy[1]==xx[1] and all(yy[i]==xx[i]+(xx[i]-cc[i])/50 for i in range(2,6))
geo1=geometry(p1['point_numerators'],p1['denominator']);geo2=geometry(b2['point_numerators'],b2['denominator'])
out={'status':'PASS_TWO_FIXED_TIMES_STRICTLY_REFUTED','pid':os.getpid(),'elapsed_seconds':perf_counter()-start,'kernel_replayed':False,'new_kernel_passes':0,'cumulative_source_new_passes':7,'hash_checks_count':len(checked_hashes),'K064_point1':K64.str(100,more=True),'K032_point1':K32first.str(100,more=True),'K032_point2':K32second.str(100,more=True),'strict_refutation_times':['4/125','8/125'],'does_not_refute_all_T':True,'global_minimum_claimed':False,'new_C6_bound':False,'point1_geometry':geo1,'point2_geometry':geo2,'sources':[{'path':str(HERE/f),'sha256':sha(HERE/f)}for f in ['gradient_results.json','point1_results.json','other_time_results.json','batch2_results.json','point2.json']],'source_sha256':sha(__file__),'python':sys.version,'python_flint':flint.__version__,'arb_bits':ctx.prec,'exit_code':0}
(HERE/'checked_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
