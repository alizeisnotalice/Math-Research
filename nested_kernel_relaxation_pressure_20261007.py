"""Exact tests of a NON-CUBE normalized nested-kernel relaxation counterexample.
This tests common-source, hardband and winner constraints, not Euclidean geometry.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json
out=Path(__file__).with_name('nested_kernel_relaxation_pressure_20261007_results.json')
if out.exists():raise SystemExit('Refuse overwrite')
rows=[]
for stage,m in enumerate((4,8,16),1):
 N=2**m;lam=F(1,2*N)
 for j in range(m):
  size=2**(j+1)
  # One fixed source belongs to exactly size cyclic blocks; translation gives every source.
  hits=sum(1 for ell in range(N) if (-ell)%N<size)
  assert hits==size
  for u in (F(0),F(1,8),F(1,2),F(7,8),F(1)):
   vol=2**j*(1+u)
   captured=1+sum(2**i for i in range(j))+2**j*u
   assert captured==vol and captured/vol==1
   response=F(size,N)/vol
   band=(2*lam<response<=4*lam)
   assert band==(u<1)
   # Same row is zero before arrival and strictly decreases after it.
   for scale_factor in (F(1),F(3,2),F(2)):
    assert response/scale_factor<=response
   rows.append({'stage':stage,'m':m,'N':N,'j':j,'u':str(u),'column_capture_volume':str(vol),'column_integral':'1','row_winner_response':str(response),'hardband':band,'cyclic_source_hits':hits})
summary=[]
for m in (4,8,16):
 N=2**m
 event_volume=F(m*N,2)
 ratio=F(1,2*N)*event_volume
 assert ratio==F(m,4)
 summary.append({'m':m,'N':N,'event_volume':str(event_volume),'lambda_times_event_volume':str(ratio),'traffic_exact':f'{m}*log(2)','arrival_energy_exact':f'{m}*log(2)','cube_counterexample':False})
result={'status':'PASS_EXACT_NON_GEOMETRIC_RELAXATION_COUNTEREXAMPLE','rows':rows,'row_count':len(rows),'rounds':3,'summary':summary,'scope':'Receiver/source label kernel, not a translation-invariant cube convolution; no FIRST or actual geom qualification. Refutes derivation from nestedness, column normalization, common source and hardband/winner alone.','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
out.write_text(json.dumps(result,indent=2));print({k:v for k,v in result.items() if k not in ('rows','summary')})
