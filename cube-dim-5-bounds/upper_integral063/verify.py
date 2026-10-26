#!/usr/bin/env python3
"""Fraction geometry/readback and independent384bit power bounds; no H/jet/log."""
from pathlib import Path
from fractions import Fraction as F
from math import prod
import gzip,json,hashlib,os,time
from flint import arb,ctx,fmpq
ctx.prec=384
def A(x):
 x=F(x);return arb(fmpq(x.numerator,x.denominator))
def upper(x):return str(x.upper().fmpq())
H=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
reg=json.loads((H/'registration.json').read_text());receipt=json.loads((H/'run_receipt.json').read_text());began=time.monotonic();checks=[]
assert sha(H/'pilot.py')==receipt['source_sha256'];assert sha(H/'registration.json')==receipt['registration_sha256']
for s in reg['sources']:assert sha(s['source'])==s['sha256']
table=json.loads(gzip.decompress(Path(reg['sources'][0]['source']).read_bytes()));vertices=table['vertices'];summary=[]
for row in reg['anchors']:
 r=json.loads((H/(row['name']+'_results.json')).read_text());assert r['fingerprint']==receipt['fingerprint'];assert r['full3504_jet_calls']==1;center=list(map(F,r['selected_center']['point']));D=r['selected_center']['D'];anchor=list(map(F,row['nominal_anchor']));nodes=[sum(int(x*D)*z for x,z in zip(center,v)) for v in vertices];assert len(set(nodes))==280;a=144*D*F(reg['eps']);nomnodes=[sum(x*z for x,z in zip(anchor,v)) for v in vertices];zeroids=[i for i,k in enumerate(nomnodes) if k==0]
 assert row['candidates_in_order'][row['selected_candidate']]==r['selected_center']
 for ci,c in enumerate(row['candidates_in_order']):
  p=list(map(F,c['point']));ns=[sum(int(x*c['D'])*z for x,z in zip(p,v)) for v in vertices];assert len(set(ns))==c['distinct'];assert min(map(abs,ns))==c['min_abs_projection']
  if ci<row['selected_candidate']:assert c['distinct']<280
 previous=None; independent=[]
 for s in r['scans']:
  h=F(s['h']);lo=[max(F(0),p-h) for p in anchor[1:]];hi=[min(F(1),p+h) for p in anchor[1:]];ul=[l-c for l,c in zip(lo,center[1:])];uh=[u-c for u,c in zip(hi,center[1:])];radii=[max(abs(l),abs(u)) for l,u in zip(ul,uh)]
  assert list(map(F,s['box_lower']))==lo and list(map(F,s['box_upper']))==hi
  assert list(map(F,s['u_lower']))==ul and list(map(F,s['u_upper']))==uh
  assert all(l<=p<=u for l,p,u in zip(lo,anchor[1:],hi));assert all(l<=c<=u for l,c,u in zip(lo,center[1:],hi));assert F(s['volume'])==prod(u-l for l,u in zip(lo,hi))
  ratios=[]
  for v,k in zip(vertices,nodes):
   low=D*sum(z*(l if z>=0 else u) for z,l,u in zip(v[1:],ul,uh));high=D*sum(z*(u if z>=0 else l) for z,l,u in zip(v[1:],ul,uh));ratios.append(max(abs(low),abs(high))**2/(a*a+k*k))
  assert max(ratios)==F(s['maximum_exact_q_squared']);assert s['guard_passed']==(max(ratios)<1)
  if not s['guard_passed']:continue
  assert max(ratios)<=F(s['max_q_rational_upper'])**2<1
  moments=[prod((u**(n+1)-l**(n+1))/((n+1)*(u-l)) for l,u,n in zip(ul,uh,c['alpha'])) for c in r['coefficients']];assert moments==list(map(F,s['exact_mean_moments']))
  # Rational endpoints independently bound the polynomial variation.
  V=sum(max(abs(F(c['enclosure']['lower'])),abs(F(c['enclosure']['upper'])))*prod(rr**n for rr,n in zip(radii,c['alpha'])) for c in r['coefficients'][1:]);vU=F(s['variation']['upper']);assert V<=vU
  rem=F(s['all_time_resolvent_remainder']['upper']);assert rem>=0
  plo=sum(min(F(c['enclosure']['lower'])*m,F(c['enclosure']['upper'])*m) for c,m in zip(r['coefficients'],moments));phi=sum(max(F(c['enclosure']['lower'])*m,F(c['enclosure']['upper'])*m) for c,m in zip(r['coefficients'],moments))
  assert abs(F(s['mean_Taylor']['lower'])-plo)<F(1,2**280) and abs(F(s['mean_Taylor']['upper'])-phi)<F(1,2**280)
  L,U=F(s['H_range_lower']),F(s['H_range_upper']);mL,mU=F(s['mean_H_lower']),F(s['mean_H_upper']);assert 0<=L<=mL<=mU<=U
  LL=max(F(0),F(r['center_H']['lower'])-V-rem);UU=F(r['center_H']['upper'])+V+rem
  # Independent outward endpoints may enlarge intermediate balls by a few ulps.
  LL=min(L,LL);UU=max(U,UU);MM=min(UU,max(mU,phi+rem))
  QQ=A(F(5,4));ch=A(s['volume'])*(A(LL)**QQ+(A(UU)**QQ-A(LL)**QQ)/A(UU-LL)*A(MM-LL))
  assert abs(F(upper(ch))-F(s['chord_integral_upper']))<F(1,2**280)
  independent.append({'h':s['h'],'range_lower':str(LL),'range_upper':str(UU),'mean_upper':str(MM),'chord_integral_upper':upper(ch)})
  assert s['positive_lower_certificate']==(L>0);assert s['remainder_le_one_percent_center']==(100*rem<=F(r['center_H']['lower']))
  assert 0<F(s['chord_integral_upper'])<=F(s['constant_range_integral_upper'])
  if previous is not None:assert rem<previous
  previous=rem
 for field,key in [('max_guard_h','guard_passed'),('max_positive_lower_h','positive_lower_certificate'),('max_one_percent_h','remainder_le_one_percent_center')]:assert r[field]==next((s['h'] for s in r['scans'] if s.get(key)),None)
 cp=r['checkpoint'];assert sha(cp['path'])==cp['sha256'];saved=json.loads(Path(cp['path']).read_text());assert saved['next_simplex']==3504 and saved['fingerprint']==r['fingerprint']
 summary.append({'name':row['name'],'nominal_distinct_projection_count':len(set(nomnodes)),'nominal_zero_projection_count':len(zeroids),'nominal_zero_witness_vertices':[vertices[i] for i in zeroids[:4]],'all_11_box_anchor_and_center_containments':True,'exact_guard_and_clipped_moments_passed':True,'remainders_strictly_decrease_on_registered_passed_scales':True,'independent_rational_envelope_and_384bit_power_readback':independent,'max_guard_h':r['max_guard_h'],'max_positive_lower_h':r['max_positive_lower_h'],'max_one_percent_h':r['max_one_percent_h']})
budget=json.loads((H/'budget.json').read_text());assert 3125*F(1,25000)**4<30**4*F(1,10000)**5
assert F(budget['conditional_cost']['upper'])<30
assert F(budget['UH_gate_with_1e6_tube']['upper'])<F(budget['UH_gate_zero_tube']['lower'])<F(1,25000)
assert F(budget['I_if_UH_equals_target_zero_tube']['lower'])>F(1,25000)
out={'status':'exact_geometry_rational_envelope_and_384bit_chord_readback_passed','pid':os.getpid(),'seconds':time.monotonic()-began,'production_pid':receipt['pid'],'production_source_sha256':sha(H/'pilot.py'),'registration_sha256':sha(H/'registration.json'),'verifier_sha256':sha(__file__),'new_H_or_jet_evaluations':0,'center_checks':summary,'conditional_cost_power_check':True,'numeric_centers_do_not_prove_global_I_or_g':True}
(H/'readback.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
