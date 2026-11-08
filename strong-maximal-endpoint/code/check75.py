import json,sys,re,math
from pathlib import Path
from decimal import Decimal
p=Path(sys.argv[1]);out=Path(sys.argv[2])
def num(s):return Decimal(re.sub(r'`[^*]*','',s).replace('*^','E'))
a=json.load(open(p/'python/transport75.json'));b=json.load(open(p/'wolfram/transport75.json'));key=lambda r:(r['N'],r['kind'],r.get('tube_power'));base={key(r):r for r in a['rows'] if r['refinement']==1};fine={key(r):r for r in a['rows'] if r['refinement']==.1};ref=max(abs(r['pushforward_ratio']/base[key(r)]['pushforward_ratio']-1) for r in fine.values());lo={key(r):r for r in b['rows'] if r['precision']==50};pe=[];ce=[]
for r in b['rows']:
 assert r['numeric']
 if r['precision']==100:
  k=key(r)
  for f in ['receiver_weighted_integral','image_integral','pushforward_ratio']:
   v=num(r[f]);pe.append(abs(v-num(lo[k][f]))/abs(v));ce.append(abs(float(v)/fine[k][f]-1))
assert b['symbolic_rank_one_identity'] and all(r['exact_pass'] for r in b['matrices'])
assert ref<1e-10 and max(pe)<Decimal('1e-12') and max(ce)<1e-10
maxlog=max(abs(r['det_log_error']) for r in a['matrices'] if r['q']>0);assert maxlog<1e-9
op=json.load(open(p/'python/optimizer75.json'));ow=json.load(open(p/'wolfram/optimizer75.json'));ob={(r['N'],s['x_offset']):(r,s) for r in op['rows'] for s in r['rows']};ol={(r['N'],r['offset']):r for r in ow['rows'] if r['precision']==60};ope=[];oce=[]
for r in ow['rows']:
 assert r['numeric']
 if r['precision']==120:
  k=(r['N'],r['offset']);j,s=ob[k]
  for f,pv in [('c_optimizer',s['c_optimizer']),('reanchored_receiver',s['reanchored_receiver']),('tensor_log_cap_upper',j['source_tensor_log_cap_upper']),('receiver_log_min_average',j['receiver_cube_log_min_average'])]:
   v=num(r[f]);ope.append(abs(v-num(ol[k][f])));oce.append(abs(float(v)-pv))
assert ow['symbolic_piecewise_integral_identity'] and max(ope)<Decimal('1e-20') and max(oce)<1e-10
cp=json.load(open(p/'python/corner75.json'));cw=json.load(open(p/'wolfram/corner75.json'));cb={r['N']:r for r in cp['rows']};cl={r['N']:r for r in cw['rows'] if r['precision']==50};cpe=[];cce=[]
for r in cw['rows']:
 assert r['numeric']
 if r['precision']==100:
  for f in ['log_Hardy_bound','log_weighted_density_center','log_weighted_density_lower','log_weighted_density_upper']:
   v=num(r[f]);cpe.append(abs(v-num(cl[r['N']][f])));cce.append(abs(float(v)-cb[r['N']][f]))
assert cw['symbolic_Orlicz_endpoint_and_convexity'] and cw['symbolic_affine_Jacobian']
assert max(cpe)<Decimal('1e-15') and max(cce)<1e-7
ep=json.load(open(p/'python/envelope75.json'));ew=json.load(open(p/'wolfram/envelope75.json'));import hashlib
assert hashlib.sha256((p/'python/envelope75.json').read_bytes()).hexdigest()==ew['candidate_SHA256']
eb={(j['N'],j['K'],s['radial_node']):s for j in ep['rows'] for s in j['rows']};el={(r['N'],r['K'],r['radial_node']):r for r in ew['rows'] if r['precision']==50};epe=[];ece=[]
for r0 in ew['rows']:
 assert r0['numeric']
 if r0['precision']!=100:continue
 k0=(r0['N'],r0['K'],r0['radial_node'])
 for fld in ['log_envelope_ratio','log_shift_comparison_margin']:
  if r0[fld] is None:assert el[k0][fld] is None and eb[k0][fld] is None;continue
  v=num(r0[fld]);epe.append(abs(v-num(el[k0][fld])));ece.append(abs(float(v)-eb[k0][fld]))
assert max(epe)<Decimal('1e-15') and max(ece)<1e-8
r={'envelope_Python_parameter_sets':len(ep['rows']),'envelope_Wolfram_nodes':len(ew['rows']),'envelope_cross_fields':len(ece),'envelope_cross_max_absolute':max(ece),'envelope_precision_max_absolute':str(max(epe)),'pass':True,'transport_Python_rows':len(a['rows']),'transport_Wolfram_rows':len(b['rows']),'transport_refinement_pairs':len(fine),'transport_refinement_max_relative':ref,'transport_precision_fields':len(pe),'transport_precision_max_relative':str(max(pe)),'transport_cross_fields':len(ce),'transport_cross_max_relative':max(ce),'Python_Jacobian_matrices':len(a['matrices']),'positive_matrix_max_logdet_error':maxlog,'Wolfram_exact_rational_determinants':len(b['matrices']),'optimizer_Python_sources':len(op['rows']),'optimizer_Wolfram_roots':len(ow['rows']),'optimizer_cross_fields':len(oce),'optimizer_cross_max_absolute':max(oce),'optimizer_precision_max_absolute':str(max(ope)),'corner_Python_rows':len(cp['rows']),'corner_Wolfram_rows':len(cw['rows']),'corner_cross_fields':len(cce),'corner_cross_max_absolute':max(cce),'corner_precision_max_absolute':str(max(cpe)),'scope':'Exact geometric/Hardy proofs are separate; numerical probe consistency is not a general endpoint or source-cost certificate.'};out.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
