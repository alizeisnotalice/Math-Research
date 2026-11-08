import json,sys,re,math,hashlib
from pathlib import Path
from decimal import Decimal
from fractions import Fraction
p=Path(sys.argv[1]);out=Path(sys.argv[2])
def num(s):return Decimal(re.sub(r'`[^*]*','',s).replace('*^','E'))
a=json.load(open(p/'python/optimize76.json'));b=json.load(open(p/'wolfram/optimize76.json'));ab={(r['N'],r['kind']):r for r in a['rows']};low={(r['N'],r['kind']):r for r in b['rows'] if r['precision']==60};pe=[];ce=[];qerr=[]
for r in b['rows']:
 assert r['numeric']
 if r['precision']!=120:continue
 key=(r['N'],r['kind']);j=ab[key];assert r['branch']==j['branch']
 for f in ['tau','q','log_bound','reference_tau_frontier']:
  if r[f] is None:assert j[f] is None and low[key][f] is None;continue
  v=num(r[f]);pe.append(abs(v-num(low[key][f])));ce.append(abs(float(v)-j[f]))
 qerr.append(abs(float(num(r['q']))/j['q']-1))
assert b['symbolic_convexity_identity'] and b['symbolic_optimal_value_identity']
assert max(pe)<Decimal('1e-20') and max(ce)<2e-7 and max(qerr)<1e-8
assert len(a['rows'])==140 and len(b['rows'])==280
assert all(r['direct_minimizer_success'] for r in a['rows'])
assert max(abs(r['direct_minimizer_gap']) for r in a['rows'])<1e-6
assert max(abs(r['q_equation_residual']) for r in a['rows'] if r['q_equation_residual'] is not None)<1e-13
cp=json.load(open(p/'python/cover76.json'));cw=json.load(open(p/'wolfram/cover76.json'));cb={(r['N'],r['K']):r for r in cp['rows']};cl={(r['N'],r['K']):r for r in cw['rows'] if r['precision']==60};cpe=[];cce=[]
for r in cw['rows']:
 assert r['numeric']
 if r['precision']!=120:continue
 key=(r['N'],r['K'])
 for f in ['maximum_grid_ratio_at_leading_frontier','leading_Hardy_charge_per_dimension']:
  v=num(r[f]);cpe.append(abs(v-num(cl[key][f])));cce.append(abs(float(v)-cb[key][f]))
assert cw['symbolic_grid_monotonicity'] and cw['symbolic_grid_limit_identity']
assert max(cpe)<Decimal('1e-20') and max(cce)<1e-10
for r in cp['exact_reuse_packet']:assert Fraction(r['summed_modular'])==r['ancestor_count']*Fraction(r['source_modular']) and r['reuse_factor']==r['ancestor_count']
pa=json.load(open(p/'python/copies76.json'));pb=json.load(open(p/'wolfram/copies76.json'));pd={(r['N'],r['copies'],r['volume_enlargement']):r for r in pa['rows']}
assert len(pa['rows'])==120 and len(pb['rows'])==120
assert pb['symbolic_score_invariance']
for r0 in pb['rows']:
 assert r0['exact_pass'];s=pd[r0['N'],r0['copies'],r0['volume_enlargement']]
 assert Fraction(s['projection_separation_margin'])==Fraction(r0['projection_separation_margin'])>0
 assert s['continuous_corner_count_lower']==s['explicit_exact_corner_count_upper']==r0['copies'] and s['source_allocation_reuse']==1
r={'exact_separated_copy_certificates':len(pa['rows']),'pass':True,'Python_optimizations':len(a['rows']),'Wolfram_optimizations':len(b['rows']),'optimization_precision_fields':len(pe),'optimization_precision_max_absolute':str(max(pe)),'optimization_cross_fields':len(ce),'optimization_cross_max_absolute':max(ce),'q_cross_max_relative':max(qerr),'direct_numerical_minimizations':len(a['rows']),'maximum_direct_minimizer_objective_gap':max(abs(r['direct_minimizer_gap']) for r in a['rows']),'Python_grid_cases':len(cp['rows']),'Wolfram_grid_cases':len(cw['rows']),'grid_precision_fields':len(cpe),'grid_precision_max_absolute':str(max(cpe)),'grid_cross_fields':len(cce),'grid_cross_max_absolute':max(cce),'exact_packet_budget_identities':len(cp['exact_reuse_packet']),'scope':'Analytic conditional corner and entropy proofs are separate; numerical reference frontiers do not certify C_N or an attained covering.'};r['source_sha256']={name:hashlib.sha256((Path(__file__).resolve().parent/name).read_bytes()).hexdigest() for name in ['optimize76.py','optimize76.wl','cover76.py','cover76.wl','copies76.py','copies76.wl','check76.py']};out.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
