import json,sys,re,hashlib,math
from pathlib import Path
from fractions import Fraction as Q
from decimal import Decimal
def dec(s):return Decimal(re.sub(r'\x60[^*]*','',s).replace('*^','E'))
p=Path(sys.argv[1]);out=Path(sys.argv[2])
src=json.loads((p/'python/columns78.json').read_text());v=json.loads((p/'wolfram/verification78.json').read_text());di=json.loads((p/'wolfram/diagnostics78.json').read_text())
assert src['case_count']==10 and src['column_count']==618 and v['pass'] and di['pass']
assert len(v['source_columns'])==618 and all(c['exact_pass'] for c in v['source_columns'])
assert len(v['masters'])==10 and len(v['triangles'])==2
for r in v['masters']+v['triangles']:assert r['exact_pass'] and r['Wolfram_master_exact'] and not r['messages']
for ci,c in enumerate(src['cases']):
 expected=[(a,s) for a,ids in enumerate(c['eligible']) for s in range(1,64) if s&sum(1<<i for i in ids)==s]
 assert [(z['anchor_index'],z['mask']) for z in c['columns']]==expected
 d=c['cover'];assert Q(d['fractional_lower'])<=Q(d['fractional_upper']) and Q(d['integer_lower'])<=Q(d['integer_upper'])
 assert Q(d['integer_lower'])<=Q(d['fractional_upper']) and Q(d['integer_upper'])>=Q(d['fractional_lower'])
 assert {(r['anchor_index'],r['mask']) for r in v['source_columns'] if r['case']==ci}==set(expected)
# The analytic half-pair dual and exact natural-corner repair avoid interpreting a sampled pricing run as a proof.
tr=json.loads((p/'python/triangle78.json').read_text());an=json.loads((p/'python/analytic78.json').read_text())['cases']
for c in tr['cases']:
 expected=[(a,s) for a,ids in enumerate(c['eligible']) for s in range(1,8) if s&sum(1<<i for i in ids)==s]
 assert [(z['anchor_index'],z['mask']) for z in c['columns']]==expected
a=an[0];lr=Q(a['integer_over_fractional_lower']);ur=Q(a['integer_over_fractional_upper'])
assert Q('1.3331609908')<lr<ur<Q('1.3331609909')
assert di['analytic_triangle']['exact_pass'] and di['analytic_triangle']['natural_corner_exact_no_gap']
assert abs(float(dec(di['analytic_triangle']['triangle_integer_fractional_ratio_lower']))-float(lr))<1e-14
assert abs(float(dec(di['analytic_triangle']['triangle_integer_fractional_ratio_upper']))-float(ur))<1e-14
assert an[1]['integer_equals_fractional_exactly']
merge=json.loads((p/'python/merge78.json').read_text());dm=di['aggregate'];assert dm['exact_pass'] and Q(merge['minimum_exact_vertex_slack'])>=0
assert Q(merge['rounding_excess_upper'])<=Q(1,10**6)
pr=dm['precision_rows'];assert [z['precision'] for z in pr]==[60,120]
merger_errors=[abs(dec(pr[0][f])-dec(pr[1][f])) for f in ['aggregate_actual_modular','weighted_actual_modular','minimum_actual_mean','maximum_upward_rounding']]
assert max(merger_errors)<Decimal('1e-30')
cp=json.loads((p/'python/cross78.json').read_text())['rows'];cw=di['cross_anchor'];assert len(cp)==14 and len(cw)==28 and all(r['numeric_pass'] for r in cw)
cm={r['N']:r for r in cp};cl={r['N']:r for r in cw if r['precision']==60};errors=[];cross={f:0. for f in ['log_minimum_actual_mean','log_aggregate_height','log_receiver_volume']};meanrel=0
for r in cw:
 if r['precision']!=120:continue
 for f in cross:
  errors.append(abs(dec(r[f])-dec(cl[r['N']][f])))
  cross[f]=max(cross[f],abs(float(dec(r[f]))-cm[r['N']][f]))
 meanrel=max(meanrel,abs(float(dec(r['log_minimum_actual_mean']))/cm[r['N']]['log_minimum_actual_mean']-1))
assert max(errors)<Decimal('1e-25') and cross['log_minimum_actual_mean']<1e-14 and cross['log_aggregate_height']<1e-13 and cross['log_receiver_volume']<1e-6 and meanrel<1e-12
r=dict(pass_all=True,complete_grid_families=10,certified_grid_source_columns=618,analytic_triangle_columns=[9,16],independent_exact_master_LPs=12,
 grid_integer_fractional_ratio_upper=max(c['cover']['certified_integer_over_fractional_upper'] for c in src['cases']),
 triangle_ratio_lower=str(lr),triangle_ratio_upper=str(ur),triangle_ratio_decimal=[float(lr),float(ur)],natural_corner_exactly_removes_gap=True,
 same_anchor_rational_merge_minimum_vertex_slack=merge['minimum_exact_vertex_slack'],merge_rounding_cost_excess_upper=merge['rounding_excess_upper'],
 aggregate_precision_difference_max=str(max(merger_errors)),cross_anchor_Python_dimensions=14,cross_anchor_Wolfram_evaluations=28,
 cross_anchor_precision_difference_max=str(max(errors)),cross_anchor_cross_software_max_absolute=cross,mean_log_max_relative_difference=meanrel,
 scope='Same-anchor consolidation is a proved countable-source theorem. All numerical covering optima refer to explicitly listed finite receiver cells and anchors. No new general C_N upper/lower or submodularity claim.')
r['source_sha256']={n:hashlib.sha256((Path(__file__).resolve().parent/n).read_bytes()).hexdigest() for n in ['columns78.py','triangle78.py','merge78.py','verify78.wl','diagnostics78.wl','check78.py','allocate77.py']}
r['input_sha256']={f:hashlib.sha256((p/'python'/f).read_bytes()).hexdigest() for f in ['columns78.json','triangle78.json','analytic78.json','merge78.json','cross78.json']}
out.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items() if not k.endswith('sha256')},indent=2))
