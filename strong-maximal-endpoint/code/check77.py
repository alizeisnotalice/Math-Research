import json,sys,re,hashlib
from pathlib import Path
from fractions import Fraction as Q
from decimal import Decimal
def dec(s):return Decimal(re.sub(r'\x60[^*]*','',s).replace('*^','E'))
p=Path(sys.argv[1]);out=Path(sys.argv[2]);base=json.loads((p/'python/allocation77.json').read_text());fine=json.loads((p/'python/allocation77_refine.json').read_text())
wb=json.loads((p/'wolfram/verification77.json').read_text());wf=json.loads((p/'wolfram/verification77_refine.json').read_text())
assert len(base['rows'])==24 and len(fine['rows'])==8 and wb['pass'] and wf['pass']
perr=[];lperr=[]
for py,wl in [(base,wb),(fine,wf)]:
 assert len(py['rows'])==len(wl['rows'])
 for a,b in zip(py['rows'],wl['rows']):
  assert (a['shape'],a['epsilon'],a['refinement'],a['assignment'])==(b['shape'],b['epsilon'],b['refinement'],b['assignment'])
  assert b['exact_pass'] and len(a['groups'])==len(b['groups'])
  for g,h in zip(a['groups'],b['groups']):
   assert h['exact'] and not h['Wolfram_LP_messages']
   perr.append(abs(dec(h['primal_60'])-dec(h['primal_120'])))
   if h['Wolfram_LP_objective'] is not None:
    assert h['Wolfram_LP_exact']
    lperr.append(abs(float(dec(h['Wolfram_LP_objective']))-g['LP_objective']))
 assert all(Q(r['allocation_lower'])<=Q(r['allocation_upper']) for r in py['rows'])
assert max(perr)<Decimal('1e-45') and len(lperr)==12 and max(lperr)<1e-7
def key(r):return r['shape'],r['epsilon'],r['refinement'],r['assignment']
bmap={key(r):r for r in base['rows']};fmap={key(r):r for r in fine['rows']}
shrinks=[]
for k,r in fmap.items():
 b=bmap[k];shrinks.append(float((Q(r['allocation_upper'])-Q(r['allocation_lower']))/(Q(b['allocation_upper'])-Q(b['allocation_lower']))))
assert max(shrinks)<.4
comparisons=[]
for shape,eps in [('notch','1/2'),('notch','1'),('checker','1'),('cross','1')]:
 g=fmap[shape,eps,2,'greedy'];m=fmap[shape,eps,2,'margin'];assert g['covered_cells']==m['covered_cells'] and g['source']==m['source']
 assert Q(g['allocation_upper'])<Q(m['allocation_lower'])
 comparisons.append(dict(shape=shape,epsilon=eps,greedy_reuse_upper=float(Q(g['allocation_upper'])/Q(g['source_modular_lower'])),margin_reuse_lower=float(Q(m['allocation_lower'])/Q(m['source_modular_upper'])),certified_conditional_cost_ratio_lower=float(Q(m['allocation_lower'])/Q(g['allocation_upper']))))
co=bmap['checker','1',1,'margin'];fi=fmap['checker','1',2,'margin']
assert [(x['anchor'],x['receiver_cells']) for x in co['groups']]==[(x['anchor'],x['receiver_cells']) for x in fi['groups']]
assert Q(fi['allocation_upper'])<Q(co['allocation_lower'])
gp=json.loads((p/'python/gap77.json').read_text());gw=json.loads((p/'wolfram/gap77.json').read_text())
assert len(gp['rows'])==56 and len(gw['rows'])==112
assert gw['symbolic_phi_second_derivative'] and gw['symbolic_perspective_Hessian_zero_determinant']
glow={(r['N'],r['scale']):r for r in gw['rows'] if r['precision']==60};gpy={(r['N'],r['scale']):r for r in gp['rows']};ge=[];gce=[]
for r in gw['rows']:
 if r['precision']!=120:continue
 k=r['N'],r['scale']
 for f in ['log_integer_optimum','log_fractional_feasible_cost','log_integrality_gap_lower','log_perspective_optimum']:
  ge.append(abs(dec(r[f])-dec(glow[k][f])));gce.append(abs(float(dec(r[f]))-gpy[k][f]))
assert max(ge)<Decimal('1e-30') and max(gce)<1e-8
r={'pass':True,'allocation_cases':32,'fixed_anchor_groups':sum(len(x['groups']) for y in [base,fine] for x in y['rows']),'Wolfram_exact_independent_LPs':12,
 'Wolfram_LP_objective_max_cross_difference':max(lperr),'primal_60_120_max_difference':str(max(perr)),
 'worst_finite_cost_interval_relative_width':max(x['relative_primal_dual_gap'] for y in [base,fine] for x in y['rows']),
 'envelope_doubling_cases':8,'maximum_interval_width_ratio_after_doubling':max(shrinks),
 'certified_grouping_comparisons':comparisons,'strict_refinement_gain_for_checker_margin':True,
 'dimension_gap_Python_cases':56,'dimension_gap_Wolfram_cases':112,'dimension_gap_precision_max_difference':str(max(ge)),'dimension_gap_cross_max_difference':max(gce),
 'scope':'Full receiver-cell certificates and exact finite fixed-anchor grouping costs at aligned refinement 2; no claim of a general C_N upper or optimal anchor grouping. The fractional gap refutes bare amplitude/threshold scaling, not legitimate weighted source-column covers.'}
r['source_sha256']={n:hashlib.sha256((Path(__file__).resolve().parent/n).read_bytes()).hexdigest() for n in ['allocate77.py','verify77.wl','gap77.py','gap77.wl','check77.py']}
out.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items() if k!='source_sha256'},indent=2))
