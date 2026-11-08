import json,sys,re,hashlib
from pathlib import Path
from fractions import Fraction as Q
from decimal import Decimal,getcontext
getcontext().prec=70
def dec(s):return Decimal(re.sub(r'\x60[^*]*','',s).replace('*^','E'))
def dq(s):q=Q(s);return Decimal(q.numerator)/Decimal(q.denominator)
p=Path(sys.argv[1]);out=Path(sys.argv[2])
g=json.loads((p/'python/grid79.json').read_text());s=json.loads((p/'python/sliver79.json').read_text());d=json.loads((p/'python/rational79.json').read_text());w=json.loads((p/'wolfram/verification79.json').read_text())
assert w['pass'] and all(w[x] for x in ['symbolic_endpoint_derivative','symbolic_slab_derivative','symbolic_smooth_primitive','symbolic_orthant_factor'])
assert len(g['cases'])==len(w['grid'])==40 and len(s['rows'])==len(w['sliver'])==66
assert all(z['exact_pass'] for part in ['grid','sliver','density','rounding'] for z in w[part])
assert sum(x['rectangles'] for x in w['grid'])==5120 and sum(x['moves'] for x in w['grid'])==1472
assert [(a['N'],a['seed']) for a in g['cases']]==[(a['N'],a['seed']) for a in w['grid']]
precision=[];maxwidth=Q(0)
for a,b in zip(s['rows'],w['sliver']):
 assert (a['N'],a['kind'],a['branch'])==(b['N'],b['kind'],b['branch'])
 assert a['direct_minimizer_success']
 nr=b['precision_rows'];assert [x['precision'] for x in nr]==[60,120]
 for field in ['log_continuous_optimum','log_grid_loss']:precision.append(abs(dec(nr[0][field])-dec(nr[1][field])))
 for z in nr:
  for field,interval in [('log_continuous_optimum','log_continuous_optimum_interval'),('log_grid_loss','log_grid_onecolumn_loss_interval')]:
   assert dq(a[interval][0])-Decimal('1e-25')<=dec(z[field])<=dq(a[interval][1])+Decimal('1e-25')
 maxwidth=max(maxwidth,Q(a['log_continuous_optimum_interval'][1])-Q(a['log_continuous_optimum_interval'][0]))
assert max(precision)<Decimal('1e-30') and float(maxwidth)<1e-10
assert max(abs(z['direct_minimizer_objective_gap']) for z in s['rows'])<1e-7
assert len(w['density'])==len(d['density_rows'])==8 and len(w['rounding'])==len(d['rounding_rows'])==14
assert [z['whole_receiver_square_success'] for z in d['density_rows']]==[z['whole_receiver_square_success'] for z in w['density']]
assert all(z['source_orthant_budget_ratio']=='1' for z in d['rounding_rows'])
assert all(z['orthant_modular_ratio']==1 for z in w['rounding'])
asym=[z for z in s['rows'] if z['N']==1048576 and z['kind']=='dimension_scale']
assert len(asym)==3 and max(abs(z['asymptotic_log_loss_residual']) for z in asym)<4e-6
from sliver79 import logi
row=next(z for z in s['rows'] if z['kind']=='fixed_dimension' and z['N']==2 and z['decimal_depth']==8)
lc=logi(Q(row['sigma'])*(1+2*Q(row['d'])))[0]-logi(Q(2))[1]-Q(row['log_candidate_cost_interval'][1])
assert lc>Q('15.8045')
r=dict(pass_all=True,exact_grid_sources=40,exact_arbitrary_rectangles=5120,exact_coordinate_moves=1472,
 certified_continuous_onecolumn_optima=66,branches={b:sum(z['branch']==b for z in s['rows']) for b in ['cap','interior','full']},
 independent_60_120_precision_rows=132,maximum_optimum_log_interval_width=str(maxwidth),maximum_optimum_log_interval_width_decimal=float(maxwidth),
 precision_difference_max=str(max(precision)),maximum_independent_direct_minimizer_gap=max(abs(z['direct_minimizer_objective_gap']) for z in s['rows']),
 F_grid_misses_slab_halo_cases=sum(z['grid_sees_no_slab_halo'] for z in s['rows']),rational_density_meshes=8,no_extra_orthant_cost_roundings=14,
 largest_dimension=1048576,largest_dimension_asymptotic_log_loss_residuals=[dict(beta=z['beta'],residual=z['asymptotic_log_loss_residual']) for z in asym],
 fixed_dimension_full_grid_fractional_loss_log_lower_certificate=str(lc),
 scope='Fixed-source step halo representation, exact continuous one-column pricing for cube/slab receiver patch, and fixed-N rational density/coverage-deficit reductions. No general sharp C_N bound or full continuous fractional-capacity optimum claimed.')
r['source_sha256']={n:hashlib.sha256((Path(__file__).resolve().parent/n).read_bytes()).hexdigest() for n in ['grid79.py','sliver79.py','rational79.py','verify79.wl','check79.py']}
r['input_sha256']={n:hashlib.sha256((p/'python'/n).read_bytes()).hexdigest() for n in ['grid79.json','sliver79.json','rational79.json']}
out.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items() if not k.endswith('sha256')},indent=2))
