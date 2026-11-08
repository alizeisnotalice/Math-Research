import json,sys,re,hashlib,math
from pathlib import Path
from fractions import Fraction as Q
from decimal import Decimal,getcontext
getcontext().prec=70
def dec(s):return Decimal(re.sub(r'\x60[^*]*','',s).replace('*^','E'))
def dq(s):q=Q(s);return Decimal(q.numerator)/Decimal(q.denominator)
p=Path(sys.argv[1]);out=Path(sys.argv[2]);a=json.loads((p/'python/frontier80.json').read_text());w=json.loads((p/'wolfram/verification80.json').read_text())
assert w['pass'] and all(w['symbolic'].values())
assert len(a['coverage'])==len(w['coverage'])==13 and sum(z['certified_nodes'] for z in a['coverage'])==8204
for x,y in zip(a['coverage'],w['coverage']):
 assert y['exact_pass'] and x['M']==y['M'] and x['certified_nodes']==y['nodes']==x['M']+1
 assert Q(x['minimum_Bernstein_coefficient'])>0 and x['maximum_subdivision_depth']==0
 assert Q(x['epsilon'])==Q(1,144*x['M']**2) and Q(x['budget_lower'])==Q(x['M'],6144)
 assert Q(x['unnormalized_margin_lower'])==Q(229,768)*Q(x['epsilon'])
 assert Q(x['discarded_layer_volume'])==Q(x['epsilon'])/24
jf=a['source_modular_interval'];jrows=w['source_modular_precision'];assert [x['precision'] for x in jrows]==[60,120]
jerr=abs(dec(jrows[0]['source_modular'])-dec(jrows[1]['source_modular']));assert jerr<Decimal('1e-50')
assert all(dq(jf[0])<=dec(r['source_modular'])<=dq(jf[1]) for r in jrows)
assert len(a['rational_step'])==len(w['rational_step'])==16 and sum(len(z['samples']) for z in a['rational_step'])==528
for x,y in zip(a['rational_step'],w['rational_step']):
 m=x['mesh'];assert y['exact_pass'] and y['mesh']==m and y['samples']==len(x['samples'])==33
 assert Q(x['source_modular_interval'][0])==Q(jf[0])-Q(12,m) and Q(x['source_modular_interval'][1])==Q(jf[1])
 assert Q(x['normalized_budget_lower'])==Q(x['budget_lower'])/Q(jf[1])
 assert all(Q(s['exact_maximal_value'])>=1+Q(68,15*m) for s in x['samples'])
last=a['rational_step'][-1];assert Q(last['normalized_budget_lower'])>Q('12.7389') and Q(last['actual_halo_ratio_upper'])<Q('11.7252')
assert len(a['high_dimension'])==len(w['high_dimension'])==14
prec=[];cross=[]
for x,y in zip(a['high_dimension'],w['high_dimension']):
 assert x['N']==y['N'] and [z['precision'] for z in y['rows']]==[60,120]
 for f in ['log_section_constant','log_truncated_cover_budget_lower']:
  prec.append(abs(dec(y['rows'][0][f])-dec(y['rows'][1][f])));cross.append(abs(float(dec(y['rows'][1][f]))-x[f]))
assert max(prec)<Decimal('1e-25') and max(cross)<1e-8
r=dict(pass_all=True,symbolic_continuum_gates=w['symbolic'],rational_anchor_cover_families=13,certified_Bernstein_cells=8204,maximum_Bernstein_subdivision_depth=0,
 rational_step_meshes=16,exact_step_source_maximal_tests=528,largest_implicit_step_mesh=last['mesh'],
 source_modular_interval=jf,source_modular_precision_difference=str(jerr),
 largest_step_normalized_column_budget_lower=last['normalized_budget_lower'],largest_step_actual_halo_ratio_upper=last['actual_halo_ratio_upper'],
 fixed_dimension_parameters=14,high_dimension_precision_evaluations=28,high_dimension_precision_difference_max=str(max(prec)),high_dimension_cross_max_difference=max(cross),
 scope='Counterexample to uniformly bounded full strict-halo fixed-corner source-column budgets, including vanishing coverage deficit. Matching epsilon^-1/2 truncated-cover rate. This is not divergence of C_N; the actual weak ratios of the rational-step counterfamily remain bounded.')
r['source_sha256']={n:hashlib.sha256((Path(__file__).resolve().parent/n).read_bytes()).hexdigest() for n in ['frontier80.py','verify80.wl','check80.py','sliver79.py']}
r['input_sha256']={'frontier80.json':hashlib.sha256((p/'python/frontier80.json').read_bytes()).hexdigest()}
out.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items() if not k.endswith('sha256')},indent=2))
