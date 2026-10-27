from pathlib import Path
from fractions import Fraction
from math import factorial
import json,hashlib,ast,os,time
os.nice(max(0,10-os.getpriority(os.PRIO_PROCESS,0)))
P=Path(__file__).resolve().parent
read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
registration=read(P/'registration.json')
for row in registration['source_pins']:
 f=Path(row['path']);assert f.stat().st_size==row['bytes'] and sha(f)==row['sha256']
Q=Path(registration['source_pins'][0]['path']).parent
saved=read(Q/'next_result.json');exit=read(Q/'next_exit_receipt.json');packet=read(Q/'next_cases.json');release=read(Q/'next_parent_release.json');started=read(Q/'next_started.json')
assert exit['result_sha256']==sha(Q/'next_result.json')
assert saved['release_sha256']==sha(Q/'next_parent_release.json')==started['release_sha256']
assert saved['run_sha256']==sha(Q/'generate_next.py')==release['run_sha256']
assert saved['case_packet_sha256']==sha(Q/'next_cases.json')==release['case_packet_sha256']
assert saved['status']==exit['status']=='COMPLETE_TWO_FIXED_EXACT_DENSITY_GENERATIONS'
assert saved['actual_counts']==exit['actual_counts']
assert saved['actual_counts']=={'density_stage_checks':7,'raw_endpoint_contributions':16,'new_k_integral_kernel_native_calls':0,'old_6D_grid_replays':0}
assert saved['new_k_integral_evaluations']==0 and saved['new_moment_upper'] is False and saved['new_C6_bound'] is False
assert release['independent_conservative_literal_endpoint_upper_bound']==34 and 'not interpret' in release['bookkeeping_limitation']
for p,h in packet['source_hashes'].items():assert sha(p)==h
cones=read(packet['cone_table_path'])['cones']
start=time.monotonic();cpu=time.process_time()
import sympy as sp
assert sp.__version__=='1.14.0'
X,Y,U,V,S,T=sp.symbols('X Y U V S T')
# Only parse arithmetic, integers and the permitted formal variables; no eval/import.
def parse(text,names):
 def convert(n):
  if isinstance(n,ast.Constant) and type(n.value) is int:return sp.Integer(n.value)
  if isinstance(n,ast.Name) and n.id in names:return names[n.id]
  if isinstance(n,ast.UnaryOp) and isinstance(n.op,(ast.UAdd,ast.USub)):
   z=convert(n.operand);return z if isinstance(n.op,ast.UAdd) else -z
  if isinstance(n,ast.BinOp):
   a,b=convert(n.left),convert(n.right)
   if isinstance(n.op,ast.Add):return a+b
   if isinstance(n.op,ast.Sub):return a-b
   if isinstance(n.op,ast.Mult):return a*b
   if isinstance(n.op,ast.Div):return a/b
   if isinstance(n.op,ast.Pow):assert b.is_Integer;return a**b
  raise ValueError(ast.dump(n))
 return convert(ast.parse(text,mode='eval').body)
def vec(v):return tuple(sp.Rational(x) for x in v)
def det(a,b):return a[0]*b[1]-a[1]*b[0]
def dot(a,b):return a[0]*b[0]+a[1]*b[1]
def direction(v):return tuple(x/sum(v) for x in v)
stage_rows=[];minimum_decay=None;final_rows=[]
assert len(saved['cases'])==len(packet['cases'])==2
for case,fixed in zip(saved['cases'],packet['cases']):
 for key in ['case_id','cone_id','body_pair_zero_based','all_six_columns']:assert case[key]==fixed[key]
 cone=cones[case['cone_id']];assert cone['id']==case['cone_id'];i,j=case['body_pair_zero_based']
 source=[list(p) for p in zip(cone['nodes_by_body'][i],cone['nodes_by_body'][j])];assert source==case['all_six_columns']
 columns=[vec(z) for z in source];assert all(c>=0 for v in columns for c in v)
 active_columns=[v for v in columns if any(v)];assert case['zero_columns_removed']==6-len(active_columns)
 # Reconstruct which unchanged columns were active at each saved stage, without invoking the generator.
 independent=next((a,b) for a in range(len(active_columns)) for b in range(a+1,len(active_columns)) if det(active_columns[a],active_columns[b]))
 order=[active_columns[k] for k in independent]+[v for k,v in enumerate(active_columns) if k not in independent]
 assert len(case['stages'])==len(order)-1
 for k,stage in enumerate(case['stages'],start=2):
  assert stage['stage']==stage['active_columns']==k
  active=order[:k];rays=sorted(set(map(direction,active)),key=lambda z:z[1]);expected_cells=list(zip(rays,rays[1:]))
  assert stage['open_sector_count']==len(stage['density_cells'])==len(expected_cells)
  assert stage['nonzero_term_count']==sum(len(c['terms']) for c in stage['density_cells'])
  total=sp.S(0);decays=[];degrees=[]
  for cell,(left,right) in zip(stage['density_cells'],expected_cells):
   a,b=vec(cell['left_ray']),vec(cell['right_ray']);assert (a,b)==(left,right)
   assert all(x>=0 for ray in [a,b] for x in ray) and sum(a)==sum(b)==1 and det(a,b)>0
   lambdas=[]
   for term in cell['terms']:
    lam=vec(term['lambda']);lambdas.append(lam);d1,d2=dot(lam,a),dot(lam,b)
    assert d1>0 and d2>0 # Each saved signed summand is absolutely integrable on its own sector at S=T=0.
    decays.extend([d1,d2]);p=sp.Poly(parse(term['polynomial'],{'X':X,'Y':Y}),X,Y,domain=sp.QQ)
    assert not p.is_zero and p.total_degree()<=k-2;degrees.append(p.total_degree())
    shifted=p.as_expr().subs({X:a[0]*U+b[0]*V,Y:a[1]*U+b[1]*V},simultaneous=True)
    coefficients=sp.Poly(sp.expand(shifted),U,V,domain=sp.QQ)
    for (m,n),coefficient in coefficients.terms():
     total+=det(a,b)*coefficient*factorial(m)*factorial(n)/(d1+dot((S,T),a))**(m+1)/(d2+dot((S,T),b))**(n+1)
   assert len(lambdas)==len(set(lambdas))
  target=sp.prod(1/(1+v[0]*S+v[1]*T) for v in active)
  numerator,denominator=sp.fraction(sp.cancel(total-target));assert numerator==0
  assert sp.cancel(parse(stage['joint_Laplace_product'],{'S':S,'T':T})-target)==0
  mass=sp.cancel(total).subs({S:0,T:0});assert mass==1
  assert stage['exact_Laplace_identity'] is True and stage['normalization']==1
  minimum_decay=min(decays) if minimum_decay is None else min(minimum_decay,*decays)
  stage_rows.append({'case_id':case['case_id'],'stage':k,'sectors':len(expected_cells),'terms':stage['nonzero_term_count'],'observed_max_polynomial_degree':max(degrees),'min_lambda_dot_boundary_ray':str(min(decays)),'rational_Laplace_difference_numerator':str(numerator),'mass':str(mass),'support_exact_prefix_column_cone':True,'each_signed_term_absolutely_integrable':True})
 final_rows.append(stage_rows[-1])
assert len(stage_rows)==7
for row in registration['source_pins']:assert sha(row['path'])==row['sha256']
for p,h in packet['source_hashes'].items():assert sha(p)==h
out={'status':'PASS_SEVEN_SAVED_DENSITY_STAGE_INDEPENDENT_EXACT_RECONSTRUCTIONS','registration_sha256':sha(P/'registration.json'),'checker_sha256':sha(__file__),'source_result_sha256':sha(Q/'next_result.json'),'source_exit_sha256':sha(Q/'next_exit_receipt.json'),'sympy_version':sp.__version__,'stage_checks':stage_rows,'final_case_summaries':final_rows,'global_min_lambda_dot_boundary_ray':str(minimum_decay),'positivity_basis':'finite signed density integrable on nonnegative support; its exact Laplace transform equals the original positive exponential pushforward law, whose positivity was independently proved in prior C147 audit; signed summands never presumed positive','reported_program_counter_endpoint_differences':16,'literal_endpoint_count_claimed':False,'independent_literal_endpoint_upper_bound':34,'m2_toy_repeats':0,'new_density_generation':0,'kernel_or_special_function_values':0,'native_calls':0,'old_6D_grid_calls':0,'new_moment_upper':False,'new_C6_bound':False,'nice':os.getpriority(os.PRIO_PROCESS,0),'wall_seconds':time.monotonic()-start,'CPU_seconds':time.process_time()-cpu}
(P/'actual_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
