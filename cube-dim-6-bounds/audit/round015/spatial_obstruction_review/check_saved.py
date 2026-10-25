#!/usr/bin/env python3
"""Independent exact audit of two saved B015 outputs; no producer imports."""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
from fractions import Fraction as F
from datetime import datetime
from decimal import Decimal,localcontext
import hashlib,itertools,json,math,os,time
start=time.monotonic()
HERE=Path(__file__).resolve().parent
B=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-upper-astra/cube-dim-6-bounds/upper/round015/spatial_minimum_search')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
read=lambda p:json.loads(p.read_text())
manifest=read(B/'manifest.json');sources={}
assert manifest['status']=='B015_COMPLETED_AND_FROZEN_FIRST_64_POINT_OBSTRUCTION'
for r in manifest['files']:
 p=B/r['path'];assert sha(p)==r['sha256'];sources[str(p)]=r['sha256']
sources[str(B/'manifest.json')]=sha(B/'manifest.json')
prep=read(B/'preparation.json');release=read(B/'release.json');state=read(B/'checkpoint.json');receipt=read(B/'run_receipt.json')
assert release['explicit_parent_release'] and release['C014_final_reviewed_by_parent']
assert release['scope']=='B015-spatial-minimum-search'
for f,digest in release['prepared_file_sha256'].items():assert sha(B/f)==digest
assert sha(Path(release['C014_final_report_path']))==release['C014_final_report_sha256']
assert sha(B/'release.json')==state['release_sha256']==receipt['release_sha256']=='01d229c70878f3a3ce4d9a4a4236c0a8d358f2a661250b25cec5a552c62e6d25'
# Older numerical files are hashed only as release inputs, not parsed/rechecked.
for r in prep['inputs']:
 assert sha(Path(r['path']))==r['sha256'];sources[r['path']]=r['sha256']
assert datetime.fromisoformat(release['time']).timestamp()<state['start_epoch']
assert state['status']=='POINT_OBSTRUCTION_TO_STRICT_64_STOP_IMMEDIATELY'
assert state['native_attempts']==state['callback_requests']==len(state['native_results'])==len(state['proposal_records'])==2
assert state['counterexample_native_index']==2 and len(state['seed_runs'])==1
assert state['all_spawned_native_processes_waited'] and state['elapsed_seconds']<1200
assert receipt['driver_exit_code']==0 and receipt['native_evaluations']==2 and not receipt['second_seed_started']
assert not receipt['active_tool_sessions'] and receipt['all_native_children_waited']
T=F(6,125);threshold=T/4096;target625=T/4000;D=10**14;gap=10**6
assert (prep['T_numerator'],prep['T_denominator'])==(6,125)
def nearest(x):
 n,r=divmod(x.numerator,x.denominator)
 return n if 2*r<x.denominator else n+1 if 2*r>x.denominator else n+(n%2)
def ceil(x):return -((-x.numerator)//x.denominator)
def floor(x):return x.numerator//x.denominator
def bounds(row):
 assert set(row)=={'lower_mantissa','upper_mantissa','common_binary_exponent'}
 a,b,e=map(int,(row['lower_mantissa'],row['upper_mantissa'],row['common_binary_exponent']))
 assert a<=b and abs(e)<=10000
 return F(a)*F(2)**e,F(b)*F(2)**e
def disp(q):
 with localcontext()as c:c.prec=60;return str(Decimal(q.numerator)/Decimal(q.denominator))
oldpoint=read(Path(prep['seeds'][0]['source_path']))
seed=[nearest(F(n,oldpoint['denominator'])*D)for n in oldpoint['candidate_numerators']]
assert seed==prep['seeds'][0]['numerators']==state['native_results'][0]['numerators']
assert any(F(n,D)!=F(o,oldpoint['denominator'])for n,o in zip(seed,oldpoint['candidate_numerators']))
center=[F(n,D)for n in seed]
lower=[max(F(1,10**8),v-F(1,20))for v in center[1:]]
upper=[min(1-F(1,10**8),v+F(1,20))for v in center[1:]]
proposal=state['proposal_records'][1]
raw=[F.from_float(float.fromhex(h))for h in proposal['proposal_binary64_hex']]
free=[min(floor(b*D),max(ceil(a*D),nearest(q*D)))for a,b,q in zip(lower,upper,raw)]
prev=D
for i in range(5):free[i]=min(free[i],prev-gap);prev=free[i]
for i in range(4,-1,-1):free[i]=max(free[i],ceil(lower[i]*D),free[i+1]+gap if i<4 else gap)
initial=[D,*free]
assert initial==proposal['initial_numerators'] and proposal['initial_denominator']==D
assert max(abs(F(n,D)-q)for n,q in zip(free,raw))<=F(1,10**10)
gaps=[initial[i]-initial[i+1]for i in range(5)]+[initial[-1]]
j=max(range(6),key=lambda i:gaps[i]);assert j==4
delta=[11**i if i!=j else 0 for i in range(6)];delta[j]=-sum(delta)
moved=[a+b for a,b in zip(gaps,delta)];point=[D]+[D-sum(moved[:i])for i in range(1,6)]
assert delta==[1,11,121,1331,-162515,161051] and sum(delta)==0 and min(moved)>=gap
assert len(proposal['repairs'])==1 and proposal['repairs'][0]['k']==1
assert proposal['repairs'][0]['gap_delta_numerators']==delta
assert point==proposal['evaluated_numerators']==state['native_results'][1]['numerators'] and point!=initial

# Independent zero-form orbit, using canonical first nonzero sign directly.
tokens=iter(map(int,Path(prep['fan_path']).read_text().split()));assert next(tokens)==21
rows=set()
for _ in range(21):
 assert next(tokens)>0
 for _ in range(6):rows.add((next(tokens),tuple(next(tokens)for _ in range(6))))
assert next(tokens,None)is None and len(rows)==18
forms=set()
for _,ray in rows:
 for perm in set(itertools.permutations(ray)):
  nz=[i for i,v in enumerate(perm)if v];g=math.gcd(*perm)
  for tail in itertools.product((-1,1),repeat=len(nz)-1):
   a=list(perm)
   for i,s in zip(nz,(1,*tail)):a[i]=s*a[i]//g
   forms.add(tuple(a))
forms=sorted(forms)
assert len(forms)==9308
assert hashlib.sha256(json.dumps(forms,separators=(',',':')).encode()).hexdigest()==prep['zero_forms_sha256_canonical_json']
assert [list(a)for a in forms if sum(v*n for v,n in zip(a,initial))==0]==[[1,-2,1,0,0,0]]
vs=[(1,*s)for s in itertools.product((-1,1),repeat=5)]
def distance(x):
 norm=sum(v*v for v in x);vals=[]
 for i,v in enumerate(vs):
  for w in vs[:i]:
   av=sum(a*b for a,b in zip(v,x));aw=sum(a*b for a,b in zip(w,x));c=sum(a*b for a,b in zip(v,w))
   vals.append(norm-(6*av*av-2*c*av*aw+6*aw*aw)/(36-c*c))
 assert len(vals)==496 and min(vals)>0
 return min(vals)
results=[]
for index,(rec,pr)in enumerate(zip(state['native_results'],state['proposal_records']),1):
 ns,den=rec['numerators'],rec['denominator'];assert den==D and ns[0]==D
 assert all(ns[i]-ns[i+1]>=gap for i in range(5)) and ns[-1]>=gap
 assert all(a<=F(n,den)<=b for a,b,n in zip(lower,upper,ns[1:]))
 assert pr['evaluated_numerators']==ns and pr['evaluated_denominator']==den
 assert rec['T']=='6/125' and rec['index']==index and rec['returncode']==0 and rec['status']=='POINT_INTERVALS_VALIDATED'
 assert rec['command']==[prep['native_binary_path'],'full',prep['fan_path'],'6','125',str(den),*map(str,ns)]
 for field,path in [('binary_sha256',prep['native_binary_path']),('source_sha256',prep['native_source_path']),('fan_sha256',prep['fan_path'])]:assert rec[field]==sha(Path(path))
 assert rec['start_epoch']>=state['start_epoch']
 assert sha(Path(rec['stdout_path']))==rec['stdout_sha256'] and sha(Path(rec['stderr_path']))==rec['stderr_sha256']
 assert Path(rec['stderr_path']).read_text()==''
 n=read(Path(rec['stdout_path'])) # Exactly two parsed native output records.
 assert n['status']=='generic_nonzero_supported' and n['precision_bits']==512
 assert n['full_group_terms']==967680 and n['computed_sign_pairs']==483840 and n['zero_node_terms']==n['confluent_terms']==0
 assert len(n['gradient'])==6 and n['unique_node_values']==16931
 dots=[sum(a*b for a,b in zip(f,ns))for f in forms];assert all(dots)
 radius=min(F(abs(d),2*den*sum(map(abs,f)))for d,f in zip(dots,forms))
 assert max(sum(abs(a)*max(map(abs,ns))for a in r)for _,r in rows)<2**63-1
 assert den*max(q for q,_ in rows)<2**63-1
 K=bounds(n['K']);p=bounds(n['p_T']);gr=[bounds(v)for v in n['gradient']]
 assert K[0]>0 and p[0]>0
 assert tuple(F(v)for v in rec['check']['K_interval'])==K
 for a,b in [K,p,*gr]:assert b-a<max(F(1),abs((a+b)/2))*F(2)**-256
 lo=5*K[0]+T*p[0];hi=5*K[1]+T*p[1]
 for xi,(a,b)in zip([F(v,den)for v in ns],gr):lo+=min(xi*a,xi*b);hi+=max(xi*a,xi*b)
 assert lo<=0<=hi
 a,b=bounds(n['Euler_x_grad_plus5K_plusTp']);assert a<=0<=b
 qp=((T*p[0]-K[1])/T**2,(T*p[1]-K[0])/T**2);assert qp[1]<0
 results.append({'index':index,'numerators':ns,'denominator':den,'K_exact_dyadic_bounds':list(map(str,K)),'strict_64_obstruction':K[1]<threshold,'strict_625_obstruction':K[1]<target625,'threshold64_minus_K_upper_exact':str(threshold-K[1]),'K_midpoint_display_only':disp(sum(K)/2),'model_cost_lower_exact':str(T/(64*K[1])),'model_cost_lower_display_only':disp(T/(64*K[1])),'dq_dT_exact_bounds':list(map(str,qp)),'dq_dT_strictly_negative':True,'Euler_exact_interval_contains_zero':True,'S2_distance_squared_exact':str(distance([F(a,den)for a in ns])),'generic_linf_radius_exact':str(radius),'raw_sha256':rec['stdout_sha256']})
assert not results[0]['strict_64_obstruction'] and not results[0]['strict_625_obstruction']
assert results[1]['strict_64_obstruction'] and results[1]['strict_625_obstruction']
K2=tuple(F(v)for v in results[1]['K_exact_dyadic_bounds'])
assert F(5574745829,500000000000000)<K2[0]<=K2[1]<F(111494916581,10**16)<threshold
assert T/(64*K2[1])>F(168169,2500)
assert results[1]['S2_distance_squared_exact']=='6020805515960314567115330493/40000000000000000000000000000'
assert results[1]['generic_linf_radius_exact']=='1/80000000000000'
assert state['native_results'][1]['start_epoch']>=state['native_results'][0]['start_epoch']+state['native_results'][0]['elapsed_seconds']
out={'status':'PASS_FIXED_BODY32_T48_ESSENTIAL_MINIMUM_OBSTRUCTION','parsed_saved_native_outputs':2,'new_kernel_calls':0,'binary_launched':False,'optimizer_launched':False,'release_pins_match':True,'first_strict_64_obstruction_stopped_at_second_call':True,'second_seed_not_started':True,'seed_rounding_is_new_exact_input':True,'second_candidate_repair_reconstructed_from_saved_binary64_proposal':True,'original_zero_candidate_never_evaluated':True,'repair_k':1,'gap_delta_numerators':delta,'zero_form_orbit_count':len(forms),'results':results,'strict_uniform_model_cost_lower':'168169/2500','fixed_T':'6/125','other_T_or_mixtures_excluded':False,'new_C6_bound':False,'source_hash_checks':sources,'pid':os.getpid(),'elapsed_seconds':time.monotonic()-start,'script_sha256':sha(Path(__file__)),'exit_code':0}
(HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items()if k not in('source_hash_checks','results')},indent=2))
print(json.dumps([{k:v for k,v in r.items()if k in('index','K_midpoint_display_only','model_cost_lower_display_only','S2_distance_squared_exact','generic_linf_radius_exact','strict_64_obstruction')}for r in results],indent=2))
