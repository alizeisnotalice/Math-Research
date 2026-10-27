import sys,os,time,json,hashlib,itertools,math,gzip
sys.set_int_max_str_digits(0);sys.dont_write_bytecode=True
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
from fractions import Fraction as F
from collections import Counter,defaultdict
P=Path(__file__).resolve().parent
rd=lambda q:json.loads(Path(q).read_text());sha=lambda q:hashlib.sha256(Path(q).read_bytes()).hexdigest()
checks=0;started=time.monotonic()
def ck(ok,label):
 global checks
 if not ok:raise AssertionError(label)
 checks+=1
def checkpoint(phase,**data):
 p=P/'verification_checkpoint.json';t=p.with_name(p.name+'.pending');t.write_text(json.dumps({'phase':phase,'pid':os.getpid(),'assertions':checks,'elapsed_seconds':time.monotonic()-started,**data},indent=2)+'\n');os.replace(t,p)
def dot(x,y):return sum(a*b for a,b in zip(x,y))
def determinant(rows):
 n=len(rows);answer=F(0)
 for perm in itertools.permutations(range(n)):
  sign=-1 if sum(perm[i]>perm[j]for i in range(n)for j in range(i+1,n))%2 else 1
  answer+=sign*math.prod(rows[i][perm[i]]for i in range(n))
 return answer
def local_series(nodes,b,m):
 # [u^r] product_(c!=b) (b-c+u)^(-multiplicity(c)), r<=m-1.
 order=m-1;coeff=[F(1)]+[F(0)]*order
 for c,mc in nodes.items():
  if c==b:continue
  delta=b-c;series=[F((-1)**k*math.comb(mc+k-1,k))*delta**(-mc-k)for k in range(order+1)]
  coeff=[sum(coeff[j]*series[k-j]for j in range(k+1))for k in range(order+1)]
 return coeff
def residue_channels(tuple_nodes):
 multiplicities=Counter(tuple_nodes);gamma=F(0);negative=[F(0)]*3
 for b,m in multiplicities.items():
  coeff=local_series(multiplicities,b,m)
  if b:gamma+=(1 if b>0 else-1)*coeff[m-1]
  else:
   for index,k in enumerate((1,3,5)):
    if k<m:negative[index]=F((-1)**((k-1)//2),k)*coeff[m-1-k]
 return(gamma,*negative)
reg=rd(P/'registration.json');ck(sha(__file__)==reg['checker_sha256'],'registered verifier hash');ck(os.getpriority(os.PRIO_PROCESS,0)>=10,'nice>=10')
A=Path(reg['author_directory']);inp=rd(A/'inputs.json');result=rd(A/'results.json');qual=rd(A/'pair_qualification.json')
for row in rd(A/'source_manifest.json')['sources']:ck(sha(row['path'])==row['sha256'],'original archived source hash')
checkpoint('RANK_AND_ALL496_GRAM_QUALIFICATION')
x=tuple(map(F,inp['point']));V=[(1,)+tail for tail in itertools.product((-1,1),repeat=5)];ck(x==(F(1),F(1),F(1),F(3,8),F(1,4),F(1,8)),'sole registered new point')
ck(qual['body_columns']==[list(v)for v in V]and qual['point']==inp['point'],'actual frame/point binding')
basis=[tuple(map(F,v))for v in qual['basis_witness']];ck(len(basis)==6 and all(tuple(v)in V for v in basis)and abs(determinant(basis))==32,'independent Leibniz rank witness')
ck([[sum(v[i]*v[j]for v in V)for j in range(6)]for i in range(6)]==[[32*int(i==j)for j in range(6)]for i in range(6)],'actual body32 tight frame')
# Signed-coordinate flips and adjacent swaps generate the full signed permutation group.
lines={tuple(v)for v in V}|{tuple(-a for a in v)for v in V}
for axis in range(6):ck(all(tuple(-a if j==axis else a for j,a in enumerate(v))in lines for v in V),'frame sign generator')
for axis in range(5):
 ck(all(tuple(v[j+1]if j==axis else v[j-1]if j==axis+1 else v[j]for j in range(6))in lines for v in V),'frame permutation generator')
minimum=None;pairno=0
for i in range(32):
 for j in range(i+1,32):
  gram=[[dot(v,w)for w in (V[i],V[j],x)]for v in(V[i],V[j],x)];D=determinant([row[:2]for row in gram[:2]]);ck(D>0,'body pair independence')
  rho=determinant(gram)/D;row=qual['pairs'][pairno]
  ck(rho>0 and row['i']==i and row['j']==j and F(row['Gram_determinant'])==D and F(row['distance_squared'])==rho,'all496 independent Gram determinant-ratio residuals')
  minimum=rho if minimum is None else min(minimum,rho);pairno+=1
ck(pairno==qual['pair_count']==496 and minimum==F(qual['minimum_squared_distance'])==F(result['minimum_pair_squared_distance'])==F(1,32)and qual['off_true_S2']is True,'true offS2 qualification')
checkpoint('EXACT_FAN_BINDING')
fan=rd(inp['fan_path']);audit=rd(inp['fan_audit_path']);ck(audit['status']=='PASS'and audit['input_sha256']==sha(inp['fan_path'])and audit['full_cells']==len(fan['cells'])==21 and F(audit['simplex_section_volume_sum'])==F(1,120),'audited complete chamber coverage binding')
cones=[]
for cell in fan['cells']:
 rays=[tuple(map(F,r))for r in cell['rays_a0_equals_one']];rates=list(map(F,cell['rates']));det=F(cell['determinant']);ck(abs(determinant(rays))==det>0,'independent ray determinant')
 ck(all(q>0 and q==sum(abs(dot(v,r))for v in V)for q,r in zip(rates,rays)),'positive original psi rates')
 ck(all(min(dot(v,r)for r in rays)>=0 or max(dot(v,r)for r in rays)<=0 for v in V),'body weak signs')
 cones.append((rays,rates,det/math.prod(rates)))
R0=F(inp['R0']);ck(46080*sum(a for _,_,a in cones)==R0==F(fan['full_exp_minus_symbol_integral'])==F(audit['full_exp_minus_symbol_integral']),'complete positive fan mass')
checkpoint('DISTINCT_ORBIT_AND_RECONSTRUCT_WEIGHTS')
ints=tuple(inp['point_integer_numerators']);den=inp['point_denominator'];ck(tuple(F(v,den)for v in ints)==x,'integer point encoding')
perms=sorted(set(itertools.permutations(ints)));stabilizer=math.prod(math.factorial(m)for m in Counter(ints).values());orbit=sorted({tuple(s*v for s,v in zip(signs,p))for p in perms for signs in itertools.product((-1,1),repeat=6)})
ck(len(perms)==120 and len(orbit)==7680 and stabilizer==6 and len(orbit)*stabilizer==46080,'independent distinct-coordinate orbit construction')
weights=defaultdict(F)
for oi,y_int in enumerate(orbit):
 y=tuple(F(v,den)for v in y_int)
 for rays,rates,alpha in cones:
  node_tuple=tuple(sorted(dot(y,r)/q for r,q in zip(rays,rates)));weights[node_tuple]+=stabilizer*alpha
 if(oi+1)%1000==0:checkpoint('DISTINCT_ORBIT_AND_RECONSTRUCT_WEIGHTS',orbit_done=oi+1,merged_so_far=len(weights))
ck(result['orbit_size']==7680 and result['stabilizer']==6 and result['compressed_orbit_cone_terms']==161280 and result['uncompressed_group_cone_terms']==967680,'saved complete group counts')
ck(len(weights)==result['merged_exact_node_tuples']==160298 and len(weights)<=161280 and sum(weights.values())==R0,'all compressed tuple weights/mass')
checkpoint('INDEPENDENT_LOCAL_POLE_COEFFICIENTS',merged_tuples=len(weights))
ck(sha(A/'exact_tuple_certificate.json.gz')==result['tuple_certificate_sha256'],'saved tuple certificate hash chain')
with gzip.open(A/'exact_tuple_certificate.json.gz','rt')as f:certificate=json.load(f)
ck(certificate['channel_order']==['Gamma','epsilon^-1','epsilon^-3','epsilon^-5']and len(certificate['tuples'])==len(weights),'complete saved tuple channel schema')
total=[F(0)]*4;zeros=repeated=0
for ni,(nodes,row)in enumerate(zip(sorted(weights),certificate['tuples'])):
 ck(tuple(map(F,row['nodes']))==nodes and F(row['weight'])==weights[nodes]>0,'independent full tuple/weight equality')
 ch=residue_channels(nodes);ck(tuple(map(F,row['channels']))==ch,'independent local negative-binomial pole coefficients')
 for j in range(4):total[j]+=weights[nodes]*ch[j]
 zeros+=0 in nodes;repeated+=len(set(nodes))<6
 if(ni+1)%5000==0:checkpoint('INDEPENDENT_LOCAL_POLE_COEFFICIENTS',tuples_done=ni+1,total_tuples=len(weights))
ck(tuple(map(F,result['complete_channels']))==tuple(total)and total[0]>0 and total[1:]==[0,0,0],'all exact Gamma/principal totals and complete cancellation')
ck(zeros==result['zero_node_tuples']==22404 and repeated==result['repeated_node_tuples']==38266,'all zero/repeated tuples retained')
checkpoint('SINGLE_PI768_NORMALIZATION_AND_EXACT_COMPARISONS')
from flint import arb,fmpq,ctx
ctx.prec=768;q=total[0]/128;ck(F(result['Q_over_pi5'])==q,'full G Gamma/128 coefficient')
v=arb(fmpq(q.numerator,q.denominator))/arb.pi()**5;ck(v.is_finite(),'finite normalization');end=[]
for z in(v.lower(),v.upper()):m,e=z.man_exp();end.append(F(int(m))*F(2)**int(e))
gl,gu=map(F,result['G_interval']);ck(0<gl<=end[0]<=end[1]<=gu,'fresh768 complete G inside saved512/384-grid interval')
gstar=F(1,63500);g0=F(1,65000);ck(result['target']=='1/63500'and result['G_upper_below_floor_target']==(gu<gstar)and result['G_lower_at_least_target']==(gl>=gstar)and result['essential_floor_target_excluded_by_offS2_continuity']==(gu<gstar),'accurate gstar point/essential-floor flags')
M6=F(rd(inp['round128_certificate_path'])['M_power6_upper']);K=F(8,9)**6*24576**5;gap=M6-K*gu**6;budget=F(9,8)**6*M6/(384**5*gu**6)
ck(F(result['saved128_M_power6_upper'])==M6 and F(result['saved128_M6_minus_threshold_at_GpointU'])==gap and F(result['saved128_certificate_budget_at_GpointU_power5'])==budget and(gap>=0)==(budget>=64**5)==result['saved128_upper_certificate_cannot_reach_below64_for_any_valid_g'],'accurate complete current-Mupper conditional budget')
ck(F(result['target_minus_G_upper_exact'])==gstar-gu and F(result['G_lower_minus_target_exact'])==gl-gstar,'exact point comparison residuals')
ck(result['full_positive_time_G']is True and result['numerical_time_truncation']is False and result['unpaid_time_tail']is False and result['actual_M_or_model_cost_or_C6_lower_bound']is False and result['global_G_floor_proved']is False and result['new_C6_bound']is False,'scalar/fulltime and limited inference scope')
out={'status':'PASS_ONE_COMPLETE_INDEPENDENT_C135_NEW_POINT','assertions':checks,'pair_count':496,'minimum_pair_distance_squared':str(minimum),'orbit_size':7680,'stabilizer':6,'raw_cone_orbit_terms':161280,'all_merged_tuples_checked':len(weights),'zero_tuples':zeros,'repeated_tuples':repeated,'method':'independent local negative-binomial pole coefficient expansion, not producer triangular DD','complete_channels':list(map(str,total)),'fresh_G_interval_768':list(map(str,end)),'saved_G_interval':result['G_interval'],'gstar_point_lower_pass':gl>=gstar,'gstar_essential_floor_excluded':gu<gstar,'g0_point_lower_pass':gl>=g0,'g0_essential_floor_excluded':gu<g0,'g0_lower_minus_target_exact':str(gl-g0),'g0_target_minus_upper_exact':str(g0-gu),'saved_Mupper_certificate_budget_obstructed':gap>=0,'actual_M_or_Gglobal_or_model_cost_or_C6_lower_bound':False,'new_points':1,'old_anchor_baseline_calls':0,'Arb_pi_normalizations':1,'precision_bits':768}
(P/'review_results.json').write_text(json.dumps(out,indent=2)+'\n');checkpoint('COMPLETE',tuples_done=len(weights));print(json.dumps({k:v for k,v in out.items()if k not in['complete_channels','fresh_G_interval_768','saved_G_interval','g0_lower_minus_target_exact','g0_target_minus_upper_exact']}),flush=True)
