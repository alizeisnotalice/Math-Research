"""Independent saved C160B audit. Standard library integers only; no producer import/root calls."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations_with_replacement,islice
from math import comb
import gzip,hashlib,json,os,subprocess,sys,time,traceback
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
P=Path(__file__).resolve().parent
TARGET=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-brascamp108/cube-dim-6-bounds/audit/upper160_saved')
S=1<<192
PASS='PASS_COMPLETE160B_SAVED_ROOTS_BOTH_LAWS_AND_NORMALIZED_M'
ENV=['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def canon(v):return json.dumps(v,sort_keys=True,separators=(',',':')).encode()
def need(b,m):
 if not b:raise ArithmeticError(m)
def save(n,v):
 with(P/n).open('x')as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
def checkpoint(n,v):
 tmp=P/(n+'.tmp')
 with tmp.open('w')as f:json.dump(v,f);f.write('\n');f.flush();os.fsync(f.fileno())
 os.replace(tmp,P/n)
def release_requirements():
 reg=read(P/'registration.json')
 return {'round':'160B_SAVED_INDEPENDENT','authorized':True,'max_science_runs':1,'registration_sha256':sha(P/'registration.json'),'inputs_sha256':reg['inputs_sha256'],'source_manifest_sha256':reg['source_manifest_sha256'],'code_sha256':reg['code_sha256'],'limits':reg['limits'],'new_C6_bound':False}
def gate():
 need(P==TARGET,'TERRA_ONLY_ROOT');reg=read(P/'registration.json');release=read(P/'parent_release.json')
 for k,v in release_requirements().items():need(release.get(k)==v,'DUAL_RELEASE:'+k)
 need(release.get('main_approved')is True and release.get('independent_review_approved')is True,'DUAL_CODE_AND_CONCRETE_GATE')
 need(sha(P/'inputs.json')==reg['inputs_sha256']and sha(P/'source_manifest.json')==reg['source_manifest_sha256'],'PREPARED_BYTES')
 for n,h in reg['code_sha256'].items():need(sha(P/n)==h,'CHECKER_CODE_IDENTITY')
 need(all(os.environ.get(k)=='1'for k in ENV)and os.getpriority(os.PRIO_PROCESS,0)>=10,'RESOURCE_ENV_NICE')
 need(str(Path(sys.executable).resolve())==read(P/'inputs.json')['python_executable'],'REGISTERED_PYTHON')
 return read(P/'inputs.json')
def pins():
 for row in read(P/'source_manifest.json')['sources']:need(sha(row['path'])==row['sha256']and Path(row['path']).stat().st_size==row['bytes'],'SOURCE_PRE_POST:'+row['path'])
def model_check(inp):
 src=Path(inp['producer_root']);d=src/'phaseB';cert=read(d/'candidate_certificate.json');tab=read(d/'integer_tables.json');candidate=read(d/'candidate.json');co=read(d/'coefficients.json');m=read(src/'inputs.json');model=m['model_input'];q=read(src/'quantization_certificate.json')
 roles={r['role']:Path(r['path'])for r in read(src/'STATIC_SOURCE_PINS.json')['sources']};old=read(roles['128_input']);prior=read(inp['source128_independent_review'])
 need(prior['status']=='PASS_COMPLETE_INDEPENDENT_SAVED_JOINT_CANDIDATE'and sha(roles['128_certificate'])in prior['source_hashes'].values(),'PRIOR128_ACTUAL_ACCEPTED_CERT_BYTES_NO_REPLAY')
 ex=read(d/'exit_receipt.json');term=read(d/'operator_terminal_receipt.json');started=read(d/'started.json');release=read(d/'parent_release.json');reg=read(d/'registration.json')
 need(cert['status']==ex['status']=='COMPLETE_ONE_NEW_NORMALIZED_FULL_H_CERTIFICATE'and sha(d/'candidate_certificate.json')==inp['producer_certificate_sha256']==ex['certificate_sha256'],'ACTUAL160B_TERMINAL_CERT')
 need(term['exit_code']==0 and term['automatic_retries']==0 and term['pid']==started['pid']and started['workers']==ex['workers_completed']==15,'ACTUAL15_WORKER_ORIGINAL_TERMINAL')
 need(sha(d/'parent_release.json')==term['release_sha256']==started['release_sha256']==cert['release_sha256']and sha(src/'run_phaseB.py')==cert['run_sha256']==release['run_sha256'],'PRODUCER_RELEASE_AND_SOURCE')
 for n,k in [('operator_science_stdout.txt','stdout_sha256'),('operator_science_stderr.txt','stderr_sha256')]:need(sha(d/n)==term[k],'PRODUCER_STDIO')
 need((d/'operator_science_stderr.txt').read_bytes()==b'','PRODUCER_NO_STDERR')
 for k,n in [('inputs_sha256','inputs.json'),('quantization_certificate_sha256','quantization_certificate.json')]:need(sha(src/n)==reg[k]==release[k],'ACTUAL_MODEL_INPUT_CHAIN')
 need(sha(d/'registration.json')==release['registration_sha256']and release['main_approved']is True and release['independent_review_approved']is True,'PRODUCER_SEPARATE_RELEASE')
 for k,n in [('candidate_sha256','candidate.json'),('coefficients_sha256','coefficients.json'),('selection_sha256','selection.json'),('jobs_sha256','jobs.json')]:need(sha(d/n)==reg[k],'FIXED_PHASEB_ARTIFACT')
 A=read(src/'phaseA/result.json');selection=read(d/'selection.json')
 need(A['candidate']==candidate==cert['candidate']and sha(src/'phaseA/result.json')==reg['phaseA_result_sha256']==selection['phaseA_result_sha256'],'ONE_ACTUAL_SELECTED_CANDIDATE')
 need(A['candidate_canonical_sha256']==selection['candidate_canonical_sha256']==hashlib.sha256(canon(candidate)).hexdigest()and selection['main_selected']is True,'CANONICAL_CANDIDATE_SELECTION')
 aexit=read(src/'phaseA/exit_receipt.json');aterm=read(src/'phaseA/operator_terminal_receipt.json');astart=read(src/'phaseA/started.json')
 need(aexit['result_sha256']==sha(src/'phaseA/result.json')and aexit['status']==A['status']and aterm['exit_code']==0 and aterm['pid']==astart['pid']and aterm['automatic_retries']==0,'ACTUAL_PHASEA_TERMINAL_NOT_FLOAT_GATE')
 need(A['optimizer_invocations']==1 and A['restarts']==0 and 0<A['objective_gradient_calls']<=150 and A['prefix_states']==735471 and A['tail_rows']==16,'SINGLE_AUTHOR_SELECTOR_NOT_REPLAYED')
 # Actual156 full accepted upper is the ONLY cap substitution, not a display or lower estimate.
 accept=read(inp['source156_acceptance']);audit156=read(TARGET.parent/'upper156_saved/certificate.json');r156=read(roles['156_result']);B=Q(r156['budgets']['finite_best']['valid_true_beta2_upper'])
 augmented=release['source156_saved_audit']
 need(augmented==cert['source156_saved_audit']and set(augmented)==set(accept)and all(augmented[k]==v for k,v in accept.items()if k!='evidence_files')and accept['main_accepted']is True and accept['independent_saved_review_passed']is True,'ACTUAL156_ORIGINAL_ACCEPTANCE_FIELDS')
 need(augmented['evidence_files'][:-1]==accept['evidence_files']and augmented['evidence_files'][-1]['sha256']=='18eca9a3a6b33c987952e42244221286efb143c038c54335a5bb02125e9713e4','ORIGINAL_EVIDENCE_WITH_ONE_ACTUAL_SEAL_ADDED')
 need(accept['result_sha256']==audit156['producer_result_sha256']==sha(roles['156_result'])==q['source_result_sha256']==cert['source156_result_sha256'],'FULL156_SOURCE_LINK')
 need(audit156['status']=='PASS_COMPLETE_SAVED156_INDEPENDENT_SCALAR768_AND_RATIONAL_CAPS'and Q(audit156['certified_true_beta2_upper'])==Q(accept['accepted_true_beta2_upper_exact'])==B==Q(r156['budgets']['finite_best']['budget_expression_interval'][1]),'EXACT_COMPLETE156_ACCEPTED_VALUE')
 need(audit156['old130_lower_le_new_BU']is True and audit156['BL_is_actual_mu2_lower']is False and Q(audit156['old130_true_complete_moment_interval'][0])<=B,'ALREADY_CERTIFIED130_LOWER_CONSISTENCY_NO_REPLAY')
 need(q['source_field_selector']==['budgets','finite_best','valid_true_beta2_upper']and Q(q['complete_source_upper_exact'])==B and q['from_display_decimal']is False,'NO_DISPLAY_OR_PARTIAL_CAP')
 D,N=int(q['D']),int(q['N']);need(D==1<<96 and 0<B.numerator*D<=N*B.denominator<B.numerator*D+B.denominator,'EXISTING_UPWARD96_WITNESS_CONSUMER_ONLY')
 cap=Q(N,D);need(Q(q['betaD_exact'])==Q(m['new_beta2_cap_exact'])==Q(cert['new_beta2_cap_exact'])==cap and m['changed_flat_index']==32,'UNCHANGED_CONSERVATIVE_NEW_BETA2')
 for key in ['model','h_values','beta0','p_intervals','Phi_table','threshold_M_power6_exact','prefix_Kmax','prefix_state_count','tail_K_values']:need(model[key]==old[key],'SAME_ACCEPTED128_MODEL:'+key)
 need(len(model['beta_by_h'])==8 and all(len(t)==6 for t in model['beta_by_h'])and model['h_values'][5]=='1/128','FULL48_FEATURE_INDEX')
 for j in range(8):
  for k in range(6):
   b=model['beta_by_h'][j][k]
   need((Q(b)==cap<Q(old['beta_by_h'][j][k]))if(j,k)==(5,1)else b==old['beta_by_h'][j][k],'ONLY_GENUINE_BETA2_CHANGED')
 phi=model['Phi_table'];need(len(phi)==33 and all(len(row)==7 and row[0]==1 and row[1]==n and row[2]==comb(n,2)and all(type(v)is int and v>=0 for v in row)for n,row in enumerate(phi)),'COMPLETE_REAL_FEATURE_TABLE')
 need(all(phi[n+1][k]>=phi[n][k]for n in range(32)for k in range(7)),'POSITIVE_MONOTONE_FEATURES')
 grid=int(candidate['denominator']);nums=candidate['numerators'];need(grid==1<<44 and len(nums)==49 and all(type(v)is int and v>=0 for v in nums)and nums[0]>0 and sum(nums)==grid,'EXACT_NEW49_SIMPLEX')
 caps=[Q(model['beta0'])]+[Q(v)for row in model['beta_by_h']for v in row];need(min(caps)>0,'POSITIVE_UPPER_CAPS');a=[Q(v,grid)/b for v,b in zip(nums,caps)]
 need(sum(v*b for v,b in zip(a,caps))==1 and co['sum_a_cap_exact']==cert['normalizing_cap_cost_exact']=='1'and cert['true_cost_upper_at_most_one']is True,'EXACT_NORMALIZED_COST_ONCE')
 need(list(map(Q,co['a_exact']))==a and list(map(Q,co['cap_vector_exact']))==caps and [Q(tab['a0'])]+[Q(v)for row in tab['a_by_h']for v in row]==a,'ACTUAL49_RAW_RATIONAL_COEFFICIENTS')
 need(co['raw_a_quantized']is False and cert['raw_a_quantized']is False and hashlib.sha256(canon(list(map(str,a)))).hexdigest()==co['canonical_a_sha256']==cert['canonical_a_sha256'],'NO_RAW_A_QUANTIZATION')
 def outward(v):return [v.numerator*S//v.denominator,-((-v.numerator*S)//v.denominator)]
 common=outward(a[0]);groups=[[outward(sum(a[1+6*j+k]*phi[n][k+1]for k in range(6)))for n in range(33)]for j in range(8)]
 need(common[0]>0 and common==list(map(int,tab['common_W_units']))and groups==[[list(map(int,z))for z in row]for row in tab['group_W_units']],'INDEPENDENT192_W_TABLES')
 ps=[];qs=[]
 for side in range(2):
  vals=[Q(row[side])*S for row in model['p_intervals']];need(len(vals)==8 and all(v.denominator==1 for v in vals),'EXACT192_CDF_GRID')
  p=list(map(int,vals));qv=[p[0]]+[p[j]-p[j-1]for j in range(1,8)]+[S-p[-1]];need(len(qv)==9 and min(qv)>0 and sum(qv)==S,'TWO_COMPLETE_NESTED_COUNT_LAWS');ps.append(p);qs.append(qv)
 need(all(l<=u for l,u in zip(ps[0],ps[1]))and ps==[list(map(int,row))for row in tab['p_units']]and qs==[list(map(int,row))for row in tab['category_units']],'WHOLE_CDF_COUPLING_NOT_INDEPENDENT_CATEGORIES')
 need(tab['scale_power2']==cert['arithmetic_bits']==192 and tab['probability_denominator_power2']==cert['probability_denominator_power2']==6144 and tab['weighted_sum_denominator_power2']==cert['weighted_sum_denominator_power2']==6336,'EXACT_DENOMINATORS')
 jobs=read(d/'jobs.json');need(len(jobs)==len(cert['chunks'])==15 and jobs[0]['start']==0 and jobs[-1]['end_exclusive']==735471,'FIXED15_COMPLETE_SCOPE')
 for i,j in enumerate(jobs):need(j['job']==i and 0<j['states']==j['end_exclusive']-j['start']<=49032 and(i==0 or jobs[i-1]['end_exclusive']==j['start']),'DISJOINT_CONTIGUOUS_FIXED_JOBS')
 need(model['prefix_Kmax']==16 and model['prefix_state_count']==735471 and model['tail_K_values']==list(range(17,33))and cert['logical_root_requests']==735487 and cert['workers']==15,'NO_SCOPE_EXPANSION')
 need(cert['tail_discarded']is False and cert['probability_renormalized']is False and cert['old_baseline_replays']==cert['optimizer_kernel_erfc_native_density_calls']==0 and cert['actual_M_lower']is False and cert['new_C6_bound']is False and cert['global_g_proved']is False,'HONEST_SCOPE_FLAGS')
 # Preserve both acknowledged phaseB guard errors; the first prepare error is only partially captured.
 need(release['known_guard_failures_preserved']==2 and read(src/'main_operator_phaseB_prepare_terminal.json')['exit_code']==1 and read(src/'main_operator_phaseB_prepare_terminal.json')['actual_new_interval_root_calls']==0,'FAILED_INDEPENDENT_CLONE_GUARD_PRESERVED_ZERO_ROOT')
 guard1=read(src/'operator_phaseB_prepare_receipt.json');guard2=read(src/'main_operator_phaseB_prepare_terminal.json')
 need(guard1['exit_code']==1 and guard1['pid']==1971 and guard2['child_pid']==2329,'BOTH_PHASEB_GUARD_FAILURES_PRESERVED')
 for receipt,pre in [(guard1,'operator_phaseB_prepare'),(guard2,'main_operator_phaseB_prepare')]:
  for ext,key in [('stdout.txt','stdout_sha256'),('stderr.txt','stderr_sha256')]:need(sha(src/(pre+'_'+ext))==receipt[key],'GUARD_ORIGINAL_STDIO_HASH')
 need(read(src/'main_operator_phaseB_prepare_v2_terminal.json')['exit_code']==0 and sha(src/'main_operator_phaseB_prepare_v2_terminal.json')==release['corrected_prepare_terminal_sha256'],'CORRECTED_EXPLICIT_PREPARE')
 return {'common':common,'groups':groups,'p':ps,'q':qs,'jobs':jobs,'producer_root':str(src),'expected_state_hash':m['expected_state_hash'],'source_tables_canonical_sha256':hashlib.sha256(canon(tab)).hexdigest()}
def rootcheck(wl,wu,vl,vu):need(0<wl<=wu and 0<vl<=vu and vl**5*wu<=S**6<=vu**5*wl,'SAVED_RECIPROCAL_FIFTH_POWER_CROSS_PRODUCTS')
def worker(jobid):
 gate();prefix='workers/job_'+str(jobid).zfill(2);model=read(P/'audit_model.json');src=Path(model['producer_root'])/'phaseB';header=read(src/('chunks/header_'+str(jobid).zfill(2)+'.json'));job=model['jobs'][jobid]
 save(prefix+'_started.json',{'pid':os.getpid(),'job':job,'audit_model_sha256':sha(P/'audit_model.json'),'root_inequalities_completed':0,'restarts':0})
 done=0;start=time.monotonic();cpu=time.process_time();err=None
 try:
  need(all(header[k]==v for k,v in job.items())and header==read(src/'candidate_certificate.json')['chunks'][jobid],'SAME_FROZEN_JOB_HEADER')
  path=src/header['root_file'];need(path.stat().st_size<100000000 and sha(path)==header['root_file_sha256'],'SAVED_ONLY_ROOT_FILE')
  with gzip.open(path,'rt',encoding='ascii')as f:records=json.load(f)
  need(len(records)==job['states']<=49032,'FIXED_ROOT_ROWS')
  powers=[[[int(q)**n for n in range(33)]for q in side]for side in model['q']];common=model['common'];groups=model['groups'];h=hashlib.sha256();counts=[0]*17;mass=[0,0];low=upper=0;first=last=None
  states=islice(combinations_with_replacement(range(17),8),job['start'],job['end_exclusive'])
  for state,pair in zip(states,records,strict=True):
   if first is None:first=state
   last=state;h.update(bytes(state));need(len(pair)==2,'EXACT_TWO_SAVED_ROOT_ENDPOINTS');vl,vu=map(int,pair);wl=common[0]+sum(groups[j][n][0]for j,n in enumerate(state));wu=common[1]+sum(groups[j][n][1]for j,n in enumerate(state));rootcheck(wl,wu,vl,vu)
   # Sequential binomials reconstruct multinomial mass independently of the producer factorial ratio.
   cats=[state[0]]+[state[j]-state[j-1]for j in range(1,8)]+[32-state[-1]];remaining=32;coef=1
   for c in cats[:-1]:coef*=comb(remaining,c);remaining-=c
   need(remaining==cats[-1],'COMPLETE32_IID_CATEGORIES');probs=[]
   for side in range(2):
    prob=coef
    for j,c in enumerate(cats):prob*=powers[side][j][c]
    probs.append(prob);mass[side]+=prob
   low+=probs[1]*vl;upper+=probs[0]*vu;counts[state[-1]]+=1;done+=1
   if done%4096==0:checkpoint(prefix+'_checkpoint.json',{'pid':os.getpid(),'verified_states':done,'last_state':list(state),'root_inequalities_completed':done,'uncheckpointed_work':'UNKNOWN until terminal receipt'})
  need(done==job['states']and header['first']==list(first)and header['last']==list(last)and h.hexdigest()==header['state_bytes_sha256'],'EXACT_STATE_ORDER_ENDPOINTS')
  need(counts==header['prefix_count_by_K']and mass==list(map(int,header['probability_numerators_by_law']))and low==int(header['lower_at_pU_numerator'])and upper==int(header['upper_at_pL_numerator']),'ALL_EXACT_CHUNK_MASSES_SUMS')
  need(header['integer_tables_sha256']==model['source_tables_canonical_sha256'],'SAME_NEW_W_TABLES')
  save(prefix+'_result.json',{'status':'PASS_FIXED_SAVED_JOB','job':job,'pid':os.getpid(),'verified_states':done,'state_bytes_sha256':h.hexdigest(),'counts_by_K':counts,'probability_numerators_by_law':list(map(str,mass)),'prefix_lower_numerator':str(low),'prefix_upper_numerator':str(upper),'new_root_or_producer_calls':0,'audit_model_sha256':sha(P/'audit_model.json')})
 except BaseException:err=traceback.format_exc()
 save(prefix+'_exit.json',{'status':'PASS_FIXED_SAVED_JOB'if err is None else'FAILED_NO_RETRY','pid':os.getpid(),'verified_states':done,'root_inequalities_completed':done,'unknown_work_on_abrupt_exit':'last checkpoint only; never infer planned as completed','failure':err,'CPU_before_receipt':time.process_time()-cpu,'wall_before_receipt':time.monotonic()-start,'automatic_retries':0})
 return 0 if err is None else 1
def collect(model,receipts):
 src=Path(model['producer_root'])/'phaseB';cert=read(src/'candidate_certificate.json');inp=read(Path(model['producer_root'])/'inputs.json');oldroles={r['role']:Path(r['path'])for r in read(Path(model['producer_root'])/'STATIC_SOURCE_PINS.json')['sources']}
 need(len(receipts)==15 and all(r['child_reaped']and r['exit_code']==0 for r in receipts),'ALL15_EXACT_ORIGINAL_CHILDREN_REAPED')
 counts=[0]*17;mass=[0,0];totalL=totalU=0;verified=0
 for i,job in enumerate(model['jobs']):
  pre='workers/job_'+str(i).zfill(2);result=read(P/(pre+'_result.json'));ex=read(P/(pre+'_exit.json'));started=read(P/(pre+'_started.json'));receipt=next(r for r in receipts if r['job']==i)
  need(result['status']==ex['status']=='PASS_FIXED_SAVED_JOB'and result['job']==job and started['job']==job and result['pid']==ex['pid']==started['pid']==receipt['pid'],'EACH15_TERMINAL_RESULT_IDENTITY')
  need(result['audit_model_sha256']==started['audit_model_sha256']==sha(P/'audit_model.json')and result['verified_states']==ex['verified_states']==job['states'],'FIXED_MODEL_NO_OMISSION')
  verified+=result['verified_states'];totalL+=int(result['prefix_lower_numerator']);totalU+=int(result['prefix_upper_numerator'])
  for k,v in enumerate(result['counts_by_K']):counts[k]+=v
  for side,v in enumerate(result['probability_numerators_by_law']):mass[side]+=int(v)
 # Geometry/order enumeration only; does not evaluate an old or new root again.
 h=hashlib.sha256();n=0
 for state in combinations_with_replacement(range(17),8):h.update(bytes(state));n+=1
 need(verified==n==cert['prefix_states']==735471 and h.hexdigest()==model['expected_state_hash']==cert['state_bytes_sha256'],'COMPLETE_735471_PREFIX_BIJECTION')
 need(counts==[comb(k+7,7)for k in range(17)]==cert['prefix_count_by_K']and mass==list(map(int,cert['prefix_probability_numerators_by_law'])),'ALL_K_COUNTS_AND_PREFIX_MASSES')
 common=model['common'];group=model['groups'];ps=model['p'];tailMass=[0,0];tailU=0;taildone=0
 need(len(cert['tail_records'])==16,'EXACT16_SAVED_TAIL_ROWS')
 for K,row in zip(range(17,33),cert['tail_records'],strict=True):
  wl,wu=common[0]+group[-1][K][0],common[1]+group[-1][K][1];need(row['K']==K and [wl,wu]==list(map(int,row['W_units'])),'POSITIVE_LAST_H_TAIL_RELAXATION')
  vl,vu=map(int,row['v_units']);rootcheck(wl,wu,vl,vu);probs=[comb(32,K)*ps[s][-1]**K*(S-ps[s][-1])**(32-K)for s in range(2)]
  need(probs==list(map(int,row['probability_numerators_by_law'])),'BOTH_BINOMIAL_K_TAIL_MASSES')
  for side,v in enumerate(probs):tailMass[side]+=v
  tailU+=probs[0]*vu;taildone+=1
 need(tailMass==list(map(int,cert['tail_probability_numerators_by_law']))and all(mass[side]+tailMass[side]==S**32 for side in range(2))and cert['both_complete_endpoint_laws_unit_mass']is True,'BOTH_COMPLETE_LAWS_EXACT_UNIT_MASS_NO_RENORMALIZATION')
 HL,HU=Q(totalL,S**33),Q(totalU+tailU,S**33)
 need(0<HL<=HU and [HL,HU]==list(map(Q,cert['H_interval']))and[HL**5,HU**5]==list(map(Q,cert['H_power5_interval'])),'FULL_H_LOWER_PREFIX_UPPER_PREFIX_PLUS_TAIL')
 U=HU**5;oldU=Q(read(oldroles['128_certificate'])['M_power6_upper']);need(U==Q(cert['M_power6_upper'])and oldU==Q(cert['old128_M_power6_upper']),'COMPLETE_NEW_OLD_M6_EXACT_VALUES')
 theta0=Q(8,9)**6*24576**5/65000**6;thetas=Q(8,9)**6*24576**5/63500**6
 need(theta0==Q(cert['threshold_g0_power6_exact'])==Q(inp['model_input']['threshold_M_power6_exact'])and thetas==Q(cert['threshold_gstar_power6_exact'])==Q(inp['threshold_gstar_power6_exact']),'ORIGINAL_ALPHA_THRESHOLDS_G0_GSTAR')
 # I is the original FULL R6 integral; this identity includes all constants, cost exactly once.
 I5=Q(3,16)**6*U/384**5;need(I5==Q(cert['full_I_upper_power5_from123'])==Q(3,2**59)*U,'WHOLE_I123_EXACT_INTERFACE')
 decisions={'improves_saved128':U<oldU,'beats_g0_M_slot':U<theta0,'beats_gstar_M_slot':U<thetas}
 need(all(cert[k]is v for k,v in decisions.items()),'STRICT_EXACT_SAVED_DECISIONS')
 for k,v in [('old128',oldU),('g0',theta0),('gstar',thetas)]:need(Q(cert['exact_margins_target_minus_upper'][k])==v-U,'EXACT_FULL_COMPARISON_MARGIN')
 for den,theta in [(65000,theta0),(63500,thetas)]:need((Q(9,8)**6*U*den**6/384**5<64**5)==(U<theta),'CONDITIONAL_U_POWER5_EQUIVALENCE')
 worker_chain=[]
 for i in range(15):
  pre='workers/job_'+str(i).zfill(2)
  worker_chain.append({'job':i,'files':{suffix:sha(P/(pre+'_'+suffix+'.json'))for suffix in ['intent','process','started','result','exit','reap']}})
 save('worker_chain.json',{'workers':worker_chain,'all_15_original_PIDs_reaped':True,'audit_model_sha256':sha(P/'audit_model.json')})
 pins()
 save('certificate.json',{'status':PASS,'producer_certificate_sha256':sha(src/'candidate_certificate.json'),'candidate_sha256':sha(src/'candidate.json'),'H_interval':list(map(str,[HL,HU])),'M_power6_upper':str(U),'full_I_upper_power5_from123':str(I5),'normalizing_true_cost_upper':'1','source156_main_acceptance_sha256':sha(read(P/'inputs.json')['source156_acceptance']),'worker_chain_sha256':sha(P/'worker_chain.json'),'audit_model_sha256':sha(P/'audit_model.json'),'decisions':decisions,'counts':{'prefix_states':verified,'tail_records':taildone,'saved_root_inequalities':verified+taildone,'workers':15,'workers_reaped':15,'complete_endpoint_laws':2},'both_complete_endpoint_laws_exact_unit_mass':True,'old_baseline_replays':0,'new_fifthroot_Arb_kernel_native_optimizer_calls':0,'actual_M_lower_bound':False,'new_C6_bound':False,'global_g_proved':False,'conditional_M_slot_only':True,'release_sha256':sha(P/'parent_release.json')})
 return decisions

def run():
 inp=gate();need(not(P/'started.json').exists()and not(P/'workers').exists(),'ONE_ATTEMPT_NO_RESTART');save('started.json',{'pid':os.getpid(),'release_sha256':sha(P/'parent_release.json'),'workers_max':15,'no_automatic_retry':True});(P/'workers').mkdir()
 start=time.monotonic();cpu=time.process_time();err=None;receipts=[];issued={};decisions=None;stage='source_model_validation'
 try:
  pins();model=model_check(inp);save('audit_model.json',model);stage='fixed15_saved_jobs'
  try:
   for job in model['jobs']:
    i=job['job'];pre='workers/job_'+str(i).zfill(2)
    save(pre+'_intent.json',{'job':job,'checker_sha256':sha(Path(__file__)),'audit_model_sha256':sha(P/'audit_model.json'),'started_observed':False})
    with(P/(pre+'_stdout.txt')).open('xb')as out,(P/(pre+'_stderr.txt')).open('xb')as stderr:
     child=subprocess.Popen([sys.executable,'-B',str(Path(__file__)),'--worker',str(i)],stdin=subprocess.DEVNULL,stdout=out,stderr=stderr,cwd=P)
    issued[child.pid]={'job':i,'pid':child.pid,'popen':child,'issued_at':time.monotonic()}
    save(pre+'_process.json',{'controller_pid':os.getpid(),'pid':child.pid,'job':job,'audit_model_sha256':sha(P/'audit_model.json')})
  finally:
   # Exact wait4 on each owned original PID. If interrupted, UNKNOWN remains; no replacement or timeout.
   for pid,row in issued.items():
    waited,status,usage=os.wait4(pid,0);need(waited==pid,'EXACT_OWNED_WORKER_REAP');rc=os.waitstatus_to_exitcode(status);row['popen'].returncode=rc
    receipt={'job':row['job'],'pid':pid,'exit_code':rc,'child_reaped':True,'wall_issue_to_reap_seconds':time.monotonic()-row['issued_at'],'user_CPU_seconds':usage.ru_utime,'system_CPU_seconds':usage.ru_stime}
    receipts.append(receipt);save('workers/job_'+str(row['job']).zfill(2)+'_reap.json',receipt)
  need(len(issued)==15,'EXACTLY15_FIXED_WORKERS_ISSUED')
  stage='full_sum_tail_comparisons';decisions=collect(model,receipts);stage='complete'
 except BaseException:err=traceback.format_exc()
 status=PASS if err is None else'FAILED_OR_INCOMPLETE_NO_RETRY'
 save('exit_receipt.json',{'status':status,'pid':os.getpid(),'stage':stage,'failure':err,'workers_issued':len(issued),'workers_reaped':len(receipts),'reaped_worker_pids':[r['pid']for r in receipts],'unreaped_issued_pids':[pid for pid in issued if pid not in {r['pid']for r in receipts}],'unissued_jobs':15-len(issued),'work_without_completed_receipt':'UNKNOWN; never inferred from planned count','certificate_sha256':sha(P/'certificate.json')if(P/'certificate.json').exists()else None,'counts':read(P/'certificate.json')['counts']if err is None else None,'decisions':decisions,'controller_CPU_before_receipt':time.process_time()-cpu,'wall_before_receipt':time.monotonic()-start,'automatic_retries':0,'new_fifthroot_Arb_kernel_native_optimizer_old_baseline_calls':0})
 print(json.dumps({'status':status,'workers_issued':len(issued),'workers_reaped':len(receipts),'decisions':decisions}),flush=True)
 return 0 if err is None else 1
if __name__=='__main__':
 if len(sys.argv)==3 and sys.argv[1]=='--worker':
  need(sys.argv[2]in list(map(str,range(15))),'FIXED_JOB_NUMBER');raise SystemExit(worker(int(sys.argv[2])))
 need(len(sys.argv)==1,'NO_NEW_PARAMETERS');raise SystemExit(run())
