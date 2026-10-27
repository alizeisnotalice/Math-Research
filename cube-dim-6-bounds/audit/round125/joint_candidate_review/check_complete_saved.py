"""One independent complete125 saved-root audit; no producer imports."""
import sys,os
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
import datetime,gzip,hashlib,json,time,traceback
from fractions import Fraction as Q
from pathlib import Path
from itertools import combinations_with_replacement
from functools import lru_cache
from math import comb
P=Path(__file__).resolve().parent
start,cpu=time.monotonic(),time.process_time()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rd(p):return json.loads(p.read_text())
def save(n,x):(P/n).write_text(json.dumps(x,indent=2)+'\n')
if any((P/n).exists() for n in ['started.json','review_results.json','run_receipt.json']):raise SystemExit('one attempt only')
reg=rd(P/'registration.json')
assert sha(Path(__file__))==reg['script_sha256']
assert sha(P/'source_manifest.json')==reg['source_manifest_sha256']
save('started.json',{'pid':os.getpid(),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'registration_sha256':sha(P/'registration.json')})
checks=0;idx=0;taildone=0;status='FIRST_FAILURE_STOP';stage='pins';failure=None
def ck(b,label):
    global checks
    assert b,label
    checks+=1
try:
    pins={}
    for row in rd(P/'source_manifest.json')['files']:
        p=Path(row['path']);ck(sha(p)==row['sha256'],'frozen source hash');ck(p.stat().st_size==row['bytes'],'frozen source bytes');pins[str(p)]=sha(p)
    A=Path(reg['author_directory']);inp=rd(A/'upstream.json');cand=rd(A/'candidate.json');tab=rd(A/'integer_tables.json');cert=rd(A/'candidate_certificate.json');out=rd(A/'results.json');opt=rd(A/'optimization_receipt.json')
    ck(out['status']=='COMPLETE_ONE_JOINT_CANDIDATE' and out['new_candidates']==1,'one complete terminal candidate')
    ck(out['candidate_sha256']==sha(A/'candidate.json') and out['certificate_sha256']==sha(A/'candidate_certificate.json') and out['optimization_sha256']==sha(A/'optimization_receipt.json'),'terminal artifact chain')
    S=1<<192;denProb=S**32;denWeighted=S**33;cross=S**6
    ck(tab['scale_power2']==192 and tab['probability_denominator_power2']==cert['probability_denominator_power2']==6144 and tab['weighted_sum_denominator_power2']==cert['weighted_sum_denominator_power2']==6336,'all integer scale schema fields')
    grid=int(cand['denominator']);nums=cand['numerators']
    ck(grid==2**44 and len(nums)==49 and min(nums)>=0 and nums[0]>0 and sum(nums)==grid,'exact49 simplex')
    beta0=Q(inp['beta0']);beta=[[Q(v) for v in row] for row in inp['beta_by_h']];phi=inp['Phi_table'];x=[Q(v,grid) for v in nums]
    a0=x[0]/beta0;aa=[[x[1+6*j+k]/beta[j][k] for k in range(6)] for j in range(8)]
    ck(Q(tab['a0'])==a0 and [[Q(v) for v in row] for row in tab['a_by_h']]==aa,'candidate coefficient table')
    ck(a0*beta0+sum(aa[j][k]*beta[j][k] for j in range(8) for k in range(6))==1,'exact EdW cap normalization')
    common=list(map(int,tab['common_W_units']));group=[[[int(v) for v in pair] for pair in row] for row in tab['group_W_units']]
    def down(q):return (q.numerator*S)//q.denominator
    def up(q):return -((-q.numerator*S)//q.denominator)
    ck(common==[down(a0),up(a0)] and common[0]>0,'common exact outward W intercept')
    ck(len(group)==8 and all(len(t)==33 for t in group),'all264 group rows')
    for j in range(8):
        for n in range(33):
            value=sum(aa[j][k]*phi[n][k+1] for k in range(6))
            ck(group[j][n]==[down(value),up(value)] and group[j][n][0]>=0,'group outward rounding')
    p=[];q=[];powers=[]
    for side in range(2):
        ps=[]
        for pair in inp['p_intervals']:
            val=Q(pair[side])*S;ck(val.denominator==1,'CDF exact192 grid');ps.append(int(val))
        qs=[ps[0]]+[ps[j]-ps[j-1] for j in range(1,8)]+[S-ps[-1]]
        ck(len(qs)==9 and min(qs)>0 and sum(qs)==S,'one honest nine-bin endpoint law')
        p.append(ps);q.append(qs);powers.append([[v**c for c in range(33)] for v in qs])
    ck(p==[[int(v) for v in row] for row in tab['p_units']] and q==[[int(v) for v in row] for row in tab['category_units']],'saved endpoint probability tables')
    ck(all(l<=u for l,u in zip(p[0],p[1])),'whole CDF componentwise coupling')
    ck(not inp['pending_new_beta_tables_used'] and not inp['MC_data_used'],'frozen124 table only')
    # Literal124 data inheritance; no old scalar roots/objectives are replayed.
    old=Path(reg['old124_directory']);oldmodel=rd(old/'model_constants.json');oldh=[rd(old/('h_'+str(j)+'.json')) for j in range(8)]
    ck(phi==oldmodel['Phi_table'],'inherited Phi table')
    for j,row in enumerate(oldh):ck(inp['beta0']==row['caps']['beta_upper'][0] and inp['beta_by_h'][j]==row['caps']['beta_upper'][1:] and inp['p_intervals'][j]==row['caps']['p_interval'],'literal frozen124 cap/CDF copy')
    oldcmp=inp['124_exact_comparisons'];ck(oldcmp['best_M_power6_upper']==oldh[4]['certificate']['M_power6_upper'] and oldcmp['family_M_power6_lowers']==[row['certificate']['fixed_cap_family_M_power6_lower'] for row in oldh],'complete124 rational comparison copies')
    @lru_cache(maxsize=65536)
    def rootcheck(wl,wu,vl,vu):
        ck(0<wl<=wu and 0<vl<=vu and vl**5*wu<=cross<=vu**5*wl,'exact state fifth-root cross product')
        return True
    def probability_pair(state):
        cats=[state[0]]+[state[j]-state[j-1] for j in range(1,8)]+[32-state[-1]]
        # Independent sequential-binomial construction of the multinomial coefficient.
        remaining=32;coefficient=1
        for c in cats[:-1]:coefficient*=comb(remaining,c);remaining-=c
        if remaining!=cats[-1]:raise AssertionError('source labels do not total32')
        left=right=coefficient
        for j,c in enumerate(cats):left*=powers[0][j][c];right*=powers[1][j][c]
        return left,right
    stage='full_saved_prefix'
    states=iter(combinations_with_replacement(range(17),8));Kcounts=[0]*17;mass=[0,0];totalL=totalU=0;wholehash=hashlib.sha256()
    ck(len(cert['chunks'])==90,'complete chunk count')
    for chunkNo,header in enumerate(cert['chunks']):
        ck(header==rd(A/('cert_chunks/header_'+str(chunkNo).zfill(3)+'.json')),'separate header file matches final certificate')
        ck(header['chunk']==chunkNo and header['start']==idx,'consecutive chunk schema')
        path=A/header['root_file'];ck(sha(path)==header['root_file_sha256'] and path.stat().st_size<100000000,'root chunk hash/size')
        with gzip.open(path,'rt',encoding='ascii') as f:records=json.load(f)
        ck(0<len(records)==header['states']<=8192,'chunk root rows')
        partL=partU=0;partMass=[0,0];partK=[0]*17;partHash=hashlib.sha256();first=last=None
        for pair in records:
            state=next(states)
            if first is None:first=state
            last=state;bs=bytes(state);wholehash.update(bs);partHash.update(bs)
            wl=common[0]+sum(group[j][n][0] for j,n in enumerate(state));wu=common[1]+sum(group[j][n][1] for j,n in enumerate(state))
            vl,vu=map(int,pair);rootcheck(wl,wu,vl,vu)
            pl,pu=probability_pair(state)
            partMass[0]+=pl;partMass[1]+=pu;mass[0]+=pl;mass[1]+=pu
            low,high=pu*vl,pl*vu;partL+=low;partU+=high;totalL+=low;totalU+=high
            Kcounts[state[-1]]+=1;partK[state[-1]]+=1;idx+=1
        ck(header['end_exclusive']==idx and header['first']==list(first) and header['last']==list(last),'full state ordering/header endpoints')
        ck(header['state_bytes_sha256']==partHash.hexdigest() and header['prefix_count_by_K']==partK,'per-chunk coverage hash/counts')
        ck(list(map(int,header['prefix_probability_numerator_by_law']))==partMass and int(header['lower_at_pU_numerator'])==partL and int(header['upper_at_pL_numerator'])==partU,'exact chunk masses and weighted sums')
        if (chunkNo+1)%10==0:
            save('progress.json',{'chunks':chunkNo+1,'verified_prefix_states':idx,'elapsed':time.monotonic()-start})
            print(json.dumps({'phase':stage,'chunks':chunkNo+1,'verified_states':idx}),flush=True)
    ck(next(states,None) is None,'no unverified prefix labels remain')
    ck(idx==cert['prefix_states']==out['prefix_states']==735471 and Kcounts==[comb(k+7,7) for k in range(17)]==cert['prefix_count_by_K'],'complete prefix bijection and per-K counts')
    ck(wholehash.hexdigest()==cert['state_bytes_sha256']==opt['state_bytes_sha256'],'independent full state-order hash')
    ck(mass==list(map(int,cert['prefix_probability_numerator_by_law'])),'complete prefix subprobabilities')
    stage='all16_saved_tails';tailU=0;tailMass=[0,0]
    ck(len(cert['tail_records'])==16,'all16 tail K values')
    for K,row in zip(range(17,33),cert['tail_records']):
        wl,wu=common[0]+group[-1][K][0],common[1]+group[-1][K][1]
        ck(row['K']==K and [wl,wu]==list(map(int,row['W_units'])),'true positive tail W lower bound')
        vl,vu=map(int,row['v_units']);rootcheck(wl,wu,vl,vu)
        probs=[comb(32,K)*p[s][-1]**K*(S-p[s][-1])**(32-K) for s in range(2)]
        ck(probs==list(map(int,row['binomial_probability_numerator_by_law'])),'both exact binomial tail probabilities')
        tailMass[0]+=probs[0];tailMass[1]+=probs[1];term=probs[0]*vu;tailU+=term
        ck(term==int(row['upper_weighted_numerator']),'positive upper-tail weighted term');taildone+=1
    ck(tailMass==list(map(int,cert['tail_probability_numerator_by_law'])) and all(mass[s]+tailMass[s]==denProb for s in range(2)),'BOTH complete laws exact unit mass')
    ck(totalL==int(cert['prefix_lower_at_pU_numerator']) and totalU==int(cert['prefix_upper_at_pL_numerator']) and tailU==int(cert['tail_upper_numerator']),'global exact positive integer sums')
    HL,HU=Q(totalL,denWeighted),Q(totalU+tailU,denWeighted)
    ck(0<HL<=HU and [HL,HU]==list(map(Q,cert['H_interval'])),'complete H interval with upper tail')
    ck(cert['root_requests']==735487 and cert['root_cache_cap']==65536 and not cert['tail_discarded'] and not cert['prefix_renormalized'] and not cert['CDF_endpoint_laws_independent_categories'],'request/scope flags')
    stage='final_rational_compares'
    upper6=HU**5;lowerObjective6=HL**5;theta=Q(inp['threshold_M_power6_exact'])
    ck(theta==Q(8,9)**6*24576**5/65000**6==Q(cert['threshold_M_power6_exact']),'complete threshold')
    ck(upper6==Q(cert['M_power6_upper']) and lowerObjective6==Q(cert['prefix_only_H_power5_lower']) and theta-upper6==Q(cert['threshold_minus_upper_exact']),'candidate exact sixth-power certificate expressions')
    ml,mu=map(Q,cert['M_upper_expression_interval']);ck(0<ml<=mu and ml**6<=upper6<=mu**6,'exact M-upper-expression root interval')
    ck(Q(cert['full_I_upper_power5_from123'])==Q(3,2**59)*upper6,'complete I upper-power transfer')
    oldbest=Q(oldcmp['best_M_power6_upper']);olddual=min(map(Q,oldcmp['family_M_power6_lowers']))
    ck(cert['beats_saved124_best_upper']==(upper6<oldbest) and cert['beats_all_eight_fixed_single_h_cap_families']==(upper6<olddual),'exact old/new comparisons')
    ck(Q(cert['old_best_upper_minus_new_upper_exact'])==oldbest-upper6 and Q(cert['old_minimum_family_lower_minus_new_upper_exact'])==olddual-upper6,'complete rational comparison residuals')
    ck(cert['status']==('STRICT_MOMENT_THRESHOLD_PASS' if upper6<theta else 'UNSEPARATED_OR_NONCOMPETITIVE_FIXED_CANDIDATE'),'future status criterion applied correctly')
    ck(not cert['actual_M_lower_bound'] and not cert['joint_family_optimum_excluded'] and not cert['global_g0_proved'] and not cert['new_C6_bound'],'candidate-only interpretation')
    ck(opt['objective_calls']==out['objective_calls']==116 and opt['optimizer_invocations']==1 and opt['restarts']==0,'one author optimization not replayed')
    rc=rd(A/'run_receipt.json');obs=rd(A/'observed_tool_receipt.json');release=rd(A/'parent_release.json')
    ck(rc['status']=='COMPLETE_ONE_JOINT_CANDIDATE' and rc['failure'] is None and rc['pid']==3475 and obs['observer_pid']==3467 and obs['exit_code']==0 and obs['child_reaped'],'terminal author exit evidence')
    ck(obs['release_sha256']==sha(A/'parent_release.json') and release['script_sha256']==sha(A/'run_prepared.py'),'author release binding')
    for path,digest in pins.items():ck(sha(Path(path))==digest,'source unchanged after complete audit')
    save('review_results.json',{'status':'PASS_COMPLETE_INDEPENDENT_SAVED_JOINT_CANDIDATE','assertions':checks,'prefix_states':idx,'tail_values':taildone,'root_requests_checked':idx+taildone,'root_cache_info':str(rootcheck.cache_info()),'both_endpoint_laws_exact_mass1':True,'all_positive_integer_sums_match':True,'candidate_M_upper_expression_interval':cert['M_upper_expression_interval'],'complete_I_upper_power5':cert['full_I_upper_power5_from123'],'threshold_pass':upper6<theta,'prefix_objective_lower_exceeds_threshold':lowerObjective6>theta,'beats_saved124_upper':upper6<oldbest,'beats_all8_fixed_single_h_family_witnesses':upper6<olddual,'actual_M_I_lower_bound':False,'joint_family_optimum_excluded':False,'globalg0_C6_unchanged':True,'source_paths_pre_post':len(pins),'source_hashes':pins,'optimizer_author_old124_MC_native_Arb_MATLAB_children':0,'candidate_verifier_configurations':1,'science_CPU_cores':1,'wall_before_result_write':time.monotonic()-start,'CPU_before_result_write':time.process_time()-cpu})
    status='COMPLETED'
except Exception:
    failure=traceback.format_exc()
finally:
    save('run_receipt.json',{'status':status,'pid':os.getpid(),'assertions':checks,'stage':stage,'prefix_states':idx,'tail_values':taildone,'failure':failure,'attempts':1,'automatic_retries':0,'children':0,'science_CPU_cores':1,'current_user_runtime_ceiling':None,'wall_before_final_receipt_write':time.monotonic()-start,'CPU_since_start_before_final_receipt_write':time.process_time()-cpu,'CPU_process_total_before_final_receipt_write':time.process_time(),'timing_scope':'one checker process before final receipt/stdout/exit; complete tool roundtrip separately observed','script_sha256':sha(Path(__file__)),'registration_sha256':sha(P/'registration.json')})
    print(json.dumps({'status':status,'assertions':checks,'prefix_states':idx,'tail_values':taildone,'pid':os.getpid()}),flush=True)
if status!='COMPLETED':raise SystemExit(1)
