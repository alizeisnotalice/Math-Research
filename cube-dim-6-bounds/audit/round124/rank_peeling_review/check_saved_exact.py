"""One E124 independent pure saved-exact configuration; no old imports."""
import sys,os
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
import datetime,hashlib,json,time,traceback
from fractions import Fraction as Q
from pathlib import Path
from math import comb,factorial,isqrt
P=Path(__file__).resolve().parent
start,cpu=time.monotonic(),time.process_time()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def save(n,x):(P/n).write_text(json.dumps(x,indent=2)+'\n')
if any((P/n).exists() for n in ['started.json','review_results.json','run_receipt.json']):raise SystemExit('one attempt only')
reg=read(P/'registration.json')
assert sha(Path(__file__))==reg['script_sha256']
assert sha(P/'source_manifest.json')==reg['source_manifest_sha256']
save('started.json',{'pid':os.getpid(),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'registration_sha256':sha(P/'registration.json')})
checks=0;status='FIRST_FAILURE_STOP';stage='pins';failure=None;roots=0;completed=0
def ck(ok,label):
    global checks
    assert ok,label
    checks+=1
def sqrt_interval(x):
    S=1<<384;a=isqrt(x.numerator*S*S//x.denominator)
    lo,hi=Q(a,S),Q(a+1,S)
    ck(lo*lo<=x<=hi*hi,'independent integer square root')
    return lo,hi
def count_basis(n,k):
    if k==0:return 1
    if k<=3:return comb(n,k) if n>=k else 0
    if n<=2**(k-2):return 0
    z=n
    for rank in range(1,k):z*=n-2**(rank-1)
    den=factorial(k)
    return -(-z//den)
def binomial(p):return [Q(comb(32,n))*p**n*(1-p)**(32-n) for n in range(33)]
def cut_bounds(h,piL,piU):
    # Even partial sums through128 are upper; remove last term for lower.
    term=Q(1);eupper=Q(1);fupper=Q(1)
    for j in range(1,129):
        term*= -h/j;eupper+=term;fupper+=term/(2*j+1)
    elower=eupper-term;flower=fupper-term/257
    ck(0<flower<=fupper and 0<elower<=eupper<1,'alternating exp/erf partial bounds')
    z=(1-eupper,1-elower)
    sqL=sqrt_interval(h/piU)[0];sqU=sqrt_interval(h/piL)[1]
    return z,(2*sqL*flower,2*sqU*fupper)
try:
    pins={}
    for row in read(P/'source_manifest.json')['files']:
        p=Path(row['path']);ck(sha(p)==row['sha256'],'source hash');ck(p.stat().st_size==row['bytes'],'source size');pins[str(p)]=sha(p)
    A=Path(reg['author_directory']);man=read(A/'manifest.json')
    ck(sha(A/'manifest.json')==reg['author_manifest_sha256'],'author final manifest')
    for row in man['whitelist']:ck(sha(A/row['path'])==row['sha256'] and (A/row['path']).stat().st_size==row['bytes'],'author complete whitelist')
    inp=read(A/'upstream.json');out=read(A/'results.json');model=read(A/'model_constants.json');cert=read(A/'certificate.json');summary=read(A/'summary.json')
    fan=read(Path(reg['fan_source']));fancheck=read(Path(reg['fan_audit']));old33=read(Path(reg['projection_source']))
    R0=Q(inp['full_fan_R0'])
    ck(R0==Q(inp['independent_fan_R0'])==Q(fan['full_exp_minus_symbol_integral'])==Q(fancheck['full_exp_minus_symbol_integral'])<Q(27,2**21),'exact full fan R0 normalization input')
    ck(fancheck['status']=='PASS' and fancheck['input_sha256']==sha(Path(reg['fan_source'])) and fancheck['source_sha256']==sha(Path(reg['fan_audit_code'])),'saved independent fan proof lineage, no rerun')
    cases={r['id']:r for r in old33['cases']}
    ck(Q(inp['R1'])==Q(cases['r1']['R_exact']),'exact singleton coefficient')
    for row in inp['R2']:ck(Q(row['R'])==Q(cases[row['id']]['R_exact']),'exact pair coefficient')
    for row in inp['R3full']:ck(Q(row['R'])==Q(cases[row['id']]['R_exact']) and row['unordered_count']==cases[row['id']]['triple_orbit_count'],'exact full triple source class')
    ck([r['unordered_count'] for r in inp['R2']]==[32*(comb(5,1)+comb(5,5))//2,32*(comb(5,2)+comb(5,4))//2,32*comb(5,3)//2] and sum(r['unordered_count'] for r in inp['R2'])==496,'complete unordered pair counts')
    ck(sum(r['unordered_count'] for r in inp['R3full'])==comb(32,3)==4960,'complete triples, no no-ones restriction')
    theta=Q(inp['threshold_M_power6_exact'])
    ck(theta==Q(8,9)**6*24576**5/65000**6==Q(read(Path(reg['source123results']))['threshold_M_power6_exact']),'complete123 threshold only, no MC input')
    ck(model['complete_threshold_M_power6']==str(theta) and inp['MC_values_used'] is False,'fixed threshold/model scope')
    pirow=inp['pi_for_saved_check'];oldpi=read(Path(reg['pi_source']))['constants']['pi']
    ck(pirow==oldpi,'complete original768bit pi')
    piL,piU=[Q(pirow[s+'_mantissa'])*Q(2)**int(pirow[s+'_binary_exponent']) for s in ['lower','upper']]
    ck(3<piL<piU<Q(22,7),'pi range')
    # Check producer metadata; no optimization or producer import.
    rc=read(A/'run_receipt.json');obs=read(A/'observed_tool_receipt.json');tool=read(A/'actual_tool_receipt.json');rel=read(A/'parent_release.json');areg=read(A/'registration.json')
    ck(out['status']==rc['status']=='COMPLETE_EIGHT_FIXED_SINGLE_H_CERTIFICATES' and rc['failure'] is None and rc['completed_h']==8,'author complete batch')
    ck(rc['pid']==tool['science_PID']==89699 and obs['observer_pid']==tool['observer_PID']==89690,'author PIDs')
    ck(obs['exit_code']==tool['exit_code']==0 and obs['child_reaped'] and obs['children_created']==1,'author full exit/reap')
    ck(rel['script_sha256']==obs['science_entry_sha256']==sha(A/'run_prepared.py') and rel['core_sha256']==sha(A/'scalar_core.py'),'released producer/core')
    ck(rel['observer_sha256']==sha(A/'observe_once.py') and rel['registration_sha256']==sha(A/'registration.json') and out['release_sha256']==obs['release_sha256']==sha(A/'parent_release.json'),'release chain')
    ck(read(A/'started.json')['pid']==89699 and read(A/'observer_started.json')['pid']==89690,'exclusive markers')
    for n,digest in areg['prepared_pins'].items():ck(sha(A/n)==digest,'prepared source pin')
    ck(cert['results_sha256']==sha(A/'results.json') and cert['model_constants_sha256']==sha(A/'model_constants.json'),'final certificate bindings')
    ck(tool['saved_only_check']['status']=='PASS_SAVED_ONLY_PURE_RATIONAL_SERIES' and tool['saved_only_check']['checker_sha256']==sha(A/'check_saved.py'),'author pure-check evidence')
    ck(tool['sealing_admin_first_attempt']['exit_code']==1 and tool['sealing_admin_first_attempt']['phase']=='document and manifest writing only; source did not parse or execute','administrative failure retained separately')
    # Independently reconstruct the fixed count lower table and B caps.
    table=[[count_basis(n,k) for k in range(7)] for n in range(33)]
    ck(table==model['Phi_table'],'full Phi table with zero/ceil branches')
    for k in range(7):ck(all(table[n][k]<=table[n+1][k] for n in range(32)),'Phi monotonicity')
    ck(all(table[n][k]==0 for k in range(4,7) for n in range(2**(k-2)+1)),'low N zero branches essential')
    rational={0:R0,1:64*Q(inp['R1']),2:4*sum(Q(r['unordered_count'])*Q(r['R']) for r in inp['R2']),3:8*sum(Q(r['unordered_count'])*Q(r['R']) for r in inp['R3full']),4:Q(comb(32,4)*5,192),6:Q(2*comb(32,6))}
    ck({str(k):str(v) for k,v in rational.items()}==model['pi_cubed_B_rational'],'exact fixed B coefficients')
    coefs={k:(v,v) for k,v in rational.items()};s5=sqrt_interval(Q(8,3));coefs[5]=tuple(Q(comb(32,5),8)*v for v in s5)
    hvalues=[Q(1,2**j) for j in range(12,4,-1)]
    ck(list(map(Q,inp['h_values']))==hvalues and len(out['completed_h'])==8,'exact eight fixed h')
    totalcalls=0;uppers=[];duals=[];rows=[]
    for idx,item in enumerate(out['completed_h']):
        stage='h_'+str(idx)
        ck(item['index']==idx and Q(item['h'])==hvalues[idx] and sha(A/item['path'])==item['sha256'],'saved h lineage')
        row=read(A/item['path']);cap,opt,c=row['caps'],row['optimization'],row['certificate'];h=Q(row['h'])
        ck(h==hvalues[idx] and opt['optimizer_restarts']==opt['random_samples']==0 and opt['objective_calls']<=2500,'fixed h, no restarts/samples')
        totalcalls+=opt['objective_calls']
        z,pcalc=cut_bounds(h,piL,piU);pL,pU=map(Q,cap['p_interval']);blo,beta=list(map(Q,cap['b_lower'])),list(map(Q,cap['beta_upper']))
        ck(0<pL<=pcalc[0]<=pcalc[1]<=pU<1,'independent true-binomial probability enclosure')
        ck(all((x*2**192).denominator==1 for x in [pL,pU]),'saved probability dyadic grid')
        for k in range(7):
            low=coefs[k][0]*z[0]**k/piU**3;up=coefs[k][1]*z[1]**k/piL**3
            ck(0<blo[k]<=low<=up<=beta[k],'independent exact mixed-moment cap enclosure')
            ck((blo[k]*2**240).denominator==(beta[k]*2**240).denominator==1,'cap grid')
        den=int(c['x_grid_denominator']);ints=c['x_grid_numerators'];x=[Q(v,den) for v in ints]
        ck(den==2**44 and len(ints)==7 and min(ints)>=0 and ints[0]>0 and sum(ints)==den,'nonnegative simplex with positive a0')
        aa=[v/b for v,b in zip(x,beta)]
        ck(list(map(Q,c['a_rational']))==aa and sum(a*b for a,b in zip(aa,beta))==1 and c['EdW_upper_exact']=='1','rigorous beta normalization')
        vl=[];vu=[];W=[]
        ck(len(c['root_certificates'])==33,'all33 states retained')
        for n,root in enumerate(c['root_certificates']):
            w=sum(a*ph for a,ph in zip(aa,table[n]));lo,hi=map(Q,root['v_interval'])
            ck(root['n']==n and Q(root['W'])==w>0 and 0<lo<=hi and lo**5<=1/w<=hi**5,'saved reciprocal fifth root')
            ck((lo*2**192).denominator==(hi*2**192).denominator==1,'root grid')
            vl.append(lo);vu.append(hi);W.append(w);roots+=1
        ck(all(W[n]<=W[n+1] for n in range(32)),'W count monotonicity')
        probsL,probsU=binomial(pL),binomial(pU)
        ck(sum(probsL)==sum(probsU)==1 and min(probsL)>0 and min(probsU)>0,'both33-state probabilities normalized')
        HL=sum(a*b for a,b in zip(probsU,vl));HU=sum(a*b for a,b in zip(probsL,vu))
        ck(0<HL<=HU and [HL,HU]==list(map(Q,c['H_interval'])),'monotonic endpoint expectation bounds')
        maxima=[]
        for n in range(33):
            candidates=[probsL[n],probsU[n]];critical=Q(n,32)
            if pL<=critical<=pU:candidates.append(Q(comb(32,n))*critical**n*(1-critical)**(32-n))
            maxima.append(max(candidates))
        terms=[pr*v**6 for pr,v in zip(maxima,vu)]
        gl=[-Q(1,5)/beta[k]*sum(terms[n]*table[n][k] for n in range(33)) for k in range(7)]
        ck(gl==list(map(Q,c['gradient_lower'])),'full probability-max gradient lower vector')
        support=Q(6,5)*HL+min(gl);gamma=min(a/b for a,b in zip(blo,beta))
        ck(support==Q(c['support_plane_lower']) and gamma==Q(c['gamma_exact_cap_rescaling_lower']) and 0<gamma<=1,'Euler support and gamma normalization')
        upper6=HU**5;dual6=gamma*max(Q(0),support)**5
        ck(upper6==Q(c['M_power6_upper']) and dual6==Q(c['fixed_cap_family_M_power6_lower']),'candidate upper and family-only lower')
        ck(theta-upper6==Q(c['threshold_minus_upper_exact']) and dual6-theta==Q(c['fixed_cap_family_lower_minus_threshold']),'exact threshold residuals')
        ml,mu=map(Q,c['M_upper_expression_interval'])
        ck(0<ml<=mu and ml**6<=upper6<=mu**6,'saved candidate sixth-root interval')
        ck(Q(c['full_I_upper_power5_from123'])==Q(3,2**59)*upper6,'complete123 I upper power')
        ck(dual6>Q(7,4)*theta and c['status']==item['status']=='FIXED_SINGLE_H_TABLE_RIGOROUSLY_OBSTRUCTED','strict fixed-family obstruction exceeds7/4 threshold')
        ck(not c['actual_M_lower_bound'] and not c['global_g0_proved'] and not c['new_C6_bound'],'scope flags')
        uppers.append(upper6);duals.append(dual6);rows.append({'h':str(h),'dual_above_7over4_theta':True,'candidate_M_upper_interval':c['M_upper_expression_interval'],'root_checks':33});completed+=1
    ck(totalcalls==out['objective_calls_total']==summary['objective_calls_total']==273,'273 author objective calls preserved, none replayed')
    ck(roots==264 and completed==8 and not out['any_strict_moment_threshold_pass'] and out['all_fixed_single_h_tables_obstructed'],'complete8h/264roots result')
    best=min(range(8),key=lambda i:uppers[i]);minimum=min(duals)
    ck(best==4 and summary['best_h']=='1/256' and summary['minimum_family_lower_witness_h']=='1/256' and duals[4]==minimum,'exact best candidate and minimum witness')
    ck(Q(summary['minimum_family_lower6_minus_7_over4_threshold_exact'])==minimum-Q(7,4)*theta>0,'full summary strict residual')
    simple=Q(cert['best_complete_moment_upper']['strict_simple_M_upper_exact'])
    ck(simple==Q(summary['strict_simple_M_upper_exact'])==Q(68475443529,10**12) and simple**6>uppers[4],'true strict simple moment upper')
    ck(cert['obstruction_scope']['all_family_lower6_exceed_7_over4_threshold'] and not cert['obstruction_scope']['actual_M_lower_bound'] and not cert['obstruction_scope']['all_h_continuum_or_multih_obstructed'],'family limits preserved')
    for p,digest in pins.items():ck(sha(Path(p))==digest,'source unchanged after exact review')
    save('review_results.json',{'status':'PASS_ONE_INDEPENDENT_PURE_SAVED_EXACT_CONFIGURATION','assertions':checks,'completed_h':8,'reciprocal_fifth_roots':roots,'candidate_sixth_root_certificates':8,'series_order':128,'integer_sqrt_bits':384,'all33_states_each_h':True,'all_beta_p_enclosures_independently_verified':True,'all_dual_gamma_witnesses_verified':True,'strict_best_M_upper':str(simple),'exact_threshold_power6':str(theta),'all_fixed_family_lower6_above_7over4_theta':True,'rows':rows,'candidate_complete_M_I_upper_valid':True,'actual_M_I_lower_bound':False,'all_h_or_multih_excluded':False,'globalg0_S_Cn_unpaid':True,'source_paths_pre_post':len(pins),'source_hashes':pins,'optimizer_author_old_function_fan_peak_tuple_Arb_MC_native_grid_MATLAB_children':0,'science_CPU_cores':1,'wall_before_result_write':time.monotonic()-start,'CPU_before_result_write':time.process_time()-cpu})
    status='COMPLETED'
except Exception:
    failure=traceback.format_exc()
finally:
    save('run_receipt.json',{'status':status,'pid':os.getpid(),'assertions':checks,'stage':stage,'completed_h':completed,'root_checks':roots,'failure':failure,'attempts':1,'automatic_retries':0,'children':0,'science_CPU_cores':1,'current_user_runtime_ceiling':None,'wall_before_final_receipt_write':time.monotonic()-start,'CPU_since_start_before_final_receipt_write':time.process_time()-cpu,'CPU_process_total_before_final_receipt_write':time.process_time(),'timing_scope':'one pure exact process before final receipt/stdout/exit; complete tool roundtrip recorded separately','script_sha256':sha(Path(__file__)),'registration_sha256':sha(P/'registration.json')})
    print(json.dumps({'status':status,'assertions':checks,'completed_h':completed,'root_checks':roots,'pid':os.getpid()}))
if status!='COMPLETED':raise SystemExit(1)
