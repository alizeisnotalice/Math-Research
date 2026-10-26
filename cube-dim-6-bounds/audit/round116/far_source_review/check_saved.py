"""One independent116 saved/exact certificate configuration; no generator."""
import sys
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
import datetime,hashlib,json,os,time,traceback
from fractions import Fraction as Q
from pathlib import Path
P=Path(__file__).resolve().parent
start,cpu=time.monotonic(),time.process_time()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def write(n,x):(P/n).write_text(json.dumps(x,indent=2)+'\n')
if any((P/n).exists() for n in ['started.json','review_results.json','run_receipt.json']):raise SystemExit('one attempt only')
reg=read(P/'registration.json')
assert sha(Path(__file__))==reg['script_sha256']
assert sha(P/'source_manifest.json')==reg['source_manifest_sha256']
write('started.json',{'pid':os.getpid(),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'registration_sha256':sha(P/'registration.json')})
count=0;stage='source';status='FIRST_FAILURE_STOP';failure=None
def ck(ok,label):
    global count
    assert ok,label
    count+=1
def findkey(obj,key,value):
    if isinstance(obj,dict):return obj.get(key)==value or any(findkey(v,key,value) for v in obj.values())
    if isinstance(obj,list):return any(findkey(v,key,value) for v in obj)
    return False
def dyad(d,side):return Q(d[side+'_mantissa'])*Q(2)**int(d[side+'_binary_exponent'])
def check_power_saved(pair,radicand,power):
    lo,hi=map(Q,pair);ck(0<=lo<=hi,'ordered saved positive power')
    ck(lo**power<=radicand<=hi**power,'exact saved power residuals')
    return lo,hi
try:
    sources=read(P/'source_manifest.json')['files'];pre={}
    for r in sources:
        p=Path(r['path']);ck(sha(p)==r['sha256'],'pinned source '+r['path']);ck(p.stat().st_size==r['bytes'],'pinned bytes');pre[str(p)]=sha(p)
    A=Path(reg['author_directory']);am=read(A/'manifest.json')
    ck(sha(A/'manifest.json')==reg['author_manifest_sha256'],'author final manifest')
    for r in am['whitelist']:ck(sha(A/r['path'])==r['sha256'],'author whitelist');ck((A/r['path']).stat().st_size==r['bytes'],'author whitelist bytes')
    inp=read(A/'upstream.json');out=read(A/'results.json');w,t=out['configurations'];cert=read(A/'certificate.json')
    old036=read(Path(reg['source036']));old033=read(Path(reg['source033']));old111=read(Path(reg['source111']))
    ck(inp['pi_saved036']==old036['configurations'][0]['pi'],'complete saved512 pi bytes')
    ck(inp['R1']==next(x['R_exact'] for x in old033['cases'] if x['id']=='r1'),'033 exact R1')
    for r in inp['R2']:
        old=next(x for x in old033['cases'] if x['id']==r['id'])
        ck(old['R_exact']==r['R_exact'] and old['Gram_determinant']==r['Gram_determinant'],'033 exact R2 provenance')
    ck(findkey(old036,'B3_pi6_exact',inp['B3_pi6_exact']),'036 exact B3 provenance')
    reference=old111['final']['comparisons']['ref108_frozen_reference']
    ck(inp['far_allowance']['lower']==reference['S_plan_lower'],'complete111 lower copied verbatim')
    ck(inp['far_allowance']['upper']==reference['S_plan_upper'],'complete111 upper copied verbatim')
    ck('Kref=1/(8pi5)' in reference['reference_assumption'] and reference['g0']=='unproved globally','111 comparison qualification')
    RL,RU=Q(inp['far_allowance']['lower']),Q(inp['far_allowance']['upper'])
    ck(0<RL<=RU,'positive complete allowance')
    targetL=RL+Q(reference['CEH_plus_middle_I_upper_certificate'])
    targetU=RU+Q(reference['CEH_plus_middle_I_lower_expression'])
    ck(0<targetL<=targetU and targetL**5<=Q(64)**5*Q(1,390000)**6<=targetU**5,'111 subtraction directions and complete target')
    # Source/continuation evidence is part of the saved certificate.
    rc0=read(A/'run_receipt.json');rc2=read(A/'run_receipt_v2.json');obs=read(A/'observed_tool_receipt_v2.json');tools=read(A/'tool_observation.json')
    ck(rc0['status']=='FAILED' and rc0['completed']==1 and '4300 digits' in rc0['failure'],'original cfg2 output failure retained')
    ck(rc0['pid']==67690 and tools['original']['exit_code']==1,'original PID and failed exit')
    ck(out['cfg1_original_checkpoint_sha256']==sha(A/'checkpoint.json') and out['original_failure_receipt_sha256']==sha(A/'run_receipt.json'),'v2 original lineage')
    ck(w==read(A/'checkpoint.json')['completed_configurations'][0],'cfg1 reused verbatim')
    ck(out['cfg_window_calls_v2']==rc2['cfg_window_calls']==0 and out['cfg_tail_calls_v2']==rc2['cfg_tail_calls']==1,'only tail continuation once')
    ck(rc2['failure'] is None and rc2['pid']==68172 and obs['scientific_child_exit_code']==0 and obs['child_reaped'],'v2 child exit/reap')
    ck(obs['observer_pid']==68164 and obs['own_children_created']==1,'observer scope')
    ck(obs['entry_sha256']==sha(A/'resume_tail_v2.py') and obs['release_sha256']==sha(A/'parent_release_v2.json'),'observer entry release pins')
    ck(tools['scientific_function_attempts']=={'cfg_window':1,'cfg_tail':2},'actual author function attempts')
    ck(not tools['first_failure_reclassified_as_success'],'failure not rewritten')
    ck(cert['results_sha256']==sha(A/'results.json') and cert['upstream_sha256']==sha(A/'upstream.json'),'certificate bindings')
    stage='saved_exact'
    piL,piU=(dyad(inp['pi_saved036'],s) for s in ['lower','upper'])
    ck(3<piL<piU<Q(22,7),'saved pi enclosure')
    rho,a,b,m=Q(1,64),Q(1,64),Q(1,32),Q(7,1760)
    ck(Q(1,32)**2>rho**2/Q(60**2*6),'alpha exceeds beta')
    ck(w['fixed_window']==['1/64','1/32'] and w['rho']=='1/64','fixed original window')
    ck([r['ordered_count'] for r in inp['R2']]==[192,480,320] and sum(r['ordered_count'] for r in inp['R2'])==992,'complete pair counts')
    B2=64*sum(Q(r['ordered_count'])*Q(r['R_exact']) for r in inp['R2'])
    ck(Q(w['B2_pi6_exact_from_saved_R2'])==B2,'saved B2 identity')
    # Verify stored endpoint certificates by clearing exact denominators.
    ds=list(map(Q,w['D_at_b_certified_lower']))
    den=[b**6,b**4*(1+1024*b*b),b*b*(1+1024*b*b)*(rho*rho+992*b*b),(1+1024*b*b)*(rho*rho+992*b*b)*(rho*rho+960*b*b)]
    coef=[Q(216,2**30),32*Q(inp['R1']),B2,Q(inp['B3_pi6_exact'])]
    for d,f,c in zip(ds,den,coef):ck(d*piU**6*f==c and d>Q(1,5),'stored pointwise-window certificate')
    prob=Q(w['source_probability_lower'])
    ck(prob==Q(23,24)**31-Q(32,180**30)>Q(15,64),'saved source probability lower')
    ck(Q(w['source_probability_minus_15_over64'])==prob-Q(15,64),'probability residual')
    mass=Q(w['source_mass_window_lower'])
    ck(mass==a*Q(15,64)/(Q(22,7)*(Q(1,4)+b*b))>m,'saved source mass lower')
    ck(Q(1,80)/piU>m,'face mass lower at a')
    # Rational majorants prove Ctube<150, without a new geometric grid.
    ck(Q(10,7)**2>2 and Q(9,4)**2>5,'root majorants for volume')
    ck(Q(22,7)**2/2*(10+10*Q(10,7)+2*Q(9,4))<150,'true tube coefficient below150')
    ck(Q(w['far_volume_lower'])==1-150*rho**4>Q(1,2) and Q(1,10)>m,'true far volume and cap lower')
    W5=Q(w['I_envelope_lower_to_fifth'])
    ck(W5==(b-a)**6*m**5/5,'same-time envelope fifth power')
    ck(Q(w['residual_to_RL_fifth'])==W5-RL**5>0 and Q(w['residual_to_RU_fifth'])==W5-RU**5>0,'both complete111 residuals strictly positive')
    Wlo,Whi=check_power_saved(w['I_envelope_lower_interval'],W5,5)
    ratioL,ratioU=check_power_saved(w['obstruction_to_RU_ratio_interval'],W5/RU**5,5)
    ck(W5>(4*RU)**5 and Wlo>Q(19,1000000),'strict coarse mechanism separation')
    ck(cert['cfg1']['I_envelope_lower_to_fifth']==w['I_envelope_lower_to_fifth'] and cert['cfg1']['I_envelope_lower_interval']==w['I_envelope_lower_interval'],'cfg1 certificate copy')
    T=Q(t['T']);pL,pU=map(Q,t['P6_interval']);tl,tu=map(Q,t['tail_norm_interval'])
    ck(T==Q(1,8) and pL*piU**6==pU*piL**6==Q(216,2**30),'saved original6D peak endpoints')
    ck(tl*5*T**5==pL and tu*5*T**5==pU,'actual integrated tail norm endpoints')
    ck(Q(t['new_over_old_exact'])==Q(135,4194304)==Q(160*216,2**30),'old095 peak exact comparison')
    tail_res=Q(t['tail_below_one_twentieth_remaining_norm_residual'])
    ck(tail_res==RL**5/20**6-tu**6>0,'tail uses less than one twentieth norm')
    sl,su=map(Q,t['short_time_norm_allocation_interval'])
    ck(0<sl<=su and (sl+tu)**6<=RL**5<=(su+tl)**6,'short norm interval and subtraction direction')
    il,iu=map(Q,t['tail_I_if_alone_interval'])
    ck(0<il<=iu and il**5<=tl**6<=tu**6<=iu**5,'tail-I alone saved powers only')
    ck(tu<Q(1372,10**9) and sl>Q(33,10**6),'strict coarse tail and remaining norm bounds')
    ck(t['norm_addition_required_before_q_power'] and not w['actual_I_far_lower_bound'],'norm interface and envelope-only obstruction')
    ck(not out['new_C6_bound'] and not out['full_far_paid'] and not t['global_g0_proved'],'unpaid global scope')
    ck(cert['cfg2']['tail_norm_interval']==t['tail_norm_interval'] and cert['cfg2']['short_time_norm_allocation_interval']==t['short_time_norm_allocation_interval'],'cfg2 certificate copy')
    for p,h in pre.items():ck(sha(Path(p))==h,'source unchanged after exact check')
    write('review_results.json',{'status':'PASS_ONE_INDEPENDENT_SAVED_EXACT_CONFIGURATION','assertions':count,'new_saved_exact_configurations':1,'source_pins_verified_pre_post':len(pre),'author_source_pins':12,'author_whitelist_files':len(am['whitelist']),'W5':str(W5),'strict_certificates':{'W5_above_both111_fifth_endpoints':True,'W_lower_exceeds_four_times_RU':True,'tail_norm_less_than_RL_power5over6_div20':True,'tail_norm_upper_less_than':'1372/1000000000','short_norm_allocation_lower_greater_than':'33/1000000'},'author_lineage':{'original_failed_completed1':True,'v2_cfg_window_calls':0,'v2_cfg_tail_calls':1,'total_author_tail_function_attempts':2,'not_single_initial_success':True},'no_actual_I_lower_bound_or_new_Cn':True,'unpaid':['actual far0<t<1/8','globalg0','fullglobalS','generalL1/allscale numerical closure'],'source_hashes':pre,'author_generators_old_functions_native_kernel_grids_MATLAB_children':0,'wall_before_result_write':time.monotonic()-start,'CPU_before_result_write':time.process_time()-cpu})
    status='COMPLETED'
except Exception:
    failure=traceback.format_exc()
finally:
    write('run_receipt.json',{'status':status,'pid':os.getpid(),'assertions':count,'stage':stage,'failure':failure,'attempts':1,'automatic_retries':0,'children':0,'current_user_runtime_ceiling':None,'wall_before_final_receipt_write':time.monotonic()-start,'CPU_since_start_before_final_receipt_write':time.process_time()-cpu,'CPU_process_total_before_final_receipt_write':time.process_time(),'timing_scope':'one process sample before final receipt/stdout/exit; complete awaited tool roundtrip recorded separately','script_sha256':sha(Path(__file__)),'registration_sha256':sha(P/'registration.json')})
    print(json.dumps({'status':status,'assertions':count,'pid':os.getpid()}))
if status!='COMPLETED':raise SystemExit(1)
