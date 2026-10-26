"""E120 one exact saved-certificate configuration, no new experiment."""
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
def save(n,x):(P/n).write_text(json.dumps(x,indent=2)+'\n')
if any((P/n).exists() for n in ['started.json','review_results.json','run_receipt.json']):raise SystemExit('one attempt only')
reg=read(P/'registration.json')
assert sha(Path(__file__))==reg['script_sha256']
assert sha(P/'source_manifest.json')==reg['source_manifest_sha256']
save('started.json',{'pid':os.getpid(),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'registration_sha256':sha(P/'registration.json')})
count=0;status='FIRST_FAILURE_STOP';stage='pins';failure=None
def ck(ok,label):
    global count
    assert ok,label
    count+=1
def power_interval(pair,x,n=5):
    lo,hi=map(Q,pair)
    ck(0<lo<=hi and lo**n<=x<=hi**n,'exact saved power endpoints')
    return lo,hi
def dyad(d,s):return Q(d[s+'_mantissa'])*Q(2)**int(d[s+'_binary_exponent'])
try:
    src=read(P/'source_manifest.json');pins={}
    for r in src['files']:
        p=Path(r['path']);ck(sha(p)==r['sha256'],'source hash');ck(p.stat().st_size==r['bytes'],'source bytes');pins[str(p)]=sha(p)
    A=Path(reg['author_directory']);am=read(A/'manifest.json')
    ck(sha(A/'manifest.json')==reg['author_manifest_sha256'],'author final seal')
    for r in am['whitelist']:ck(sha(A/r['path'])==r['sha256'] and (A/r['path']).stat().st_size==r['bytes'],'full author whitelist')
    inp=read(A/'upstream.json');out=read(A/'results.json');g,b=out['configurations'];cert=read(A/'certificate.json')
    oldpi=read(Path(reg['pi_source']))['configurations'][0]['pi']
    ck(inp['pi_saved036']==oldpi,'complete saved512 pi copy')
    piL,piU=(dyad(oldpi,s) for s in ['lower','upper'])
    ck(3<piL<piU<Q(22,7),'saved pi range')
    ref=read(Path(reg['source111']))['final']['comparisons']['ref108_frozen_reference']
    ck(inp['far_allowance']['lower']==ref['S_plan_lower'] and inp['far_allowance']['upper']==ref['S_plan_upper'],'complete111 endpoints verbatim')
    ck('Kref=1/(8pi5)' in ref['reference_assumption'] and ref['g0']=='unproved globally','correct111 comparison')
    RL,RU=Q(ref['S_plan_lower']),Q(ref['S_plan_upper'])
    TL=RL+Q(ref['CEH_plus_middle_I_upper_certificate']);TU=RU+Q(ref['CEH_plus_middle_I_lower_expression'])
    ck(0<RL<=RU and 0<TL<=TU,'allowance and target order')
    ck(TL**5<=Q(64)**5*Q(1,390000)**6<=TU**5,'complete111 target subtraction direction')
    ck(cert['proof_sha256']==sha(A/'proof.md') and cert['results_sha256']==sha(A/'results.json') and cert['upstream_sha256']==sha(A/'upstream.json'),'certificate core bindings')
    rc=read(A/'run_receipt.json');obs=read(A/'observed_tool_receipt.json');tool=read(A/'actual_tool_receipt.json');rel=read(A/'parent_release.json');areg=read(A/'registration.json')
    ck(out['status']==rc['status']=='COMPLETE_TWO_FIXED_CONFIGURATIONS' and rc['completed_configurations']==2 and rc['failure'] is None,'producer completed two')
    ck(rc['pid']==tool['science_PID']==71808 and obs['observer_pid']==tool['observer_PID']==71800,'producer PIDs')
    ck(obs['exit_code']==tool['science_tool_exit_code']==0 and obs['child_reaped'] and obs['children_created']==1,'observer final child exit')
    ck(rc['children']==0 and rc['native_kernel_or_old_program_calls']==0,'science no child/old calls')
    ck(obs['science_entry_sha256']==rel['script_sha256']==sha(A/'run_prepared.py'),'producer source binding')
    ck(rel['observer_sha256']==sha(A/'observe_once.py') and rel['registration_sha256']==sha(A/'registration.json'),'observer and registration release')
    ck(out['release_sha256']==obs['release_sha256']==sha(A/'parent_release.json'),'release chain')
    ck(tool['unique_science_run_count']==1 and tool['unique_scientific_configurations']==2,'single producer invocation')
    ck(read(A/'started.json')['pid']==71808 and read(A/'observer_started.json')['pid']==71800,'exclusive-start lineage')
    ck(read(A/'checkpoint.json')['configurations']==out['configurations'],'complete checkpoint consistency')
    for n,h in areg['prepared_pins'].items():ck(sha(A/n)==h,'producer prepared pins')
    ck(tool['saved_only_check']['status']=='PASS_SAVED_ONLY' and tool['saved_only_check']['exit_code']==0 and tool['saved_only_check']['checker_sha256']==sha(A/'check_saved.py'),'author saved evidence')
    stage='saved_exact'
    q=Q(inp['q']);R,rho=Q(inp['R']),Q(inp['rho'])
    ck(q==Q(6,5) and Q(5,2)*q==3 and q-1==Q(1,5),'critical Jensen powers')
    ck(6*2*32==384 and Q(3,16)**6*64/Q(384)**5==Q(3,2**53),'critical constant collapse pi^3 coefficient')
    ck(R==256 and rho==Q(1,64) and rho*R==4,'single fixed geometry')
    ck(32*31//2==g['all_unordered_actual_pairplanes']==496,'all actual unordered pairs')
    tail=Q(g['tail_coefficient_over_pi']);tube=Q(g['tube_coefficient_over_pi_squared'])
    ck(tail==Q(384)/R==Q(3,2),'coordinate marginal tail coefficient')
    ck(tube==Q(496,2)*(rho*R)**4*Q(14,3)**2/Q(32)**4==Q(1519,1152),'rank4 ball-union coefficient')
    P0=Q(g['P_far_lower_coarse']);Pref=Q(g['P_far_lower_refined'])
    ck(P0==1-tail/3-tube/9==Q(3665,10368)>0,'coarse true far probability')
    ck(Pref==1-tail/piL-tube/piL**2 and P0<Pref<1,'same fixed R refined probability lower')
    ck(g['original_time_first'] and g['original_full_law_projection_only'],'original law qualifiers')
    # Constants from the analytic p1(0) lower/upper witnesses, no peak evaluation.
    ck(Q(120,2**36)==Q(15,2**33),'radial Fourier lower coefficient pi^-3')
    ck(Q(3,2**53)*Q(15,2**33)==Q(45,2**86),'Jensen lower coefficient')
    ck(Q(3,2**53)*Q(216,2**30)==Q(81,2**80),'complete actual far upper coefficient pi^-3')
    W5=Q(b['J_far_lower_fifth_exact']);Wref5=Q(b['J_far_refined_lower_fifth'])
    ck(W5==Q(45,2**86)*P0**5 and Wref5==Q(45,2**86)*Pref**5>W5,'saved exact Jensen lower fifth powers')
    ck(Q(b['residual_vs_RL_fifth'])==W5-RL**5>0 and Q(b['residual_vs_RU_fifth'])==W5-RU**5>0,'both complete111 comparisons')
    wlo,whi=power_interval(b['J_far_lower_interval'],W5)
    reflo,refhi=power_interval(b['J_far_refined_lower_interval'],Wref5)
    rlo,rhi=power_interval(b['intrinsic_lower_over_RU_interval'],W5/RU**5)
    ck(wlo>Q(5,10**6) and rlo>Q(104,100),'strict coarse obstruction margin')
    lower5,upper5=Q(81,2**80)/piU**3,Q(81,2**80)/piL**3
    ck(Q(b['accepted108_complete_actual_far_upper_fifth_input'])==upper5,'true upper coefficient uses pi lower')
    ulo,uhi=map(Q,b['accepted108_complete_actual_far_I_upper_interval'])
    ck(0<ulo<=uhi and ulo**5<=lower5<=upper5<=uhi**5,'complete coarse upper coefficient saved enclosure')
    ck(RU<ulo and uhi<Q(185,10**7),'coarse upper coefficient noncompetitive and finite')
    ck(cert['true_far_event']['probability_lower_coarse']==g['P_far_lower_coarse'] and cert['true_far_event']['probability_lower_refined']==g['P_far_lower_refined'],'event certificate copies')
    ck(cert['intrinsic_obstruction']['J_far_lower_fifth_exact']==b['J_far_lower_fifth_exact'] and cert['intrinsic_obstruction']['J_far_lower_interval']==b['J_far_lower_interval'],'obstruction certificate copies')
    ck(cert['complete_coarse_actual_far_upper']['interval']==b['accepted108_complete_actual_far_I_upper_interval'],'coarse bound certificate copy')
    ck(not b['actual_I_far_lower_bound'] and not b['actual_I_far_competitive'] and not out['full_far_competitive'],'J lower is not actual I lower')
    ck(not out['new_C6_bound'] and not b['global_g0_proved'] and inp['g0_proved'] is False,'global unpaid scope')
    for p,h in pins.items():ck(sha(Path(p))==h,'source unchanged after exact review')
    save('review_results.json',{'status':'PASS_ONE_INDEPENDENT_SAVED_EXACT_CONFIGURATION','assertions':count,'saved_configurations':1,'source_paths_verified_pre_post':len(pins),'author_upstream_sources':17,'author_whitelist_files':len(am['whitelist']),'exact_W5':str(W5),'exact_probability_lower':str(P0),'critical_J_fifth_coefficient':'3*pi^3*a/2^53','strict_residuals_vs_both_complete111_endpoints':True,'scope':{'actual_I_F_le_J_F':True,'J_F_lower_obstructs_only_specified_Jensen_family':True,'actual_I_F_lower_claim':False,'complete_coarse_actual_far_upper_valid':True,'competitive_far_globalg0_fullS_unpaid':True,'new_Cn':False},'author_evidence':{'science_runs':1,'completed_configs':2,'science_PID':71808,'observer_PID':71800,'child_exit_code':0,'complete_child_wall':obs['complete_science_wall_seconds'],'complete_child_CPU':obs['complete_science_CPU_seconds'],'observer_own_final_write_exit_excluded':True},'author_generators_old_functions_MC_native_kernel_grids_MATLAB_children':0,'source_hashes':pins,'wall_before_result_write':time.monotonic()-start,'CPU_before_result_write':time.process_time()-cpu})
    status='COMPLETED'
except Exception:
    failure=traceback.format_exc()
finally:
    save('run_receipt.json',{'status':status,'pid':os.getpid(),'assertions':count,'stage':stage,'failure':failure,'attempts':1,'automatic_retries':0,'children':0,'current_user_runtime_ceiling':None,'wall_before_final_receipt_write':time.monotonic()-start,'CPU_since_start_before_final_receipt_write':time.process_time()-cpu,'CPU_process_total_before_final_receipt_write':time.process_time(),'timing_scope':'one Python process before receipt/stdout/exit; complete tool roundtrip recorded separately','script_sha256':sha(Path(__file__)),'registration_sha256':sha(P/'registration.json')})
    print(json.dumps({'status':status,'assertions':count,'pid':os.getpid()}))
if status!='COMPLETED':raise SystemExit(1)
