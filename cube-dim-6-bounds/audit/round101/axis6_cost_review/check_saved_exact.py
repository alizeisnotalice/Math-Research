#!/usr/bin/env python3
"""Independent axes6 upper-family-cost audit; no author/old kernel execution."""
import sys
sys.dont_write_bytecode=True
import signal
signal.signal(signal.SIGALRM,signal.SIG_DFL)
signal.setitimer(signal.ITIMER_REAL,179)
import resource
resource.setrlimit(resource.RLIMIT_CPU,(15,15))
import os,json,hashlib,time,ast,difflib,subprocess
from pathlib import Path
from fractions import Fraction as Q
from math import comb
from datetime import datetime,timezone
START=time.monotonic();CPU=time.process_time()
OUT=Path(__file__).resolve().parent
AUTHOR=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-localmax-audit/cube-dim-6-bounds/upper/round101/axis6_cost_preflight')
HASHES={}

def need(ok,message):
    if not ok:raise ValueError(message)

def sha(p):
    p=Path(p);value=hashlib.sha256(p.read_bytes()).hexdigest()
    old=HASHES.setdefault(str(p),value);need(old==value,'Changing original '+str(p));return value

def read(p):sha(p);return json.loads(Path(p).read_text())

def save(name,data):
    with (OUT/name).open('x') as f:
        json.dump(data,f,indent=2,sort_keys=True);f.write('\n');f.flush();os.fsync(f.fileno())

def poly_mul(a,b):
    z=[Q(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):z[i+j]+=x*y
    return z

def check_fuse(text,seconds):
    tree=ast.parse(text);fuses=[];alarm_handlers=[]
    for node in ast.walk(tree):
        if not isinstance(node,ast.Call) or not isinstance(node.func,ast.Attribute):continue
        obj=ast.unparse(node.func.value);name=node.func.attr
        if obj=='signal' and name=='signal' and ast.unparse(node.args[0])=='signal.SIGALRM':
            alarm_handlers.append(ast.unparse(node.args[1]))
        if obj=='signal' and name=='setitimer':
            need(ast.unparse(node.args[0])=='signal.ITIMER_REAL','Only real total fuse')
            fuses.append(ast.literal_eval(node.args[1]))
        if obj=='signal' and name=='alarm':need(False,'No alarm cancellation/replacement allowed')
    need(alarm_handlers==['signal.SIG_DFL'] and fuses==[seconds],'Default hard fuse retained in entire source')
    return tree

def cpu_limits(tree):
    return sorted((n.lineno,ast.literal_eval(n.args[1])) for n in ast.walk(tree)
        if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and
        ast.unparse(n.func)=='resource.setrlimit' and ast.unparse(n.args[0])=='resource.RLIMIT_CPU')

def main():
    reg=read(OUT/'registration.json')
    need(sha(AUTHOR/'manifest.json')==reg['author_final_manifest_sha256'],'Final101 manifest identity')
    manifest=read(AUTHOR/'manifest.json')
    for item in manifest['artifacts']:
        p=AUTHOR/item['path'];need(p.parent==AUTHOR and p.stat().st_size==item['bytes'] and sha(p)==item['sha256'],'Scientific artifact binding')
    source_reg=read(AUTHOR/'registration.json');execution=read(AUTHOR/'execution_registration.json')
    need(len(source_reg['source_pins'])==22 and len(execution['execution_pins'])==9,'22source and9 execution pins')
    for p,digest in source_reg['source_pins'].items():need(sha(p)==digest,'Original source pin')
    for name,digest in execution['execution_pins'].items():need(sha(AUTHOR/name)==digest,'Original execution pin')
    need({x['path']:x['sha256'] for x in manifest['readonly_source_pins']}==source_reg['source_pins'] and
         all(x['matches'] is True for x in manifest['readonly_source_pins']),'Manifest source set complete')
    inputs=read(AUTHOR/'inputs.json');result=read(AUTHOR/'results.json')
    need(sha(AUTHOR/'inputs.json')==source_reg['inputs_sha256'],'Inputs manifest hash')
    need(result['self_sha256']==sha(AUTHOR/'check_exact.py'),'Original result/math source binding')
    history_path=next(p for p in source_reg['source_pins'] if p.endswith('/axis_frame_cost_probe.json'))
    history=read(history_path);history_ref=read(AUTHOR/'historical_axis_probe_reference.json')
    need(history['kind']=='numerical_candidate_only' and history['g_at_cube_corner_exact']==inputs['g_exact_inherited_corner_formula']=='3/(512*pi^5)',
         'Historical formula handoff only, not a new corner evaluation')
    need(history_ref['copied_fields']['kind']=='numerical_candidate_only' and history_ref['program_executed_now'] is False,'Historical numerical pilot not proof')
    pi_path=next(p for p in source_reg['source_pins'] if p.endswith('/audit/round098/direct_green_review/review_results.json'))
    inherited=read(pi_path);pi_lo,pi_hi=map(Q,inputs['pi_enclosure_inherited_E098'])
    need(list(map(Q,inherited['pi_interval']))==[pi_lo,pi_hi],'Exact inherited pi certificate bytes')
    P,Qpi=map(Q,inputs['pi_coarse_bounds']);need(P==Q(157,50)<pi_lo<pi_hi<Qpi==Q(22,7),'Coarse pi directions')
    q=Q(inputs['q']);need(q==Q(6,5) and 5/q==Q(25,6) and 1-(q-1)*5/q==Q(1,6),'Actual short-time norm exponents')
    need(Q(1,1)/(1+Q(1,6))==Q(6,7) and Q(25,6)>4 and P**11>8**6 and P>3,'Both full-time norm endpoint gates')
    # Only inherited Gamma recurrence coefficients are audited, no p_t or
    # G(corner), quadrature, Sobol point or old axis program is evaluated.
    beta_pi=Q(1,2)*Q(15,8)*Q(3,4)/120
    need(beta_pi==Q(3,512),'Inherited Beta/pi normalization')
    face_factor=6*64;radial_factor=1-Q(1,6);mass_factor=face_factor*radial_factor
    need(face_factor==384 and mass_factor==320 and mass_factor/64==5,'Actual full-space positive-part mass factors')
    # Fifth power of the spatial cutoff factor is5^5/[s(1-s)^5].
    # Verify derivative of its denominator factorizes as(1-s)^4(1-6s).
    denominator=poly_mul([0,1],[Q((-1)**k*comb(5,k)) for k in range(6)])
    derivative=[(k+1)*denominator[k+1] for k in range(len(denominator)-1)]
    factorization=poly_mul([Q((-1)**k*comb(4,k)) for k in range(5)],[1,-6])
    need(derivative==factorization,'Unique all-spatial-cutoff optimizer sign')
    s=Q(1,6);need(5**5/(s*(1-s)**5)==6**6,'True optimal family coefficient')
    need(-1/(5*s)+1/(1-s)==0 and 5*(q-1)==1,'Cost minimum atc=g/6 and correct homogeneity')
    # Critical independent exact residuals use the correlated exact g to cancel
    # pi before its lower bound. They never infer true I from an upper budget.
    a,b=map(Q,inputs['fixed_mean_window']);need((a,b)==(Q(1,4),Q(1,2)) and 2**2>3,'Fixed atan>pi/3 window')
    argument=(1+b*b)/(1+a*a);need(argument==Q(20,17),'Correct positive-time mean primitive')
    x=argument-1;log_lower=x-x*x/2;need(log_lower==Q(93,578)>0,'Strict analytic log lower coefficient')
    #1-(1-x)(1+x)=x^2 verifies the derivative residual sign.
    need(poly_mul([1,-1],[1,1])==[1,0,-1],'Analytic log polynomial residual')
    lower=read(AUTHOR/'results.json')['true_lower_certificate'];upper=result['true_upper_certificate']
    mean_lower=log_lower/(486*Qpi);Rlower=Q(512,243)*P**4*log_lower
    lower_residual=Rlower**6-Q(66)**5
    need(Q(lower['mean_A_lower'])==mean_lower and Q(lower['six_B_over_g_lower'])==Rlower and
         Q(lower['power_residual_R6_minus66power5'])==lower_residual>0,'True mean/Jensen family cost exceeds66')
    Bupper=Q(3,448)+Q(1,3645);Rupper=8*Qpi**5;upper_residual=Q(12000)**5-Rupper**6
    need(Q(upper['B_rational_upper'])==Bupper<Q(1,128) and Q(upper['six_B_over_g_upper'])==Rupper and
         Q(upper['power_residual_12000power5_minusR6'])==upper_residual>0,'Whole-face/all-time true cost upper')
    I_lo,I_hi=map(Q,inputs['fixed_face_I_targets'])
    I_lower_residual=mean_lower**6-I_lo**5;I_upper_residual=I_hi**5-Q(1,128)**6
    need((I_lo,I_hi)==(Q(1,100000),Q(1,300)) and min(I_lower_residual,I_upper_residual)>0,'Complete actual I interval')
    need(result['actual_model_cost_strict_interval']==['66','12000'] and result['actual_face_I_strict_interval']==['1/100000','1/300'] and
         result['actual_C6_lower_bound_claimed'] is False and result['new_C6_improvement'] is False and
         result['arbitrary_axes_kernel_or_anisotropic_weights_excluded'] is False,'Upper-family obstruction only')
    # Exact limited guard repair and direct original receipts. Neither guard nor
    # scientific child is imported/executed by this reviewer.
    old=(AUTHOR/'run_guard_once.py').read_text();new=(AUTHOR/'run_guard_v2.py').read_text()
    transformed=old.replace('resource.setrlimit(resource.RLIMIT_CPU,(2,2))','resource.setrlimit(resource.RLIMIT_CPU,(2,30))',1)
    transformed=transformed.replace("run_receipt.json","run_receipt_v2.json")
    needle='stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,start_new_session=True)\n'
    need(transformed.count(needle)==1,'Unique owned child launch')
    transformed=transformed.replace(needle,needle+'resource.setrlimit(resource.RLIMIT_CPU,(2,2))\n')
    need(transformed==new,'Only CPU-inheritance and receipt-name repair')
    diff=''.join(difflib.unified_diff(old.splitlines(keepends=True),new.splitlines(keepends=True),
        fromfile='run_guard_once.py',tofile='run_guard_v2.py'))
    need(diff==(AUTHOR/'guard_v2_diff.patch').read_text(),'Exact guard repair diff')
    oldtree=check_fuse(old,177);newtree=check_fuse(new,177)
    child=(AUTHOR/'check_exact.py').read_text();childtree=check_fuse(child,150)
    need([x[1] for x in cpu_limits(oldtree)]==[(2,2)] and [x[1] for x in cpu_limits(newtree)]==[(2,30),(2,2)] and
         [x[1] for x in cpu_limits(childtree)]==[(28,28)],'Correct hardCPU inheritance scope')
    popens=[n for n in ast.walk(newtree) if isinstance(n,ast.Call) and ast.unparse(n.func)=='subprocess.Popen']
    need(len(popens)==1 and cpu_limits(newtree)[0][0]<popens[0].lineno<cpu_limits(newtree)[1][0] and
         any(k.arg=='start_new_session' and ast.literal_eval(k.value) is True for k in popens[0].keywords),
         'Owned PG and parentCPU lowering after spawn')
    need('os._exit(' in new and 'communicate(timeout=160)' in new and 'communicate(timeout=5)' in new,'Bounded wait/cleanup/final exit source')
    continuation=read(AUTHOR/'continuation_registration.json');release=read(AUTHOR/'continuation_release.json')
    need(release['status']=='PARENT_STATIC_REVIEW_RELEASED_ONE_CORRECTED_GUARD_ATTEMPT' and
         release['primary_checked_source_pins']==22 and release['primary_checked_execution_pins']==9 and
         release['primary_readonly_check_chunk']=='74bf96' and release['permitted_corrected_attempts']==1 and
         release['upper_task_only'] is True and release['no_lower_bound_work'] is True,'Actual primary release')
    need(sha(AUTHOR/'run_guard_v2.py')==release['new_guard_sha256']==continuation['new_guard_sha256'] and
         sha(AUTHOR/'check_exact.py')==release['unchanged_math_script_sha256']==continuation['original_math_program_unchanged'],
         'Released exact guard/math bytes')
    for name,digest in continuation['preserved_failed_attempt'].items():
        if name.endswith(('.py','.json')):need(sha(AUTHOR/name)==digest,'Preserved first failed attempt bytes')
    first=read(AUTHOR/'run_receipt.json');second=read(AUTHOR/'run_receipt_v2.json')
    firsttool=read(AUTHOR/'tool_execution_receipt.json');secondtool=read(AUTHOR/'tool_execution_receipt_v2.json')
    need(first['exit_code']==firsttool['exit_code']==1 and first['child_reaped'] is True and first['stdout']=='' and
         'ValueError: not allowed to raise maximum limit' in first['stderr'] and 'line 14' in first['stderr'],
         'First setup failure before mathematical try')
    need(manifest['first_failed_entry_frozen']['mathematical_configurations']==0 and manifest['first_failed_entry_frozen']['mathematical_assertions']==0,
         'No math negative evidence from CPU setup error')
    need(second['exit_code']==secondtool['exit_code']==0 and second['child_reaped'] is True and second['stderr']=='' and
         second['hard_fuse_still_armed'] is True and second['automatic_retries']==0,'Unique corrected attempt complete')
    printed=json.loads(second['stdout']);need(printed['status']==result['status'] and printed['assertions']==result['assertions']==77 and
         len(result['configurations'])==3 and all(c['pass'] is True for c in result['configurations']),'Actual3-config saved result/receipt binding')
    need(secondtool['authorized_corrected_attempt']==1 and secondtool['session_id'] is None and
         Q(str(secondtool['complete_tool_roundtrip_wall_seconds']))==Q(91,1000)<180 and
         'Date.now' in secondtool['roundtrip_clock'],'Independent complete tool-roundtrip scope')
    known=[first['guard_pid'],first['owned_child_pid'],second['guard_pid'],second['owned_child_pid']]
    process=subprocess.run(['/bin/ps','-p',','.join(map(str,known)),'-o','pid=,stat=,command='],capture_output=True,text=True,timeout=5)
    need(process.returncode in (0,1) and process.stdout.strip()=='','A listed original task PID still present')
    for p,digest in list(HASHES.items()):need(hashlib.sha256(Path(p).read_bytes()).hexdigest()==digest,'Original changed during audit')
    review=dict(status='PASS_INDEPENDENT_ACTUAL_AXIS6_UPPER_FAMILY_COST_OBSTRUCTION',author_manifest_sha256=sha(AUTHOR/'manifest.json'),
        source_hashes=HASHES,source_pin_count=22,execution_pin_count=9,
        beta_pi_coefficient=str(beta_pi),corner_numeric_evaluations=0,global_cube_floor='3/(512*pi^5) via analytic monotonicity and corner continuity',
        full_cube_boundary_factor=384,radial_positive_part_mass_factor=320,K_mass_factor=5,
        optimal_c_over_g='1/6',minimum_spatial_factor_fifth_power=str(6**6),
        true_cost_strict_interval=['66','12000'],true_entire_face_I_strict_interval=['1/100000','1/300'],
        correlated_ratio_lower=str(Rlower),cost_lower_exact_positive_residual=str(lower_residual),
        ratio_upper=str(Rupper),cost_upper_exact_positive_residual=str(upper_residual),
        I_lower_exact_positive_residual=str(I_lower_residual),I_upper_exact_positive_residual=str(I_upper_residual),
        actual_operator_C6_lower_claim=False,upper_construction_family_only=True,
        arbitrary_axes_kernel_or_weight_families_excluded=False,general_L1_all_common_radii_transfer_scope='analytic proof review with axes C=pi^-6 and norm-comparable symbol',
        first_setup_error_preserved_zero_math=True,guard_v2_exact_limited_diff=True,original_guard_or_author_math_executions=0,
        original_outer_default_hard_fuse_seconds=177,original_child_default_hard_fuse_seconds=150,
        fuses_statically_retained_through_last_receipt_and_exit=True,actual_timeout_observed=False,
        original_corrected_complete_tool_roundtrip_seconds=0.091,
        original_corrected_guard_timer_before_final_receipt=second['wall_seconds_before_final_receipt_write'],
        first_raw_tool_wait_field=firsttool['exec_command_complete_wall_seconds'],first_raw_field_not_established_complete_process_wall=True,
        all_listed_original_PIDs_absent=known,ps_exit_code=process.returncode,signals_sent=0,all_descendants_reaped_claim=False,
        old_Sobol_65536_replays=0,old_96_time_replays=0,old_corner_numeric_parameter_replays=0,
        author_three_configuration_replays=0,kernel_calls=0,MATLAB_calls=0,source_writes=0,new_C6_bound=False)
    save('review_results.json',review)
    save('run_receipt.json',dict(status='PASS',pid=os.getpid(),exit_code=0,utc=datetime.now(timezone.utc).isoformat(),
        wall_seconds_including_result_save=time.monotonic()-START,cpu_seconds=time.process_time()-CPU,
        default_hard_SIGALRM_seconds=179,cpu_hard_limit_seconds=15,alarm_kept_through_receipt_and_exit=True,
        saved_exact_group=1,automatic_retries=0,original_guard_or_scientific_executions=0,result_sha256=sha(OUT/'review_results.json')))
    print(json.dumps(dict(status=review['status'],source_pins=22,execution_pins=9,strict_family_cost=['66','12000'],
        true_I=['1/100000','1/300'],original_guard_runs=0,actual_C6_lower_claim=False)))

if __name__=='__main__':
    try:main()
    except BaseException as error:
        save('failure_receipt.json',dict(status='STOP_NO_RETRY',error=repr(error),pid=os.getpid(),
            wall_seconds=time.monotonic()-START,cpu_seconds=time.process_time()-CPU,author_runs=0,old_kernel_replays=0))
        raise
