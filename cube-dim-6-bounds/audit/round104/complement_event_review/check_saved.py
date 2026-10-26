#!/usr/bin/env python3
"""A104 source/saved-budget audit; MAIN STATIC RELEASE REQUIRED before use."""
import signal
signal.signal(signal.SIGALRM,signal.SIG_DFL)
signal.setitimer(signal.ITIMER_REAL,179)
import resource
resource.setrlimit(resource.RLIMIT_CPU,(20,20))
import sys
sys.dont_write_bytecode=True
import os,json,hashlib,time,ast
from pathlib import Path
from fractions import Fraction as Q
from math import comb
from datetime import datetime,timezone
START=time.monotonic();CPU=time.process_time()
OUT=Path(__file__).resolve().parent
AUTHOR=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-upper098/cube-dim-6-bounds/upper/round104/complement_event_budget')
BASE=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds')
HASHES={}

def need(ok,message):
    if not ok:raise ValueError(message)
def sha(p):
    p=Path(p);value=hashlib.sha256(p.read_bytes()).hexdigest()
    old=HASHES.setdefault(str(p),value);need(old==value,'Changing input '+str(p));return value
def read(p):sha(p);return json.loads(Path(p).read_text())
def save(name,value):
    with (OUT/name).open('x') as f:
        json.dump(value,f,indent=2,sort_keys=True);f.write('\n');f.flush();os.fsync(f.fileno())
def I(raw):
    need(set(raw)=={'lower_mantissa','upper_mantissa','common_binary_exponent'},'Exact dyadic object required')
    q=Q(2)**int(raw['common_binary_exponent']);a=int(raw['lower_mantissa'])*q;b=int(raw['upper_mantissa'])*q
    need(a<=b,'Dyadic endpoint order');return a,b
def add(a,b):return a[0]+b[0],a[1]+b[1]
def sub(a,b):return a[0]-b[1],a[1]-b[0]
def mul(a,b):
    values=[x*y for x in a for y in b];return min(values),max(values)
def scale(a,n):return mul(a,(Q(n),Q(n)))
def div(a,b):need(b[0]>0,'Positive denominator');return mul(a,(1/b[1],1/b[0]))
def power(a,n):need(a[0]>=0,'Nonnegative power base');return a[0]**n,a[1]**n
def overlaps(a,b):return max(a[0],b[0])<=min(a[1],b[1])
def bind(a,b,message):need(overlaps(a,b),message)
def text(a):return list(map(str,a))

def main():
    # These checks precede any certificate arithmetic. The primary writes this
    # release only after approving the exact scripts and registered inputs.
    approval=read(OUT/'primary_release.json');registration=read(OUT/'registration.json')
    need(approval['authorized'] is True and approval['round']=='104-independent-review' and
         approval['review_registration_sha256']==sha(OUT/'registration.json') and
         approval['script_sha256'][Path(__file__).name]==sha(Path(__file__)), 'No exact main static release')
    for name,digest in registration['input_sha256'].items():need(sha(AUTHOR/name)==digest,'Registered original input')
    need(sha(AUTHOR/'manifest.json')==approval['author_final_manifest_sha256']==registration['author_final_manifest_sha256'],'Author final seal')
    manifest=read(AUTHOR/'manifest.json')
    rows=manifest.get('whitelist',manifest.get('files',manifest.get('artifacts',[])));need(bool(rows),'No final artifact whitelist')
    for row in rows:
        name=row.get('name',row.get('path'));p=AUTHOR/name
        need(p.resolve().is_relative_to(AUTHOR.resolve()) and '..' not in Path(name).parts,'Manifest path escaped author scope')
        need(sha(p)==row['sha256'],'Final author artifact hash')
        if 'bytes' in row:need(p.stat().st_size==row['bytes'],'Final author bytes')
    sources=read(AUTHOR/'source_manifest.json')
    for row in sources['sources']:need(sha(row['path'])==row['sha256'],'Readonly source pin')
    data=read(AUTHOR/'results.json');receipt=read(AUTHOR/'run_receipt.json');release=read(AUTHOR/'parent_release.json')
    old=read(AUTHOR/'upstream_constants.json')
    need(data['status']=='PASS_FOUR_CONFIGURATIONS_SOURCE_PROOF_REVIEW_SEPARATE' and
         len(data['configurations'])==receipt['completed_configurations']==4,'Actual producer completion')
    need(receipt['pid']==42765,'Original actual PID binding')
    need(receipt['status']=='COMPLETED' and receipt['failure'] is None and receipt['retries']==0 and receipt['children_created']==0 and
         receipt['SIGALRM_handler_is_default'] is True and receipt['fuse_not_cancelled'] is True and
         Q(str(receipt['fuse_remaining_before_final_receipt_write']))>0,'Actual no-retry single-process/default-fuse receipt')
    need(release['authorized'] is True and release['script_sha256']==data['source_code_sha256']==sha(AUTHOR/'run_prepared.py') and
         release['registration_sha256']==data['registration_sha256']==sha(AUTHOR/'registration.json') and
         receipt['release_sha256']==sha(AUTHOR/'parent_release.json'),'Unique original release binding')
    need(data['upstream_constants_sha256']==sha(AUTHOR/'upstream_constants.json') and
         data['new_native_kernel_peak_grid_replays']==data['children']==0 and data['Cn_improved'] is False,'Original scope flags')
    tree=ast.parse((AUTHOR/'run_prepared.py').read_text());handlers=[];timers=[]
    for node in ast.walk(tree):
        if not isinstance(node,ast.Call):continue
        name=ast.unparse(node.func)
        if name=='signal.signal' and ast.unparse(node.args[0])=='signal.SIGALRM':handlers.append(ast.unparse(node.args[1]))
        if name=='signal.setitimer':timers.append((ast.unparse(node.args[0]),ast.literal_eval(node.args[1])))
        need(name not in ('signal.alarm','os.system','os.fork','subprocess.Popen','subprocess.run'),'Producer must retain its only default fuse/no child')
    need(handlers==['signal.SIG_DFL'] and timers==[('signal.ITIMER_REAL',177.0)],'Full producer source has only retained default177 fuse')
    tool=read(AUTHOR/'execution_evidence/actual_tool_receipt.json');gate=read(AUTHOR/'execution_evidence/release_receipt.json')
    need(tool['tool_result']['exit_code']==0 and tool['tool_result']['chunk_id']=='d5b280' and
         tool['scientific_invocations']==1 and tool['native_calls']==0,'Original main tool execution')
    toolout=json.loads(tool['tool_result']['output'])
    need(toolout['status']=='COMPLETED' and toolout['completed_configurations']==4 and toolout['pid']==42765,
         'Original completion output binding')
    need(0<Q(str(tool['full_tool_roundtrip_seconds']))<180 and gate['only_hash_and_AST_static_checks'] is True and
         gate['science_calls']==0 and gate['release_sha256']==sha(AUTHOR/'parent_release.json'),'Original tool timing/release scope')
    cfg={x['configuration']:x for x in data['configurations']};c1,c2,c3,c4=(cfg[k] for k in (1,2,3,4))
    # Analytic Hamming/coarea combinatorics, not a replay of a pair/fan baseline.
    partner=[comb(5,1)+comb(5,5),comb(5,2)+comb(5,4),comb(5,3)]
    full=[32*x for x in partner];removed=[2*x for x in partner];counts=[a-b for a,b in zip(full,removed)]
    need(partner==[6,15,10] and counts==c1['ordered_pair_counts']==[180,450,300] and sum(counts)==31*30==930,'All930 ordered pairs')
    need([6**2-dot**2 for dot in (4,2,0)]==c1['J2_squared']==[20,32,36],'True two-column coarea squares')
    need(c1['Crem_missing_peak_loss_and_third_deletion_used'] is False and c1['actual_full_symbol_coercivity_reused_from098']=='32/sqrt3',
         'Restored actual full quotient, no omission/triple loss')
    need(Q(512)*3/32**3==Q(3,64),'Restored D coefficient normalization')
    need(Q(7,40)-Q(9,256)>Q(1,64),'Uniform positive lower for every k pair distance')
    need([x['k'] for x in c1['all_pair_gaps']]==list(range(1,6)) and
         all(I(x['all930_pair_distance_lower'])[0]>Q(1,64) for x in c1['all_pair_gaps']),'All recorded pair-distance enclosures strictly outside every other tube')
    # Bind actual old scalar certificates only; no051/074 function/grid run.
    source100=read(AUTHOR.parents[1]/'round100/time_first_budget_refinement/results.json')
    source062=read(BASE/'upper/round062/primary_joint_account/actual_run/results.json')
    source074=read(BASE/'upper/round074/middle_long_budget/results.json')
    need(old['norm100_core']==source100['configurations'][2]['new_full_G_norm'] and
         c2['large_branch_norm']==source100['configurations'][2]['Pnorm'],'Identical old100 core/P object')
    need(Q(old['paid_CEH_I_upper'])==Q(source062['paid_C_E_H_I_upper']),'Paid CEH is actual full-G^q allocation')
    need(Q(old['long074_rational_upper'])==Q(source074['strict_long_norm_rational_upper']) and
         Q(old['short051_rational_upper'])==Q(source074['comparison']['old_short_upper_read_only']),
         'Identical old short/long middle budgets')
    paid=Q(old['paid_CEH_I_upper']);g0=Q(1,65000);target5=64**5*(g0/6)**6
    P,N,B=(I(c2[k]) for k in ('large_branch_norm','full_G_core_norm','complement_norm'))
    bind(N,add(P,B),'Core full-time spatial Minkowski')
    bind(power(I(c2['I_one_core']),5),power(N,6),'Core norm-to-I relation')
    bind(I(c2['I_ten_cores']),scale(I(c2['I_one_core']),10),'Ten disjoint core allocations')
    ratio=I(c2['ratio_to_saved100_same_core_norm']);bind(ratio,div(N,(I(old['norm100_core'])[1],)*2),'Same-core norm ratio denominator is saved upper')
    need(Q('0.1925')<ratio[0]<=ratio[1]<Q('0.1927')<1,'Local core budget improvement only')
    core_plan=add((paid,paid),I(c2['I_ten_cores']));bind(core_plan,I(c2['paid_CEH_plus_core_plan_I']),'CEH plus cores')
    need(core_plan[1]**5<target5 and I(c2['S_paid_CEH_plus_core_other_unpaid_zero'])[0]>0 and
         c2['operator_breakthrough'] is False and 'UNPAID' in c2['global_I_status'],'Positive core remainder is partial, not full S')
    rows=c3['cases'];need([x['k'] for x in rows]==list(range(1,6)) and
         [x['multiplicity'] for x in rows]==[comb(5,k-1) for k in range(1,6)] and sum(x['multiplicity'] for x in rows)==31,'All31 middle domains')
    total=(Q(0),Q(0))
    for row in rows:
        norm=I(row['full_norm_one']);bind(norm,add(I(row['Pnorm']),I(row['Bnorm'])),'Every type pays P+B')
        bind(power(I(row['weighted_I']),5),scale(power(norm,6),Q(row['multiplicity'])**5),'Multiplicity before summing actual same-domain I budgets')
        total=add(total,I(row['weighted_I']))
    newI=I(c3['new_middle_I_upper']);oldI=I(c3['old051short_plus074long_same_middle_I_upper_expression'])
    bind(total,newI,'Complete middle positive sum')
    old_norm=Q(old['short051_rational_upper'])+Q(old['long074_rational_upper'])
    need(power(oldI,5)[0]<=old_norm**6<=power(oldI,5)[1],'Old comparison is same full G short+long norm before q')
    all_ratio=I(c3['new_over_old_middle_expression']);bind(all_ratio,div(newI,oldI),'Same full-middle expression ratio')
    need(Q('4.74')<all_ratio[0]<=all_ratio[1]<Q('4.76') and newI[0]>oldI[1],'Full middle expression is worse')
    plan=add((paid,paid),newI);bind(plan,I(c3['paid_CEH_plus_new_middle_I_plan']),'New partial face I plan')
    need(plan[0]**5>target5 and I(c3['S_new_plan_far_optimistically_zero'])[1]<0 and
         c3['actual_I_or_C6_lower_bound'] is False and c3['g0_status']=='unproved globally' and 'UNPAID' in c3['far_status'],
         'Negative S expression only; far/g unpaid')
    need(c4['three_nonparallel_normal_contributions_sum_zero'] and c4['no_two_normal_contributions_cancel'],'Original source-witness flags')
    for p,digest in list(HASHES.items()):need(hashlib.sha256(Path(p).read_bytes()).hexdigest()==digest,'Postcheck evidence mutation')
    review=dict(status='PASS_SAVED_BINDINGS_RESTORED_PAIR_ACCOUNT_AND_LOCAL_GLOBAL_SCOPE',source_hashes=HASHES,
        author_final_manifest_sha256=sha(AUTHOR/'manifest.json'),ordered_pair_counts=counts,coarea_J_squared=[20,32,36],ordered_pairs=930,
        source_restoration_precedes_peak_and_time=True,omission_Crem_and_third_deletion_not_used=True,
        core_ratio=text(ratio),core_partial_I_plan=text(core_plan),core_remaining_allowance_is_not_global_S=True,
        full_middle_multiplicities=[1,5,10,10,5],full_middle_ratio=text(all_ratio),full_middle_I_allocation=text(newI),
        full_middle_same_G_source_events_and_all_time_paid=True,far_unpaid=True,global_g_unproved=True,complete_I_or_S_proved=False,
        positive_part_I_interface='full original G^q on the entire face is required by mass of(G-c)+; local budgets and CEH only do not close it',
        actual_I_model_cost_C6_lower_claim=False,new_C6_bound=False,
        original_main_tool_chunk='d5b280',original_main_PID=42765,original_main_exit_code=0,
        original_tool_roundtrip_seconds=tool['full_tool_roundtrip_seconds'],
        original_internal_before_final_receipt_seconds=receipt['elapsed_before_final_receipt_write'],
        original_default177_fuse_retained=True,original_single_process_no_children=True,
        original_startup_before_timer_not_independently_hard_capped_by_entry=True,
        author_generator_runs=0,old_configurations_or_functions_replayed=0,grids=0,physical_kernel_calls=0,MATLAB=0,source_writes=0)
    save('saved_results.json',review)
    save('saved_run_receipt.json',dict(status='PASS',pid=os.getpid(),exit_code=0,utc=datetime.now(timezone.utc).isoformat(),
        wall_seconds_including_result_save=time.monotonic()-START,cpu_seconds=time.process_time()-CPU,
        default_hard_SIGALRM_seconds=179,CPU_hard_limit_seconds=20,alarm_kept_through_receipt_and_exit=True,
        saved_exact_group=1,automatic_retries=0,author_generator_runs=0,result_sha256=sha(OUT/'saved_results.json')))
    print(json.dumps(dict(status=review['status'],ordered_pairs=930,core_ratio_interval=text(ratio),full_middle_ratio_interval=text(all_ratio),
        complete_g_I_S=False,author_generator_runs=0)))

if __name__=='__main__':
    try:main()
    except BaseException as error:
        save('saved_failure_receipt.json',dict(status='STOP_NO_RETRY',error=repr(error),pid=os.getpid(),
            wall_seconds=time.monotonic()-START,cpu_seconds=time.process_time()-CPU,author_generator_runs=0))
        raise
