#!/usr/bin/env python3
"""A100 saved endpoints/source audit, not a four-configuration replay."""
import sys
sys.dont_write_bytecode=True
import signal
signal.signal(signal.SIGALRM,signal.SIG_DFL)
signal.setitimer(signal.ITIMER_REAL,179.5)
import resource
resource.setrlimit(resource.RLIMIT_CPU,(20,20))
import os,json,hashlib,time,subprocess
from pathlib import Path
from fractions import Fraction as Q
from datetime import datetime,timezone
START=time.monotonic();CPU=time.process_time()
OUT=Path(__file__).resolve().parent
AUTHOR=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-upper098/cube-dim-6-bounds/upper/round100/time_first_budget_refinement')
BASE=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds')
HASHES={}

def need(ok,message):
    if not ok:raise ValueError(message)

def sha(p):
    p=Path(p);value=hashlib.sha256(p.read_bytes()).hexdigest()
    old=HASHES.setdefault(str(p),value);need(old==value,'Changing source '+str(p));return value

def read(p):sha(p);return json.loads(Path(p).read_text())

def save(name,data):
    with (OUT/name).open('x') as f:
        json.dump(data,f,indent=2,sort_keys=True);f.write('\n');f.flush();os.fsync(f.fileno())

def interval(raw):
    need(set(raw)=={'lower_mantissa','upper_mantissa','common_binary_exponent'},'Dyadic shape')
    a,b,e=(int(raw[k]) for k in ('lower_mantissa','upper_mantissa','common_binary_exponent'))
    need(a<=b and abs(e)<10000,'Dyadic order/exponent');return a*Q(2)**e,b*Q(2)**e

def add(a,b):return a[0]+b[0],a[1]+b[1]
def sub(a,b):return a[0]-b[1],a[1]-b[0]
def mul(a,b):
    vals=[x*y for x in a for y in b];return min(vals),max(vals)
def div(a,b):
    need(b[0]>0,'Positive denominator needed');return mul(a,(1/b[1],1/b[0]))
def scale(a,q):return mul(a,(Q(q),Q(q)))
def power(a,n):
    need(a[0]>0,'Positive power base');return a[0]**n,a[1]**n
def overlaps(a,b):return max(a[0],b[0])<=min(a[1],b[1])
def binding(a,b,message):need(overlaps(a,b),message)
def text(a):return list(map(str,a))

def det3(g):
    return g[0][0]*(g[1][1]*g[2][2]-g[1][2]*g[2][1])-g[0][1]*(g[1][0]*g[2][2]-g[1][2]*g[2][0])+g[0][2]*(g[1][0]*g[2][1]-g[1][1]*g[2][0])

def main():
    reg=read(OUT/'registration.json')
    need(sha(AUTHOR/'manifest.json')==reg['author_final_manifest_sha256'],'Final author seal')
    for item in read(AUTHOR/'manifest.json')['files']:
        p=AUTHOR/item['name'];need(p.parent==AUTHOR and p.stat().st_size==item['bytes'] and sha(p)==item['sha256'],'Author artifact binding')
    source=read(AUTHOR/'source_manifest.json')
    for item in source['source_records']:
        p=Path(item['path']);need(p.stat().st_size==item['bytes'] and sha(p)==item['sha256'],'Frozen source binding')
    data=read(AUTHOR/'results.json');certificate=read(AUTHOR/'certificate.json');receipt=read(AUTHOR/'run_receipt.json')
    need(data['status']=='PASS_FOUR_FINITE_CONFIGURATIONS_NO_GLOBAL_BOUND' and len(data['configurations'])==4,'Original result status')
    need(data['native_G_H_p_calls']==data['old_program_replays']==data['time_grid']==data['matlab_calls']==0 and data['Cn_improved'] is False,'Original no-replay/no-bound scope')
    need(receipt['status']=='COMPLETED' and receipt['exit_code']==0 and receipt['reaped'] is True and receipt['retries']==0 and receipt['stderr']=='','Original actual completion')
    need(receipt['child_pid']==data['pid'] and certificate['results_sha256']==sha(AUTHOR/'results.json') and
         certificate['receipt_sha256']==sha(AUTHOR/'run_receipt.json'),'Original execution binding')
    configs={x['configuration']:x for x in data['configurations']};c1,c2,c3,c4=(configs[k] for k in (1,2,3,4))
    old=read(AUTHOR/'upstream_constants.json')
    old098=read(AUTHOR.parents[1]/'round098/direct_green_geometry/results.json')
    need(old['B098_coarse']==old098['configurations'][2]['B0']==c1['old_B0_source_dyadic'],'Original098 same B0 bytes')
    need(old['I098_tiny_Q']==old098['configurations'][3]['local_I_upper'],'Original098 tinyQ bytes')
    peaks=read(BASE/'upper/round033/exact_projection_peaks/results.json')
    cases={x['id']:x for x in peaks['cases']};labels=['4110','3210','2220','3111','2211']
    counts=[435,1740,435,580,1305];full=[480,1920,480,640,1440];dets=[64,96,128,160,192]
    need(old['excluded_body_unordered_triple_counts']==counts and old['Gram_determinants']==dets,'Exact original class data')
    weighted=Q(0)
    for label,n,nfull,det in zip(labels,counts,full,dets):
        case=cases['r3_'+label];need(case['r']==case['dimension']==3 and case['Gram_determinant']==det,'Original033 quotient dimension/determinant')
        need(Q(old['projection_R'][label])==Q(case['R_exact']),'Original033 exact R matching')
        rows=case['normal_rows'];gram=[[sum(x*y for x,y in zip(a,b)) for b in rows] for a in rows]
        need(det3(gram)==det and 32*n==29*nfull,'True coarea and omission class count')
        weighted+=n*Q(case['R_exact'])
    need(6*sum(counts)==26970 and 24576*weighted==Q(44743079,2970)==Q(old['source_Crem_times_pi6']),
         'Complete ordered triple coefficient arithmetic')
    Bfull=interval(c1['B_middle_full']);B0=interval(old['B098_coarse'])
    need(Bfull[1]<B0[0],'New uniform middle expression strictly smaller')
    binding(interval(c1['ratio_to_saved098_B0']),div(Bfull,(B0[0],B0[0])),'Uniform expression ratio')
    d3={k:interval(v) for k,v in c3.items() if isinstance(v,dict) and 'common_binary_exponent' in v}
    P,B,N=(d3[k] for k in ('Pnorm','joint_Bnorm','new_full_G_norm'))
    binding(N,add(P,B),'Minkowski sum endpoints')
    binding(power(d3['new_I_one'],5),power(N,6),'q6/5 norm-to-integral relationship')
    binding(d3['new_I_ten'],scale(d3['new_I_one'],10),'Ten same-domain upper allocations')
    g0=Q(1,65000);gate=Q(4,48750);target=d3['conditional_whole_I_target'];ten=d3['new_I_ten']
    need(64**5*(g0/6)**6==gate**6,'Exact conditional B gate')
    need(power(target,5)[0]<=gate**6<=power(target,5)[1],'Target I exact algebraic endpoint enclosure')
    derivedS=sub(target,ten)
    binding(d3['S_ten_optimistic'],derivedS,'Negative term uses lower endpoint for upper difference')
    need(derivedS[1]<0 and d3['S_ten_optimistic'][1]<0 and ten[0]>target[1], 'Constructed budget fails robust endpoint comparison')
    binding(d3['U_ten_optimistic'],scale(div(ten,target),64),'Conditional U_plan relationship')
    old_same=d3['old098_saved_B0_same_domain_norm']
    need(N[1]<old_same[0],'Same R budget improvement')
    binding(d3['new_over_098_saved_B0'],div(N,old_same),'Same R expression ratio')
    need(Q('0.0717')<d3['new_over_098_saved_B0'][0] and d3['new_over_098_saved_B0'][1]<Q('0.0718'),'7.17 percent display scope')
    binding(d3['background_fraction_of_norm'],div(B,N),'Complement budget fraction')
    need(Q('0.9995')<d3['background_fraction_of_norm'][0] and d3['background_fraction_of_norm'][1]<1,'99.95 percent of allocated norm')
    gamma=d3['Crem_factor_needed_even_if_other_regions_zero']
    necessary_norm=add(P,mul(gamma,B))
    need(power(necessary_norm,6)[0]<=gate**6/10**5<=power(necessary_norm,6)[1],'Gamma gate equation with saved P/B')
    need(Q('0.4301')<gamma[0] and gamma[1]<Q('0.4302')<1,'Necessary coefficient factor only under fixed allocation')
    newer,older=interval(c4['new_lawful_Crem_after_omission']),interval(c4['old_exact_peak_Crem_after_omission'])
    ratio=interval(c4['new_over_old']);binding(ratio,div(newer,older),'Lawful moment / exact-peak ratio')
    need(newer[0]>older[1] and Q('2.70')<ratio[0] and ratio[1]<Q('2.71') and c4['automatic_replacement_improves'] is False,
         'New moment coefficient is worse, not a budget improvement')
    need(c3['actual_I_or_C6_lower_bound'] is False and c3['g0_global_status']=='UNPAID' and c3['complementary_face_budget']=='UNPAID',
         'Optimistic budget failure is not actual model or C6 lower bound')
    outcomes=certificate['outcomes']
    for key,value in [('B_middle_full',c1['B_middle_full']),('new_local_norm',c3['new_full_G_norm']),('new_I_one',c3['new_I_one']),
                      ('new_I_ten',c3['new_I_ten']),('S_ten_plan',c3['S_ten_optimistic']),('U_ten_plan',c3['U_ten_optimistic']),
                      ('coercivity_replacement_ratio',c4['new_over_old'])]:need(outcomes[key]==value,'Certificate exact outcome binding')
    appendix=read(AUTHOR/'tool_runtime_addendum.json')
    need(sha(AUTHOR/'tool_runtime_addendum.json')=='ae3fdb3378ca2b61818755a0643e09433b510919ea3c7f7a792c617fe6ecc213','Main supplied administrative addendum binding')
    # Runtime semantics are reported from the entire sealed guard, not merely
    # from cancellation in the inner scientific entry. No outer guard is assumed.
    guard=(AUTHOR/'run_once.py').read_text()
    need(guard.index('signal.alarm(0)')<guard.index("receipt.write_text"),'Last receipt written after sealed guard alarm cancellation')
    known=[receipt['guard_pid'],receipt['child_pid']]
    process=subprocess.run(['/bin/ps','-p',','.join(map(str,known)),'-o','pid=,stat=,command='],capture_output=True,text=True,timeout=5)
    need(process.returncode in (0,1) and process.stdout.strip()=='','Original owned100 PID still present')
    for p,digest in list(HASHES.items()):need(hashlib.sha256(Path(p).read_bytes()).hexdigest()==digest,'Evidence changed during audit')
    result=dict(status='PASS_SAVED_BINDINGS_AND_CONSTRUCTED_BUDGET_COMPARISONS',source_hashes=HASHES,
        original_author_configurations_not_replayed=4,exact_original033_R_matching=True,ordered_triples=26970,
        exact_Crem_pi6=str(24576*weighted),same_R_ratio=text(d3['new_over_098_saved_B0']),
        complement_fraction_of_allocated_norm=text(d3['background_fraction_of_norm']),
        ten_region_constructed_I=text(ten),conditional_I_target=text(target),derived_S_plan=text(derivedS),
        necessary_Crem_factor_in_fixed_budget=text(gamma),new_moment_over_old_exact_peak=text(ratio),
        all_reported_ratios_are_upper_expression_comparisons=True,actual_I_or_C6_lower_claim=False,
        source_event_total_payment_requires_reviewed043_074_and098_proofs=True,
        author_actual_guard_measurement=receipt['guard_seconds_after_child_cleanup_before_final_receipt_write'],
        author_actual_guard_measurement_scope='after child cleanup, before final receipt write; not exact full process wall',
        author_full_cleanup_inclusive_hard_runtime_policy_certified=False,future_same_unfixed_generator_approved=False,
        complete_author_process_wall='unavailable',actual_timeout_evidence=False,administrative_addendum_sha256=sha(AUTHOR/'tool_runtime_addendum.json'),
        original_PIDs_observed_absent=known,ps_exit_code=process.returncode,ps_stdout=process.stdout,signals_sent=0,
        author_program_imports_or_executions=0,old062_formula_evaluations=0,old043074_replays=0,grids=0,kernel_calls=0,MATLAB=0,source_writes=0)
    save('saved_results.json',result)
    save('saved_run_receipt.json',dict(status='PASS',pid=os.getpid(),exit_code=0,utc=datetime.now(timezone.utc).isoformat(),
        wall_seconds_including_result_save=time.monotonic()-START,cpu_seconds=time.process_time()-CPU,
        hard_default_alarm_seconds=179.5,cpu_hard_limit_seconds=20,alarm_kept_through_receipt_and_exit=True,
        saved_exact_interface_group=1,automatic_retries=0,result_sha256=sha(OUT/'saved_results.json')))
    print(json.dumps(dict(status=result['status'],ordered_triples=26970,source_files_bound=len(HASHES),
        negative_S_endpoint_direction_verified=True,author_runtime_hard_policy_certified=False)))

if __name__=='__main__':
    try:main()
    except BaseException as error:
        save('saved_failure_receipt.json',dict(status='STOP_NO_RETRY',error=repr(error),pid=os.getpid(),
            wall_seconds=time.monotonic()-START,cpu_seconds=time.process_time()-CPU,author_program_runs=0))
        raise
