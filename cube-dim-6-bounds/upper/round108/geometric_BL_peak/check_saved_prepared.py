"""108 limited SAVED checker, PREPARED ONLY; independent proof review separate."""
import signal
signal.signal(signal.SIGALRM,signal.SIG_DFL)
signal.setitimer(signal.ITIMER_REAL,177.0)
import sys
sys.dont_write_bytecode=True
import datetime,hashlib,json,os,resource,time,traceback
from pathlib import Path
from fractions import Fraction as Q
resource.setrlimit(resource.RLIMIT_CPU,(25,25))
P=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
rpath=P/'checker_release.json'
if not rpath.exists():raise SystemExit('108 SAVED CHECK PREPARED_NOT_RELEASED')
rel=json.loads(rpath.read_text())
if not(rel.get('authorized')is True and rel.get('round')=='108_saved'
       and rel.get('script_sha256')==sha(Path(__file__))
       and rel.get('registration_sha256')==sha(P/'checker_registration.json')):
    raise SystemExit('108 saved checker release mismatch')
if any((P/n).exists()for n in ['checker_started.json','checked_saved.json','checker_receipt.json']):
    raise SystemExit('108 saved checker one attempt only')
start=time.monotonic();startutc=datetime.datetime.now(datetime.timezone.utc).isoformat()
(P/'checker_started.json').write_text(json.dumps({'pid':os.getpid(),'started_utc':startutc,
    'release_sha256':sha(rpath),'default177_fuse_active':True},indent=2)+'\n')
checks=[];status='FIRST_FAILURE_STOP';failure=None
def ends(d):
    e=int(d['common_binary_exponent'])
    l=Q(int(d['lower_mantissa']))*Q(2)**e
    u=Q(int(d['upper_mantissa']))*Q(2)**e
    assert l<=u
    return l,u
def overlap(a,b):return max(a[0],b[0])<=min(a[1],b[1])
try:
    cr=json.loads((P/'checker_registration.json').read_text())
    assert sha(P/'saved_check_bindings.json')==cr['bindings_sha256']
    bind=json.loads((P/'saved_check_bindings.json').read_text())
    for n,h in bind['files'].items():assert sha(P/n)==h,(n,'byte mismatch')
    checks.append('all pinned final local bytes unchanged')
    manifest=json.loads((P/'source_manifest.json').read_text())
    snap=json.loads((P/'administrative_source_snapshots.json').read_text())
    for row in snap['files']:
        assert sha(P/row['snapshot'])==row['asof_sha256']
        old=next(x for x in manifest['files']if x['path']==row['original_path'])
        assert old['sha256']==row['asof_sha256']
    checks.append('historical mutable admin snapshots match original source manifest hashes')
    result=json.loads((P/'results.json').read_text());receipt=json.loads((P/'run_receipt.json').read_text())
    release=json.loads((P/'parent_release.json').read_text());tool=json.loads((P/'actual_tool_receipt.json').read_text())
    assert result['status']=='PASS_THREE_CONFIGURATIONS_PROOF_AUDIT_SEPARATE'
    assert len(result['configurations'])==3 and receipt['completed_configurations']==3
    assert receipt['status']=='COMPLETED' and receipt['failure']is None
    assert receipt['children_created']==0 and receipt['retries']==0
    assert receipt['release_sha256']==sha(P/'parent_release.json')
    assert receipt['SIGALRM_handler_is_default']and receipt['fuse_never_cancelled']
    assert release['authorized']is True and release['max_configurations']==3
    assert tool['tool_chunk_id']=='5bdf65'and tool['exit_code']==0
    assert receipt['pid']==tool['observed_PID']==50992
    assert result['source_code_sha256']==sha(P/'run_prepared.py')==release['script_sha256']
    assert result['registration_sha256']==sha(P/'registration.json')==release['registration_sha256']
    assert result['source_manifest_sha256']==sha(P/'source_manifest.json')
    assert result['saved104_inputs_sha256']==sha(P/'saved104_inputs.json')
    checks.append('unique three-configuration producer release, PID, exit and retained default fuse')
    a,b,c=result['configurations']
    assert a['zero_projection_witness_histogram']=={'0':4,'2':12,'4':12,'6':4}
    assert a['zero_projection_witness_EZ']=='3'and a['zero_projection_witness_EZ2']=='12'
    assert len(a['nonzero_quadratic_moment_coefficients'])==36
    assert c['dimensions'][2]['peak_coefficient_including_pi']==b['Q3_new']
    assert Q(c['full6_pi_free_coefficient'])==Q(216,2**30)
    assert Q(c['full6_ratio_to_saved001_det160_peak'])==Q(135,4194304)
    assert Q(c['dimensions'][5]['exact_square_of_pi_free_peak_coefficient'])==Q(216,2**30)**2
    checks.append('saved finite-sign witness flags and exact r3/full6 constants')
    ll,lu=ends(b['ratio_to098']);sl,su=ends(b['conditional_sharpKhinchine_ratio_to098'])
    fl,fu=ends(b['integral_scaling_lambda_power6over5'])
    assert 0<ll<=lu<sl<=su<1
    assert sl**2<=Q(8,27)<=su**2
    assert fl>0 and fl**5<=lu**6 and fu**5>=ll**6
    inp=json.loads((P/'saved104_inputs.json').read_text())
    for src,dst,factor in [('core_full_norm','one_core_full_G_norm',(ll,lu)),
            ('core_ten_integral','ten_core_I_allocation',(fl,fu)),
            ('complete_middle_integral','complete31_middle_I_allocation',(fl,fu))]:
        lo,hi=ends(inp[src]); assert lo>0
        assert overlap(ends(b[dst]),(lo*factor[0],hi*factor[1]))
    paid=Q(inp['paid_CEH_I_upper']);target_power5=Q(64)**5*Q(1,390000)**6
    core_lo,core_hi=ends(inp['core_ten_integral'])
    middle_lo,middle_hi=ends(inp['complete_middle_integral'])
    safe_core_upper=paid+core_hi*fu
    safe_middle_lower=paid+middle_lo*fl
    assert safe_core_upper>0 and safe_core_upper**5<target_power5
    assert safe_middle_lower>0 and safe_middle_lower**5>target_power5
    assert ends(b['conditional_core_remaining_allowance'])[0]>0
    assert ends(b['conditional_complete_middle_S_far_zero'])[1]<0
    assert b['core_allowance_positive']and b['complete_middle_allocation_S_negative']
    checks.append('exact saved rational-power/rescaling consistency and conditional core/middle signs')
    assert result['Cn_improved']is False and b['Cn_improvement']is False
    assert result['old_generator_native_kernel_grid_calls']==0 and result['children']==0
    checks.append('scope preserves no-C6 conclusion and no old scientific source replay')
    output={'status':'PASS_LIMITED_SAVED_BINDINGS_AND_EXACT_CONSISTENCY',
        'checks':checks,'scientific_proof_review_separate':True,
        'scope':'saved byte, execution and exact algebra consistency only; no theorem/geometry/log/radial regeneration',
        'new_scientific_configurations':0,'generator_native_kernel_grid_calls':0,
        '104_application_dependency':'K104 independent scientific review still separate',
        'global_g_fullI_far_S_unpaid':True}
    (P/'checked_saved.json').write_text(json.dumps(output,indent=2)+'\n')
    status='COMPLETED'
except Exception:
    failure=traceback.format_exc()
finally:
    receipt={'round':'108_saved','status':status,'pid':os.getpid(),'started_utc':startutc,
        'receipt_constructed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'elapsed_before_final_receipt_write':time.monotonic()-start,
        'CPU_process_total_before_final_receipt_write':time.process_time(),
        'default_SIGALRM_still_active':signal.getsignal(signal.SIGALRM)==signal.SIG_DFL,
        'remaining_fuse_seconds':signal.getitimer(signal.ITIMER_REAL)[0],
        'fuse_never_cancelled':True,'fuse_scope':'continues through receipt/print/normal exit',
        'children':0,'retries':0,'checks_completed':checks,'failure':failure}
    (P/'checker_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':status,'pid':os.getpid(),'checks':len(checks)}))
if status!='COMPLETED':raise SystemExit(1)
