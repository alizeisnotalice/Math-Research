"""108 PREPARED; three fixed pure exact/Arb configurations after parent release.
Default SIGALRM177 remains active through final receipt, print and normal exit.
CPU25; no child, old generator, native, MATLAB, kernel, fan or grid call.
"""
import signal
signal.signal(signal.SIGALRM, signal.SIG_DFL)
signal.setitimer(signal.ITIMER_REAL, 177.0)

import sys
sys.dont_write_bytecode = True
import datetime, hashlib, itertools, json, os, resource, time, traceback
from fractions import Fraction as Q
from pathlib import Path
resource.setrlimit(resource.RLIMIT_CPU, (25, 25))
P = Path(__file__).resolve().parent
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
release_path = P/'parent_release.json'
if not release_path.exists(): raise SystemExit('108 PREPARED_NOT_RELEASED: no mathematics')
release = json.loads(release_path.read_text())
if not (release.get('authorized') is True and release.get('round') == '108'
        and release.get('script_sha256') == sha(Path(__file__))
        and release.get('registration_sha256') == sha(P/'registration.json')):
    raise SystemExit('108 parent release mismatch: no mathematics')
if any((P/n).exists() for n in ['started.json','results.json','run_receipt.json','checkpoints.json']):
    raise SystemExit('108 one attempt only: existing run artifact')
started=time.monotonic(); cpu_start=time.process_time()
started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()
(P/'started.json').write_text(json.dumps({'round':'108','status':'RELEASED_STARTED',
    'pid':os.getpid(),'started_utc':started_utc,'SIGALRM_handler':'SIG_DFL',
    'fuse177_active':True,'CPU_seconds':25,'children':0,
    'release_sha256':sha(release_path)},indent=2)+'\n')
completed=[]

def done(row):
    completed.append(row)
    (P/'checkpoints.json').write_text(json.dumps({'round':'108',
        'completed_configurations':completed,'remaining_fuse_seconds':signal.getitimer(signal.ITIMER_REAL)[0],
        'SIGALRM_handler_is_default':signal.getsignal(signal.SIGALRM)==signal.SIG_DFL},indent=2)+'\n')

def run_math():
    from flint import arb, fmpq, ctx
    ctx.prec=384
    reg=json.loads((P/'registration.json').read_text())
    assert sha(P/'source_manifest.json') == reg['source_manifest_sha256']
    assert sha(P/'saved104_inputs.json') == reg['saved104_inputs_sha256']
    data=json.loads((P/'saved104_inputs.json').read_text())
    def A(x):
        x=Q(x); return arb(fmpq(x.numerator,x.denominator))
    def endpoints(d):
        exponent=int(d['common_binary_exponent'])
        return tuple(Q(int(d[side+'_mantissa']))*Q(2)**exponent for side in ['lower','upper'])
    def bounds(lo,hi=None):
        if hi is None: hi=lo
        lm,le=lo.lower().man_exp(); um,ue=hi.upper().man_exp()
        le=int(le); ue=int(ue); e=min(le,ue)
        return {'lower_mantissa':str(int(lm)*2**(le-e)),
                'upper_mantissa':str(int(um)*2**(ue-e)),
                'common_binary_exponent':str(e)}
    def scaled_saved(d,factor):
        lo,hi=endpoints(d); assert 0 < lo <= hi
        return A(lo)*factor,A(hi)*factor
    pi=arb.pi(); q=Q(6,5)

    # 1. Exact generic quadratic sign polynomial identity, not a sample fit.
    signs=[(1,)+s for s in itertools.product([-1,1],repeat=5)]
    assert len(signs)==32
    gram=[[sum(v[i]*v[j] for v in signs) for j in range(6)] for i in range(6)]
    assert gram==[[32*int(i==j) for j in range(6)] for i in range(6)]
    entries=[(i,i) for i in range(6)]+[(i,j) for i in range(6) for j in range(i+1,6)]
    identity=[]
    for a,(i,j) in enumerate(entries):
        for b in range(a,len(entries)):
            k,l=entries[b]
            coeff=Q(sum((1 if i==j else 2*v[i]*v[j]) *
                        (1 if k==l else 2*v[k]*v[l]) for v in signs),32)
            if a!=b: coeff*=2
            expected=(Q(1) if a==b and i==j else Q(4) if a==b else
                      Q(2) if i==j and k==l else Q(0))
            assert coeff==expected
            if coeff: identity.append({'entries':[[i,j],[k,l]],'coefficient':str(coeff)})
    # Exact orthogonal rank-three projection with four zero representatives.
    proj=[[Q(0) for _ in range(6)] for _ in range(6)]
    for b in range(0,6,2):
        proj[b][b]=proj[b+1][b+1]=Q(1,2)
        proj[b][b+1]=proj[b+1][b]=Q(-1,2)
    assert all(sum(proj[i][k]*proj[k][j] for k in range(6))==proj[i][j]
               for i in range(6) for j in range(6))
    assert sum(proj[i][i] for i in range(6))==3
    vals=[sum(Q(v[i])*proj[i][j]*v[j] for i in range(6) for j in range(6)) for v in signs]
    assert Q(sum(vals),32)==3 and Q(sum(z*z for z in vals),32)==12
    histogram={str(z):vals.count(z) for z in sorted(set(vals))}
    assert histogram=={'0':4,'2':12,'4':12,'6':4}
    done({'configuration':1,'kind':'exact finite-sign projection moment polynomial',
        'body32_frame_gram':gram,'nonzero_quadratic_moment_coefficients':identity,
        'general_identity':'EZ2=(trP)^2+4sum_i<j Pij^2; rank-r projector =>r(r+2)-2sumPii^2',
        'zero_projection_witness_histogram':histogram,'zero_projection_witness_EZ':'3',
        'zero_projection_witness_EZ2':'12','scope':'exact algebra, universal analytic proof separate'})

    # 2. Only algebraically rescale complete saved104 dyadics; no old function.
    lam=8/(3*A(3).sqrt()*pi)
    sharp=A(Q(8,27)).sqrt()
    assert lam.upper()<sharp.lower() and (lam/sharp-2*A(2).sqrt()/pi).contains(0)
    assert lam.lower()>0 and lam.upper()<A(Q(1,2)).lower()
    qold=1/(pi*pi*(32/A(3).sqrt())**3)
    qnew=1/(4096*pi**3)
    assert (qnew-lam*qold).contains(0)
    factor=lam**A(q)
    core_norm=scaled_saved(data['core_full_norm'],lam)
    core_I=scaled_saved(data['core_ten_integral'],factor)
    middle_I=scaled_saved(data['complete_middle_integral'],factor)
    paid=A(Q(data['paid_CEH_I_upper']))
    target=64*A(Q(1,65000)/6)**A(q)
    core_plan=(paid+core_I[0],paid+core_I[1])
    core_allowance=(target-paid-core_I[1],target-paid-core_I[0])
    middle_plan=(paid+middle_I[0],paid+middle_I[1])
    middle_S=(target-paid-middle_I[1],target-paid-middle_I[0])
    oldlo,oldhi=endpoints(data['old051short074long_middle_expression'])
    ratio_old=(middle_I[0]/A(oldhi),middle_I[1]/A(oldlo))
    actual_cost=(middle_plan[0]/A(Q(1,65000)/6)**A(q),
                 middle_plan[1]/A(Q(1,65000)/6)**A(q))
    done({'configuration':2,'kind':'universal3D coefficient and saved104 algebraic transfer',
        'Q3_new':bounds(qnew),'ratio_to098':bounds(lam),
        'conditional_sharpKhinchine_ratio_to098':bounds(sharp),
        'ratio_to_conditional_sharpKhinchine':bounds(lam/sharp),
        'integral_scaling_lambda_power6over5':bounds(factor),
        'one_core_full_G_norm':bounds(*core_norm),'ten_core_I_allocation':bounds(*core_I),
        'paidCEH_plus_ten_core_plan':bounds(*core_plan),
        'conditional_core_remaining_allowance':bounds(*core_allowance),
        'complete31_middle_I_allocation':bounds(*middle_I),
        'new_over_old051short074long_middle_expression':bounds(*ratio_old),
        'paidCEH_plus_complete_middle_plan':bounds(*middle_plan),
        'conditional_complete_middle_S_far_zero':bounds(*middle_S),
        'conditional_complete_middle_U_far_zero':bounds(*actual_cost),
        'core_allowance_positive':core_allowance[0].lower()>0,
        'complete_middle_allocation_S_negative':middle_S[1].upper()<0,
        'g0_status':'1/65000 unproved globally','far_status':'UNPAID; zero only in diagnostic',
        'inputs':'full exact saved104 dyadic endpoints, no rounded decimal',
        'old_generator_invocations':0,'Cn_improvement':False})

    # 3. Exact general-r coefficients and the direct full6 special case.
    rows=[]
    for r in range(1,7):
        M=Q(2)+Q(2*r,3)
        coeff_squared=M**r/Q(32)**(2*r)
        coefficient=A(M)**A(Q(r,2))/32**r
        rows.append({'dimension':r,'M_r':str(M),
            'exact_square_of_pi_free_peak_coefficient':str(coeff_squared),
            'peak_coefficient_including_pi':bounds(coefficient/pi**r)})
    c6=Q(216,2**30)
    assert c6*c6==Q(rows[5]['exact_square_of_pi_free_peak_coefficient'])
    assert c6/(Q(1,160))==Q(135,4194304)
    assert Q(4)**3/Q(32)**6==Q(1,4096)**2
    done({'configuration':3,'kind':'r1to6 peak constants and direct full6 frame',
        'dimensions':rows,'full6_pi_free_coefficient':str(c6),
        'full6_ratio_to_saved001_det160_peak':'135/4194304',
        'full6_peak_including_pi':bounds(A(c6)/pi**6),
        'scope':'q_E,t(x)<=q_E,t(0)<=C_r t^-r; no isotropic spatial-decay shape',
        'kernel_or_peak_numerical_evaluations':0})
    assert len(completed)==reg['configuration_cap']==3
    return {'round':'108','status':'PASS_THREE_CONFIGURATIONS_PROOF_AUDIT_SEPARATE',
        'precision_bits':384,'configurations':completed,'source_code_sha256':sha(Path(__file__)),
        'registration_sha256':sha(P/'registration.json'),
        'source_manifest_sha256':sha(P/'source_manifest.json'),
        'saved104_inputs_sha256':sha(P/'saved104_inputs.json'),
        'old_generator_native_kernel_grid_calls':0,'children':0,'Cn_improved':False}

status='FIRST_FAILURE_STOP'; failure=None
try:
    result=run_math()
    (P/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    status='COMPLETED'
except Exception:
    failure=traceback.format_exc()
finally:
    receipt={'round':'108','status':status,'pid':os.getpid(),
        'completed_configurations':len(completed),'started_utc':started_utc,
        'receipt_constructed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'elapsed_before_final_receipt_write':time.monotonic()-started,
        'CPU_since_start_marker_before_final_receipt_write':time.process_time()-cpu_start,
        'CPU_process_total_before_final_receipt_write':time.process_time(),
        'SIGALRM_handler_is_default':signal.getsignal(signal.SIGALRM)==signal.SIG_DFL,
        'fuse_never_cancelled':True,'remaining_fuse_seconds':signal.getitimer(signal.ITIMER_REAL)[0],
        'fuse_scope':'continues through final receipt, print and normal interpreter exit',
        'children_created':0,'cleanup_required':'none; one process','retries':0,
        'failure':failure,'release_sha256':sha(release_path)}
    (P/'run_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'round':'108','status':status,'completed_configurations':len(completed),
        'pid':os.getpid(),'default_fuse_active_at_receipt':receipt['SIGALRM_handler_is_default']}))
if status!='COMPLETED': raise SystemExit(1)
