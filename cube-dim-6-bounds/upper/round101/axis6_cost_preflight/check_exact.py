#!/usr/bin/env python3
"""One new pure-exact axes6 cost certificate; no old point or grid replay."""
import signal
signal.signal(signal.SIGALRM,signal.SIG_DFL)
signal.setitimer(signal.ITIMER_REAL,150)
from pathlib import Path
from fractions import Fraction as Q
import hashlib
import json
import resource
import sys
import time

resource.setrlimit(resource.RLIMIT_CPU,(28,28))
signal.signal(signal.SIGXCPU,signal.SIG_DFL)
CPU,WALL=time.process_time(),time.monotonic()
H=Path(__file__).resolve().parent
checks=0
out={'status':'RUNNING','configurations':[],
     'old_parameter_points_recomputed':0,'old_grids_replayed':0,
     'body_native_calls':0,'p_G_kernel_calls':0,'MATLAB_calls':0,
     'automatic_retries':0,'hard_fuse_not_cancelled':True}
def ck(ok,msg):
    global checks
    checks+=1
    if not ok:raise AssertionError(msg)
    if time.process_time()-CPU>=26:raise TimeoutError('CPU prelimit')
def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def done(i,name):
    ck(len(out['configurations'])+1==i<=4,'configuration cap/order')
    out['configurations'].append({'id':i,'name':name,'pass':True})
try:
    reg=load(H/'registration.json');inp=load(H/'inputs.json');bind=load(H/'execution_registration.json')
    for p,expected in reg['source_pins'].items():ck(sha(p)==expected,'source pin '+p)
    ck(sha(H/'inputs.json')==reg['inputs_sha256'],'input pin')
    for name,expected in bind['execution_pins'].items():ck(sha(H/name)==expected,'execution pin '+name)
    old=load(next(p for p in reg['source_pins'] if p.endswith('/upper/theory/axis_frame_cost_probe.json')))
    ck(old['kind']=='numerical_candidate_only','old pilot status')
    ck(old['g_at_cube_corner_exact']==inp['g_exact_inherited_corner_formula']=='3/(512*pi^5)','inherited exact corner formula handoff; no point evaluation')
    pi_l,pi_u=map(Q,inp['pi_enclosure_inherited_E098'])
    P,Qp=map(Q,inp['pi_coarse_bounds'])
    ck(P==Q(157,50)<pi_l<pi_u<Qp==Q(22,7),'read-only inherited pi certificate')
    q=Q(inp['q']);ck(q==Q(6,5),'q')

    # 1. Normalization algebra and model/source gates; no kernel evaluation.
    # The known Beta coefficient is audited by Gamma half-integer recurrence,
    # not by reevaluating the old G(ones) point or its numerical pilot.
    beta_coefficient=Q(1,2)*Q(15,8)*Q(3,4)/120
    ck(beta_coefficient==Q(3,512),'inherited Beta/pi coefficient identity')
    J_factor=6*64;mass_factor=J_factor*Q(5,6)
    ck(J_factor==384 and mass_factor==320 and mass_factor/64==5,'full boundary/radial/source cost factors')
    s=Q(1,6)
    ck(-1/(5*s)+1/(1-s)==0,'optimal spatial truncation ratio')
    ck((5/(1-s))**5/s==6**6,'minimum coefficient equals 6^(6/5)')
    ck(1<q and 6-1==5 and q*5==6,'order-one dimension/homogeneity')
    done(1,'inherited_floor_formula_handoff_and_new_full_space_cost_normalization')

    # 2. New fixed-window true face-mean lower certificate.
    a,b=map(Q,inp['fixed_mean_window'])
    ck((a,b)==(Q(1,4),Q(1,2)) and 2**2>3,'all-window atan lower branch')
    ratio=(1+b*b)/(1+a*a)
    ck(ratio==Q(20,17),'exact fixed-window primitive argument')
    x=ratio-1;log_lower=x-x*x/2
    ck(log_lower==Q(93,578)>0,'strict analytic logarithm lower polynomial')
    A_lower=log_lower/(486*Qp)
    R_lower=Q(512,243)*P**4*log_lower
    lower_cost=Q(inp['fixed_lower_cost_target'])
    ck(lower_cost==66 and R_lower**6>lower_cost**5,'true ideal model cost strictly exceeds66')
    I_lower=Q(inp['fixed_face_I_targets'][0])
    ck(A_lower**6>I_lower**5,'true entire-face I lower')
    out['true_lower_certificate']={'window':[str(a),str(b)],'mean_A_lower':str(A_lower),
      'log_20_over17_lower':str(log_lower),'six_B_over_g_lower':str(R_lower),
      'cost_lower':'66','power_residual_R6_minus66power5':str(R_lower**6-lower_cost**5),
      'I_lower':str(I_lower),'uses_actual_mean_not_upper_budget':True}
    done(2,'fixed_window_actual_mean_and_Jensen_true_cost_lower')

    # 3. New analytic whole-time/whole-face upper certificate.
    ck(P>3 and P**11>8**6,'coarse pi fractional-power gate')
    ck(25>24,'two-to-minus25over6 is below one-sixteenth')
    B_upper=Q(3,448)+Q(1,3645)
    ck(B_upper<Q(1,128),'complete short_plus_infinite_tail norm upper')
    R_upper=8*Qp**5
    upper_cost=Q(inp['fixed_upper_cost_target'])
    ck(upper_cost==12000 and R_upper**6<upper_cost**5,'true ideal model cost below12000')
    I_upper=Q(inp['fixed_face_I_targets'][1])
    ck(Q(1,128)**6<I_upper**5,'true entire-face I upper')
    out['true_upper_certificate']={'time_split':'(0,1] and [1,infinity)',
      'B_rational_upper':str(B_upper),'B_simplified_upper':'1/128',
      'six_B_over_g_upper':str(R_upper),'cost_upper':'12000',
      'power_residual_12000power5_minusR6':str(upper_cost**5-R_upper**6),'I_upper':str(I_upper)}
    done(3,'analytic_whole_face_all_time_Minkowski_and_complete_cost_interval')
    for p,expected in reg['source_pins'].items():ck(sha(p)==expected,'post-run source immutability')
    ck(len(out['configurations'])==3,'actual configuration count')
    out.update(status='STRICT_ACTUAL_AXES6_GREEN_POSITIVE_PART_FAMILY_NONCOMPETITIVE',
      exact_global_cube_floor='3/(512*pi^5)',actual_model_cost_strict_interval=['66','12000'],
      actual_face_I_strict_interval=inp['fixed_face_I_targets'],
      all_spatial_truncation_c_excluded_from_improving64=True,
      actual_C6_lower_bound_claimed=False,new_C6_improvement=False,
      arbitrary_axes_kernel_or_anisotropic_weights_excluded=False,
      proof_scope='actual equal-axis Green-positive-part family only; full-space source and all-common-radius/general-L1 transfer checked in theorem.md',
      inherited_pi_source='E098 sealed exact rational certificate, read only')
except Exception as exc:
    out.update(status='FAILED_FIRST_FAILURE_NO_RETRY',error=repr(exc))
finally:
    out.update(assertions=checks,CPU_seconds_before_save=time.process_time()-CPU,
      wall_seconds_before_save=time.monotonic()-WALL,CPU_hard_seconds=28,
      hard_wall_fuse_seconds=150,hard_fuse_disposition='SIG_DFL; never cancelled',
      self_sha256=sha(Path(__file__)))
    with (H/'results.json').open('x') as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
    print(json.dumps({k:out[k] for k in ['status','assertions','CPU_seconds_before_save','wall_seconds_before_save']}),flush=True)
# The default hard timer remains armed through the output and normal exit.
sys.exit(0 if out['status']=='STRICT_ACTUAL_AXES6_GREEN_POSITIVE_PART_FAMILY_NONCOMPETITIVE' else 1)
