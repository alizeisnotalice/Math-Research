"""One authorized pure-rational analytic lower-window run; no F/L point replay."""
import os,sys,time,signal,resource,json,hashlib,traceback
from pathlib import Path
from fractions import Fraction as Q
from math import comb
sys.set_int_max_str_digits(0)
resource.setrlimit(resource.RLIMIT_CPU,(30,31))
def halt(sig,frame):raise TimeoutError('CPU30 or wall175 stop')
signal.signal(signal.SIGXCPU,halt);signal.signal(signal.SIGALRM,halt);signal.alarm(175)
START=time.monotonic();CPU=time.process_time();H=Path(__file__).resolve().parent
BASE=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/upper')
rd=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(n,d):(H/n).write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n')
def ep(v):return Q(v['lower_mantissa'])*Q(2)**v['lower_binary_exponent'],Q(v['upper_mantissa'])*Q(2)**v['upper_binary_exponent']
def pair(a,b):return {'lower':str(a),'upper':str(b),'display_only':[float(a),float(b)]}
if (H/'lower_receipt.json').exists():raise RuntimeError('One attempt only')
reg=rd(H/'supplement_registration.json');par=rd(H/'supplement_parameters.json')
scale=10**par['root_decimal_scale'];downscale=10**par['coefficient_downward_decimal_scale']
receipt={'status':'RUNNING','pid':os.getpid(),'attempted_configurations':0,'completed_configurations':0,
         'new_envelope_endpoint_functions':0,'old_F_L_recalculations':0,'kernel_calls':0,'automatic_retries':0}
out={'status':'INCOMPLETE','configurations':[],'root_certificates':[],'f_coefficients':{}}
def record():
    receipt.update(elapsed_seconds=time.monotonic()-START,cpu_seconds=time.process_time()-CPU)
    save('lower_receipt.json',receipt);save('lower_results.json',out)
def begin(i):
    assert i==receipt['attempted_configurations']+2<=4
    receipt['attempted_configurations']+=1;record()
def done(i,name):
    out['configurations'].append({'id':i,'name':name,'pass':True})
    receipt['completed_configurations']+=1;record()
def pins():
    assert sha(H/'supplement_parameters.json')==reg['parameters_sha256']
    assert sha(H/'supplement_theorem.md')==reg['theorem_sha256']
    for p,h in reg['source_sha256'].items():assert sha(p)==h,p
    for n,h in rd(H/'lower_execution_binding.json').items():assert sha(H/n)==h,n
def root(x,n,label):
    assert x>=0
    target=(x.numerator*scale**n)//x.denominator
    lo=0;hi=1<<((target.bit_length()+n-1)//n+1)
    while hi-lo>1:
        mid=(lo+hi)//2
        if mid**n<=target:lo=mid
        else:hi=mid
    a,b=Q(lo,scale),Q(lo+1,scale)
    assert a**n<=x<b**n
    out['root_certificates'].append({'label':label,'degree':n,'radicand':str(x),'lower':str(a),'upper':str(b)})
    return a,b
def down(x):
    z=Q((x*downscale).__floor__(),downscale)
    assert 0<z<=x<z+Q(1,downscale)
    return z

try:
    record();pins()
    old=rd(H.parent.parent/'round074/middle_long_budget/results.json')
    old51=rd(BASE/'round051/all_unsigned_middles/results.json')
    old43=rd(BASE/'round043/two_parameter_middle/results.json')['configurations'][2]
    old49=rd(BASE/'round049/middle_joint_marginal/results.json')['configurations'][2]
    epsilon,b=Q(par['epsilon']),Q(par['b']);assert b==Q(1,64)
    begin(2)
    pi_lo,pi_hi=ep(old51['shared_constants']['pi']);assert pi_lo>0
    Klo=Q(1000,729)*ep(old51['shared_constants']['Crem'])[0]*ep(old51['shared_constants']['cauchy_average_factor'])[0]
    containment=[];f={};V={}
    for k in range(1,6):
        source=old51['cases'][str(k)];geom=old43 if k==3 else source;limits=old49 if k==3 else source
        Llo,Lhi=ep(limits['z_lower']);Ulo,Uhi=ep(limits['z_upper'])
        assert Lhi<=Q(1,4) and Ulo>=Q(3,4)
        containment.append({'k':k,'z_lower_upper_endpoint':str(Lhi),'z_upper_lower_endpoint':str(Ulo),'contains_common_interval':True})
        h={n:ep(geom[n])[1] for n in ['hv0','hv1','hw0','hw1','g']}
        assert all(x>0 for x in h.values())
        raw=Klo/(h['g']**2+928*b*b)*(1/((h['hv0']**2+992*b*b)*(h['hv1']**2+960*b*b))+1/((h['hw0']**2+992*b*b)*(h['hw1']**2+960*b*b)))
        f[k]=down(raw);V[k]=down(ep(old['cases'][str(k)]['V_middle'])[0])
        out['f_coefficients'][str(k)]={'uniform_F_lower':str(f[k]),'V_lower':str(V[k]),'n':comb(5,k-1),'h_upper_inputs':{n:str(v) for n,v in h.items()},'positive_P_dropped':True,'downward_fraction_check':True}
    cmax=1024*b*b;derivative_max=2*(Q(9,16)+cmax-1)
    denominator_max=(Q(25,16)+cmax)*(Q(9,16)+cmax)
    assert cmax==Q(1,4) and derivative_max==Q(-3,8) and denominator_max==Q(377,256)
    cJ=down(Q(64)/(pi_hi*pi_hi*denominator_max))
    cW=min(5*cJ,down(1/(b*pi_hi*(1+cmax))))
    assert Q(15,64)<=Q(1,4)<Q(3,4)<=Q(49,64)
    out.update(common_interval_checks=containment,K_lower=str(Klo),denominator_max=str(denominator_max),derivative_upper=str(derivative_max),cJ_lower=str(cJ),cW_lower=str(cW),pi_upper=str(pi_hi))
    done(2,'all_saved_coordinate_intervals_and_uniform_positive_lower_coefficients')
    begin(3)
    fifth={k:root(f[k],5,f'f{k}-fifth')[0] for k in f}
    qpower={k:root(f[k]**6,5,f'f{k}-six-fifths')[0] for k in f}
    A=root(sum(comb(5,k-1)*V[k]*qpower[k] for k in f)**5,6,'volume-constant-lower')[0]
    J=root((cJ*sum(comb(5,k-1)*fifth[k] for k in f))**5,6,'joint-t-power-coefficient')[0]
    W=root((cW*max(fifth.values()))**5,6,'whole-t-power-coefficient')[0]
    C=min(J,W);xlo,xhi=root((A/C)**3,5,'power-clipping-switch')
    assert epsilon==Q(1,1000) and b==Q(1,64)
    a83=Q(1,10**8);b83=Q(1,65536)
    assert a83**3==epsilon**8 and b83**3==b**8
    if xhi<=epsilon:
        ilo=ihi=A*(b-epsilon);case='constant_on_whole_window'
    elif xlo>=b:
        ilo=ihi=Q(3,8)*C*(b83-a83);case='power_on_whole_window'
    elif epsilon<xlo<=xhi<b:
        ilo=A*b-Q(5,8)*A*xhi-Q(3,8)*C*a83
        ihi=A*b-Q(5,8)*A*xlo-Q(3,8)*C*a83
        case='interior_switch'
    else:raise ArithmeticError('UNRESOLVED switch placement; stop')
    assert 0<ilo<=ihi
    out.update(branch_lower_coefficients={'volume_constant':str(A),'joint_t5over3':str(J),'whole_t5over3':str(W),'minimum_t5over3':str(C)},
               lower_function='min(A,C*t^(5/3)) on [1/1000,1/64], zero elsewhere',
               switch=pair(xlo,xhi),integration_case=case,analytic_lower_integral=pair(ilo,ihi))
    done(3,'three_branch_minimum_and_exact_full_domain_lower_integral')
    begin(4)
    ledger=rd(BASE/'round062/primary_joint_account/actual_run/results.json')
    allowance=Q(ledger['sufficient_remaining_D_F_B_allowance']);short=Q(old['comparison']['old_short_upper_read_only'])
    paid=Q(ledger['paid_C_E_H_I_upper']);g0=Q(par['g0']);R=Q(old['strict_long_norm_rational_upper'])
    Rslab=Q(old['slab_exact_upper_sum']);tail=ep(old['infinite_tail'])[1]
    assert 0<ilo<Rslab<Rslab+tail<R
    I_lower=paid+root((short+ilo)**6,5,'new-envelope-plan-I-lower')[0]
    I_upper=paid+root((short+R)**6,5,'old-rectangle-plan-I-upper')[1]
    Tlo,Thi=root((g0/6)**6,5,'conditional-I-target-unit');Tlo*=64;Thi*=64
    factorlo,factorhi=root((6/g0)**6,5,'conditional-U-factor')
    Blo=root(I_lower**5,6,'envelope-plan-B-lower')[0]
    Bhi=root(I_upper**5,6,'old-rectangle-plan-B-upper')[1]
    Slo,Shi=Tlo-I_upper,Thi-I_lower
    Ulo,Uhi=factorlo*I_lower,factorhi*I_upper
    status='STRICT_FIXED_ENVELOPE_NONCOMPETITIVE' if ilo>allowance and Shi<0 else 'UNRESOLVED_STOP'
    out.update(status=status,conditional_g0=str(g0),global_g0_certified=False,
        same_E_integral_enclosure=pair(ilo,R),saved_D_F_B_allowance=str(allowance),
        envelope_lower_minus_allowance=str(ilo-allowance),
        analytic_lower_over_allowance=pair(ilo/allowance,ihi/allowance),
        conditional_plan_I=pair(I_lower,I_upper),conditional_plan_B=pair(Blo,Bhi),
        conditional_I_target=pair(Tlo,Thi),conditional_plan_U=pair(Ulo,Uhi),conditional_plan_S_margin=pair(Slo,Shi),
        remaining_irrecoverable_plan_I_deficit=str(I_lower-Thi),
        maximum_86_rectangle_error=str(Rslab-ilo),
        maximum_full_time_integration_gain=str(R-ilo),
        old_infinite_tail_upper=str(tail),
        maximum_recoverable_S_gain=str(I_upper-I_lower),
        maximum_recoverable_total_B_gain=str(Bhi-Blo),
        far_assumed_zero_only_for_optimistic_plan_diagnosis=True,
        true_norm_lower_bound=False,actual_model_cost_interval=False,actual_C6_bound=False,
        legacy_cost_results_S_fields_mean='U; metadata-label correction only, no repeated numerical evaluation',
        scope='Only the fixed074 envelope E with this conservative regional accounting; not true p/G norm, R_c, model, or C6')
    done(4,'strict_gap_and_correct_cost_U_margin_S_accounting')
    pins();receipt.update(status='COMPLETED_STOPPED',exit_code=0);record()
    receipt['results_sha256']=sha(H/'lower_results.json');save('lower_receipt.json',receipt)
    print(json.dumps({'status':status,'pid':os.getpid(),'CPU':receipt['cpu_seconds'],'wall':receipt['elapsed_seconds'],
        'J0_display':out['analytic_lower_integral']['display_only'],
        'ratio_display':out['analytic_lower_over_allowance']['display_only'],
        'U_display':out['conditional_plan_U']['display_only'],'S_margin_display':out['conditional_plan_S_margin']['display_only'],
        'clipping_switch':out['switch']['display_only']}))
except BaseException as ex:
    receipt.update(status='FAILED_STOP_NO_RETRY',exit_code=1,error=repr(ex),traceback=traceback.format_exc());record();raise
finally:signal.alarm(0)
