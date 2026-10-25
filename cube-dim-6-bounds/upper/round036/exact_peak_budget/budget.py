#!/usr/bin/env python3
"""B036: twelve registered Arb512 configurations, blocked until F035/release."""
import sys
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
import hashlib,json,os,signal,time,traceback
from fractions import Fraction as Q
from pathlib import Path

HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def save(name,obj):
    p=HERE/(name+'.tmp');p.write_text(json.dumps(obj,indent=2,ensure_ascii=False)+'\n');os.replace(p,HERE/name)
def bounds(row):
    return(Q(int(row['lower_mantissa']))*Q(2)**int(row['lower_binary_exponent']),
           Q(int(row['upper_mantissa']))*Q(2)**int(row['upper_binary_exponent']))

def run():
    assert not(HERE/'run_receipt.json').exists(),'one batch already started'
    assert(HERE/'release.json').exists(),'NOT RELEASED: zero scalar configurations'
    reg=read(HERE/'registration.json');release=read(HERE/'release.json')
    assert release['explicit_release'] and release['F035_PASS_or_main_release']
    assert sha(HERE/'registration.json')==release['registration_sha256']
    for path,digest in release['gate_source_sha256'].items():assert sha(path)==digest
    for path,digest in reg['source_sha256'].items():assert sha(path)==digest
    assert sha(HERE/'budget.py')==reg['code_sha256']['budget.py']
    peaks=read(reg['peaks_path']);cases={row['id']:row for row in peaks['cases']}
    assert cases['r1']['R_exact']=='67320683/3251298631680'
    expected=[('r3_4110',480,64),('r3_3210',1920,96),('r3_2220',480,128),('r3_3111',640,160),('r3_2211',1440,192)]
    for name,n,d in expected:
        assert cases[name]['triple_orbit_count']==n and cases[name]['Gram_determinant']==d
    assert reg['delta']=='1/16' and reg['epsilon']=='1/1000'
    assert reg['g_target']=='1/65000' and reg['H_norm_goal']=='1/14000'
    from flint import arb,fmpq,ctx
    ctx.prec=512
    def B(value):
        value=Q(value);return arb(fmpq(value.numerator,value.denominator))
    def enclosed(value):
        assert value.is_finite()
        l,e=value.lower().man_exp();u,f=value.upper().man_exp()
        return dict(lower_mantissa=str(l),lower_binary_exponent=int(e),upper_mantissa=str(u),upper_binary_exponent=int(f),display_only=str(value))
    tic=time.perf_counter()
    receipt=dict(status='RUNNING',pid=os.getpid(),start_epoch=time.time(),max_wall_seconds=1200,
                 registration_sha256=sha(HERE/'registration.json'),release_sha256=sha(HERE/'release.json'),
                 completed_configurations=0,attempted_configurations=0,spatial_kernel_calls=0,grid_points=0,
                 native_G_H_K_calls=0,matlab_calls=0)
    out=dict(status='INCOMPLETE',precision_bits=512,configurations=[],new_I_certificate=False,new_C6_bound=False)
    save('run_receipt.json',receipt)
    def timeout(signum,frame):raise TimeoutError('registered 20 minute limit')
    signal.signal(signal.SIGALRM,timeout);signal.alarm(1200)
    def begin(number,name):
        assert number==receipt['attempted_configurations']+1<=12
        receipt['attempted_configurations']=number;receipt['current_configuration']=name;save('run_receipt.json',receipt)
    def done(number,name,data):
        out['configurations'].append(dict(index=number,name=name,**data));receipt['completed_configurations']=number
        save('results.json',out);save('run_receipt.json',receipt)
    code=1
    try:
        begin(1,'rank1_true_peak_and_face_constant')
        pi=arb.pi();R1=Q(cases['r1']['R_exact']);Q5=B(6).sqrt()*B(R1)/(2*pi)**5
        Cf=32*B(R1)/pi**6;A=Cf**B(Q(1,6))*pi**B(Q(-5,6));K=B(Q(6,7))*A
        assert (1024*Q5/(pi*B(6).sqrt())-Cf).contains(0)
        cf_ratio=7776*R1;assert 0<cf_ratio<1
        done(1,'rank1_true_peak_and_face_constant',dict(R1=str(R1),pi=enclosed(pi),Q5=enclosed(Q5),Cf=enclosed(Cf),A=enclosed(A),K=enclosed(K),new_over_old_Cf_exact=str(cf_ratio)))
        weighted=Q(0);terms=[]
        for i,(name,n,d)in enumerate(expected,2):
            begin(i,name)
            r=Q(cases[name]['R_exact']);kappa=B(r)/(2*pi)**3;peak=B(d).sqrt()*kappa
            weighted+=n*r;terms.append(6*n*kappa)
            done(i,name,dict(R3=str(r),unordered_count=n,ordered_count=6*n,Gram_determinant=d,
                             Q3=enclosed(peak),Q3_over_sqrt_Gram=enclosed(kappa),ordered_contribution=enclosed(terms[-1])))
        begin(7,'all_ordered_triples_B3')
        coef=24576*weighted;B3=B(coef)/pi**6
        assert ((32/pi)**3*sum(terms,arb(0))-B3).contains(0)
        done(7,'all_ordered_triples_B3',dict(unordered_total=4960,ordered_total=29760,weighted_R3_sum=str(weighted),B3_pi6_exact=str(coef),B3=enclosed(B3)))
        begin(8,'single_delta_new_geometry')
        delta=B(Q(1,16));C=pi**2/2*(10+10*B(2).sqrt()+2*B(5).sqrt());V=C*delta**4
        assert V>0 and V<1
        oldrow=read(reg['geometry_results_path'])['rows'][0];assert oldrow['delta']=='1/16'
        vlo,vhi=bounds(enclosed(V));oldlo,oldhi=Q(oldrow['new_volume']['lower_exact']),Q(oldrow['new_volume']['upper_exact'])
        assert max(vlo,oldlo)<=min(vhi,oldhi)
        done(8,'single_delta_new_geometry',dict(delta='1/16',Cnew=enclosed(C),volume_upper=enclosed(V),overlaps_frozen_exact_034_endpoints=True))
        begin(9,'single_epsilon_whole_face_short_tail')
        eps=B(Q(1,1000));tail=K*eps**B(Q(7,6))
        done(9,'single_epsilon_whole_face_short_tail',dict(epsilon='1/1000',tail_upper_budget=enclosed(tail)))
        begin(10,'single_delta_epsilon_off_S2_short_tail')
        off=B(Q(6,11))*B3**B(Q(1,6))*pi**B(Q(-5,6))*delta**B(Q(-2,3))*eps**B(Q(11,6))
        done(10,'single_delta_epsilon_off_S2_short_tail',dict(delta='1/16',epsilon='1/1000',off_S2_tail_budget=enclosed(off),off_over_whole_budget=enclosed(off/tail)))
        begin(11,'single_new_conditional_Green_cost')
        g=B(Q(1,65000));h=B(Q(1,14000));target=B(Q(4,48750));goal=h+tail
        Igoal=goal**B(Q(6,5));cost=(6*goal/g)**B(Q(6,5))
        blo,bhi=bounds(enclosed(goal));clo,chi=bounds(enclosed(cost))
        competitive=(bhi<Q(4,48750) and chi<64)
        done(11,'single_new_conditional_Green_cost',dict(g_target='1/65000',H_norm_goal='1/14000',B_critical_exact='4/48750',B_goal=enclosed(goal),I_goal=enclosed(Igoal),conditional_cost=enclosed(cost),strictly_below64_if_both_global_premises_proved=competitive,both_global_premises_proved=False))
        begin(12,'single_tube_two_envelope_budget_and_stop')
        e=(Cf*pi*V)**B(Q(1,5));R=1/(3*e**3)-1024/e+32768*(pi/2-(32*e).atan())
        assert R>0
        short=K*e**B(Q(7,6));long=V**B(Q(5,6))*Cf*R
        U=short+long;L=short/(1+1024*e**2)+long
        assert U>L>0
        ul,uu=bounds(enclosed(U));ll,lu=bounds(enclosed(L))
        failure=ll>Q(4,48750)
        lower_ratio=ll/Q(4,48750);simple=Q((lower_ratio*10**6).__floor__(),10**6)
        done(12,'single_tube_two_envelope_budget_and_stop',dict(delta='1/16',derived_epsilon_selection_rule='(Cf*pi*V)^(1/5)',epsilon_star=enclosed(e),R=enclosed(R),allocated_upper_budget=enclosed(U),best_two_envelope_lower_budget=enclosed(L),upper_over_new_target=enclosed(U/target),lower_over_new_target=enclosed(L/target),safe_lower_method_ratio=str(simple),two_envelope_method_fails_new_target=failure,old_cap_formula_stopped=True,actual_tube_G_norm_lower_bound=False))
        out['status']=('EXACT_PEAKS_IMPROVE_CONSTANTS_BUT_SINGLE_TUBE_FORMULA_FAILS'if failure else'SINGLE_TUBE_DIAGNOSTIC_COMPLETE_NO_GLOBAL_CERTIFICATE')
        out['conditional_cost_below64']=competitive;out['ordered_weights_complete']=True
        out['new_global_g_lower_bound']=False;out['new_global_H_norm_upper_bound']=False
        out['new_C6_bound']=False;save('results.json',out)
        receipt['status']='COMPLETED';code=0
    except BaseException as exc:
        receipt['status']='STOPPED_FAILURE_NO_RETRY';receipt['error']=repr(exc);receipt['traceback']=traceback.format_exc()
    finally:
        signal.alarm(0);receipt['elapsed_seconds']=time.perf_counter()-tic;receipt['exit_code']=code
        if(HERE/'results.json').exists():receipt['results_sha256']=sha(HERE/'results.json')
        save('run_receipt.json',receipt)
    print(json.dumps({k:v for k,v in receipt.items()if k!='traceback'}),flush=True)
    return code

if __name__=='__main__':
    assert sys.argv[1:]==['--run']
    raise SystemExit(run())
