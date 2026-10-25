#!/usr/bin/env python3
"""Six fixed new configurations for the integrated two-parameter mechanism."""
import sys
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
import json,hashlib,time,datetime,os,traceback
import sympy as sp
from flint import arb,fmpq,ctx
ctx.prec=512
HERE=Path(__file__).resolve().parent
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(name,obj):(HERE/name).write_text(json.dumps(obj,indent=2)+'\n')
def A(x):
    x=Q(x);return arb(fmpq(x.numerator,x.denominator))
def enclose(x):
    assert x.is_finite()
    a,e=x.lower().man_exp();b,f=x.upper().man_exp()
    lo=Q(int(a))*Q(2)**int(e);hi=Q(int(b))*Q(2)**int(f)
    assert 0<=hi-lo<=max(Q(1),abs((lo+hi)/2))*Q(2)**-450
    return dict(lower_mantissa=str(a),lower_binary_exponent=int(e),upper_mantissa=str(b),upper_binary_exponent=int(f),display_only=str(x))
def run():
    assert not(HERE/'run_receipt.json').exists(),'one batch only'
    reg=read(HERE/'registration.json')
    for p,h in reg['source_sha256'].items():assert sha(p)==h,p
    for f,h in reg['code_sha256'].items():assert sha(HERE/f)==h,f
    deadline=datetime.datetime.fromisoformat(reg['deadline_utc']).timestamp()
    assert time.time()<deadline
    receipt=dict(status='RUNNING',pid=os.getpid(),started_epoch=time.time(),attempted_configurations=0,completed_configurations=0,registration_sha256=sha(HERE/'registration.json'),native_calls=0,spatial_grid_points=0,matlab_calls=0)
    result=dict(status='INCOMPLETE',precision_bits=512,configurations=[],new_I_certificate=False,new_C6_bound=False)
    save('run_receipt.json',receipt);tic=time.monotonic();code=1
    def begin(i):
        assert i==receipt['attempted_configurations']+1<=6 and time.time()<deadline
        receipt['attempted_configurations']=i;save('run_receipt.json',receipt)
    def done(i,data):
        result['configurations'].append(dict(id=i,**data));receipt['completed_configurations']=i
        save('results.json',result);save('run_receipt.json',receipt)
    try:
        begin(1)
        V=[(1,)+s for s in product((-1,1),repeat=5)]
        values=[sum(v)for v in V];assert sum(map(abs,values))==60
        assert Q(60,6)==10 and 1-Q(1,10)==Q(9,10)
        done(1,dict(name='symbol_comparison_for_one_omitted_parameter',body_directions=32,psi_ones=60,psi_over_l1_lower=10,missing_symbol_fraction='9/10',dimension3_peak_loss='1000/729',new_projection_peak_integrals=0))
        begin(2)
        v=sp.ones(6,1);w=sp.Matrix([1,1,1,-1,-1,-1])
        b=sp.Matrix([1,1,1,1,1,-1]);c=sp.Matrix([1,1,1,-1,-1,1]);B=sp.Matrix.hstack(b,c)
        x=sp.Matrix([1,1,1,sp.Rational(1,4),sp.Rational(1,4),sp.Rational(1,4)]);a=sp.Rational(1,4);y=x-a*v
        P=sp.eye(6)-B*(B.T*B).inv()*B.T
        assert y==sp.Rational(3,8)*(b+c) and P*y==sp.zeros(6,1) and (P*x).dot(P*x)==sp.Rational(1,6)
        assert sp.Matrix.hstack(v,b,c).rank()==3
        J2sq=sp.Matrix.hstack(v,w).T*sp.Matrix.hstack(v,w);assert J2sq.det()==36
        L,s,gamma,aa,lam=sp.symbols('L s gamma aa lam',positive=True)
        normz=L**2*(aa**2+gamma**2)-lam*s**2
        diff=sp.expand(L**2*(aa**2+(s+gamma)**2)-(normz+(lam+L**2)*s**2+2*L**2*s*gamma))
        assert diff==0
        done(2,dict(name='true_zero_and_cauchy_resolvent_identity',spatial_u='1/4',retained_parameter_a='1/4',other_pair=[list(map(int,b)),list(map(int,c))],residual_coefficients=['3/8','3/8'],residual_pair_distance_squared='0',original_x_pair_distance_squared='1/6',selected_pair_Gram_determinant=36,selected_pair_coarea_J2=6,resolvent_denominator_identity_symbolic='0',uniform_fiber_distance_claim_refuted=True))
        begin(3)
        pi=arb.pi();rho=A(Q(1,64));eps=A(Q(1,1000));eta=A(3).sqrt()*rho
        hv0=A(3).sqrt()/2-eta;hv1=A(Q(3,32)).sqrt()-eta
        hw0=15*A(3).sqrt()/16-eta;hw1=A(Q(3,128)).sqrt()-eta
        g=A(Q(1,32)).sqrt()-eta
        assert all(x>0 for x in [hv0,hv1,hw0,hw1,g])
        F=1+A(Q(6,928)).sqrt()
        done(3,dict(name='single_middle_tube_geometry',rho='1/64',epsilon='1/1000',middle_interval=['1/4','3/4'],thresholds=['1/2','1/16'],pi=enclose(pi),eta=enclose(eta),hv0=enclose(hv0),hv1=enclose(hv1),hw0=enclose(hw0),hw1=enclose(hw1),g=enclose(g),cauchy_average_factor=enclose(F),deleted_frame_Gram_bounds=[26,32],true_time_coefficients=[992,960,928]))
        begin(4)
        peaks=read(HERE.parents[1]/'round033/exact_projection_peaks/results.json');pc={r['id']:r for r in peaks['cases']}
        R2=Q(pc['r2_k3']['R_exact']);assert R2==Q(2051,47038464)
        C2=64*A(R2)/pi**6
        good=A(Q(2,3))*C2**A(Q(1,6))*pi**A(Q(-5,6))*eps**A(Q(3,2))
        done(4,dict(name='both_specified_parameters_large',R2_exact=str(R2),C2=enclose(C2),big_branch_short_norm_upper=enclose(good),quotient_dimension=4,coarea_J2=6,coefficient_C2_times_pi6=str(64*R2)))
        begin(5)
        old=read(HERE.parents[1]/'round041/one_large_jump_cap/results.json')
        coef=Q(old['configurations'][1]['Crem_times_pi6_exact']);assert coef==Q(44743079,2970)
        Crem=A(coef)/pi**6;loss=A(Q(1000,729))
        common=loss*Crem*F/(g*g)
        Mv=common/(hv0*hv0*hv1*hv1);Mw=common/(hw0*hw0*hw1*hw1);M=Mv+Mw
        done(5,dict(name='both_complete_complement_bounds',Crem_exact_numerator_over_pi6=str(coef),Crem_source_reused=True,ordered_triples_per_retained_direction=26970,missing_peak_loss_exact='1000/729',integrated_distance_coefficient=enclose(common),Mv=enclose(Mv),Mw=enclose(Mw),M_combined=enclose(M),third_distance_uniformly_positive_assumed=False,actual_Cauchy_density_integrated=True,omitted_parameter_events=0,Q2_used=False))
        begin(6)
        bad=A(Q(6,11))*M**A(Q(1,6))*pi**A(Q(-5,6))*eps**A(Q(11,6))
        total=good+bad;target=A(Q(4,48750));enc=enclose(total)
        upper=Q(int(enc['upper_mantissa']))*Q(2)**enc['upper_binary_exponent']
        simple=Q((upper*10**18).__ceil__(),10**18);assert upper<simple
        done(6,dict(name='one_middle_short_potential_fully_paid',combined_actual_complement_short_norm_upper=enclose(bad),total_short_norm_upper=enc,simple_rational_upper=str(simple),conditional_Bcrit='4/48750',total_over_conditional_Bcrit=enclose(total/target),strictly_below_Bcrit=bool(total<target),complement_fraction_of_upper_budget=enclose(bad/total),single_actual_complement_L1_budget=True,whole_face_or_full_time_claim=False))
        result['status']='COMPLETE_SIX_CONFIGURATIONS_TRUE_SHORT_MIDDLE_UPPER';save('results.json',result);receipt['status']='COMPLETED';code=0
    except BaseException as e:
        receipt['status']='STOPPED_FAILURE_NO_RETRY';receipt['error']=repr(e);receipt['traceback']=traceback.format_exc()
    finally:
        receipt['exit_code']=code;receipt['elapsed_seconds']=time.monotonic()-tic;receipt['finished_epoch']=time.time();receipt['no_active_children']=True
        if(HERE/'results.json').exists():receipt['results_sha256']=sha(HERE/'results.json')
        save('run_receipt.json',receipt)
    print(json.dumps(dict(status=receipt['status'],configurations=receipt['completed_configurations'],elapsed_seconds=receipt['elapsed_seconds'],error=receipt.get('error'))))
    return code
if __name__=='__main__':
    assert sys.argv[1:]==['--run'];raise SystemExit(run())
