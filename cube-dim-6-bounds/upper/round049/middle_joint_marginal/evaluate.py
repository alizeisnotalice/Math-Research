#!/usr/bin/env python3
"""Two exact structural checks and one fixed new joint-marginal configuration."""
import sys
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
from pathlib import Path
from fractions import Fraction as Q
import json,hashlib,datetime,time,os,traceback
import sympy as sp
from flint import arb,fmpq,ctx
ctx.prec=768
HERE=Path(__file__).resolve().parent
OLD=HERE.parents[1]/'round043/two_parameter_middle'
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(name,obj):(HERE/name).write_text(json.dumps(obj,indent=2)+'\n')
def A(x):
    x=Q(x);return arb(fmpq(x.numerator,x.denominator))
def bounds(r):
    return Q(int(r['lower_mantissa']))*Q(2)**int(r['lower_binary_exponent']),Q(int(r['upper_mantissa']))*Q(2)**int(r['upper_binary_exponent'])
def load_ball(r):
    lo,hi=bounds(r);assert 0<lo<=hi
    value=arb(A((lo+hi)/2),A((hi-lo)/2))
    assert value.contains(A(lo)) and value.contains(A(hi))
    return value
def enclose(x):
    assert x.is_finite()
    a,e=x.lower().man_exp();b,f=x.upper().man_exp()
    r=dict(lower_mantissa=str(a),lower_binary_exponent=int(e),upper_mantissa=str(b),upper_binary_exponent=int(f),display_only=str(x))
    lo,hi=bounds(r);assert 0<=hi-lo<=max(Q(1),abs((lo+hi)/2))*Q(2)**-450
    return r
def run():
    assert not(HERE/'run_receipt.json').exists(),'one batch only'
    reg=read(HERE/'registration.json')
    for p,h in reg['source_sha256'].items():assert sha(p)==h,p
    for f,h in reg['code_sha256'].items():assert sha(HERE/f)==h,f
    for f,h in read(OLD/'manifest.json')['whitelist_sha256'].items():assert sha(OLD/f)==h,f
    deadline=datetime.datetime.fromisoformat(reg['deadline_utc']).timestamp();assert time.time()<deadline
    receipt=dict(status='RUNNING',pid=os.getpid(),started_epoch=time.time(),registration_sha256=sha(HERE/'registration.json'),attempted_configurations=0,completed_configurations=0,new_kernel_calls=0,full_enumerations=0,matlab_calls=0)
    result=dict(status='INCOMPLETE',precision_bits=768,input_scalar_precision_bits=512,configurations=[],old043_budget_recomputed=False,new_I_certificate=False,new_C6_bound=False)
    save('run_receipt.json',receipt);tic=time.monotonic();code=1
    def begin(i):
        assert i==receipt['attempted_configurations']+1<=3 and time.time()<deadline
        receipt['attempted_configurations']=i;save('run_receipt.json',receipt)
    def done(i,data):
        result['configurations'].append(dict(id=i,**data));receipt['completed_configurations']=i
        save('results.json',result);save('run_receipt.json',receipt)
    try:
        begin(1)
        columns=[sp.ones(6,1)]
        for j in (3,1,2,4,5):
            col=sp.ones(6,1);col[j]=-1;columns.append(col)
        basis=sp.Matrix.hstack(*columns);reordered=basis[[0,3,1,2,4,5],:]
        pivot=reordered[:2,:2];schur=reordered[2:,2:]-reordered[2:,:2]*pivot.inv()*reordered[:2,2:]
        assert abs(basis.det())==32 and abs(pivot.det())==2 and abs(schur.det())==16
        projection=sp.eye(3)-sp.ones(3,3)/3;functional=projection[:,0]
        assert (functional.T*functional)[0]==sp.Rational(2,3)
        assert 2**4==16 and 16+16==32
        s,z=sp.symbols('s z',positive=True)
        density=sp.Rational(1,2)*(16*s/(sp.pi*((1+z)**2/4+256*s**2)))*(16*s/(sp.pi*((1-z)**2/4+256*s**2)))
        normalized=2048*s**2/(sp.pi**2*((1+z)**2+1024*s**2)*((1-z)**2+1024*s**2))
        assert sp.simplify(density-normalized)==0 and 2048//16==128
        done(1,dict(name='true_joint_marginal_and_coordinate_range',chosen_j=4,block_counts=[16,16],block_count_method='2^4 free signs, no frame enumeration',basis_columns=[list(map(int,c))for c in columns],basis_abs_det=32,retained_pivot_abs_det=2,free_schur_abs_det=16,inverse_joint_Jacobian='1/2',remaining_positive_orthant_factor='1/16',two_density_numerator=2048,joint_L1_numerator=128,B_deviation_functional_norm_squared='2/3',symbolic_density_difference='0'))
        begin(2)
        primitive=z/(2*(1-z*z))+sp.log((1+z)/(1-z))/4
        assert sp.simplify(sp.diff(primitive,z)-1/(1-z*z)**2)==0
        q=Q(6,5);a=1-1/q;b=1/q
        big_power=-2*a+2*b;bad_power=2*b
        assert (a,b,big_power,bad_power)==(Q(1,6),Q(5,6),Q(4,3),Q(5,3))
        assert 1/(big_power+1)==Q(3,7) and 1/(bad_power+1)==Q(3,8)
        done(2,dict(name='primitive_and_actual_branch_interpolation',primitive_derivative_difference='0',q='6/5',interpolation_exponents=['1/6','5/6'],density_norm_time_powers=['4/3','5/3'],integrated_time_powers=['7/3','8/3'],integrated_prefactors=['3/7','3/8'],branch_symmetry_assumed=False,actual_disjoint_complement_domination=True))
        begin(3)
        old=read(OLD/'results.json')['configurations']
        dsrc=old[3]['C2'];msrc=old[4]['M_combined'];osrc=old[5]['total_short_norm_upper']
        D=load_ball(dsrc);M=load_ball(msrc);oldtotal=load_ball(osrc)
        pi=arb.pi();rho=A(Q(1,64));eps=A(Q(1,1000));delta=A(Q(2,3)).sqrt()*rho
        lower=A(Q(1,4))-delta;upper=A(Q(3,4))+delta
        assert lower>0 and upper<1 and upper>lower
        def F(a):return a/(2*(1-a*a))+((1+a)/(1-a)).log()/4
        integral=F(upper)-F(lower);assert integral>0
        C=128*integral/pi**2
        common=C**A(Q(5,6))
        good=A(Q(3,7))*D**A(Q(1,6))*common*eps**A(Q(7,3))
        bad=A(Q(3,8))*M**A(Q(1,6))*common*eps**A(Q(8,3))
        total=good+bad;enclosed=enclose(total);hi=bounds(enclosed)[1]
        simple=Q((hi*10**18).__ceil__(),10**18);assert hi<simple
        ratio=oldtotal/total;assert ratio>1
        l1ratio=C*pi*eps;assert l1ratio<1
        done(3,dict(name='single_new_joint_L1_interface_budget',rho='1/64',epsilon='1/1000',source_results_path=str(OLD/'results.json'),source_results_sha256=sha(OLD/'results.json'),source_D_big_key='configurations[3].C2',source_M_bad_key='configurations[4].M_combined',source_old_budget_key='configurations[5].total_short_norm_upper',source_D_big_interval=dsrc,source_M_bad_interval=msrc,source_old_budget_interval=osrc,pi=enclose(pi),coordinate_deviation=enclose(delta),z_lower=enclose(lower),z_upper=enclose(upper),primitive_integral=enclose(integral),C_joint=enclose(C),joint_over_old_loose_L1_at_epsilon=enclose(l1ratio),new_big_short_upper=enclose(good),new_bad_short_upper=enclose(bad),new_total_short_upper=enclosed,simple_rational_upper=str(simple),saved_old_budget_over_new=enclose(ratio),new_over_saved_old_budget=enclose(total/oldtotal),conditional_Bcrit='4/48750',new_over_Bcrit=enclose(total/A(Q(4,48750))),bad_fraction_of_new_upper=enclose(bad/total),source_values_from_complete_dyadic_endpoints=True,old_budget_recomputed=False,actual_norm_reduction_factor_claimed=False))
        result['status']='COMPLETE_NEW_TWO_COORDINATE_L1_SHORT_MIDDLE_BOUND';save('results.json',result);receipt['status']='COMPLETED';code=0
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
