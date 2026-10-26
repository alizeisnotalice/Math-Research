#!/usr/bin/env python3
"""One fixed high-endpoint joint radial-angular budget, no kernel or M integral."""
import sys
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
from pathlib import Path
from fractions import Fraction as Q
from math import comb
import json,hashlib,datetime,time,os,traceback
import sympy as sp
from flint import arb,fmpq,ctx
ctx.prec=768
H=Path(__file__).resolve().parent
OLD=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-audit/cube-dim-6-bounds/upper/round037/pair_resolved_green_envelope')
def rd(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(n,o):(H/n).write_text(json.dumps(o,indent=2)+'\n')
def A(x):
    x=Q(x);return arb(fmpq(x.numerator,x.denominator))
def bounds(r):
    if 'common_binary_exponent'in r:
        e=int(r['common_binary_exponent']);return Q(int(r['lower_mantissa']))*Q(2)**e,Q(int(r['upper_mantissa']))*Q(2)**e
    return Q(int(r['lower_mantissa']))*Q(2)**int(r['lower_binary_exponent']),Q(int(r['upper_mantissa']))*Q(2)**int(r['upper_binary_exponent'])
def load(r):
    l,u=bounds(r);assert 0<l<=u
    z=arb(A((l+u)/2),A((u-l)/2));assert z.contains(A(l))and z.contains(A(u));return z
def enc(x):
    assert x.is_finite()and x>0
    a,e=x.lower().man_exp();b,f=x.upper().man_exp()
    row=dict(lower_mantissa=str(a),lower_binary_exponent=int(e),upper_mantissa=str(b),upper_binary_exponent=int(f),display_only=str(x))
    l,u=bounds(row);assert 0<=u-l<max(Q(1),u)*Q(2)**-450
    return row
def run():
    assert not(H/'run_receipt.json').exists(),'one production batch only'
    reg=rd(H/'registration.json');binding=rd(H/'execution_binding.json')
    assert binding['registration_sha256']==sha(H/'registration.json')
    for f,x in binding['code_sha256'].items():assert sha(H/f)==x,f
    for f,x in reg['source_sha256'].items():assert sha(f)==x,f
    for f,x in rd(OLD/'manifest.json')['whitelist_sha256'].items():assert sha(OLD/f)==x,f
    dep=rd(OLD/'audit_dependency.json');assert dep['status']=='F035_PASS_DEPENDENCY_DISCHARGED'
    for f,x in dep['sha256'].items():assert sha(f)==x,f
    deadline=datetime.datetime.fromisoformat(reg['deadline_utc']).timestamp();assert time.time()<deadline
    receipt=dict(status='RUNNING',pid=os.getpid(),started_epoch=time.time(),deadline_epoch=deadline,registration_sha256=sha(H/'registration.json'),binding_sha256=sha(H/'execution_binding.json'),attempted_configurations=0,completed_configurations=0,old_M_recomputations=0,new_kernel_calls=0,old_cap_reruns=0,full_enumerations=0,grid_calls=0,matlab_calls=0)
    out=dict(status='INCOMPLETE',precision_bits=768,configurations=[],source_M_intervals={},new_global_I=False,new_C6_bound=False)
    save('run_receipt.json',receipt);tic=time.monotonic();code=1
    def begin(i):
        assert i==receipt['attempted_configurations']+1<=9 and time.time()<deadline
        receipt['attempted_configurations']=i;save('run_receipt.json',receipt)
    def done(i,**data):
        out['configurations'].append(dict(id=i,**data));receipt['completed_configurations']=i
        save('results.json',out);save('run_receipt.json',receipt)
    try:
        begin(1);n=[comb(5,1)+comb(5,5),comb(5,2)+comb(5,4),comb(5,3)];N=[192,480,320]
        assert n==[6,15,10]and sum(n)==31 and sum(N)==992
        bg=[x-2*y for x,y in zip(N,n)];assert bg==[180,450,300]and bg==[30*y for y in n]
        rhoq=Q(1,64);aq=Q(47,64);R2=Q(385,4096);r02=1+5*aq**2;D2=Q(10,3)*aq**2
        assert R2==rhoq**2+Q(3,32)and 31*rhoq**4<R2**2 and 31>1
        assert 17**2>2 and 1-rhoq>Q(3,4)
        old=rd(OLD/'results.json');raw={str(x['k']):x['M']for x in old['configurations'][:3]};assert set(raw)=={'1','2','3'}
        out['source_M_intervals']=raw
        done(1,name='incident_ones_and_all_background_exact_counts',incident_planes=31,incident_ordered=62,ones_first_counts=n,reverse_counts=n,background_counts=bg,total_ordered=992,no_full_enumeration=True)
        begin(2);G=sp.eye(5)-sp.ones(5,5)/6;assert G.det()==sp.Rational(1,6)
        assert max(k*(6-k)for k in range(1,6))==9
        assert min(h*(6-h)for h in range(1,6))==5
        done(2,name='projected_face_Jacobian_and_radius',Gram_determinant='1/6',Jacobian='sqrt(6)',symmetric_cone_fraction='1/6',R_squared=str(R2),r0_squared=str(r02),D_squared=str(D2),g_squared=str(2*aq**2))
        begin(3);c=sp.symbols('c',nonnegative=True)
        cdf=1-sp.Rational(3,2)*c+sp.Rational(1,2)*c**3
        assert sp.expand(cdf-(1-c)**2*(c+2)/2)==0
        assert sp.expand(2*(1+c)**2-(c+2)-(3*c+2*c**2))==0
        assert 1-Q(6,5)/4==Q(7,10)
        done(3,name='angular_CDF_and_layercake_coefficients',CDF_factor='(1-c)^2(c+2)/2',nonnegative_difference='3c+2c^2',cover_count=31,rearrangement_factor='10/7',rearrangement_exponent='7/10')
        begin(4);rr,ss,RR=sp.symbols('rho s0 R',positive=True);b=sp.Rational(7,5)
        p3=(ss**b-rr**b)/b+ss**(2*b)*(ss**(-b)-RR**(-b))/b
        closed3=(2*ss**b-rr**b-ss**(2*b)/RR**b)/b
        p1=(ss**sp.Rational(19,5)-rr**sp.Rational(19,5))/sp.Rational(19,5)+ss**sp.Rational(14,5)*(RR-ss)
        closed1=ss**sp.Rational(14,5)*RR-sp.Rational(14,19)*ss**sp.Rational(19,5)-sp.Rational(5,19)*rr**sp.Rational(19,5)
        pv=(ss**5-rr**5)/5+ss**4*(RR-ss);closedv=ss**4*(RR-sp.Rational(4,5)*ss)-rr**5/5
        residuals=[sp.simplify(p3-closed3),sp.simplify(p1-closed1),sp.simplify(pv-closedv)];assert residuals==[0,0,0]
        done(4,name='fixed_radial_primitive_identities',symbolic_residuals=[str(x)for x in residuals],joint_singularity_exponent='s^(-3q)',no_line_distance_rho_collapse=True)
        begin(5);pi=arb.pi();rho=A(rhoq);a=A(aq);R=A(R2).sqrt();s0=A(31).root(4)*rho;sqrt6=A(6).sqrt();sqrt2=A(2).sqrt();g=sqrt2*a
        assert rho<s0<R and g>0
        M={int(k):load(v)for k,v in raw.items()};W=sum((z*M[k]for k,z in enumerate(n,1)),arb(0))
        rb=rho**A(Q(7,5));sb=s0**A(Q(7,5));Rb=R**A(Q(7,5));sphere=8*pi*pi/3;K=sphere/sqrt6*A(Q(10,7))
        out['constants']={k:enc(v)for k,v in dict(pi=pi,sqrt6=sqrt6,sqrt2=sqrt2,R=R,s0=s0,rho_power_7_over5=rb,s0_power_7_over5=sb,R_power_7_over5=Rb,W=W,sphere_area=sphere,K=K,g=g).items()}
        out['rational_constants']=dict(rho=str(rhoq),a=str(aq),R_squared=str(R2),r0_squared=str(r02),D_squared=str(D2))
        done(5,name='fixed_constants_and_complete_source_interval_loading',source_M_precision_bits=old['precision_bits'])
        begin(6);J3=K*A(Q(5,7))*(2*sb-rb-sb*sb/Rb);assert J3>0
        out['J3']=enc(J3);done(6,name='joint_radial_angular_J3')
        begin(7);J1=K*(sb*sb*R-A(Q(14,19))*sb*sb*s0-A(Q(5,19))*rb*rb*rho);assert J1>0
        out['J1']=enc(J1);done(7,name='reverse_order_J1')
        begin(8);V=sphere/sqrt6*(31*rho**4*(R-A(Q(4,5))*s0)-rho**5/5);assert V>0
        out['V']=enc(V);done(8,name='union_angular_cover_volume')
        begin(9);coef=[W/A(r02),W/A(r02*D2),30*W/(A(r02*D2)*g)];ints=[J3,J1,V];names=['ones_first_joint','reverse_order','all930_background']
        vals=[c*j**A(Q(5,6))for c,j in zip(coef,ints)];total=sum(vals,arb(0));e=enc(total);lo,hi=bounds(e);simple=Q((hi*10**18).__ceil__()+1,10**18);target=Q(4,48750)
        out['branches']={name:dict(coefficient=enc(c),integral_upper=enc(j),norm_upper=enc(v))for name,c,j,v in zip(names,coef,ints,vals)}
        outcome='FIXED_UPPER_FORMULA_EXCEEDS_TARGET'if lo>target else('FIXED_UPPER_FORMULA_BELOW_TARGET'if hi<target else'UNRESOLVED_COMPARISON')
        out['union']=dict(full_time_high_norm_upper=e,simple_rational_upper=str(simple),conditional_whole_face_target=str(target),formula_over_target=enc(total/A(target)),comparison=outcome,all992_terms_paid=True,all31_plane_overlap_retained=True,cap_subtracted_by_s_ge_rho=True,low_mid_subtractions_safely_relaxed=True,actual_G_lower_bound=False,original_Rc_excluded=False,whole_face_paid=False,new_C6_bound=False)
        done(9,name='complete_same_model_full_time_high_budget',comparison=outcome)
        out['status']='COMPLETE_FIXED_HIGH_JOINT_RADIAL_ANGULAR_UPPER';save('results.json',out);receipt['status']='COMPLETED';code=0
    except BaseException as e:
        receipt['status']='STOPPED_FAILURE_NO_RETRY';receipt['error']=repr(e);receipt['traceback']=traceback.format_exc()
    finally:
        receipt['exit_code']=code;receipt['elapsed_seconds']=time.monotonic()-tic;receipt['finished_epoch']=time.time();receipt['no_active_children']=True
        if(H/'results.json').exists():receipt['results_sha256']=sha(H/'results.json')
        save('run_receipt.json',receipt)
    print(json.dumps(dict(status=receipt['status'],configurations=receipt['completed_configurations'],seconds=receipt['elapsed_seconds'],error=receipt.get('error'))));return code
if __name__=='__main__':
    assert sys.argv[1:]==['--run'];raise SystemExit(run())
