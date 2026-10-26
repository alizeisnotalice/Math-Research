#!/usr/bin/env python3
"""Fixed new high-event budget; nine producer configurations, no old scalar rerun."""
import sys
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
from pathlib import Path
from fractions import Fraction as Q
from math import comb
import json,hashlib,time,os,traceback
import sympy as sp
from flint import arb,fmpq,ctx
ctx.prec=768
H=Path(__file__).resolve().parent
def rd(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(n,o):(H/n).write_text(json.dumps(o,indent=2)+'\n')
def A(x):
    x=Q(x);return arb(fmpq(x.numerator,x.denominator))
def bounds(r):return Q(int(r['lower_mantissa']))*Q(2)**int(r['lower_binary_exponent']),Q(int(r['upper_mantissa']))*Q(2)**int(r['upper_binary_exponent'])
def enc(x):
    assert x.is_finite()and x>0
    a,e=x.lower().man_exp();b,f=x.upper().man_exp()
    r=dict(lower_mantissa=str(a),lower_binary_exponent=int(e),upper_mantissa=str(b),upper_binary_exponent=int(f),display_only=str(x));l,u=bounds(r)
    assert 0<=u-l<max(Q(1),u)*Q(2)**-600
    return r
def run():
    assert not(H/'run_receipt.json').exists(),'one production batch only'
    reg=rd(H/'registration.json');release=rd(H/'main_release.json');binding=rd(H/'execution_binding.json')
    assert release['status']=='MAIN_RELEASED_FIXED_062_AFTER_THEORY_READ'
    assert release['registration_sha256']==binding['registration_sha256']==sha(H/'registration.json')
    assert release['theorem_sha256']==sha(H/'theorem.md')and binding['main_release_sha256']==sha(H/'main_release.json')
    for f,x in binding['code_sha256'].items():assert sha(H/f)==x,f
    for f,x in reg['source_sha256'].items():assert sha(f)==x,f
    assert reg['fixed_threshold_b']==release['fixed_b']=='1/4'and reg['fixed_complement_coordinate_lower_d']==release['fixed_d']=='31/64'
    guard=rd(H/'producer_process_receipt.json');assert guard['status']=='RUNNING'and guard['owned_pid']==os.getpid()and guard['wall_timeout_seconds']==180
    receipt=dict(status='RUNNING',pid=os.getpid(),started_epoch=time.time(),registration_sha256=sha(H/'registration.json'),main_release_sha256=sha(H/'main_release.json'),binding_sha256=sha(H/'execution_binding.json'),attempted_configurations=0,completed_configurations=0,program_timeout_seconds=180,old_M_recomputations=0,old_cap_budget_reruns=0,old060_reruns=0,new_kernel_calls=0,full_enumerations=0,grid_calls=0,matlab_calls=0)
    out=dict(status='INCOMPLETE',precision_bits=768,configurations=[],all_positive_time=True,new_global_I=False,new_C6_bound=False)
    save('run_receipt.json',receipt);tic=time.monotonic();code=1
    def begin(i):
        assert i==receipt['attempted_configurations']+1<=9 and time.monotonic()-tic<180
        receipt['attempted_configurations']=i;save('run_receipt.json',receipt)
    def done(i,**data):
        out['configurations'].append(dict(id=i,**data));receipt['completed_configurations']=i
        save('results.json',out);save('run_receipt.json',receipt)
    try:
        begin(1);rhoq=Q(1,64);bq=Q(1,4);aq=Q(47,64);dq=Q(31,64)
        assert aq-bq==dq>0 and aq==Q(3,4)-rhoq and 17**2>2
        assert min(h*(6-h)for h in range(1,6))==5
        assert [32*(31-k)for k in range(3)]==[992,960,928]
        r1=Q(reg['fixed_R1']);crem_exact=Q(44743079,2970)
        a041=Path(next(p for p in reg['source_sha256']if p.endswith('/round041/one_large_jump_cap/results.json')));old=rd(a041)
        assert any(c.get('R1_exact')==str(r1)for c in old['configurations'])
        assert any(c.get('Crem_times_pi6_exact')==str(crem_exact)for c in old['configurations'])
        done(1,name='new_positive_event_geometry_and_frozen_exact_coefficients',b=str(bq),d=str(dq),distance_squares=['6*d^2','10*d^2/3','2*d^2'],all_remaining_pairs_mixed_by_two_body_line_lemma=True,old_geometry_enumeration=False,R1_exact=str(r1),Crem_times_pi6_exact=str(crem_exact))
        begin(2);e,bb,C,pp,VV=sp.symbols('epsilon b C pi V',positive=True)
        tail=1/(3*bb**2*e**3)-1/(bb**4*e)+sp.atan(bb/e)/bb**5
        res_tail=sp.simplify(sp.diff(tail,e)+e**-4/(bb**2+e**2));assert res_tail==0
        partial=bb**-2*e**-4-bb**-4*e**-2+bb**-4/(bb**2+e**2)
        assert sp.factor(partial-e**-4/(bb**2+e**2))==0
        loose=sp.Rational(6,7)*C**sp.Rational(1,6)*pp**sp.Rational(-5,6)*e**sp.Rational(7,6)+VV**sp.Rational(5,6)*C/(3*e**3)
        res_opt=sp.simplify(sp.diff(loose,e).subs(C,e**5/(pp*VV)));assert res_opt==0
        done(2,name='general_b_complete_tail_and_unique_loose_optimizer',tail_derivative_residual=str(res_tail),loose_optimizer_residual=str(res_opt),epsilon_identity='epsilon^5=C*pi*V',exact_tail_optimized=False)
        begin(3);aa,bb,cc=sp.symbols('a b c',positive=True)
        I2=lambda x,y:sp.pi/(2*x*y*(x+y))
        J=sp.pi*(aa+bb+cc)/(2*aa*bb*cc*(aa+bb)*(aa+cc)*(bb+cc))
        res=sp.factor((I2(aa,cc)-I2(bb,cc))/(bb**2-aa**2)-J);assert res==0
        repeated=sp.pi*(2*aa+cc)/(4*aa**3*cc*(aa+cc)**2)
        assert sp.factor(J.subs(bb,aa)-repeated)==0
        assert sp.factor(repeated.subs(cc,aa)-3*sp.pi/(16*aa**5))==0
        done(3,name='three_positive_quadratics_full_time_identity',partial_fraction_residual=str(res),double_root_residual='0',triple_root_residual='0',positive_root_continuity=True)
        begin(4);ledger=[]
        for k in range(1,6):
            assert (sp.eye(k-1)-sp.ones(k-1,k-1)/k).det()==sp.Rational(1,k)
            ledger.append(dict(k=k,m=6-k,n=comb(5,k-1),Gram_determinant=str(Q(1,k)),cone_fraction=str(Q(1,k)),mean_interval_length='1/4'))
        assert [x['n']for x in ledger]==[1,5,10,10,5]and sum(x['n']for x in ledger)==31
        done(4,name='true_high_volume_Jacobians_and_overlapping_union',ledger=ledger,no_disjointness_assumption=True,restricted_P_symmetry_assumed=False)
        begin(5);pi=arb.pi();rho=A(rhoq);b=A(bq);d=A(dq)
        terms=[comb(5,k-1)*A(Q(6-k,k)).sqrt()for k in range(1,6)];shape=sum(terms,arb(0));V=pi*pi/8*rho**4*shape
        D=A(r1)/(32*pi**6);Cbig=D/(b*b);Crem=A(crem_exact)/pi**6;Vp=V**A(Q(5,6))
        out['constants']={k:enc(v)for k,v in dict(pi=pi,shape_sum=shape,V=V,V_power_5_over6=Vp,D=D,Cbig=Cbig,Crem=Crem).items()}
        out['volume_terms']=[dict(k=k,weighted_sqrt=enc(v))for k,v in enumerate(terms,1)]
        done(5,name='fixed_new_high_union_volume_and_event_constants')
        begin(6);eps=(Cbig*pi*V)**A(Q(1,5));Pshort=A(Q(6,7))*Cbig**A(Q(1,6))*pi**A(Q(-5,6))*eps**A(Q(7,6))
        assert eps>0 and eps/b<A(Q(1,2))
        out['epsilon']=enc(eps);out['Pshort']=enc(Pshort);out['epsilon_over_b']=enc(eps/b)
        done(6,name='one_analytic_epsilon_and_P_short',epsilon_scan=False)
        begin(7);T=1/(3*b*b*eps**3)-1/(b**4*eps)+(b/eps).atan()/b**5;looseT=1/(3*b*b*eps**3)
        assert 0<T<looseT
        Plong=Vp*D*T;out['Ptail_integral']=enc(T);out['Plong']=enc(Plong)
        done(7,name='complete_P_tail_to_infinity',tail_positive=True,strictly_below_loose_tail=True)
        begin(8);Avals=[992,960,928];cvals=[A(6).sqrt()*d,A(Q(10,3)).sqrt()*d,A(2).sqrt()*d];alphas=[c/A(a).sqrt()for c,a in zip(cvals,Avals)]
        x,y,z=alphas;JB=pi*(x+y+z)/(2*x*y*z*(x+y)*(x+z)*(y+z)*A(992*960*928));Bfull=Vp*Crem*JB
        out['B_quadratics']=dict(c=[enc(x)for x in cvals],alpha=[enc(x)for x in alphas],time_coefficients=Avals)
        out['B_time_integral']=enc(JB);out['Bfull']=enc(Bfull)
        done(8,name='complete_complement_B_all_positive_time',time_axis='0<t<infinity',all26970_ordered_triples_paid=True)
        begin(9);total=Pshort+Plong+Bfull;e=enc(total);lo,hi=bounds(e);simple=Q((hi*10**18).__ceil__()+1,10**18);target=Q(4,48750)
        outcome='UPPER_FORMULA_BELOW_CONDITIONAL_TARGET'if hi<target else('UPPER_FORMULA_EXCEEDS_CONDITIONAL_TARGET'if lo>target else'UNRESOLVED_COMPARISON')
        out['union']=dict(full_time_high_norm_upper=e,simple_rational_upper=str(simple),conditional_whole_face_target=str(target),formula_over_target=enc(total/A(target)),comparison=outcome,all_parameter_events_paid=True,all_positive_time_paid=True,one_whole_union_volume=True,cap_low_mid_exclusions_safely_relaxed=True,whole_face_paid=False,global_g_proved=False,new_C6_bound=False)
        done(9,name='complete_new_high_event_budget',comparison=outcome)
        out['status']='COMPLETE_NEW_HIGH_POSITIVE_EVENT_FULL_TIME_UPPER';save('results.json',out);receipt['status']='COMPLETED';code=0
    except BaseException as ex:
        receipt['status']='STOPPED_FAILURE_NO_RETRY';receipt['error']=repr(ex);receipt['traceback']=traceback.format_exc()
    finally:
        receipt['exit_code']=code;receipt['elapsed_seconds']=time.monotonic()-tic;receipt['finished_epoch']=time.time();receipt['no_active_children']=True
        if(H/'results.json').exists():receipt['results_sha256']=sha(H/'results.json')
        save('run_receipt.json',receipt)
    print(json.dumps(dict(status=receipt['status'],configurations=receipt['completed_configurations'],seconds=receipt['elapsed_seconds'],error=receipt.get('error'))));return code
if __name__=='__main__':
    assert sys.argv[1:]==['--run'];raise SystemExit(run())
