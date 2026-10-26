#!/usr/bin/env python3
"""Four new types, one frozen type, and one positive union budget; no kernels."""
import sys
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
from pathlib import Path
from fractions import Fraction as Q
from math import comb
import hashlib,json,time,datetime,os,traceback
import sympy as sp
from flint import arb,fmpq,ctx
ctx.prec=768
HERE=Path(__file__).resolve().parent;UPPER=HERE.parents[1]
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(name,obj):(HERE/name).write_text(json.dumps(obj,indent=2)+'\n')
def A(x):
    x=Q(x);return arb(fmpq(x.numerator,x.denominator))
def bounds(r):return Q(int(r['lower_mantissa']))*Q(2)**int(r['lower_binary_exponent']),Q(int(r['upper_mantissa']))*Q(2)**int(r['upper_binary_exponent'])
def load_ball(r):
    lo,hi=bounds(r);assert 0<lo<=hi
    x=arb(A((lo+hi)/2),A((hi-lo)/2));assert x.contains(A(lo)) and x.contains(A(hi));return x
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
    for rd in ['round043/two_parameter_middle','round049/middle_joint_marginal']:
        olddir=UPPER/rd
        for f,h in read(olddir/'manifest.json')['whitelist_sha256'].items():assert sha(olddir/f)==h,f
    deadline=datetime.datetime.fromisoformat(reg['deadline_utc']).timestamp();assert time.time()<deadline
    receipt=dict(status='RUNNING',pid=os.getpid(),started_epoch=time.time(),registration_sha256=sha(HERE/'registration.json'),attempted_configurations=0,completed_configurations=0,new_kernel_calls=0,full_enumerations=0,matlab_calls=0)
    result=dict(status='INCOMPLETE',precision_bits=768,configurations=[],cases={},new_I_certificate=False,new_C6_bound=False)
    save('run_receipt.json',receipt);tic=time.monotonic();code=1
    def begin(i):
        assert i==receipt['attempted_configurations']+1<=7 and time.time()<deadline
        receipt['attempted_configurations']=i;save('run_receipt.json',receipt)
    def done(i,data):
        result['configurations'].append(dict(id=i,**data));receipt['completed_configurations']=i
        save('results.json',result);save('run_receipt.json',receipt)
    try:
        peaks=read(UPPER/'round033/exact_projection_peaks/results.json');pc={r['id']:r for r in peaks['cases']}
        begin(1);geometry=[];p,q=sp.symbols('p q',positive=True);rhoq=Q(1,64)
        for k in range(1,6):
            m=6-k;P=sp.eye(k)-sp.ones(k,k)/k
            assert (P[:,0].T*P[:,0])[0]==sp.Rational(k-1,k)
            line=sp.simplify(k*p*p+m*q*q-(k*p+m*q)**2/6-sp.Rational(k*m,6)*(p-q)**2);assert line==0
            pairdet=36-(k-m)**2;assert pairdet==4*k*m==pc[f'r2_k{min(k,m)}']['Gram_determinant']
            gates={'c0v':Q(1,2)-rhoq,'c0w':Q(15,16)-rhoq,'c1v_squared_reduced':Q(m,96)-rhoq**2,'c1w_squared_reduced':Q(m,384)-rhoq**2,'g_squared':Q(1,32)-k*rhoq**2}
            assert all(v>0 for v in gates.values())
            geometry.append(dict(k=k,m=m,A_coordinate_functional_squared=str(Q(k-1,k)),line_identity='0',pair_Gram=pairdet,pair_type=min(k,m),positive_gate_margins={a:str(b)for a,b in gates.items()}))
        done(1,dict(name='all_five_true_geometry_and_positivity_gates',cases=geometry,no_spatial_enumeration=True,all_signed_lines_and_planes_handled_in_theorem=True))
        begin(2);joint=[]
        for k in range(1,6):
            m=6-k;P=sp.eye(m)-sp.ones(m,m)/m
            assert (P[:,0].T*P[:,0])[0]==sp.Rational(m-1,m)
            joint.append(dict(k=k,m=m,B_coordinate_functional_squared=str(Q(m-1,m)),n=comb(5,k-1)))
        assert [x['n']for x in joint]==[1,5,10,10,5] and sum(x['n']for x in joint)==31
        r=sp.Rational(1,64)
        assert sp.simplify(sp.sqrt(sp.Rational(9,6))/4-sp.sqrt(sp.Rational(3,32)))==0
        assert sp.simplify(sp.sqrt(sp.Rational(9,6))/8-sp.sqrt(sp.Rational(3,128)))==0
        old043=read(UPPER/'round043/two_parameter_middle/results.json')['configurations']
        old049=read(UPPER/'round049/middle_joint_marginal/results.json')['configurations'][2]
        assert old043[2]['rho']==old049['rho']=='1/64' and old043[2]['epsilon']==old049['epsilon']=='1/1000'
        assert old049['source_D_big_interval']==old043[3]['C2'] and old049['source_M_bad_interval']==old043[4]['M_combined']
        result['cases']['3']=dict(k=3,m=3,n=10,provenance='FROZEN_049_043_READ_ONLY_NO_RECOMPUTATION',D_big=old043[3]['C2'],M_combined=old043[4]['M_combined'],C_joint=old049['C_joint'],short_big_upper=old049['new_big_short_upper'],short_bad_upper=old049['new_bad_short_upper'],short_total_upper=old049['new_total_short_upper'],source043_sha256=sha(UPPER/'round043/two_parameter_middle/results.json'),source049_sha256=sha(UPPER/'round049/middle_joint_marginal/results.json'))
        done(2,dict(name='joint_range_weights_and_frozen_k3_identity',cases=joint,positive_slice_factor='1/16',joint_numerator=128,singleton_B_deviation_exact_zero=True,weights=[1,5,10,10,5],total_domains=31,k3_symbolic_identity=True,k3_new_scalar_evaluations=0))
        pi=None
        for index,k in enumerate([1,2,4,5],3):
            begin(index);m=6-k
            if pi is None:
                pi=arb.pi();rho=A(Q(1,64));eps=A(Q(1,1000));F=1+A(Q(6,928)).sqrt();Crem=A(Q(44743079,2970))/pi**6;loss=A(Q(1000,729))
                result['shared_constants']=dict(pi=enclose(pi),cauchy_average_factor=enclose(F),Crem=enclose(Crem),Crem_times_pi6_exact='44743079/2970',missing_peak_loss_exact='1000/729',rho='1/64',epsilon='1/1000')
            eta=A(k).sqrt()*rho
            hv0=A(k).sqrt()*(A(Q(1,2))-rho);hv1=A(Q(k*m,6)).sqrt()/4-eta
            hw0=A(k).sqrt()*(A(Q(15,16))-rho);hw1=A(Q(k*m,6)).sqrt()/8-eta
            g=A(Q(1,32)).sqrt()-eta;assert all(x>0 for x in [hv0,hv1,hw0,hw1,g])
            dev=A(Q(m-1,m)).sqrt()*rho;L=A(Q(1,4))-dev;U=A(Q(3,4))+dev
            assert dev>=0 and L>0 and U<1
            def primitive(z):return z/(2*(1-z*z))+((1+z)/(1-z)).log()/4
            C=128*(primitive(U)-primitive(L))/pi**2
            R2=Q(pc[f'r2_k{min(k,m)}']['R_exact']);D=64*A(R2)/pi**6
            common=loss*Crem*F/(g*g);Mv=common/(hv0*hv0*hv1*hv1);Mw=common/(hw0*hw0*hw1*hw1);M=Mv+Mw
            good=A(Q(3,7))*D**A(Q(1,6))*C**A(Q(5,6))*eps**A(Q(7,3))
            bad=A(Q(3,8))*M**A(Q(1,6))*C**A(Q(5,6))*eps**A(Q(8,3));total=good+bad
            enc=enclose(total);hi=bounds(enc)[1];simple=Q((hi*10**20).__ceil__()+1,10**20)
            row=dict(k=k,m=m,n=comb(5,k-1),provenance='NEW_FIXED_TYPE_CONFIGURATION',pair_type=min(k,m),R2_exact=str(R2),eta=enclose(eta),hv0=enclose(hv0),hv1=enclose(hv1),hw0=enclose(hw0),hw1=enclose(hw1),g=enclose(g),joint_coordinate_deviation=enclose(dev),z_lower=enclose(L),z_upper=enclose(U),C_joint=enclose(C),D_big=enclose(D),Mv=enclose(Mv),Mw=enclose(Mw),M_combined=enclose(M),short_big_upper=enclose(good),short_bad_upper=enclose(bad),short_total_upper=enc,simple_rational_case_upper=str(simple),all_complementary_events_paid=True)
            result['cases'][str(k)]=row;done(index,dict(name=f'new_type_k{k}_m{m}',case_key=str(k)))
        begin(7);terms=[];summed=arb(0)
        for k in range(1,6):
            row=result['cases'][str(k)];n=comb(5,k-1);uk=load_ball(row['short_total_upper']);term=n*uk**A(Q(6,5));summed+=term
            terms.append(dict(k=k,n=n,weighted_power=enclose(term),case_budget_provenance=row['provenance']))
        union=summed**A(Q(5,6));enc=enclose(union);hi=bounds(enc)[1];simple=Q((hi*10**18).__ceil__()+1,10**18)
        result['union']=dict(domain_count=31,weights=[1,5,10,10,5],q='6/5',weighted_terms=terms,q_power_sum_upper=enclose(summed),short_union_norm_upper=enc,simple_rational_union_upper=str(simple),conditional_Bcrit='4/48750',union_over_Bcrit=enclose(union/A(Q(4,48750))),overlap_permitted=True,domains_claimed_disjoint=False,spatial_endpoints_paid=False,long_time_paid=False)
        done(7,dict(name='all31_domains_positive_union_budget',domain_count=31,k3_new_scalar_evaluations=0))
        result['status']='COMPLETE_ALL31_UNSIGNED_MIDDLES_SHORT_UPPER';save('results.json',result);receipt['status']='COMPLETED';code=0
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
