#!/usr/bin/env python3
"""Five fixed full-time endpoint formulas; reads old M without recomputing it."""
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

def enclose(x):
    assert x.is_finite()
    a,e=x.lower().man_exp();b,f=x.upper().man_exp()
    r=dict(lower_mantissa=str(a),lower_binary_exponent=int(e),upper_mantissa=str(b),upper_binary_exponent=int(f),display_only=str(x))
    l,u=bounds(r);assert 0<=u-l<max(Q(1),abs((l+u)/2))*Q(2)**-450
    return r

def run():
    assert not(H/'run_receipt.json').exists(),'one batch only'
    reg=rd(H/'registration.json')
    for p,h in reg['source_sha256'].items():assert sha(p)==h,p
    for p,h in reg['code_sha256'].items():assert sha(H/p)==h,p
    for f,h in rd(OLD/'manifest.json')['whitelist_sha256'].items():assert sha(OLD/f)==h,f
    dep=rd(OLD/'audit_dependency.json');assert dep['status']=='F035_PASS_DEPENDENCY_DISCHARGED'
    for p,h in dep['sha256'].items():assert sha(p)==h,p
    deadline=datetime.datetime.fromisoformat(reg['deadline_utc']).timestamp();assert time.time()<deadline
    receipt=dict(status='RUNNING',pid=os.getpid(),started_epoch=time.time(),registration_sha256=sha(H/'registration.json'),attempted_configurations=0,completed_configurations=0,new_kernel_calls=0,old_M_recomputations=0,full_enumerations=0,matlab_calls=0)
    result=dict(status='INCOMPLETE',precision_bits=768,configurations=[],cases={},new_I_certificate=False,new_C6_bound=False)
    save('run_receipt.json',receipt);tic=time.monotonic();code=1
    def begin(i):
        assert i==receipt['attempted_configurations']+1<=7 and time.time()<deadline
        receipt['attempted_configurations']=i;save('run_receipt.json',receipt)
    def done(i,data):
        result['configurations'].append(dict(id=i,**data));receipt['completed_configurations']=i
        save('results.json',result);save('run_receipt.json',receipt)
    try:
        begin(1);N=[192,480,320];rhoq=Q(1,64);gates=[]
        for k in range(1,6):
            m=6-k;j=min(k,m);c=2**m;bg=[n-(c if ell==j else 0)for ell,n in enumerate(N,1)]
            assert all(x>=0 for x in bg)and sum(bg)==992-c
            matrix=sp.eye(k-1)-sp.ones(k-1,k-1)/k
            assert matrix.det()==sp.Rational(1,k)
            # Positivity after factoring sqrt(k), and the pair gap square.
            assert Q(3*m,32)>rhoq**2 and Q(9,32)>k*rhoq**2
            assert sp.Rational(4)-sp.Rational(6,5)==sp.Rational(14,5)
            gates.append(dict(k=k,m=m,pair_class=j,incident_planes=2**(m-1),incident_ordered_pairs=c,background_by_class=bg,A_normal_Gram_determinant=str(Q(1,k)),normal_dimension=4,cone_fraction=str(Q(1,k)),line_gap_squared_margin=str(Q(3*m,32)-rhoq**2),pair_gap_squared_margin=str(Q(9,32)-k*rhoq**2),n=comb(5,k-1)))
        assert sum(x['n']for x in gates)==31
        old=rd(OLD/'results.json');raw={str(x['k']):x['M']for x in old['configurations'][:3]}
        assert set(raw)=={'1','2','3'}
        result['source_M_intervals']=raw
        result['source_dependency_status']=dep['status']
        done(1,dict(name='exact_incidence_counts_Jacobian_gates_and_positive_accounting',cases=gates,all_ordered_pairs=992,no_direction_enumeration=True))
        pi=None
        for k in range(1,6):
            begin(k+1);m=6-k;j=min(k,m);c=2**m;n=comb(5,k-1)
            if pi is None:
                pi=arb.pi();rho=A(rhoq);M={int(j):load(x)for j,x in raw.items()}
                result['shared_constants']=dict(pi=enclose(pi),rho='1/64',q='6/5',mean_length='1/4',ordered_class_counts=N,source_M_precision_bits=old['precision_bits'])
            eta=A(k).sqrt()*rho;r0=A(k).sqrt()*(1-rho)
            d0=A(Q(3,4))*A(Q(k*m,6)).sqrt()-eta;g=A(Q(9,32)).sqrt()-eta
            assert r0>0 and d0>0 and g>0
            denom=r0*r0*d0*d0
            incident=c*M[j]/denom
            bgcounts=[z-(c if ell==j else 0)for ell,z in enumerate(N,1)]
            numerator=sum((a*M[ell]for ell,a in enumerate(bgcounts,1)),arb(0))
            background=numerator/(denom*g)
            measurefactor=A(Q(1,4))*A(Q(m,k)).sqrt()
            V=measurefactor*pi*pi/2*rho**4
            J=measurefactor*2*pi*pi/A(Q(14,5))*rho**A(Q(14,5))
            singular=incident*J**A(Q(5,6));regular=background*V**A(Q(5,6));total=singular+regular
            enc=enclose(total);hi=bounds(enc)[1];simple=Q((hi*10**18).__ceil__()+1,10**18)
            row=dict(k=k,m=m,n=n,pair_class=j,incident_planes=2**(m-1),incident_ordered_pairs=c,background_by_class=bgcounts,background_ordered_pairs=sum(bgcounts),eta=enclose(eta),r0=enclose(r0),d0=enclose(d0),g=enclose(g),incident_coefficient=enclose(incident),background_coefficient=enclose(background),face_measure_factor=enclose(measurefactor),V_upper=enclose(V),J_upper=enclose(J),singular_norm_upper=enclose(singular),background_norm_upper=enclose(regular),total_norm_upper=enc,simple_rational_upper=str(simple),all_pairs_paid=True,time_range='all t>0')
            result['cases'][str(k)]=row;done(k+1,dict(name=f'fixed_full_time_endpoint_type_k{k}_m{m}',case_key=str(k)))
        begin(7);S=arb(0);terms=[]
        for k in range(1,6):
            row=result['cases'][str(k)];term=row['n']*load(row['total_norm_upper'])**A(Q(6,5));S+=term
            terms.append(dict(k=k,n=row['n'],weighted_power=enclose(term)))
        union=S**A(Q(5,6));enc=enclose(union);hi=bounds(enc)[1];lo=bounds(enc)[0];simple=Q((hi*10**18).__ceil__()+1,10**18)
        target=Q(4,48750)
        outcome='UPPER_FORMULA_BELOW_CONDITIONAL_WHOLE_TARGET'if hi<target else('UPPER_FORMULA_EXCEEDS_CONDITIONAL_WHOLE_TARGET'if lo>target else'UNRESOLVED_COMPARISON')
        result['union']=dict(domain_count=31,weights=[1,5,10,10,5],q='6/5',weighted_terms=terms,q_power_sum=enclose(S),full_time_union_norm_upper=enc,simple_rational_union_upper=str(simple),conditional_Bcrit=str(target),union_over_Bcrit=enclose(union/A(target)),budget_comparison=outcome,actual_G_lower_bound=False,overlap_permitted=True,domains_claimed_disjoint=False,high_endpoint_paid=False,middle_long_time_paid=False,whole_face_paid=False)
        done(7,dict(name='all31_low_endpoint_full_time_positive_union',domain_count=31,budget_comparison=outcome))
        result['status']='COMPLETE_ALL31_LOW_ENDPOINTS_FULL_TIME_UPPER';save('results.json',result);receipt['status']='COMPLETED';code=0
    except BaseException as e:
        receipt['status']='STOPPED_FAILURE_NO_RETRY';receipt['error']=repr(e);receipt['traceback']=traceback.format_exc()
    finally:
        receipt['exit_code']=code;receipt['elapsed_seconds']=time.monotonic()-tic;receipt['finished_epoch']=time.time();receipt['no_active_children']=True
        if(H/'results.json').exists():receipt['results_sha256']=sha(H/'results.json')
        save('run_receipt.json',receipt)
    print(json.dumps(dict(status=receipt['status'],configurations=receipt['completed_configurations'],elapsed_seconds=receipt['elapsed_seconds'],error=receipt.get('error'))));return code
if __name__=='__main__':
    assert sys.argv[1:]==['--run'];raise SystemExit(run())
