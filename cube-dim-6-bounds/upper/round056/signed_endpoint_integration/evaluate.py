#!/usr/bin/env python3
"""Gated fixed signed endpoint sharpening; old complete formulas remain inputs."""
import sys
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
from pathlib import Path
from fractions import Fraction as Q
from math import comb
import json,hashlib,datetime,time,os,traceback
from flint import arb,fmpq,ctx
ctx.prec=768
H=Path(__file__).resolve().parent;OLD=H.parents[1]/'round054/lower_endpoint_green_budget'
def rd(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(n,o):(H/n).write_text(json.dumps(o,indent=2)+'\n')
def A(x):
    x=Q(x);return arb(fmpq(x.numerator,x.denominator))
def bounds(r):return Q(int(r['lower_mantissa']))*Q(2)**int(r['lower_binary_exponent']),Q(int(r['upper_mantissa']))*Q(2)**int(r['upper_binary_exponent'])
def load(r):
    l,u=bounds(r);assert 0<l<=u
    x=arb(A((l+u)/2),A((u-l)/2));assert x.contains(A(l))and x.contains(A(u));return x

def enc(x):
    assert x.is_finite()
    a,e=x.lower().man_exp();b,f=x.upper().man_exp()
    r=dict(lower_mantissa=str(a),lower_binary_exponent=int(e),upper_mantissa=str(b),upper_binary_exponent=int(f),display_only=str(x));l,u=bounds(r)
    assert 0<=u-l<max(Q(1),abs((u+l)/2))*Q(2)**-450;return r

def run():
    assert not(H/'run_receipt.json').exists(),'one actual producer batch only'
    reg=rd(H/'continuation_registration.json');legacy=rd(H/'registration.json');release=rd(H/'gate_release.json')
    assert release['final_F054_PASS']is True,'F054 gate remains closed'
    assert release['parent_final_release']is True,'parent final release remains closed'
    assert reg['predecessor']['registration_sha256']==sha(H/'registration.json')
    assert reg['fixed_plan']==legacy['fixed_plan'] and reg['max_configurations']==10
    assert release['continuation_registration_sha256']==sha(H/'continuation_registration.json')
    for p,h in reg['source_sha256'].items():assert sha(p)==h,p
    for p,h in release['audit_source_sha256'].items():assert sha(p)==h,p
    for p,h in release['production_code_sha256'].items():assert sha(H/p)==h,p
    assert set(release['production_code_sha256'])=={'theorem.md','evaluate.py','check_saved.py'}
    assert reg['theorem_sha256']==sha(H/'theorem.md')
    for f,h in rd(OLD/'manifest.json')['whitelist_sha256'].items():assert sha(OLD/f)==h,f
    assert release['bound_A054_manifest_sha256']==sha(OLD/'manifest.json')
    assert reg['scientific_execution']['max_wall_seconds']==1200
    started=time.time();deadline=started+1200
    receipt=dict(status='RUNNING',pid=os.getpid(),started_epoch=started,scientific_deadline_epoch=deadline,scientific_budget_seconds=1200,registration_sha256=sha(H/'continuation_registration.json'),predecessor_registration_sha256=sha(H/'registration.json'),gate_release_sha256=sha(H/'gate_release.json'),attempted_configurations=0,completed_configurations=0,new_kernel_calls=0,old_M_recomputations=0,old_full_budget_recomputations=0,k5_new_evaluations=0,full_enumerations=0,matlab_calls=0)
    result=dict(status='INCOMPLETE',precision_bits=768,configurations=[],R={},cases={},new_C6_bound=False)
    save('run_receipt.json',receipt);tic=time.monotonic();code=1
    def begin(i):
        assert i==receipt['attempted_configurations']+1<=9 and time.time()<deadline
        receipt['attempted_configurations']=i;save('run_receipt.json',receipt)
    def done(i,data):
        result['configurations'].append(dict(id=i,**data));receipt['completed_configurations']=i;save('results.json',result);save('run_receipt.json',receipt)
    try:
        begin(1);q=Q(6,5)
        linear=(4-q)*q/((q-1)*(5-q));power=(4-q)/(4*(q-1))
        assert linear==Q(84,19)and power==Q(7,2)and linear-power==Q(35,38)and linear-q*power==Q(21,95)
        assert Q(1,64)<Q(1,4)
        done(1,dict(name='exact_radial_ratio_coefficients_and_positivity',q=str(q),linear=str(linear),power=str(power),R_at1=str(linear-power),derivative_lower=str(linear-q*power),normal_dimension=4,positive_B_bound='distance^2 >= (h/m)||b||^2 >= h*u^2',clipping_enlarged_only_after_positive_pointwise_bound=True))
        begin(2);old=rd(OLD/'results.json');ledger=[]
        assert old['status']=='COMPLETE_ALL31_LOW_ENDPOINTS_FULL_TIME_UPPER'
        for k in range(1,6):
            m=6-k;N={}
            for h in range(1,m//2+1):
                z=comb(m,h)
                if 2*h==m:assert z%2==0;z//=2
                N[str(h)]=z
            assert 1+sum(N.values())==2**(m-1)
            assert old['cases'][str(k)]['incident_ordered_pairs']==2**m
            ledger.append(dict(k=k,m=m,N_mh=N,total_planes=2**(m-1),ordered_incident=2**m,n=comb(5,k-1)))
        assert [x['n']for x in ledger]==[1,5,10,10,5]
        result['source_054_intervals']={k:{f:row[f]for f in ['singular_norm_upper','background_norm_upper','total_norm_upper']}for k,row in old['cases'].items()}
        result['cases']['5']=dict(k=5,m=1,n=5,provenance='FROZEN054_COMPLETE_K5_READ_ONLY_NO_REEVALUATION',singular_norm_upper=old['cases']['5']['singular_norm_upper'],background_norm_upper=old['cases']['5']['background_norm_upper'],total_norm_upper=old['cases']['5']['total_norm_upper'])
        done(2,dict(name='mixed_sign_projective_counts_and_complete054_reuse',ledger=ledger,k5_new_evaluations=0,source054_sha256=sha(OLD/'results.json')))
        ratios={}
        for h in [1,2]:
            begin(h+2);alpha=A(Q(1,16))/A(h).sqrt();R=A(Q(84,19))*alpha-A(Q(7,2))*alpha**A(Q(6,5))
            assert 0<R<1
            Rp=R**A(Q(5,6));ratios[h]=Rp
            result['R'][str(h)]=dict(h=h,alpha=enc(alpha),R=enc(R),R_power_5_over6=enc(Rp))
            done(h+2,dict(name=f'fixed_ratio_h{h}',h=h))
        for k in range(1,5):
            begin(k+4);m=6-k;row=old['cases'][str(k)];counts=ledger[k-1]['N_mh']
            factor=(1+sum((n*ratios[int(h)]for h,n in counts.items()),arb(0)))/2**(m-1)
            assert 0<factor<1
            oldsing=load(row['singular_norm_upper']);bg=load(row['background_norm_upper']);sing=factor*oldsing;total=sing+bg
            assert sing<oldsing and total<load(row['total_norm_upper'])
            e=enc(total);hi=bounds(e)[1];simple=Q((hi*10**18).__ceil__()+1,10**18)
            result['cases'][str(k)]=dict(k=k,m=m,n=comb(5,k-1),provenance='NEW_SIGNED_INTEGRAL_FACTOR_WITH_FROZEN054_BACKGROUND',N_mh=counts,incident_norm_factor=enc(factor),singular_norm_upper=enc(sing),background_norm_upper=row['background_norm_upper'],total_norm_upper=e,simple_rational_upper=str(simple),new_over_old=enc(total/load(row['total_norm_upper'])),background_recomputed=False,all_incident_pairs_and_background_paid=True)
            done(k+4,dict(name=f'new_signed_endpoint_type_k{k}_m{m}',case_key=str(k)))
        begin(9);S=arb(0);terms=[]
        for k in range(1,6):
            row=result['cases'][str(k)];n=comb(5,k-1);term=n*load(row['total_norm_upper'])**A(Q(6,5));S+=term;terms.append(dict(k=k,n=n,weighted_power=enc(term)))
        union=S**A(Q(5,6));olde=old['union']['full_time_union_norm_upper'];oldunion=load(olde);assert union<oldunion
        e=enc(union);hi=bounds(e)[1];simple=Q((hi*10**18).__ceil__()+1,10**18)
        result['union']=dict(domain_count=31,weights=[1,5,10,10,5],q='6/5',weighted_terms=terms,q_power_sum=enc(S),new_full_time_union_norm_upper=e,simple_rational_union_upper=str(simple),source_old_union_interval=olde,new_over_old=enc(union/oldunion),conditional_Bcrit='4/48750',new_over_Bcrit=enc(union/A(Q(4,48750))),overlap_permitted=True,domains_claimed_disjoint=False,k5_new_evaluations=0,backgrounds_recomputed=False,whole_face_paid=False,new_C6_bound=False)
        done(9,dict(name='all31_same_domain_signed_sharpened_union',domain_count=31))
        result['status']='COMPLETE_SIGNED_ENDPOINT_FULL_TIME_SHARPENING';save('results.json',result);receipt['status']='COMPLETED';code=0
    except BaseException as e:
        receipt['status']='STOPPED_FAILURE_NO_RETRY';receipt['error']=repr(e);receipt['traceback']=traceback.format_exc()
    finally:
        receipt['exit_code']=code;receipt['elapsed_seconds']=time.monotonic()-tic;receipt['finished_epoch']=time.time();receipt['no_active_children']=True
        if(H/'results.json').exists():receipt['results_sha256']=sha(H/'results.json')
        save('run_receipt.json',receipt)
    print(json.dumps(dict(status=receipt['status'],configurations=receipt['completed_configurations'],elapsed_seconds=receipt['elapsed_seconds'],error=receipt.get('error'))));return code
if __name__=='__main__':
    assert sys.argv[1:]==['--run'];raise SystemExit(run())
