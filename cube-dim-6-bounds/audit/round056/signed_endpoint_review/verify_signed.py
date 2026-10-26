#!/usr/bin/env python3
"""One pre-registered independent Fraction/integer-root batch, ten configurations.
No author imports, Arb, old integral, old complete budget, or k5 evaluation.
"""
import sys
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
from pathlib import Path
from fractions import Fraction as Q
from math import comb
import hashlib,json,time,datetime,os,traceback
H=Path(__file__).resolve().parent
SRC=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-curvature-astra/cube-dim-6-bounds/upper/round056/signed_endpoint_integration')
OLD=SRC.parents[1]/'round054/lower_endpoint_green_budget'
SCALE=10**70
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text())
roots=[]
def pair(v):
    if 'common_binary_exponent' in v:
        z=Q(2)**int(v['common_binary_exponent'])
        return Q(int(v['lower_mantissa']))*z,Q(int(v['upper_mantissa']))*z
    return (Q(int(v['lower_mantissa']))*Q(2)**int(v['lower_binary_exponent']),
            Q(int(v['upper_mantissa']))*Q(2)**int(v['upper_binary_exponent']))
def root(x,n,label):
    x=Q(x);assert x>0 and n>=2
    integer=x.numerator*SCALE**n//x.denominator
    lo=0;hi=1<<((integer.bit_length()+n-1)//n)
    while lo+1<hi:
        mid=(lo+hi)//2
        if mid**n<=integer:lo=mid
        else:hi=mid
    l=Q(lo,SCALE);u=Q(lo+1,SCALE)
    assert l**n<=x<u**n
    roots.append({'label':label,'n':n,'input':str(x),'lower':str(l),'upper':str(u),'exact_power_bracketing':True})
    return l,u
def enc(v):return [str(x) for x in v]
def contains(v,author):
    a,b=pair(author);assert v[0]<=a<=b<=v[1]
def run():
    reg=read(H/'registration.json')
    def guard():
        assert time.time()-reg['started_unix']<reg['maximum_scientific_wall_seconds']
        for p,h in reg['source_sha256'].items():assert sha(p)==h,p
    guard();old=read(OLD/'results.json');author=read(SRC/'results.json')
    out={'status':'RUNNING','configurations':[],'source_pin_count':len(reg['source_sha256']),'ratios':{},'cases':{}}
    q=Q(6,5);rho=Q(1,64);L=Q(1,4)
    # 1: exact algebraic coefficients and combinatorial identities only.
    A=q*(4-q)/((q-1)*(5-q));B=(4-q)/(4*(q-1))
    assert (A,B)==(Q(84,19),Q(7,2))
    assert A-B==Q(35,38) and A-B*q==Q(21,95)>0
    multiplicity={}
    for m in range(1,6):
        row={h:comb(m,h)//(2 if 2*h==m else 1) for h in range(1,m//2+1)}
        assert 1+sum(row.values())==2**(m-1)
        assert 2*(1+sum(row.values()))==2**m
        multiplicity[m]=row
    assert multiplicity=={1:{},2:{1:1},3:{1:3},4:{1:4,2:3},5:{1:5,2:10}}
    weights=[comb(5,k-1) for k in range(1,6)];assert weights==[1,5,10,10,5] and sum(weights)==31
    out['configurations'].append({'id':1,'name':'exact integration ratio coefficients and complete mixed-sign ledgers','ratio_coefficients':[str(A),str(B)],'R_at_one':str(A-B),'derivative_lower_on_unit_interval':str(A-B*q),'multiplicity':multiplicity,'weights':weights,'quadrature_calls':0})
    # 2--3: signed integral ratios, all endpoints rational and outward.
    R={};P={}
    for h in [1,2]:
        guard()
        c=(Q(1),Q(1)) if h==1 else root(Q(h),2,'h2:sqrt_h')
        alpha=(rho/L/c[1],rho/L/c[0]);assert 0<alpha[0]<=alpha[1]<1
        tL=root(alpha[0]**6,5,f'h{h}:alpha_q_lower')[0]
        tU=root(alpha[1]**6,5,f'h{h}:alpha_q_upper')[1]
        R[h]=(A*alpha[0]-B*tU,A*alpha[1]-B*tL)
        assert 0<R[h][0]<=R[h][1]<Q(35,38)<1
        P[h]=(root(R[h][0]**5,6,f'h{h}:R_5over6_lower')[0],root(R[h][1]**5,6,f'h{h}:R_5over6_upper')[1])
        out['ratios'][str(h)]={'sqrt_h':enc(c),'alpha':enc(alpha),'R':enc(R[h]),'R_power_5_over6':enc(P[h])}
        out['configurations'].append({'id':h+1,'name':f'independent h{h} exact radial-mean ratio'})
    # 4--7: improve only k1--4, carrying full054 singular/background intervals.
    U={};F={};singular={};source={}
    for k in range(1,5):
        guard();m=6-k;counts=multiplicity[m]
        F[k]=((1+sum((n*P[h][0] for h,n in counts.items()),Q(0)))/2**(m-1),
              (1+sum((n*P[h][1] for h,n in counts.items()),Q(0)))/2**(m-1))
        assert 0<F[k][0]<=F[k][1]<1
        oc=old['cases'][str(k)];s=pair(oc['singular_norm_upper']);b=pair(oc['background_norm_upper'])
        assert 0<s[0]<=s[1] and 0<b[0]<=b[1]
        source[str(k)]={v:oc[v] for v in ['singular_norm_upper','background_norm_upper','total_norm_upper']}
        singular[k]=(F[k][0]*s[0],F[k][1]*s[1])
        U[k]=(singular[k][0]+b[0],singular[k][1]+b[1])
        assert U[k][1]<pair(oc['total_norm_upper'])[0]
        out['cases'][str(k)]={'k':k,'m':m,'multiplicity':counts,'factor':enc(F[k]),'new_singular':enc(singular[k]),'background_complete_interval':oc['background_norm_upper'],'new_total':enc(U[k]),'strictly_below_old_formula_lower':True,'display_only_total_upper':float(U[k][1])}
        out['configurations'].append({'id':k+3,'name':f'fixed k{k} mixed-sign norm with frozen054 inputs'})
    # 8: identity-only preservation, including saved k5 q-power. No k5 powers evaluated.
    for k in range(1,6):
        oc=old['cases'][str(k)]
        source[str(k)]={v:oc[v] for v in ['singular_norm_upper','background_norm_upper','total_norm_upper']}
        assert author['source_054_intervals'][str(k)]==source[str(k)]
        assert author['cases'][str(k)]['background_norm_upper']==oc['background_norm_upper']
    for key in source['5']:assert author['cases']['5'][key]==source['5'][key]
    old5=old['union']['weighted_terms'][4]
    assert old5['k']==5 and old5['n']==5 and author['union']['weighted_terms'][4]==old5
    old_union=old['union']['full_time_union_norm_upper'];old_bounds=pair(old_union)
    assert author['union']['source_old_union_interval']==old_union
    out['source_054_complete_intervals']=source
    out['k5_complete_saved_case']=source['5'];out['k5_saved_weighted_q_power']=old5
    out['configurations'].append({'id':8,'name':'complete background/k5/old-union identity checks','k5_new_evaluations':0,'background_recomputations':0,'old_union_complete_interval':old_union})
    # 9: same31 union. The saved weighted k5 term enters directly.
    terms=[];lo=Q(0);hi=Q(0)
    for k in range(1,5):
        term=(weights[k-1]*root(U[k][0]**6,5,f'union:k{k}:q_lower')[0],
              weights[k-1]*root(U[k][1]**6,5,f'union:k{k}:q_upper')[1])
        lo+=term[0];hi+=term[1];terms.append({'k':k,'weight':weights[k-1],'interval':enc(term)})
    copied=pair(old5['weighted_power']);lo+=copied[0];hi+=copied[1]
    terms.append({'k':5,'weight':5,'interval':enc(copied),'copied_without_reevaluation':True})
    final=(root(lo**5,6,'union:five_sixths_lower')[0],root(hi**5,6,'union:five_sixths_upper')[1])
    claim=Q(author['union']['simple_rational_union_upper'])
    assert final[1]<claim<old_bounds[0]
    ratio=(final[0]/old_bounds[1],final[1]/old_bounds[0])
    out['union']={'weights':weights,'domain_count':31,'weighted_terms':terms,'q_power_sum':enc((lo,hi)),'independent_interval':enc(final),'claimed_strict_rational_upper':str(claim),'strict_margin_to_claim':str(claim-final[1]),'old_formula_complete_interval':old_union,'strict_margin_to_old_formula_lower':str(old_bounds[0]-final[1]),'new_over_old_formula':enc(ratio),'display_only_upper':float(final[1]),'display_only_ratio_upper':float(ratio[1]),'overlaps_allowed':True,'all_positive_time':True,'whole_face_paid':False}
    out['configurations'].append({'id':9,'name':'same31 positive q-union and strict old054-formula comparison'})
    # 10: compare independent outward intervals to all corresponding frozen056 intervals.
    for h in [1,2]:
        ar=author['R'][str(h)];row=out['ratios'][str(h)]
        for key in ['alpha','R','R_power_5_over6']:contains(tuple(map(Q,row[key])),ar[key])
    for k in range(1,5):
        ar=author['cases'][str(k)]
        assert ar['N_mh']=={str(h):n for h,n in multiplicity[6-k].items()}
        contains(F[k],ar['incident_norm_factor']);contains(singular[k],ar['singular_norm_upper']);contains(U[k],ar['total_norm_upper'])
        assert U[k][1]<Q(ar['simple_rational_upper'])
    contains(final,author['union']['new_full_time_union_norm_upper']);contains(ratio,author['union']['new_over_old'])
    assert author['union']['weights']==weights and author['union']['domain_count']==31
    assert author['union']['k5_new_evaluations']==0 and not author['union']['backgrounds_recomputed']
    out['configurations'].append({'id':10,'name':'independent interval containment and frozen certificate comparison','all_corresponding_author_intervals_contained':True})
    guard();assert len(out['configurations'])==10
    out.update(status='PASS_TEN_INDEPENDENT_EXACT_CONFIGURATIONS',integer_root_enclosures=roots,old_M_integral_recomputations=0,old_complete054_budget_recomputations=0,k5_new_evaluations=0,background_recomputations=0,native_or_kernel_calls=0,grid_calls=0,Arb_calls=0,Matlab_calls=0,author_code_runs=0,new_C6_bound=False)
    return out
if __name__=='__main__':
    assert not(H/'results.json').exists() and not(H/'failure_receipt.json').exists(),'one batch; no retry'
    tic=time.monotonic()
    try:
        out=run();(H/'results.json').write_text(json.dumps(out,indent=2)+'\n')
        receipt={'status':'PASS','pid':os.getpid(),'exit_code':0,'seconds':time.monotonic()-tic,'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_sha256':sha(__file__),'registration_sha256':sha(H/'registration.json'),'theory_before_results_sha256':sha(H/'theory_before_results.md'),'results_sha256':sha(H/'results.json'),'python':sys.version,'python_executable':sys.executable,'python_executable_sha256':sha(sys.executable),'configurations':len(out['configurations']),'integer_root_enclosures':len(roots),'retry_count':0,'child_processes':[]}
        (H/'run_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
        print(json.dumps({'status':out['status'],'pid':os.getpid(),'seconds':receipt['seconds'],'configurations':10,'roots':len(roots),'union_upper_display_only':out['union']['display_only_upper'],'new_over_old_formula_upper_display_only':out['union']['display_only_ratio_upper']}))
    except BaseException as e:
        (H/'failure_receipt.json').write_text(json.dumps({'status':'STOP_FAILURE','pid':os.getpid(),'seconds':time.monotonic()-tic,'exception':repr(e),'traceback':traceback.format_exc(),'source_sha256':sha(__file__)},indent=2)+'\n')
        raise
