#!/usr/bin/env python3
"""Eight independent fixed configurations; Fraction/integers only.
No author imports, no Arb, no old M integral, no direction enumeration.
"""
import sys
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
from pathlib import Path
from fractions import Fraction as Q
from math import comb
import hashlib,json,time,datetime,os,traceback
H=Path(__file__).resolve().parent
SRC=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-curvature-astra/cube-dim-6-bounds/upper/round054/lower_endpoint_green_budget')
OLD=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-audit/cube-dim-6-bounds/upper/round037/pair_resolved_green_envelope')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
SCALE=10**60
root_ledger=[]
def read(p):return json.loads(Path(p).read_text())
def bounds(v):
    if 'common_binary_exponent' in v:
        e=int(v['common_binary_exponent']);return Q(int(v['lower_mantissa']))*Q(2)**e,Q(int(v['upper_mantissa']))*Q(2)**e
    return Q(int(v['lower_mantissa']))*Q(2)**int(v['lower_binary_exponent']),Q(int(v['upper_mantissa']))*Q(2)**int(v['upper_binary_exponent'])
def nth_root_box(x,n,label):
    x=Q(x);assert x>=0 and n>=2
    integer=(x.numerator*SCALE**n)//x.denominator
    low=0;high=1<<((integer.bit_length()+n-1)//n)
    while low+1<high:
        mid=(low+high)//2
        if mid**n<=integer:low=mid
        else:high=mid
    l=Q(low,SCALE);u=Q(low+1,SCALE)
    assert l**n<=x<u**n
    root_ledger.append({'label':label,'n':n,'lower':str(l),'upper':str(u),'exact_power_bracketing':True,
                        'input_fraction_sha256':hashlib.sha256(str(x).encode()).hexdigest()})
    return l,u

def determinant(matrix):
    a=[list(map(Q,row))for row in matrix]; n=len(a);out=Q(1)
    for k in range(n):
        pivot=next(i for i in range(k,n)if a[i][k])
        if pivot!=k:a[k],a[pivot]=a[pivot],a[k];out=-out
        q=a[k][k];out*=q
        for i in range(k+1,n):
            v=a[i][k]/q
            for j in range(k+1,n):a[i][j]-=v*a[k][j]
    return out

def atan_box(reciprocal,even_count):
    assert reciprocal>1 and even_count%2==0
    partial=sum((Q((-1)**n,(2*n+1)*reciprocal**(2*n+1))for n in range(even_count)),Q(0))
    return partial,partial+Q(1,(2*even_count+1)*reciprocal**(2*even_count+1))

def run():
    reg=read(H/'registration.json');plan=read(H/'execution_plan.json')
    assert sha(H/'registration.json')==plan['registered_input_sha256']
    pins=dict(reg['source_sha256']);pins.update(plan['additional_dependency_pins'])
    def guard():
        assert time.time()-reg['started_unix']<900
        for p,h in pins.items():assert sha(p)==h,p
    guard();out={'status':'RUNNING','configurations':[],'cases':{},'source_pin_count':len(pins)}
    rho=Q(1,64);N=[192,480,320];weights=[comb(5,k-1)for k in range(1,6)]
    author=read(SRC/'results.json');old=read(OLD/'results.json')
    # 1: exact structural data; no enumeration of physical pairs.
    structures=[]
    for k in range(1,6):
        m=6-k;j=min(k,m);c=2**m;bg=[v-(c if i==j else 0)for i,v in enumerate(N,1)]
        assert all(v>=0 for v in bg)and sum(bg)+c==992
        mat=[[Q(int(a==b))-Q(1,k)for b in range(k-1)]for a in range(k-1)]
        assert determinant(mat)==Q(1,k)
        assert Q(3*m,32)>rho**2 and Q(9,32)>k*rho**2
        ar=author['cases'][str(k)]
        assert ar['incident_planes']==2**(m-1) and ar['incident_ordered_pairs']==c
        assert ar['background_by_class']==bg and ar['n']==weights[k-1] and ar['pair_class']==j
        structures.append({'k':k,'m':m,'class':j,'incident_planes':2**(m-1),'ordered_incident':c,
                           'background':bg,'normal_Gram_det':str(Q(1,k)),'normal_dimension':4,
                           'maximal_cone_fraction':str(Q(1,k)),'weight':weights[k-1],
                           'line_positivity_margin':str(Q(3*m,32)-rho**2),'plane_positivity_margin':str(Q(9,32)-k*rho**2)})
    assert weights==[1,5,10,10,5] and sum(weights)==31 and sum(N)==992
    out['configurations'].append({'id':1,'name':'exact five-type structure and Gram checks','rows':structures})
    # 2: independent pi enclosure, exact alternating-series errors.
    l5,u5=atan_box(5,100);l239,u239=atan_box(239,30)
    piL=16*l5-4*u239;piU=16*u5-4*l239
    assert Q(3)<piL<piU<Q(22,7) and piU-piL<Q(1,10**135)
    out['configurations'].append({'id':2,'name':'fixed rational Machin pi bracket','identity':'pi=16 atan(1/5)-4 atan(1/239)',
                                  'atan_1_5_terms':100,'atan_1_239_terms':30,'pi_lower':str(piL),'pi_upper':str(piU),
                                  'width_less_than':'1/10^135'})
    raw={str(c['k']):c['M']for c in old['configurations'][:3]}
    assert set(raw)=={'1','2','3'} and raw==author['source_M_intervals']
    endpoints={j:bounds(v)for j,v in raw.items()};assert all(0<l<=u for l,u in endpoints.values())
    M={int(j):u for j,(l,u)in endpoints.items()};out['source_M_complete_intervals']=raw
    U={}
    # 3--7: independent upper-endpoint arithmetic with exact integer roots.
    for k in range(1,6):
        guard();m=6-k;j=min(k,m);c=2**m
        bg=[v-(c if i==j else 0)for i,v in enumerate(N,1)]
        etaU=nth_root_box(k*rho**2,2,f'k{k}:eta')[1]
        dL=nth_root_box(Q(3*k*m,32),2,f'k{k}:line_baseline')[0]-etaU
        gL=nth_root_box(Q(9,32),2,f'k{k}:plane_baseline')[0]-etaU
        assert dL>0 and gL>0
        r2=k*(1-rho)**2
        fU=nth_root_box(Q(m,16*k),2,f'k{k}:face_factor')[1]
        I=c*M[j]/(r2*dL**2)
        B=sum((v*M[i]for i,v in enumerate(bg,1)),Q(0))/(r2*dL**2*gL)
        V=fU*piU**2/2*rho**4
        # J^5 identity removes fractional radius powers before root bounding.
        J5=(fU*Q(5,7)*piU**2)**5*rho**14
        singular=nth_root_box(I**6*J5,6,f'k{k}:singular_norm')[1]
        regular=nth_root_box(B**6*V**5,6,f'k{k}:background_norm')[1]
        U[k]=singular+regular
        simple=Q(author['cases'][str(k)]['simple_rational_upper'])
        assert U[k]<simple
        row={'k':k,'m':m,'weight':weights[k-1],'exact_r2':str(r2),'d_lower':str(dL),'g_lower':str(gL),
             'face_factor_upper':str(fU),'singular_norm_upper':str(singular),'background_norm_upper':str(regular),
             'total_norm_upper':str(U[k]),'author_simple_rational_upper':str(simple),'independent_upper_strictly_below_author_simple':True,
             'display_only_total':float(U[k]),'source_M_integrals_recomputed':False}
        out['cases'][str(k)]=row;out['configurations'].append({'id':k+2,'name':f'fixed type k{k}', 'case_key':str(k)})
    # 8: overlap-safe q-union, all operations positive with integer powers.
    terms=[];S=Q(0)
    for k in range(1,6):
        v=nth_root_box(U[k]**6,5,f'union:k{k}:six_fifths')[1]*weights[k-1]
        terms.append({'k':k,'weight':weights[k-1],'weighted_q_power_upper':str(v)});S+=v
    final=nth_root_box(S**5,6,'union:five_sixths')[1]
    target=Q(10293038485441,5*10**17);conditional=Q(4,48750)
    assert target==Q(author['union']['simple_rational_union_upper'])
    assert final<target<Q(251,1000)*conditional
    out['union']={'terms':terms,'independent_rational_upper':str(final),'claimed_strict_rational_upper':str(target),
                  'strict_margin':str(target-final),'display_only_upper':float(final),'conditional_target':str(conditional),
                  'ratio_upper_display_only':float(final/conditional),'relative_upper':'251/1000','weights':weights,
                  'domain_count':31,'overlaps_allowed':True,'whole_face_paid':False,'all_time':True}
    out['configurations'].append({'id':8,'name':'positive weighted q-union and fixed exact threshold'})
    guard();assert len(out['configurations'])==8 and len(root_ledger)==36
    out.update(status='PASS_EIGHT_INDEPENDENT_EXACT_CONFIGURATIONS',integer_root_enclosures=root_ledger,
               old_M_integral_recomputations=0,full_pair_enumerations=0,native_or_kernel_calls=0,Arb_calls=0,
               author_writer_runs=0,Matlab_calls=0,all_time_lower_endpoint_union_only=True,new_C6_bound=False)
    return out

if __name__=='__main__':
    assert not(H/'results.json').exists()and not(H/'failure_receipt.json').exists(),'one batch, no automatic retry'
    tic=time.monotonic()
    try:
        out=run();(H/'results.json').write_text(json.dumps(out,indent=2)+'\n')
        receipt={'status':'PASS','pid':os.getpid(),'exit_code':0,'seconds':time.monotonic()-tic,
                 'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
                 'source_sha256':sha(__file__),'registration_sha256':sha(H/'registration.json'),
                 'execution_plan_sha256':sha(H/'execution_plan.json'),'results_sha256':sha(H/'results.json'),
                 'python':sys.version,'python_executable':sys.executable,'python_executable_sha256':sha(sys.executable),
                 'configurations':8,'integer_root_enclosures':len(root_ledger),'all_new_kernels':0,'old_M_recomputations':0,
                 'author_writer_runs':0,'retry_count':0,'child_processes':[]}
        (H/'run_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
        print(json.dumps({'status':out['status'],'pid':receipt['pid'],'seconds':receipt['seconds'],'configurations':8,
                          'union_upper_display_only':out['union']['display_only_upper'],'rational_threshold':out['union']['claimed_strict_rational_upper']}))
    except BaseException as exc:
        (H/'failure_receipt.json').write_text(json.dumps({'status':'STOP_FAILURE','pid':os.getpid(),'seconds':time.monotonic()-tic,
                  'exception':repr(exc),'traceback':traceback.format_exc(),'source_sha256':sha(__file__)},indent=2)+'\n')
        raise
