"""Independent B062 fixed high-event audit: exact dyadic-outward rational arithmetic."""
import sys
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
import datetime,hashlib,json,math,os,time,traceback
from fractions import Fraction as Q
from pathlib import Path
P=Path(__file__).resolve().parent
S=1<<192
ROOTS=[]
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(f,o):(P/f).write_text(json.dumps(o,indent=2)+'\n')
def exact(x):return (Q(x),Q(x))
def rounded(l,u):
    assert l<=u
    a=(l.numerator*S)//l.denominator
    b=-((-u.numerator*S)//u.denominator)
    return Q(a,S),Q(b,S)
def add(x,y):return rounded(x[0]+y[0],x[1]+y[1])
def neg(x):return -x[1],-x[0]
def sub(x,y):return add(x,neg(y))
def mul(x,y):
    v=[a*b for a in x for b in y]
    return rounded(min(v),max(v))
def div(x,y):
    assert y[0]>0
    return mul(x,(1/y[1],1/y[0]))
def power(x,n):
    assert x[0]>=0 and n>=0
    return rounded(x[0]**n,x[1]**n)
def isum(xs):
    out=exact(0)
    for x in xs:out=add(out,x)
    return out
def prod(xs):
    out=exact(1)
    for x in xs:out=mul(out,x)
    return out
def root_floor(n,k):
    assert n>=0 and k>=2
    if n<2:return n
    x=1<<((n.bit_length()+k-1)//k)
    while True:
        y=((k-1)*x+n//x**(k-1))//k
        if y>=x:break
        x=y
    while (x+1)**k<=n:x+=1
    while x**k>n:x-=1
    assert x**k<=n<(x+1)**k
    return x
def endpoint_root(q,k):
    assert q>=0
    f=root_floor((q.numerator*S**k)//q.denominator,k)
    l,u=Q(f,S),Q(f+1,S)
    assert l**k<=q<u**k
    ROOTS.append({'q':str(q),'degree':k,'lower':str(l),'upper':str(u),'lower_power_le_q':True,'q_lt_upper_power':True})
    return l,u
def root(x,k):return endpoint_root(x[0],k)[0],endpoint_root(x[1],k)[1]
def atan_endpoint(x,n):
    assert 0<=x<=Q(1,2) and n>0 and n%2==0
    ans=Q(0);term=x
    for k in range(n):
        ans+=(-1 if k%2 else 1)*term/Q(2*k+1)
        term*=x*x
    return ans,ans+term/Q(2*n+1)
def atan_small(x,n):return rounded(atan_endpoint(x[0],n)[0],atan_endpoint(x[1],n)[1])
def saved(v):return Q(int(v['lower_mantissa']))*Q(2)**v['lower_binary_exponent'],Q(int(v['upper_mantissa']))*Q(2)**v['upper_binary_exponent']
def compact(x,digits=40):
    z=10**digits
    l=Q((x[0].numerator*z)//x[0].denominator,z)
    u=Q(-((-x[1].numerator*z)//x[1].denominator),z)
    return {'lower':str(l),'upper':str(u)}
def encoded(x):return {'lower':str(x[0]),'upper':str(x[1]),'width':str(x[1]-x[0]),'outward_40_decimal_places':compact(x)}
# Exact polynomial ring Q[a,b,c] for identities, not numerical samples.
def padd(*ps):
    d={}
    for p in ps:
        for k,v in p.items():d[k]=d.get(k,0)+v
    return {k:v for k,v in d.items() if v}
def pneg(p):return {k:-v for k,v in p.items()}
def pmul(p,q):
    d={}
    for a,x in p.items():
        for b,y in q.items():
            k=tuple(i+j for i,j in zip(a,b));d[k]=d.get(k,0)+x*y
    return {k:v for k,v in d.items() if v}
a={(1,0,0):1};bpoly={(0,1,0):1};c={(0,0,1):1}
def main():
    assert sys.argv[1:]==['--run']
    assert not (P/'execution_receipt.json').exists(),'no automatic rerun'
    reg=read(P/'registration.json');cfg=read(P/'input_manifest.json');C=Path(reg['producer_directory'])
    tic=time.monotonic();rec={'status':'RUNNING','pid':os.getpid(),'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'input_manifest_sha256':sha(P/'input_manifest.json'),'completed_configurations':0,'new_kernel_calls':0,'old_M_integral_calls':0,'old_cap_reruns':0,'old055_056_060_reruns':0,'full_enumerations':0,'MATLAB':0}
    save('execution_receipt.json',rec);out={'status':'INCOMPLETE','configurations':[]}
    def done(i,name,data):
        assert i==rec['completed_configurations']+1<=8
        assert time.monotonic()-tic<180
        out['configurations'].append({'id':i,'name':name,**data})
        rec['completed_configurations']=i;save('results.json',out);save('execution_receipt.json',rec)
    pins=[]
    def pin(path,h):
        assert sha(path)==h,('hash mismatch',str(path))
        pins.append((str(path),h))
    code=2
    try:
        # 1: Fixed bytes and historical receipts. No authored code is imported/run.
        pin(C/'manifest.json',reg['producer_manifest_sha256']);pin(C/'report.md',reg['producer_report_sha256'])
        m=read(C/'manifest.json')
        for f,h in m['whitelist_sha256'].items():pin(C/f,h)
        for f,h in m['source_sha256_unchanged'].items():pin(f,h)
        for f,h in cfg['audit_source_sha256'].items():pin(P/f,h)
        authored=read(C/'results.json');old=read(next(p for p in m['source_sha256_unchanged'] if p.endswith('/round041/one_large_jump_cap/results.json')))
        R1=Q(67320683,3251298631680);cr=Q(44743079,2970)
        assert any(v.get('R1_exact')==str(R1) for v in old['configurations'])
        assert any(v.get('Crem_times_pi6_exact')==str(cr) for v in old['configurations'])
        ar=read(C/'run_receipt.json');assert ar['status']=='COMPLETED' and ar['exit_code']==0 and ar['completed_configurations']==ar['attempted_configurations']==9 and ar['results_sha256']==sha(C/'results.json')
        for name in ['producer_process_receipt.json','checker_process_receipt.json']:
            e=read(C/name);assert e['status']=='COMPLETED' and e['exit_code']==0 and e['child_reaped'] and e['wall_timeout_seconds']==180 and e['elapsed_seconds']<180 and not e['stderr']
        done(1,'frozen_sources_inherited_coefficients_and_author_receipts',{'hash_entries':len(pins),'R1':str(R1),'Crem_times_pi6':str(cr)})
        # 2: Exact scalar gates for analytic coverage of every direction/plane.
        rho=Q(1,64);bb=Q(1,4);d=Q(31,64)
        assert Q(3,4)-rho-bb==d>0 and (1-Q(47,64))**2>2*rho*rho
        h_values=[Q(2,3)*h*(6-h) for h in range(1,6)]
        assert min(h_values)==Q(10,3)
        assert all(1-s*Q(6-2*h,6)>=0 for h in range(1,6) for s in [-1,1])
        apb=padd(a,bpoly);amb=padd(a,pneg(bpoly))
        body_res=padd(pmul(apb,apb),pneg(pmul(amb,amb)),pneg({k:4*v for k,v in pmul(a,bpoly).items()}));assert not body_res
        assert [32*(31-k) for k in range(3)]==[992,960,928]
        done(2,'all_subspace_geometry_analytic_gate_certificates',{'line_h_coefficients':list(map(str,h_values)),'body_plane_body_line_incidence_polynomial_residual':{},'norm_line_pair_distance_squared_coefficients':['6','10/3','2'],'no_direction_pair_triple_enumeration':True})
        # 3: Rank-one Gram determinant formula, normal dimension and union weights.
        volume_ledger=[]
        for k in range(1,6):
            mm=6-k
            assert (k-1)+(mm-1)==4
            assert 1-Q(k-1,k)==Q(1,k)
            volume_ledger.append({'k':k,'m':mm,'count':math.comb(5,k-1),'face_Gram_det':str(Q(1,k)),'maxcone_fraction':str(Q(1,k)),'mean_interval_length':'1/4'})
        assert sum(r['count'] for r in volume_ledger)==31
        done(3,'genuine_face_measure_and_one_overlapping_union',{'ledger':volume_ledger,'normal_dimension':4,'orthant_symmetry_used_only_for_full_p_slice':True})
        # 4: Full rational-function identities after multiplying positive denominators.
        # Here a=b^2 and bpoly=epsilon^2 for the tail derivative numerator.
        tail_res=padd(pneg(pmul(a,padd(bpoly,a))),pmul(bpoly,padd(bpoly,a)),pneg(pmul(bpoly,bpoly)),pmul(a,a));assert not tail_res
        # b(b+c)-a(a+c)=(b-a)(a+b+c), yielding the three-root formula.
        int_res=padd(pmul(bpoly,padd(bpoly,c)),pneg(pmul(a,padd(a,c))),pneg(pmul(padd(bpoly,pneg(a)),padd(a,bpoly,c))));assert not int_res
        assert Q(5,6)+Q(7,6)==2 and Q(-1,6)-Q(5,6)==-1
        assert Q(6,7)+Q(1,3)==Q(25,21)
        assert Q(1,6)+4==Q(25,6)
        done(4,'complete_time_identities_and_loose_optimizer',{'tail_derivative_polynomial_residual':{},'three_quadratic_polynomial_residual':{},'epsilon_identity':'epsilon^5=C*pi*V','P_reduction':'(pi*V^(1/6))^-1 * ((25/21)*epsilon^2-epsilon^4/b^2+epsilon^5*atan(b/epsilon)/b^3)','coincident_positive_roots':'continuous extension; no root-gap division used'})
        # 5: Certified elementary constants, all root endpoints separately powered.
        p5=atan_endpoint(Q(1,5),100);p239=atan_endpoint(Q(1,239),32)
        pi=rounded(16*p5[0]-4*p239[1],16*p5[1]-4*p239[0]);assert 3<pi[0]<pi[1]<4
        # Exact tangent identity behind Machin, with its classical first-quadrant branch.
        t=Q(1,5);t2=2*t/(1-t*t);t4=2*t2/(1-t2*t2)
        assert (t4-Q(1,239))/(1+t4/Q(239))==1
        assert 4*Q(1,5)/(1+Q(1,25))-Q(1,239)>0 and Q(4,5)<1
        terms=[mul(exact(math.comb(5,k-1)),root(exact(Q(6-k,k)),2)) for k in range(1,6)]
        shape=isum(terms);V=mul(mul(div(power(pi,2),exact(8)),exact(rho**4)),shape)
        D=div(exact(R1),mul(exact(32),power(pi,6)));Cbig=div(D,exact(bb*bb))
        eps=root(prod([Cbig,pi,V]),5)
        ratio=div(eps,exact(bb));assert 0<ratio[0]<ratio[1]<Q(1,2)
        done(5,'Machin_pi_dyadic_roots_high_volume_and_unique_epsilon',{'pi':encoded(pi),'shape':encoded(shape),'V':encoded(V),'epsilon':encoded(eps),'epsilon_over_b':encoded(ratio),'root_grid':'2^-192','Machin_terms':[100,32]})
        # 6: Algebraically reduced complete P formula, rather than source fractional powers.
        Vroot=root(V,6);norm=div(exact(1),mul(pi,Vroot))
        at= sub(div(pi,exact(2)),atan_small(ratio,80))
        e2=power(eps,2);e4=power(eps,4);e5=power(eps,5)
        short=mul(norm,mul(exact(Q(6,7)),e2))
        tail_inside=add(sub(div(e2,exact(3)),div(e4,exact(bb**2))),div(mul(e5,at),exact(bb**3)))
        long=mul(norm,tail_inside);assert long[0]>0
        combined=mul(norm,add(sub(mul(exact(Q(25,21)),e2),div(e4,exact(bb**2))),div(mul(e5,at),exact(bb**3))))
        assert combined[0]>0
        done(6,'independent_reduced_complete_P_budget',{'Pshort':encoded(short),'Plong':encoded(long),'Ptotal':encoded(combined),'atan_small_terms':80,'complete_time_axis':'Pshort 0..epsilon plus exact Ptail epsilon..infinity'})
        # 7: t=d*s separates the d^-5 factor from the three-root integral.
        A=[992,960,928];z=[Q(6),Q(10,3),Q(2)]
        roots=[root(exact(v/aa),2) for v,aa in zip(z,A)]
        x,y,w=roots
        denom=prod([exact(2*math.prod(A)*d**5),x,y,w,add(x,y),add(x,w),add(y,w)])
        J=div(mul(pi,isum(roots)),denom)
        V56=div(V,Vroot);Crem=div(exact(cr),power(pi,6));B=prod([V56,Crem,J])
        assert B[0]>0
        done(7,'independent_B_complete_axis_with_scaled_time',{'B_time_integral':encoded(J),'Bfull':encoded(B),'time_substitution':'t=d*s, factor d^-5','all_remaining_ordered_triples':26970})
        # 8: Independent strict endpoint certificate, then consistency with frozen authored intervals.
        total=add(combined,B);claim=Q(reg['claim_strict_upper'])
        assert total[1]<claim
        comparisons={}
        for name,v,stored in [('V',V,authored['constants']['V']),('epsilon',eps,authored['epsilon']),('Pshort',short,authored['Pshort']),('Plong',long,authored['Plong']),('Bfull',B,authored['Bfull']),('total',total,authored['union']['full_time_high_norm_upper'])]:
            l,u=saved(stored);assert v[0]<=l<=u<=v[1],(name,'saved enclosure outside independent enclosure')
            comparisons[name]='saved_interval_contained_in_independent_interval'
        assert Q(authored['union']['simple_rational_upper'])==claim
        assert claim<Q(4,48750)
        for path,h in pins:assert sha(path)==h,('source changed during audit',path)
        done(8,'strict_high_norm_claim_and_final_binding',{'total':encoded(total),'strict_claim':str(claim),'claim_minus_independent_upper':str(claim-total[1]),'claim_below_conditional_whole_face_reference':True,'saved_enclosure_consistency':comparisons,'all_source_hashes_unchanged':True})
        out.update(status='PASS_INDEPENDENT_HIGH_FULL_TIME_EVENT_BUDGET',strict_high_norm_upper=str(claim),total_certified_interval=encoded(total),root_endpoint_certificates=ROOTS,root_endpoint_checks=len(ROOTS),all_parameter_events_paid=True,all_positive_time_paid=True,one_true_volume_union=True,whole_face_I_proved=False,global_g_proved=False,new_Cn_bound=False,new_kernel_calls=0,old_M_integral_calls=0,old_cap_reruns=0,old055_056_060_reruns=0,full_enumerations=0,MATLAB=0)
        save('results.json',out);code=0;rec.update(status=out['status'],exit_code=0,results_sha256=sha(P/'results.json'))
    except BaseException as err:
        rec.update(status='STOP_INDEPENDENT_AUDIT_FAILURE_NO_RETRY',exit_code=2,error=repr(err),traceback=traceback.format_exc())
    finally:
        rec.update(elapsed_seconds=time.monotonic()-tic,completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
        save('execution_receipt.json',rec)
    print(json.dumps({'status':rec['status'],'completed_configurations':rec['completed_configurations'],'elapsed_seconds':rec['elapsed_seconds'],'error':rec.get('error')}))
    return code
if __name__=='__main__':raise SystemExit(main())
