"""Ten fixed A069 checks; 384 tiny axis6 tuple-field terms, no body32 kernel."""
import signal
import time
START=time.monotonic()
def timeout(signum,frame):raise TimeoutError('180-second A069 program wall limit')
signal.signal(signal.SIGALRM,timeout)
signal.setitimer(signal.ITIMER_REAL,180)

from pathlib import Path
from fractions import Fraction as F
from itertools import product
from collections import Counter
from math import comb
import hashlib,json,os,sys
import sympy as s
from derivative_contract import (qualify_off_S3,local_tail_contract,bind_derivative_request,
                                  formal_dd_channels,formal_cone_field,CHANNELS)

P=Path(__file__).resolve().parent
if (P/'execution_receipt.json').exists() or (P/'results.json').exists():
    raise RuntimeError('Fixed single-attempt batch cannot be replayed')
receipt={'status':'RUNNING_FIXED_TEN_CONFIGURATIONS','pid':os.getpid(),'program_wall_timeout_seconds':180,
         'completed_configurations':0,'tiny_axis_tuple_field_terms_completed':0,'body32_full_kernel_calls':0,
         'old_baseline_reruns':0,'automatic_retries':0}
checks=[]
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text())
def checkpoint():
    receipt['wall_seconds']=time.monotonic()-START
    (P/'execution_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    (P/'checkpoint.json').write_text(json.dumps({'completed_checks':checks,'receipt':receipt},indent=2)+'\n')
def done(row):
    checks.append(row);receipt['completed_configurations']=len(checks);checkpoint()

try:
    checkpoint();pins=read(P/'input_manifest.json')
    for path,digest in pins['sha256'].items():assert sha(path)==digest,path
    axis=[[int(i==j) for i in range(6)] for j in range(6)]
    point=['1','2','3','4','0','0']
    qualification=qualify_off_S3(axis,point)
    assert qualification['status']=='OFF_TRUE_S3' and qualification['min_squared_distance']=='1'
    assert qualification['triple_count']==20
    tail=local_tail_contract(qualification,'1/8','1/100')
    assert tail['uniform_distance_delta']=='3/8'
    fields=[(),(0,),(4,),(0,0),(4,4),(4,5)]
    request=bind_derivative_request(qualification,{'frame_sha256':qualification['frame_sha256'],
         'complete_fan_proved':True,'multiplicities_included_once':True,
         'proof_reference':'axis6 complete fan: all64 orthants, disjoint interiors, q=1 and determinant magnitude1'},fields)
    done({'configuration':1,'pass':True,'qualification':qualification,'local_tail_contract':tail,'derivative_request_contract':request})

    redundant=axis+[axis[0][:],[1,1,0,0,0,0]]
    guard=[]
    for p,expected in [(point,'OFF_TRUE_S3'),(['1','1','3','4','0','0'],'ON_TRUE_S3_UNSUPPORTED')]:
        q=qualify_off_S3(redundant,p)
        assert q['status']==expected and q['triple_count']==56
        assert any(row['rank']<3 for row in q['triples'])
        rejected=False
        if expected!='OFF_TRUE_S3':
            try:local_tail_contract(q,'0','1/100')
            except ValueError:rejected=True
            assert rejected
        guard.append({'point':p,'status':q['status'],'triple_count':q['triple_count'],
                      'rank_distribution':{str(r):sum(row['rank']==r for row in q['triples']) for r in (1,2,3)},
                      'minimum_squared_distance':q['min_squared_distance'],'strict_derivative_gate_rejected':rejected})
    done({'configuration':2,'pass':True,'frame_columns':redundant,'guard_cases':guard,'total_exact_triple_residuals':112})

    epsilon=s.symbols('epsilon',positive=True);z=s.symbols('z',real=True)
    for cfg,nodes in [(3,[F(0)]*6+[F(1,3)]),(4,[F(0)]*8)]:
        channels,weights=formal_dd_channels(nodes);count=Counter(nodes);n=len(nodes)-1
        conditions=[(b,k) for b,m in sorted(count.items()) for k in range(m)]
        M=s.Matrix([[s.binomial(j,k)*s.Rational(b.numerator,b.denominator)**(j-k) if j>=k else 0
                     for j in range(n+1)] for b,k in conditions])
        inverse_weights=list(M.inv()[n,:])
        assert [F(str(v)) for v in inverse_weights]==[w for b,k,w in weights]
        actual=0
        for w,(b,k) in zip(inverse_weights,conditions):
            jet=s.diff(s.atan(z/epsilon),z,k).subs(z,s.Rational(b.numerator,b.denominator))/s.factorial(k)
            actual+=w*s.series(jet,epsilon,0,1).removeO()
        expected=s.pi*s.Rational(channels['Gamma'].numerator,channels['Gamma'].denominator)/2
        for order in (1,3,5,7):
            value=channels['L'+str(order)]
            expected+=s.Rational(value.numerator,value.denominator)/epsilon**order
        assert s.simplify(actual-expected)==0
        assert all(sum(w*(F(comb(degree,k))*b**(degree-k) if degree>=k else 0)
                       for b,k,w in weights)==int(degree==n) for degree in range(n+1))
        if cfg==4:assert channels['L7']==F(-1,7) and channels['Gamma']==0
        done({'configuration':cfg,'pass':True,'formal_only':True,'not_a_derivative_at_an_ineligible_point':True,
              'nodes':list(map(str,nodes)),'divided_difference_order':n,'channels':{k:str(v) for k,v in channels.items()},
              'actual_Laurent_expression':str(s.simplify(actual)),'normalized_Hermite_oracle_agrees':True,
              'exact_residue_weights':[[str(b),k,str(w)] for b,k,w in weights]})

    # Independent density route: exact time PF, no cone or DD code.
    t=s.symbols('t',positive=True);xx=s.symbols('x0:6',real=True)
    density=t**6/s.prod(t*t+x*x for x in xx)
    location=dict(zip(xx,map(s.Rational,point)))
    def direct_integral_over_pi(field):
        expression=density
        for k in field:expression=s.diff(expression,xx[k])
        expression=s.factor(expression.subs(location))
        partial=s.apart(expression,t)
        value=s.Rational(0);terms=[]
        for term in s.Add.make_args(partial):
            if term==0:continue
            numerator,denominator=s.fraction(term)
            assert not numerator.has(t)
            coefficient,factors=s.factor_list(denominator,t)
            assert len(factors)==1
            quadratic,power=factors[0];poly=s.Poly(quadratic,t)
            assert poly.degree()==2 and poly.nth(2)==1 and poly.nth(1)==0 and poly.nth(0)>0
            a=s.sqrt(poly.nth(0));assert a.is_Rational and a>0
            integral_over_pi=s.Rational(numerator,coefficient)*s.binomial(2*power-2,power-1)/(2**(2*power-1)*a**(2*power-1))
            value+=integral_over_pi
            terms.append({'PF_term':str(term),'integral_over_pi':str(integral_over_pi)})
        return F(str(s.cancel(value))),str(expression),str(partial),terms

    names=['G','D1G','D5G','D11G','D55G','D56G']
    for cfg,(field,name) in enumerate(zip(fields,names),5):
        direct,integrand,partial,terms=direct_integral_over_pi(field)
        total={key:F(0) for key in CHANNELS};tuple_count=0;zero_count=0;repeated_count=0
        for signs in product((-1,1),repeat=6):
            base=[F(sign)*F(x) for sign,x in zip(signs,point)]
            J=[[F(signs[i] if i==j else 0) for j in range(6)] for i in range(6)]
            cone=formal_cone_field(base,J,field)
            for key in CHANNELS:total[key]+=cone[key]
            tuple_count+=1;zero_count+=any(b==0 for b in base);repeated_count+=len(set(base))<6
            receipt['tiny_axis_tuple_field_terms_completed']+=1
        assert tuple_count==64 and zero_count==64 and repeated_count==64
        assert all(total['L'+str(order)]==0 for order in (1,3,5,7))
        assert total['Gamma']/128==direct,(name,total,direct)
        done({'configuration':cfg,'pass':True,'field':name,'physical_coordinate_indices_zero_based':list(field),
              'frame':'axis6, not body32','strict_off_S3_qualification_configuration':1,
              'complete_orthant_terms':64,'zero_base_terms':zero_count,'repeated_base_terms':repeated_count,
              'complete_constant_and_principal_channels':{k:str(v) for k,v in total.items()},
              'direct_density_derivative_without_pi_minus6':integrand,'direct_time_partial_fractions':partial,
              'direct_PF_integrals_over_pi':terms,'direct_integral_over_pi':str(direct),
              'Green_field_exact':str(direct)+'/pi^5','complete_finite_part_matches_direct_integral':True})

    assert len(checks)==10 and receipt['tiny_axis_tuple_field_terms_completed']==384
    for path,digest in pins['sha256'].items():assert sha(path)==digest,path
    result={'status':'PASS_TEN_CONFIGURATIONS_C2_OFF_TRUE_S3_AND_SIX_TINY_AXIS_FIELDS','checks':checks,
            'program_wall_timeout_seconds':180,'pid':os.getpid(),'wall_seconds':time.monotonic()-START,
            'python':sys.version,'sympy':s.__version__,'configurations':10,
            'tiny_complete_frame':'axis6','tiny_complete_fields':6,'tiny_orthant_terms_each_field':64,
            'tiny_tuple_field_contributions':384,'body32_full_kernel_calls':0,'body32_all_triple_qualification_runs':0,
            'old_baseline_reference_reruns':0,'matlab_calls':0,'automatic_retries':0,
            'source_hashes_unchanged':True,'G_C2_off_true_S3_proved_analytically':True,
            'claim_about_actual_body32_derivatives_on_S3':False,'whole_face_continuous_I_certified':False,
            'global_g_certificate':False,'new_C6_bound':False}
    with (P/'results.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    receipt.update(status='PASS_FIXED_BATCH_STOPPED',results_sha256=sha(P/'results.json'))
    print(json.dumps({k:result[k] for k in ['status','configurations','tiny_tuple_field_contributions','body32_full_kernel_calls','wall_seconds']},indent=2))
except BaseException as error:
    receipt.update(status='FAILED_OR_TIMEOUT_STOP_NO_RETRY',error=repr(error));raise
finally:
    signal.setitimer(signal.ITIMER_REAL,0);checkpoint()
