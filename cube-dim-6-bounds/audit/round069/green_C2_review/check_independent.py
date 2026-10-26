"""Seven preregistered exact checks. No full-frame kernel or old reference replay."""
import signal,time
START=time.monotonic()
def timeout(signum,frame): raise TimeoutError('B069 single-attempt180-second limit')
signal.signal(signal.SIGALRM,timeout)
signal.setitimer(signal.ITIMER_REAL,180)
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from math import comb,lcm
import hashlib,json,os,sys,importlib.util
import sympy as s

HERE=Path(__file__).resolve().parent
SRC=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-literature/cube-dim-6-bounds/upper/round069/green_C2_off_S3')
sys.dont_write_bytecode=True
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:json.loads(Path(p).read_text())
if (HERE/'execution_receipt.json').exists() or (HERE/'results.json').exists():
    raise RuntimeError('Single-attempt batch already has a receipt; no replay')
reg=read(HERE/'registration.json')
receipt={'status':'RUNNING','pid':os.getpid(),'configurations_completed':0,
         'program_wall_timeout_seconds':180,'automatic_retries':0,
         'body32_kernel_calls':0,'body32_qualification_runs':0,
         'full_fan_or_orthant_sums':0,'old_reference_replays':0}
checks=[]
def save(name,obj): (HERE/name).write_text(json.dumps(obj,indent=2,ensure_ascii=False)+'\n')
def checkpoint():
    receipt['wall_seconds']=time.monotonic()-START
    save('execution_receipt.json',receipt)
    save('checkpoint.json',{'checks':checks,'receipt':receipt})
def done(row):
    row['pass']=True;checks.append(row);receipt['configurations_completed']=len(checks);checkpoint()
def verify_sources():
    assert sha(SRC/'manifest.json')==reg['producer_manifest_sha256']
    final=read(SRC/'manifest.json');pre=read(SRC/'input_manifest.json')
    for name,digest in final['artifact_sha256'].items(): assert sha(SRC/name)==digest,name
    for path,digest in pre['sha256'].items(): assert sha(path)==digest,path
    return {'manifest_sha256':sha(SRC/'manifest.json'),'artifact_pins':final['artifact_sha256'],
            'pre_execution_pins':pre['sha256']}

# Independent Fraction Gaussian elimination, not SymPy projection inversion.
def rank(rows):
    A=[list(map(F,row)) for row in rows]
    if not A:return 0
    r=0
    for j in range(len(A[0])):
        pivot=next((i for i in range(r,len(A)) if A[i][j]),None)
        if pivot is None:continue
        A[r],A[pivot]=A[pivot],A[r];d=A[r][j];A[r]=[v/d for v in A[r]]
        for i in range(r+1,len(A)):
            q=A[i][j]
            A[i]=[u-q*v for u,v in zip(A[i],A[r])]
        r+=1
        if r==len(A):break
    return r
def det(rows):
    A=[list(map(F,row)) for row in rows];n=len(A);answer=F(1)
    for j in range(n):
        pivot=next((i for i in range(j,n) if A[i][j]),None)
        if pivot is None:return F(0)
        if pivot!=j:A[j],A[pivot]=A[pivot],A[j];answer=-answer
        d=A[j][j];answer*=d
        for i in range(j+1,n):
            q=A[i][j]/d
            for k in range(j,n):A[i][k]-=q*A[j][k]
    return answer
def dot(a,b):return sum((x*y for x,y in zip(a,b)),F(0))
def gram(columns):return [[dot(a,b) for b in columns] for a in columns]
def independent_indices(columns,indices):
    selected=[]
    for j in indices:
        if rank([columns[k] for k in selected+[j]])>len(selected):selected.append(j)
    return selected
def residual(columns,x,indices):
    selected=independent_indices(columns,indices)
    U=[columns[j] for j in selected]
    D=det(gram(U));assert D>0
    return det(gram(U+[x]))/D,selected
def all_residuals(columns,x,k):
    return [{'indices':list(J),'residual':str(residual(columns,x,J)[0]),
             'rank':len(residual(columns,x,J)[1])} for J in combinations(range(len(columns)),k)]
def point(coeff):return [sum(F(coeff[j])*basis[j][i] for j in range(6)) for i in range(6)]
def invert(rows):
    n=len(rows);A=[[F(v) for v in row]+[F(i==j) for j in range(n)] for i,row in enumerate(rows)]
    for j in range(n):
        i=next(i for i in range(j,n) if A[i][j]);A[j],A[i]=A[i],A[j]
        d=A[j][j];A[j]=[v/d for v in A[j]]
        for i in range(n):
            if i!=j:
                q=A[i][j];A[i]=[u-q*v for u,v in zip(A[i],A[j])]
    return [row[n:] for row in A]
def reject(call):
    try:call()
    except (ValueError,TypeError):return True
    raise AssertionError('An unsupported contract input was accepted')

try:
    checkpoint();sources=verify_sources()
    spec=importlib.util.spec_from_file_location('frozen_A069_contract',SRC/'derivative_contract.py')
    api=importlib.util.module_from_spec(spec);spec.loader.exec_module(api)
    y=s.symbols('y',real=True);den=(1+y*y)**2
    score=-2*y/(1+y*y);logsecond=2*(y*y-1)/den
    polys=[s.expand(den*(1-score**2)),s.expand(s.cancel(den*(2-logsecond))),
           s.expand(s.cancel(den*(2+logsecond)))]
    assert s.simplify(polys[0]-(y*y-1)**2)==0
    assert polys[1]==2*y**4+2*y**2+4 and polys[2]==2*y**4+6*y**2
    done({'id':1,'score_unit_bound_polynomial':str(s.factor(polys[0])),
          'log_second_two_sided_nonnegative_polynomials':list(map(str,polys[1:])),
          'interpretation':'y=z/t; fixed-t positive kernel relative derivative bounds, no indicator differentiation'})

    basis=[[F(v) for v in col] for col in reg['weighted_fixture_basis_columns']]
    columns=basis+[[2*v for v in basis[0]],
                   [a+b for a,b in zip(basis[0],basis[1])],
                   [a+b for a,b in zip(basis[1],basis[2])]]
    x=point([1,2,4,8,16,32]);oracle=all_residuals(columns,x,3)
    q=api.qualify_off_S3(columns,x)
    assert q['status']=='OFF_TRUE_S3' and len(oracle)==84
    for expected,actual in zip(oracle,q['triples']):
        assert expected['indices']==actual['triple']
        assert expected['residual']==actual['squared_residual'] and expected['rank']==actual['rank']
    rho=min(F(row['residual']) for row in oracle);assert rho>0 and str(rho)==q['min_squared_distance']
    assert any(row['rank']<3 for row in oracle)
    done({'id':2,'columns':[[str(v) for v in c] for c in columns],'point':list(map(str,x)),
          'independent_augmented_Gram_triples':oracle,'producer_agrees':True,
          'min_squared_distance':str(rho),'frame_sha256':q['frame_sha256']})

    xs3=point([1,0,0,1,1,0]);pairs=all_residuals(columns,xs3,2);triples=all_residuals(columns,xs3,3)
    pairmin=min(F(v['residual']) for v in pairs);triplemin=min(F(v['residual']) for v in triples)
    assert len(pairs)==36 and len(triples)==84 and pairmin>0 and triplemin==0
    qs3=api.qualify_off_S3(columns,xs3);assert qs3['status']=='ON_TRUE_S3_UNSUPPORTED'
    assert reject(lambda:api.local_tail_contract(qs3,'0','1/7'))
    assert reject(lambda:api.bind_derivative_request(qs3,{},[[]]))
    done({'id':3,'point':list(map(str,xs3)),'all_pair_count':36,'all_triple_count':84,
          'off_S2_min_squared_distance':str(pairmin),'on_S3_min_squared_distance':str(triplemin),
          'zero_residual_triples':[v for v in triples if F(v['residual'])==0],
          'derivative_tail_and_request_rejected':True,'actual_differentiability_not_decided':True})

    order=[0,1,6,7,3,4];chosen=independent_indices(columns,order)
    assert rank([columns[j] for j in order[:4]])==2 and len(chosen)==4
    xlarge=[sum(columns[j][i] for j in order) for i in range(6)]
    distances=all_residuals(columns,xlarge,3);assert min(F(v['residual']) for v in distances)>0
    d=lcm(*(v.denominator for col in columns for v in col));quad=[]
    for J in combinations(range(9),4):
        D=det(gram([columns[j] for j in J]));scaled=D*d**8
        assert D>=0 and scaled.denominator==1 and (D==0 or scaled>=1)
        quad.append({'indices':list(J),'rank':rank([columns[j] for j in J]),'Gram_det':str(D),'integer_scaled_Gram_det':str(scaled)})
    assert len(quad)==126 and any(v['Gram_det']=='0' for v in quad)
    done({'id':4,'large_order':order,'first_four_rank':2,'rank_greedy_independent_four':chosen,
          'point':list(map(str,xlarge)),'min_true_S3_squared_distance':str(min(F(v['residual']) for v in distances)),
          'common_denominator':d,'all126_quadruple_Gram_records':quad,
          'independent_quadruple_count':sum(v['Gram_det']!='0' for v in quad)})

    # A different Fraction inversion reconstructs all fixed-basis coordinate constants.
    N=[[basis[j][i] for j in range(6)] for i in range(6)];inv=invert(N)
    assert [[str(v) for v in row] for row in inv]==q['basis_inverse']
    C=sum(v*v for row in inv for v in row);c=1/(1+C);S=sum(abs(v) for col in columns for v in col)
    delta0=min(rho,F(1))/2;radius=delta0/3;delta=delta0-radius;eps=F(1,7)
    R=delta/(2*S);M=F(comb(9,4)*d**4)/(64*c*c*R**8)
    c1=[sum(abs(inv[i][j]) for i in range(6)) for j in range(6)]
    c2=[[c1[j]*c1[k]+2*sum(abs(inv[i][j]*inv[i][k]) for i in range(6)) for k in range(6)] for j in range(6)]
    tail=api.local_tail_contract(q,radius,eps);half=api.local_tail_contract(q,radius,eps/2)
    assert F(tail['Mbar_delta'])==M and F(tail['uniform_distance_delta'])==delta
    assert list(map(F,tail['C1_coordinate']))==c1
    assert [[F(v) for v in row] for row in tail['C2_coordinate']]==c2
    assert F(tail['value_tail_upper'])==M*eps**3/3==8*F(half['value_tail_upper'])
    for j in range(6):
        assert F(tail['gradient_tail_absolute'][j])==M*c1[j]*eps**2/2==4*F(half['gradient_tail_absolute'][j])
        for k in range(6):assert F(tail['Hessian_tail_absolute'][j][k])==M*c2[j][k]*eps==2*F(half['Hessian_tail_absolute'][j][k])
    bad_tail=[('radius_at_margin',delta0,eps),('negative_radius',F(-1),eps),
              ('zero_epsilon',F(0),F(0)),('negative_epsilon',F(0),F(-1))]
    for name,r,e in bad_tail:assert reject(lambda r=r,e=e:api.local_tail_contract(q,r,e))
    fan={'frame_sha256':q['frame_sha256'],'complete_fan_proved':True,'multiplicities_included_once':True,
         'proof_reference':'Synthetic API boundary fixture only; no fan is computed or certified by this test.'}
    request=api.bind_derivative_request(q,fan,[[],[0],[0,5]])
    assert request['fan_proof_independently_authenticated_by_this_function'] is False
    bad_fan=[dict(fan,frame_sha256='0'*64),dict(fan,complete_fan_proved=False),
             dict(fan,multiplicities_included_once=False),dict(fan,proof_reference='')]
    for item in bad_fan:assert reject(lambda item=item:api.bind_derivative_request(q,item,[[]]))
    for field in ([6],[-1],[True],[0,1,2]):assert reject(lambda field=field:api.bind_derivative_request(q,fan,[field]))
    assert reject(lambda:api.qualify_off_S3([[0]*6]+columns,x))
    assert reject(lambda:api.qualify_off_S3([basis[0]]*6,x))
    assert reject(lambda:api.qualify_off_S3(columns,[0.5]+x[1:]))
    done({'id':5,'local_tail_contract':tail,'epsilon_halving_ratios':[8,4,2],
          'bad_tail_rejections':[v[0] for v in bad_tail],'fan_boundary_rejections':4,'field_boundary_rejections':4,
          'zero_column_rank_float_rejections':3,'external_fan_authentication_remains_conditional':True})

    # DD_n(z^k)=h_(k-n), with h complete homogeneous. Work in labelled variables
    # before any collisions; this is independent of the producer residue code.
    def h(vars,n):
        if n<0:return s.Integer(0)
        coeff=[s.Integer(1)]+[s.Integer(0)]*n
        for v in vars:
            new=[sum(coeff[j]*v**(k-j) for j in range(k+1)) for k in range(n+1)]
            coeff=new
        return s.expand(coeff[n])
    def ddpoly(nodes,n):return h(nodes,8-n)+3*h(nodes,7-n)-2*h(nodes,5-n)
    labels=s.symbols('b0:6');values=list(map(s.Rational,['0','0','1/2','1/2','-2','3']))
    dd=ddpoly(labels,5);sub=dict(zip(labels,values));rows=[]
    for field in ((0,),(0,0),(0,1),(2,3)):
        direct=dd
        for i in field:direct=s.diff(direct,labels[i])
        direct=s.expand(direct).subs(sub)
        copied=values+[values[i] for i in field]
        factor=2 if len(field)==2 and field[0]==field[1] else 1
        other=s.expand(factor*ddpoly(copied,5+len(field)))
        assert direct==other
        rows.append({'labels':list(field),'equal_values_distinct_labels':len(field)==2 and field[0]!=field[1] and values[field[0]]==values[field[1]],
                     'diagonal_factor':factor,'exact_derivative':str(direct)})
    assert rows[1]['exact_derivative']!='0' and rows[2]['exact_derivative']!='0'
    done({'id':6,'base_nodes':list(map(str,values)),'polynomial':'z^8+3*z^7-2*z^5','labelled_checks':rows})

    nodes=list(map(F,['-2','-2','0','0','0','0','1/2','1/2']))
    unique=sorted(set(nodes));conditions=[(b,k) for b in unique for k in range(nodes.count(b))]
    H=s.Matrix([[s.binomial(j,k)*s.Rational(b.numerator,b.denominator)**(j-k) if j>=k else 0 for j in range(8)] for b,k in conditions])
    weights=list(H.inv()[7,:]);epsilon=s.symbols('epsilon',positive=True);z=s.symbols('z',real=True)
    actual=0
    for w,(b,k) in zip(weights,conditions):
        jet=s.diff(s.atan(z/epsilon),z,k).subs(z,s.Rational(b.numerator,b.denominator))/s.factorial(k)
        actual+=w*s.series(jet,epsilon,0,1).removeO()
    channels,_=api.formal_dd_channels(nodes)
    expected=s.pi*s.Rational(str(channels['Gamma']))/2+sum(s.Rational(str(channels['L'+str(k)]))/epsilon**k for k in (1,3,5,7))
    assert s.simplify(actual-expected)==0 and any(channels['L'+str(k)] for k in (1,3,5,7))
    done({'id':7,'formal_only':True,'nodes':list(map(str,nodes)),
          'Hermite_weights':list(map(str,weights)),'actual_atan_Laurent_through_constant':str(s.simplify(actual)),
          'producer_channels':{k:str(v) for k,v in channels.items()},
          'nonzero_principal_channels_retained':True,'no_Green_derivative_claim_for_single_tuple':True})

    assert len(checks)==7;assert sources==verify_sources()
    result={'status':'PASS_SEVEN_NEW_EXACT_CHECKS','checks':checks,'scientific_configurations':7,
            'pid':os.getpid(),'elapsed_seconds':time.monotonic()-START,'python':sys.version,'sympy':s.__version__,
            'source_pins':sources,'registration_sha256':sha(HERE/'registration.json'),
            'checker_sha256':sha(__file__),'body32_kernel_calls':0,'body32_qualification_runs':0,
            'old_reference_replays':0,'orthant_sums':0,'new_C6_bound':False}
    save('results.json',result);receipt.update(status='PASS_SEVEN_NEW_EXACT_CHECKS_STOPPED',results_sha256=sha(HERE/'results.json'))
    print(json.dumps({'status':result['status'],'scientific_configurations':7,'pid':result['pid'],
                      'elapsed_seconds':result['elapsed_seconds'],'body32_kernel_calls':0,'old_reference_replays':0}))
except BaseException as error:
    receipt.update(status='FAILED_STOP_NO_AUTOMATIC_RETRY',error=repr(error));raise
finally:
    signal.setitimer(signal.ITIMER_REAL,0);checkpoint()
