#!/usr/bin/env python3
"""F052 independent four-configuration audit; two REAL 64-term tiny-frame references.
No author code import/execution, no production body32/152 frame or native engine.
DD coefficients use Fraction negative-binomial convolution, not matrix inverse,
SymPy residue differentiation, or a triangular divided-difference evaluator.
"""
import sys
sys.dont_write_bytecode = True
import json, hashlib, os, time, datetime, traceback
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter
from itertools import product, combinations
from math import comb, factorial
import sympy as s
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

def check_pins(reg):
    assert time.time()-reg['actual_started_unix']<reg['max_actual_seconds']
    for p,h in reg['source_sha256'].items(): assert sha(p)==h,p

def coefficients(nodes):
    """Residues of f(z)/prod(z-node); all repeated and zero nodes retained."""
    assert len(nodes)==6
    mult=Counter(map(Q,nodes)); answer=[Q(0)]*4
    for b,m in sorted(mult.items()):
        poly=[Q(1)]+[Q(0)]*(m-1)
        for c,n in sorted(mult.items()):
            if c==b: continue
            delta=b-c; assert delta
            factor=[Q((-1)**r*comb(n+r-1,r),1)/delta**(n+r) for r in range(m)]
            poly=[sum((poly[k]*factor[r-k] for k in range(r+1)),Q(0)) for r in range(m)]
        # Normalized jet k of atan(b/eps): constants / singular channels.
        for k in range(m):
            weight=poly[m-1-k]
            if b and k==0: answer[0]+=weight*(1 if b>0 else -1)
            if not b and k in (1,3,5): answer[(k+1)//2]+=weight*Q((-1)**((k-1)//2),k)
    return answer

def axis_reference(point, expected, gate):
    rows=[]; total=[Q(0)]*4
    for signs in product((-1,1),repeat=6):
        nodes=[Q(sg)*Q(x) for sg,x in zip(signs,point)]
        values=coefficients(nodes)
        total=[a+b for a,b in zip(total,values)]
        rows.append({'signs':list(signs),'nodes':list(map(str,nodes)),'alpha':'1','Gamma_L1_L3_L5':list(map(str,values))})
    assert len(rows)==64 and total==list(map(Q,expected)),total
    distances=[]
    for i,j in combinations(range(6),2):
        d=sum((Q(x)**2 for k,x in enumerate(point) if k not in (i,j)),Q(0))
        assert d>=0
        distances.append({'pair_1based':[i+1,j+1],'distance_squared':str(d)})
    rho=min(Q(r['distance_squared']) for r in distances)
    assert (rho>0)==(gate=='OFF_S2_FINITE')
    return {'point':list(map(str,point)),'frame':'six standard coordinate columns',
            'complete_fan':'64 signed coordinate orthants, boundaries null, alpha=q=1',
            'full_reference_functional_terms':64,'terms_skipped':0,
            'nodes_equal_zero_are_retained':True,'Gamma_L1_L3_L5':list(map(str,total)),
            'true_S2_all_15_pair_distances':distances,'minimum_distance_squared':str(rho),
            'eligibility':gate,'term_ledger':rows}

def run():
    reg=json.loads((HERE/'input_manifest.json').read_text()); check_pins(reg)
    e=s.Symbol('epsilon',positive=True); b=s.Symbol('b',positive=True); z=s.Symbol('z',real=True)
    t=s.Symbol('t',positive=True)
    cases=[]
    # CONFIGURATION 1: one symbolic jet family, including order zero.
    jets=[]
    for k in range(6):
        jet=s.diff(s.atan(z/e),z,k)/factorial(k)
        at_zero=s.simplify(jet.subs(z,0))
        expected=0 if k%2==0 else s.Rational((-1)**((k-1)//2),k)/e**k
        assert s.simplify(at_zero-expected)==0
        limit=s.limit(jet.subs(z,b),e,0,dir='+')
        assert s.simplify(limit-(s.pi/2 if k==0 else 0))==0
        jets.append({'k':k,'normalized_zero_jet':str(at_zero),'nonzero_positive_b_constant':str(limit)})
    assert s.simplify(s.atan(-z/e)+s.atan(z/e))==0
    assert s.simplify((2*s.pi)**-6*(s.pi/2)-1/(128*s.pi**5))==0
    cases.append({'id':1,'type':'one symbolic normalized-jet family k=0..5 plus normalization','jets':jets,
                  'negative_b_by_oddness':True,'normalization':'1/(128*pi**5)'})
    # CONFIGURATION 2: a new repeated-node scalar, disjoint from author's 12 tuples.
    nodes=[0,0,0,0,0,3]
    value=coefficients(nodes); assert value==[Q(1,243),Q(-1,81),Q(1,27),Q(0)]
    f=(s.atan(3/e)-3/e+9/e**3)/243
    # Exact Taylor-remainder DD: f(3) minus degree-four jet at zero, divided by 3^5.
    laurent=s.series((s.pi/2-s.atan(e/3)-3/e+9/e**3)/243,e,0,1).removeO()
    assert s.simplify(laurent-(s.pi/486-1/(81*e)+1/(27*e**3)))==0
    # Resolvent phase: real product for nodes 0^5,3 equals 1/[e^4(e^2+9)].
    assert s.simplify(s.diff(f,e)+1/(e**4*(e**2+9)))==0
    assert coefficients([-x for x in nodes])==value
    cases.append({'id':2,'nodes':nodes,'Gamma_L1_L3_L5':list(map(str,value)),
                  'independent_Taylor_remainder_expression':str(f),'Laurent_part':str(laurent),
                  'negative_time_derivative':'1/(epsilon**4*(epsilon**2+9))',
                  'global_node_negation_invariance':True,'spatial_G_claim':False})
    # CONFIGURATION 3: actual complete tiny analytic frame, off its true S2.
    finite=axis_reference([1,2,3,0,0,0],[Q(16,15),0,0,0],'OFF_S2_FINITE')
    direct=1/((t*t+1)*(t*t+4)*(t*t+9))
    pieces=s.Rational(1,24)/(t*t+1)-s.Rational(1,15)/(t*t+4)+s.Rational(1,40)/(t*t+9)
    assert s.cancel(direct-pieces)==0
    integral=s.pi/2*(s.Rational(1,24)-s.Rational(1,15)/2+s.Rational(1,40)/3)
    assert s.simplify(integral-s.pi/120)==0
    assert s.simplify(s.Rational(finite['Gamma_L1_L3_L5'][0])/(128*s.pi**5)-integral/s.pi**6)==0
    # Orthant sum identity follows from product_j[(t-ix_j)^-1+(t+ix_j)^-1].
    orthant_product=s.prod(2*t/(t*t+s.Rational(x)**2) for x in [1,2,3,0,0,0])
    assert s.cancel(orthant_product/64-direct)==0
    finite.update({'id':3,'density_times_pi6':str(direct),'exact_direct_G':'1/(120*pi**5)',
                   'reference_finite_part_is_G':True,'density_complete_orthant_identity':True})
    cases.append(finite)
    # CONFIGURATION 4: actual complete tiny frame on S2; finite part NOT Green.
    singular=axis_reference([1,1,0,0,0,0],[-96,64,0,0],'S2_REJECT_INFINITE_GREEN')
    d=1/(t*t*(t*t+1)**2)
    assert s.cancel(d-(1/t**2-1/(t*t+1)-1/(t*t+1)**2))==0
    H=1/e-3*s.pi/4+s.Rational(3,2)*s.atan(e)+e/(2*(1+e*e))
    assert s.simplify(s.diff(H,e)+d.subs(t,e))==0
    assert s.limit(H,e,s.oo)==0 and s.limit(H,e,0,dir='+')==s.oo
    expected=s.pi*s.Rational(singular['Gamma_L1_L3_L5'][0])/128+s.Rational(singular['Gamma_L1_L3_L5'][1])/64/e
    assert s.simplify(s.series(H,e,0,1).removeO()-expected)==0
    assert s.simplify(s.Rational(singular['Gamma_L1_L3_L5'][0])/(128*s.pi**5)+s.Rational(3,4)/s.pi**5)==0
    singular.update({'id':4,'density_times_pi6':str(d),'H_epsilon_times_pi6':str(H),
                     'actual_G':'positive_infinity','constant_finite_part':'-3/(4*pi**5)',
                     'reference_finite_part_is_G':False,'hard_rejection_required':True})
    cases.append(singular)
    check_pins(reg)
    return {'status':'PASS_FOUR_INDEPENDENT_CONFIGURATIONS','cases':cases,
            'actual_new_configurations':4,'tiny_reference_full_objects':2,
            'tiny_reference_full_terms':128,'author_twelve_case_rerun':False,
            'production_body32_or_152_full_kernels':0,'native_calls':0,
            'gradient_or_Hessian_extension':False,'new_C6_bound':False}

if __name__=='__main__':
    assert not (HERE/'results.json').exists() and not (HERE/'failure_receipt.json').exists(),'one actual batch only'
    start=time.monotonic()
    try:
        result=run()
        (HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n')
        receipt={'status':'PASS','pid':os.getpid(),'exit_code':0,'seconds':time.monotonic()-start,
                 'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
                 'source_sha256':sha(__file__),'input_manifest_sha256':sha(HERE/'input_manifest.json'),
                 'results_sha256':sha(HERE/'results.json'),'python':sys.version,
                 'python_executable':sys.executable,'python_executable_sha256':sha(sys.executable),
                 'sympy_version':s.__version__,'sympy_module_sha256':sha(s.__file__),
                 'configurations':4,'tiny_complete_references':2,'tiny_complete_terms':128,
                 'production_full_kernels':0,'native_calls':0,'retry_count':0,'child_processes':[]}
        (HERE/'run_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
        print(json.dumps({k:receipt[k] for k in ('status','pid','seconds','configurations','tiny_complete_references','tiny_complete_terms','production_full_kernels')}))
    except BaseException as exc:
        failure={'status':'STOP_FAILURE','pid':os.getpid(),'seconds':time.monotonic()-start,
                 'exception':repr(exc),'traceback':traceback.format_exc(),'source_sha256':sha(__file__)}
        (HERE/'failure_receipt.json').write_text(json.dumps(failure,indent=2)+'\n')
        raise
