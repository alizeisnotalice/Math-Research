#!/usr/bin/env python3
"""Two new fixed interface certificates; no patterns, densities or body32 data."""
from fractions import Fraction as F
from itertools import combinations
from math import prod
from pathlib import Path
import importlib.util
import json
import sys
import traceback

ROOT=Path(__file__).resolve().parent
OUTPUT={'schema':'B086-minimal-frame-results-v1','status':'STARTING',
        'configurations_started':0,'configurations_completed':0,'configurations':[],
        'body32_constant_instantiations':0,'density_or_kernel_calls':0,
        'scalar_patterns_enumerated':0,'automatic_retries':0}


def serial(value):
    if isinstance(value,F): return str(value)
    if isinstance(value,dict): return {k:serial(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)): return list(map(serial,value))
    return value


def save():
    (ROOT/'results.json').write_text(json.dumps(serial(OUTPUT),indent=2)+'\n')


def require(ok,message):
    if not ok: raise AssertionError(message)


def mm(a,b):
    return [[sum((a[i][k]*b[k][j] for k in range(len(b))),F(0))
             for j in range(len(b[0]))] for i in range(len(a))]


def transpose(a): return list(map(list,zip(*a)))


def add(a,b,scale=F(1)):
    return [[x+scale*y for x,y in zip(ra,rb)] for ra,rb in zip(a,b)]


def main():
    reg=json.loads((ROOT/'registration.json').read_text())
    require(len(reg['configurations'])==2,'Exactly two registered configurations')
    spec=importlib.util.spec_from_file_location('sealed_rank4_aggregate',reg['aggregate_path'])
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    OUTPUT.update(status='RUNNING',active_configuration=reg['configurations'][0]['id'],
                  configurations_started=1)
    save()
    base=reg['fixed_frame']
    H=[[F(v) for v in row] for row in base['unweighted_shear_rows']]
    weights=list(map(F,base['positive_column_weights']))
    require(all(w>0 for w in weights),'Positive spectral weights remain physical column factors')
    I=[[F(i==j) for j in range(6)] for i in range(6)]
    U=add(H,I,F(-1))
    power=I
    Hinv=I
    for k in range(1,6):
        power=mm(power,U)
        Hinv=add(Hinv,power,F((-1)**k))
    require(mm(power,U)==[[F(0)]*6 for _ in range(6)],'Strictly upper triangular nilpotence')
    require(mm(H,Hinv)==I,'Independent finite geometric-series inverse')
    A=[[H[i][j]*weights[j] for j in range(6)] for i in range(6)]
    Ainv=[[Hinv[i][j]/weights[i] for j in range(6)] for i in range(6)]
    det=prod(weights)
    require(mm(A,Ainv)==I and det==F(21,8),'Physical weighted inverse and determinant')
    R,T=F(reg['R']),F(reg['T'])
    # Independent closed form: every one of the 15 four-subsets is independent;
    # each original basis row is tail in C(5,3)=10 and completion in C(5,4)=5.
    a1,a2=4/R,20/R**2
    b1,b2=1+4/R,2+8/R+12/R**2
    dt=(a2+3*a1*a1)*T**3/3+a1*b1*T*T
    dc=(b2+b1*b1)*T+2*a1*b1*T*T
    require(dt==F(4401,512) and dc==F(729,16),'New exact cutoff/time constants')
    scale=1/(R**8*det)
    gram=mm(transpose(Ainv),Ainv)
    expected=[[scale*(10*dt+5*dc)*x for x in row] for row in gram]
    lam=sum((abs(x) for row in Ainv for x in row),F(0))
    hs=[F(1),max(a1,b1),max(a2,b2,max(a1,b1)**2)]
    expected_K=[15*scale*hs[k]*lam**k for k in range(3)]

    cert=module.certificate(A,R,T)
    require(cert['independent_four_count']==15,'Complete m=6 independent-four count')
    records=cert['basis_records']
    require([r['J'] for r in records]==list(combinations(range(6),4)),'Complete lexicographic J list')
    for record in records:
        J,K=record['J'],record['K']
        require(K==tuple(i for i in range(6) if i not in J),'Unique completion is ordered complement')
        require(abs(record['determinant'])==det,'Determinant preserves physical weight product')
        require(record['inverse']==[Ainv[i] for i in J+K],'Solved coordinates use J then K order')
    require(cert['pi6_M4']==expected,'Matrix equals independent 10-tail/5-completion closed form')
    require(cert['pi6_K012']==expected_K,'Scalar aggregate constants retain all15 J contributions')
    require(cert['pi6_trace_M4']==sum(expected[i][i] for i in range(6)),'Trace has correct normalization')
    OUTPUT['configurations'].append({'id':reg['configurations'][0]['id'],'status':'PASS',
        'columns':6,'independent_four_subsets':15,'determinant':det,
        'tail_count_per_physical_column':10,'completion_count_per_physical_column':5,
        'd_tail':dt,'d_completion':dc,'pi6_M4':expected,
        'pi6_K012':expected_K,'pi6_trace_M4':cert['pi6_trace_M4'],
        'reference_inverse_method':'Finite nilpotent geometric series, independently of the sealed Gaussian elimination',
        'scope':'New weighted nonaxis m=6 frame only; no full Green or body32 constant.'})
    OUTPUT['configurations_completed']=1
    OUTPUT['configurations_started']=2
    OUTPUT['active_configuration']=reg['configurations'][1]['id']
    save()

    zero_index=reg['configurations'][1]['insert_zero_column_at_zero_based_index']
    require(zero_index==2,'Fixed zero-column insertion index')
    augmented=[row[:zero_index]+[F(0)]+row[zero_index:] for row in A]
    zcert=module.certificate(augmented,R,T)
    require(zcert['independent_four_count']==15,'Zero insertion cannot create independent subsets')
    for record in zcert['basis_records']:
        require(zero_index not in record['J']+record['K'],'Zero column remains unsolved')
    require(zcert['pi6_M4']==expected and zcert['pi6_K012']==expected_K,
            'Zero-column insertion leaves the exact aggregate certificate unchanged')
    OUTPUT['configurations'].append({'id':reg['configurations'][1]['id'],'status':'PASS',
        'columns':7,'zero_column_index':zero_index,'independent_four_subsets':15,
        'matrix_and_scalar_certificate_equal_to_previous':True,
        'scope':'One fixed zero insertion. No pattern enumeration and no density evaluation.'})
    OUTPUT.update(status='PASS_BOTH_NEW_CONFIGURATIONS',configurations_completed=2,
                  active_configuration=None,python=sys.version,
                  dependencies=['Python standard library Fraction','Source-pinned aggregate.py'])
    save()
    print(json.dumps({'status':OUTPUT['status'],'new_configurations':2,
                      'independent_four_records_checked_total':30,'kernel_calls':0}))


if __name__=='__main__':
    try:
        main()
    except BaseException as exc:
        OUTPUT.update(status='FAILED_FIRST_STOP',error=repr(exc),traceback=traceback.format_exc())
        save()
        raise
