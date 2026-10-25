#!/usr/bin/env python3
"""One-candidate Fraction/geometry audit; imports no author code, calls no union."""
from decimal import Decimal, localcontext
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
from math import prod
import hashlib
import json
import os
import sys
import time

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
MAIN = HERE.parents[2]
DROOT = HERE.parents[2]/'lower/round008/spacing_chamber'
EROOT = HERE.parents[2]/'audit/lower008/parametric_union'
DEN = 10**120
ZERO = (Q(0),)*36


def field(d):
    result = list(ZERO)
    for key,value in d.items():
        i,j = map(int,key.split(','))
        assert 0 <= i < 6 and 0 <= j < 6 and type(value) in (int,str)
        result[6*i+j] = Q(value)
    return tuple(result)


def add(a,b):
    return tuple(x+y for x,y in zip(a,b))


def scale(a,c):
    return tuple(c*x for x in a)


def sub(a,b):
    return add(a,scale(b,-1))


def monomial(i,j,c=Q(1)):
    out = list(ZERO)
    out[6*(i%6)+j%6] = c*2**(i//6)*3**(j//6)
    return tuple(out)


def root_interval(m):
    # Own isolator, no author rootfloor/floor_root imports.
    lo,hi = DEN,2*DEN
    target = m*DEN**6
    assert lo**6 <= target < hi**6
    while hi-lo > 1:
        mid = (lo+hi)//2
        if mid**6 <= target:
            lo = mid
        else:
            hi = mid
    assert lo**6 <= target < (lo+1)**6
    return Q(lo,DEN),Q(lo+1,DEN)


UL,UH = root_interval(2)
VL,VH = root_interval(3)


@lru_cache(None)
def bounds(a):
    lo=hi=Q(0)
    for k,c in enumerate(a):
        if c:
            i,j = divmod(k,6)
            low,high=UL**i*VL**j,UH**i*VH**j
            lo += c*(low if c>0 else high)
            hi += c*(high if c>0 else low)
    return lo,hi


def sign(a):
    if a == ZERO:
        return 0
    lo,hi = bounds(a)
    if lo>0:
        return 1
    if hi<0:
        return -1
    raise AssertionError('unresolved sign; no tolerance or rounding equality')


def radius(ell):
    mass=prod(ell)
    m=mass
    p2=p3=0
    while m%2==0:
        p2+=1
        m//=2
    while m%3==0:
        p3+=1
        m//=3
    assert m==1
    # (2r)^6=2^p2 3^p3=mass, with positive roots.
    return monomial(p2,p3,Q(1,2))


def decimal(q):
    with localcontext() as ctx:
        ctx.prec=90
        return str(Decimal(q.numerator)/Decimal(q.denominator))


def main():
    started=time.monotonic()
    paths=[MAIN/'audit/lower008/primary/inward_smallstep.py',
           MAIN/'audit/lower008/primary/inward_smallstep_results.json',
           MAIN/'rounds/008/frozen_lower_before_candidate.json',
           MAIN/'coordination/phase_20261025/frozen_thresholds.json',
           MAIN/'lower/round001_rectangular_certificate.py',
           MAIN/'lower/round001_transfer_proof.md',
           DROOT/'chamber.py',DROOT/'check_interface_geometry.py',
           DROOT/'chamber_interface.json',DROOT/'chamber_results.json',DROOT/'geometry_results.json',
           EROOT/'parametric_union.py',EROOT/'candidate_results.json',EROOT/'interface_proof.md',
           HERE/'inward_smallstep.py',HERE/'inward_smallstep_results.json']
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    primary=json.loads(paths[1].read_text())
    replay=json.loads((HERE/'inward_smallstep_results.json').read_text())
    stable=lambda d:{k:v for k,v in d.items() if k not in ('elapsed_seconds','code_sha256')}
    assert stable(primary)==stable(replay)
    original=paths[0].read_text()
    assert (HERE/'inward_smallstep.py').read_text()==original.replace(
        'ROOT=Path(__file__).resolve().parents[3]',f'ROOT=Path("{MAIN}")')
    R=Q(replay['exact_rational_lower'])
    target=Q(2067397333267,10**12)
    finite=replay['finite_L1_witness']
    K=finite['receiver_index_radius']
    theta=(R+target)/(2*R)
    eta=Q(finite['bump_halfwidth'])
    assert K==2**47 and finite['source_index_radius']==K+2
    assert 0<theta<1 and theta==Q(finite['threshold'])
    assert eta==Q(5,10**36)
    factor=Q(2*K+1,2*K+5)**6
    finite_ratio=theta*R*factor
    phase=json.loads((MAIN/'coordination/phase_20261025/frozen_thresholds.json').read_text())
    frozen=json.loads((MAIN/'rounds/008/frozen_lower_before_candidate.json').read_text())
    phase_upper=Q(phase['lower']['exact_interval'][1])
    old_upper=Q(frozen['lower']['exact_interval'][1])
    assert frozen['lower']['id']=='R6-007-spacing2-smallstep'
    assert finite_ratio>target>old_upper>phase_upper
    packet=json.loads((DROOT/'chamber_interface.json').read_text())
    d=json.loads((DROOT/'chamber_results.json').read_text())
    e=json.loads((EROOT/'candidate_results.json').read_text())
    p=[field(c) for c in e['ratio_polynomial_coefficients']]
    assert len(p)==2 and e['degree']==1
    value=add(p[0],p[1])
    assert value==field(e['candidate_ratio_coefficients'])==field(d['candidate']['coefficients'])
    assert p[0]==field(d['current_Rideal_coefficients'])
    assert Q(e['previous_complete_Rideal_upper'])==old_upper
    assert Q(packet['normalization_factor'])==64 and Q(e['normalization_factor'])==64
    assert sign(p[1])>0
    ideal_lo,ideal_hi=bounds(value)
    e_lo,e_hi=map(Q,e['candidate_ratio_interval'])
    d_lo,d_hi=map(Q,d['candidate']['exact_interval'])
    assert e_lo<=ideal_lo<=ideal_hi<=e_hi
    assert d_lo<=ideal_lo<=ideal_hi<=d_hi
    assert ideal_lo>R>target and ideal_lo>old_upper and ideal_lo>phase_upper

    # Independently verify all exported source slabs and all order gaps;
    # no clipping-selection assumptions and no union routine are imported.
    ep={name:(field(x['c0']),field(x['c1'])) for name,x in packet['endpoints'].items()}
    evaluate=lambda name,t:add(ep[name][0],scale(ep[name][1],t))
    s0,s1=field(packet['spacing2']['c0']),field(packet['spacing2']['c1'])
    sstar=add(monomial(2,2,Q(1,3)),monomial(2,4,Q(1,6)))
    assert add(s0,s1)==sstar
    assert len(packet['boxes'])==len(packet['box_lengths'])==127
    source_checks=cell_checks=0
    max_reflected_source_index=0
    for t in (Q(0),Q(1)):
        spacing=(add(monomial(0,0),monomial(4,0,Q(1,2))),add(s0,scale(s1,t)),
                 *((add(monomial(0,0,Q(1,2)),monomial(1,0,Q(1,2))),)*4))
        for box,ell in zip(packet['boxes'],packet['box_lengths']):
            assert all(type(l) is int and 1<=l<=4 for l in ell)
            r=radius(ell)
            for ids,l,a in zip(box,ell,spacing):
                lo,hi=(evaluate(name,t) for name in ids)
                assert sign(lo)>=0 and sign(sub(hi,lo))>0 and sign(sub(scale(a,Q(1,2)),hi))>=0
                cell_checks+=3
                q=-(l//2) if l%2 else 1-l//2
                max_reflected_source_index=max(max_reflected_source_index,abs(q),abs(q+l-1))
                for x in (lo,hi):
                    for index in (q,q+l-1):
                        delta=sub(x,scale(a,index))
                        assert sign(add(r,delta))>=0 and sign(sub(r,delta))>=0
                        source_checks+=1
    assert max_reflected_source_index==2
    gaps=ties=0
    for groups in packet['order_certificate']:
        for group in groups:
            assert all(ep[name]==ep[group[0]] for name in group)
        for low,high in zip(groups,groups[1:]):
            signs=[sign(sub(evaluate(high[0],t),evaluate(low[0],t))) for t in (Q(0),Q(1))]
            assert min(signs)>=0 and max(signs)>0
            gaps+=1
            ties+=0 in signs
    assert source_checks==6096 and gaps==120 and ties==1
    result={
      'status':'PASS', 'scope':'one lower008 candidate; independent full-chain audit, no search',
      'pid':os.getpid(),'root_denominator':str(DEN),'root_intervals':{
          'u':[str(UL),str(UH)],'v':[str(VL),str(VH)]},
      'primary_replay_count':1,'primary_replay_exit_code':0,
      'primary_replay_results_match_except_timing_and_source_hash':True,
      'union_calls_in_this_independent_checker':0,'author_imports':0,
      'candidate_boxes':127,'inward_retained_boxes':96,
      'source_endpoint_checks':source_checks,'cell_face_and_width_checks':cell_checks,
      'order_gaps_checked':gaps,'boundary_ties':ties,
      'maximum_reflected_source_index':max_reflected_source_index,
      'K':K,'source_index_radius':K+2,'eta':str(eta),'theta':str(theta),
      'theta_strictly_between_zero_and_one':True,
      'boundary_factor':str(factor),'finite_ratio':str(finite_ratio),
      'finite_ratio_decimal':decimal(finite_ratio),
      'finite_ratio_minus_target':str(finite_ratio-target),
      'finite_ratio_minus_round008_frozen_complete_007_upper':str(finite_ratio-old_upper),
      'finite_ratio_minus_phase_frozen_complete_upper':str(finite_ratio-phase_upper),
      'all_36_candidate_coefficients_agree_D_E':True,
      'polynomial_degree':1,'polynomial_constant_all_36_match_007':True,
      'ideal_interval_120_digits':[str(ideal_lo),str(ideal_hi)],
      'ideal_interval_width':str(ideal_hi-ideal_lo),
      'ideal_lower_decimal':decimal(ideal_lo),
      'tighter_interval_inside_D_and_E_intervals':True,
      'claim':'C6 > 2067397333267/1000000000000 for Lebesgue centered axis-parallel common-radius cubes on general L1, all positive scales',
      'new_mechanism':False,'no_MATLAB':True,'no_agents_spawned':True,
      'input_sha256':hashes,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      'elapsed_seconds':time.monotonic()-started,
      'actual_process_status':'finished; no further lower calculation authorized by this audit'}
    (HERE/'independent_review_results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('status','pid','primary_replay_count','union_calls_in_this_independent_checker',
        'source_endpoint_checks','order_gaps_checked','K','eta','finite_ratio_decimal','ideal_lower_decimal','elapsed_seconds','claim')},indent=2))


if __name__=='__main__':
    main()
