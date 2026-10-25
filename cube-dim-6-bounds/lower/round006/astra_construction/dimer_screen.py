#!/usr/bin/env python3
"""Registered joint gap/mass screen. Floating discovery, never a C6 proof."""
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import ceil, floor, prod
from pathlib import Path
import hashlib
import json
import os
import time

HERE = Path(__file__).resolve().parent
MAIN = Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds')
OLD = Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-lower-clean/cube-dim-6-bounds/lower/round006')
RAW_CAP, NODE_CAP, WALL_CAP = 5000, 1000000, 1200
NODE_COUNT = 0


def base_spacing():
    u, v = 2**(1/6), 3**(1/6)
    return (1+u**4/2, (u**3+u**4*v**3)/4, *((1+u)/2 for _ in range(4)))


def budget(started):
    if time.monotonic()-started > WALL_CAP:
        raise TimeoutError('global 1200-second wall cap')
    if NODE_COUNT > NODE_CAP:
        raise RuntimeError('global million union-state cap')


def make_boxes(axis, e, d, started):
    order = (axis,)+tuple(j for j in range(6) if j != axis)
    a = tuple(base_spacing()[j] for j in order)
    scale = (1+e)**(1/6)
    b = (a[0]*scale/(1+e),)+tuple(x*scale for x in a[1:])
    def position(k):
        return b[0]*(k+e*(k%2))
    raw = []
    for ell in product(range(1,9), *[range(1,5) for _ in range(5)]):
        budget(started)
        tail_mass = prod(ell[1:])
        max_radius = (tail_mass*(ell[0]+d))**(1/6)/2
        # p(q) <= b*(q+e), p(q+l-1) >= b*(q+l-1).
        # Necessary meeting conditions therefore put q inside this superset.
        qlo = floor(-max_radius/b[0]-e)-1
        qhi = ceil(2+max_radius/b[0]-(ell[0]-1))+1
        for q in range(qlo,qhi+1):
            active_mass = ell[0]+(d*(-1 if q%2 else 1) if ell[0]%2 else 0)
            mass = tail_mass*active_mass
            radius = mass**(1/6)/2
            low, high = position(q), position(q+ell[0]-1)
            coords = [max(0., high-radius), min(2*b[0],low+radius)]
            if coords[1] <= coords[0]:
                continue
            source_starts = [q]
            source_endpoints = [(low,high)]
            for aj,l in zip(b[1:],ell[1:]):
                start = -(l//2) if l%2 else 1-l//2
                lo,hi = start*aj,(start+l-1)*aj
                c,h = max(0.,hi-radius), min(aj/2,lo+radius)
                if h <= c:
                    break
                coords.extend((c,h))
                source_starts.append(start)
                source_endpoints.append((lo,hi))
            else:
                for j,(lo,hi) in enumerate(source_endpoints):
                    for x in coords[2*j:2*j+2]:
                        if max(abs(x-lo),abs(x-hi)) > radius+4e-12:
                            raise AssertionError('endpoint/common-radius containment failed')
                raw.append({'lengths_reordered':ell,'starts_reordered':source_starts,
                            'mass_float':mass,'radius_float':radius,'box':tuple(coords)})
                if len(raw) > RAW_CAP:
                    raise RuntimeError('5000 raw-box cap')
    return raw,order,b


def prune(raw):
    unique = {tuple(x['box']):x for x in raw}
    # No tolerance-assisted containment: a smaller retained union would
    # still be a discovery lower witness, but avoid masking narrow features.
    return [item for b,item in unique.items()
            if not any(c != b and all(c[2*j] <= b[2*j] and b[2*j+1] <= c[2*j+1]
                                     for j in range(6)) for c in unique)]


def union_volume(boxes, started):
    @lru_cache(None)
    def rec(X):
        global NODE_COUNT
        NODE_COUNT += 1
        budget(started)
        if not X:
            return 0.
        if len(X[0]) == 2:
            ints=sorted(X)
            lo,hi=ints[0]
            value=0.
            for c,d in ints[1:]:
                if c > hi:
                    value += hi-lo
                    lo,hi=c,d
                elif d > hi:
                    hi=d
            return value+hi-lo
        cuts=sorted({x for b in X for x in b[:2]})
        value=0.
        for lo,hi in zip(cuts,cuts[1:]):
            child=tuple(sorted({b[2:] for b in X if b[0] <= lo and hi <= b[1]}))
            value += (hi-lo)*rec(child)
        return value
    return rec(tuple(sorted(tuple(x['box']) for x in boxes)))


def evaluate(axis,e,d,started):
    tick=time.monotonic()
    initial_nodes=NODE_COUNT
    raw,order,b=make_boxes(axis,float(e),float(d),started)
    kept=prune(raw)
    volume=union_volume(kept,started)
    result={'axis_one_based':axis+1,'epsilon':str(e),'delta':str(d),
            'coordinate_order_zero_based':order,'spacings_float':b,
            'raw_boxes':len(raw),'maximal_boxes':len(kept),'half_cell_volume_float':volume,
            'ratio_float':16*volume,'union_states':NODE_COUNT-initial_nodes,
            'elapsed_seconds':time.monotonic()-tick}
    return result,raw,kept


def main():
    started=time.monotonic()
    inputs=[MAIN/'results/bounds_ledger.json',
            MAIN/'coordination/phase_20261025/frozen_thresholds.json',
            MAIN/'coordination/phase_20261025/protocol.md',
            MAIN/'coordination/phase_20261025/roles.json',
            MAIN/'lower/round005/strict_extension_certificate.md',
            OLD/'correlated/report.md',OLD/'staggered/report.md',
            HERE/'registration.md',Path(__file__)]
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    ledger=json.loads(inputs[0].read_bytes())
    frozen=json.loads(inputs[1].read_bytes())
    threshold=Fraction(ledger['lower']['exact_interval'][1])
    frozen_threshold=Fraction(frozen['lower']['exact_interval'][1])
    cases=[(axis,Fraction(e),Fraction(d)) for axis in (0,1,2)
           for e in ('1/8','1/3','2/3','7/8') for d in ('0','1/5')]
    assert len(cases)==24
    out={'status':'running','scope':'numeric discovery only, no C6 claim',
         'pid':os.getpid(),'input_sha256':hashes,
         'current_complete_upper':str(threshold),'frozen_complete_upper':str(frozen_threshold),
         'normalization':'period mass=2, full reflected union=32|U_H|, ratio=16|U_H|',
         'limits':{'candidates':24,'raw_boxes_each':RAW_CAP,'union_states_global':NODE_CAP,
                   'wall_seconds_global':WALL_CAP},
         'cases':[{'axis':a+1,'epsilon':str(e),'delta':str(d)} for a,e,d in cases],
         'results':[]}
    def save():
        out['elapsed_seconds']=time.monotonic()-started
        out['union_states_global']=NODE_COUNT
        (HERE/'dimer_results.json').write_text(json.dumps(out,indent=2)+'\n')
    save()
    try:
        cal,_,_=evaluate(0,Fraction(0),Fraction(0),started)
        cal['difference_from_current_upper_float']=cal['ratio_float']-float(threshold)
        out['calibration']=cal
        save()
        if abs(cal['difference_from_current_upper_float'])>2e-9:
            raise RuntimeError('calibration mismatch')
        for axis,e,d in cases:
            result,raw,kept=evaluate(axis,e,d,started)
            result['margin_over_complete_upper_float']=result['ratio_float']-float(threshold)
            result['promotion_trigger']=result['margin_over_complete_upper_float']>1e-9
            out['results'].append(result)
            print(json.dumps(result),flush=True)
            save()
            if result['promotion_trigger']:
                (HERE/'promoted_candidate_boxes.json').write_text(json.dumps(
                    {'candidate':result,'raw':raw,'maximal':kept},indent=2)+'\n')
                out['status']='promotion_stop_requires_independent_exact_certification'
                break
        else:
            out['status']='completed_no_promotion'
    except (RuntimeError,TimeoutError) as exc:
        out['status']='stop'
        out['stop_reason']=str(exc)
    out['best']=max(out['results'],key=lambda x:x['ratio_float'],default=None)
    out['actual_process_status']='finished; no child process; no MATLAB'
    save()
    print(json.dumps({'status':out['status'],'best':out['best'],
                      'elapsed_seconds':out['elapsed_seconds'],
                      'union_states':NODE_COUNT}),flush=True)


if __name__=='__main__':
    main()
