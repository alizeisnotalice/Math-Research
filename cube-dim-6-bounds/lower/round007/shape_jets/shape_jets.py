#!/usr/bin/env python3
"""Exact right first variations of the full registered finite witness union."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from pathlib import Path
import hashlib
import importlib.util
import json
import os
import sys
import time

sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
MAIN=ROOT
ARITHMETIC=HERE.parents[1]/'round001_rectangular_certificate.py'
spec=importlib.util.spec_from_file_location('exact_shape_field',ARITHMETIC)
A=importlib.util.module_from_spec(spec)
spec.loader.exec_module(A)
F=A.F
ZERO,ONE=F(0),F(1)
STARTED=time.monotonic()
NODES=0
NODE_CAP=1000000
WALL_CAP=1200


def budget():
    if NODES>NODE_CAP:
        raise RuntimeError('global million union-state cap')
    if time.monotonic()-STARTED>WALL_CAP:
        raise TimeoutError('global 1200-second wall cap')


class J:
    __slots__=('c','d')
    def __init__(self,c=0,d=0):
        self.c=c if isinstance(c,F) else F(c)
        self.d=d if isinstance(d,F) else F(d)
    def __hash__(self):
        return hash((self.c,self.d))
    def __eq__(self,o):
        return self.c==o.c and self.d==o.d
    def __lt__(self,o):
        return self.c<o.c or (self.c==o.c and self.d<o.d)
    def __add__(self,o):
        return J(self.c+o.c,self.d+o.d)
    def __sub__(self,o):
        return J(self.c-o.c,self.d-o.d)
    def __mul__(self,o):
        c=ZERO if self.c==ZERO or o.c==ZERO else self.c*o.c
        d1=ZERO if self.c==ZERO or o.d==ZERO else self.c*o.d
        d2=ZERO if self.d==ZERO or o.c==ZERO else self.d*o.c
        return J(c,d1+d2)


JZERO=J()


def volume(boxes):
    """Slab union on integer endpoint ranks; lengths are exact dual numbers."""
    if not boxes:
        return JZERO
    dimension=len(boxes[0])//2
    points=[sorted({b[2*j+k] for b in boxes for k in (0,1)})
            for j in range(dimension)]
    ranks=[{v:i for i,v in enumerate(p)} for p in points]
    integer_boxes=tuple(sorted({tuple(ranks[j][b[2*j+k]]
                                    for j in range(dimension) for k in (0,1))
                                for b in boxes}))
    @lru_cache(None)
    def width(level,lo,hi):
        return points[level][hi]-points[level][lo]
    @lru_cache(None)
    def rec(B):
        global NODES
        NODES+=1
        budget()
        if not B:
            return JZERO
        level=dimension-len(B[0])//2
        if len(B[0])==2:
            ints=sorted(B)
            lo,hi=ints[0]
            value=JZERO
            for c,d in ints[1:]:
                if c>hi:
                    value=value+width(level,lo,hi)
                    lo,hi=c,d
                elif d>hi:
                    hi=d
            return value+width(level,lo,hi)
        cuts=sorted({x for b in B for x in b[:2]})
        value=JZERO
        for lo,hi in zip(cuts,cuts[1:]):
            child=tuple(sorted({b[2:] for b in B if b[0]<=lo and hi<=b[1]}))
            value=value+width(level,lo,hi)*rec(child)
        return value
    answer=rec(integer_boxes)
    return answer


def self_checks():
    cases=[
        ('separating_touch',[(J(0),J(1)),(J(1,1),J(2))],J(2,-1)),
        ('overlapping_touch',[(J(0),J(1,1)),(J(1),J(2))],J(2,0)),
        ('one_new_thin_direction',[(J(1),J(1,2),J(0),J(3))],J(0,6)),
        ('two_new_thin_directions',[(J(1),J(1,2),J(0),J(0,3))],J(0,0))]
    output=[]
    for name,boxes,expected in cases:
        got=volume(boxes)
        assert got==expected,(name,got.c.a,got.d.a)
        output.append({'name':name,'status':'PASS',
                       'constant':str(got.c.a[0]),'derivative':str(got.d.a[0])})
    return output


def spacings():
    u,v=A.u,A.v
    return (ONE+u**4/2,(u**3+u**4*v**3)/4,*((ONE+u)/2 for _ in range(4)))


def coefficients(x):
    return {f'{i},{j}':str(x.a[6*i+j]) for i in range(6) for j in range(6)
            if x.a[6*i+j]}


def baseline_family(a):
    amin=min(a)
    assert ONE<amin
    excluded=[]
    for ell in range(5,9):
        excess=(amin*(ell-1))**6-F(ell*4**5)
        lo,hi=A.bounds(excess.a)
        assert lo>0
        excluded.append({'active_length':ell,'strict_span_excess_interval':[str(lo),str(hi)]})
    # Permutation putting original active axis 3 first.
    order=(2,0,1,3,4,5)
    b=tuple(a[j] for j in order)
    raw=[]
    rejected_tail=0
    zero_boxes=0
    for ell in product(range(1,5),repeat=6):
        budget()
        mass=1
        for l in ell:
            mass*=l
        radius=A.root(mass)/2
        starts=[None]
        tails=[]
        for aj,l in zip(b[1:],ell[1:]):
            q=-(l//2) if l%2 else 1-l//2
            lo,hi=max(ZERO,(q+l-1)*aj-radius),min(aj/2,q*aj+radius)
            if hi<lo:
                rejected_tail+=1
                break
            starts.append(q)
            tails.extend((lo,hi))
        else:
            for q in range(-4,7):
                lo,hi=max(ZERO,(q+ell[0]-1)*b[0]-radius),min(2*b[0],q*b[0]+radius)
                if hi<lo:
                    continue
                box=(lo,hi,*tails)
                zerodims=sum(box[2*j]==box[2*j+1] for j in range(6))
                zero_boxes+=bool(zerodims)
                raw.append({'ell':ell,'q':q,'starts':(q,*starts[1:]),'radius':radius,
                            'zero_baseline_dimensions':zerodims})
    return raw,order,b,{'raw_boxes_including_baseline_degenerate':len(raw),
                        'baseline_degenerate_boxes':zero_boxes,
                        'active_lengths_5_to_8_strict_exclusion':excluded,
                        'tail_strictly_empty_length_vectors':rejected_tail,
                        'active_start_superset':[-4,6]}


def endpoint(c0=ZERO,c1=ZERO,c2=ZERO,rad=ZERO,kappa=Q(0)):
    """A polynomial plus rad*(1+kappa*t)^(1/6), with canonical germ label."""
    if kappa==0:
        c0=c0+rad
        rad=ZERO
    jet=J(c0+rad,c1+rad*(kappa/6))
    key=(c0.a,c1.a,c2.a,rad.a,str(kappa))
    return jet,frozenset((key,))


def choose(x,y,take_min):
    if x[0]==y[0]:
        return x[0],x[1]|y[1]
    if (x[0]<y[0])==take_min:
        return x
    return y


def direction_boxes(raw,order,b,h_original,D,E):
    h=tuple(h_original[j] for j in order)
    germ_groups=[{} for _ in range(6)]
    boxes=[]
    stats={'closed_to_first_order':0,'newborn_kept':0,'two_or_more_thin_dimensions':0,
           'negative_jet_width_discarded':0}
    for record in raw:
        ell,q,starts,r=record['ell'],record['q'],record['starts'],record['radius']
        kappa=Q(D*(-1 if q%2 else 1),ell[0]) if ell[0]%2 else Q(0)
        coordinates=[]
        invalid=False
        thin=0
        for j,(aj,l,qj,hj) in enumerate(zip(b,ell,starts,h)):
            last=qj+l-1
            low_c0,high_c0=aj*qj,aj*last
            low_c1,high_c1=aj*qj*hj,aj*last*hj
            low_c2=high_c2=ZERO
            if j==0:
                low_c1=low_c1+aj*E*(qj%2)
                high_c1=high_c1+aj*E*(last%2)
                low_c2=aj*hj*E*(qj%2)
                high_c2=aj*hj*E*(last%2)
            lo=choose(endpoint(),endpoint(high_c0,high_c1,high_c2,-r,kappa),False)
            cap=2*aj if j==0 else aj/2
            hi=choose(endpoint(cap,cap*hj),endpoint(low_c0,low_c1,low_c2,r,kappa),True)
            for jet,keys in (lo,hi):
                germ_groups[j].setdefault(jet,set()).update(keys)
            if hi[0]<lo[0]:
                stats['negative_jet_width_discarded']+=1
                invalid=True
                break
            if hi[0]==lo[0]:
                stats['closed_to_first_order']+=1
                invalid=True
                break
            if hi[0].c==lo[0].c:
                thin+=1
            coordinates.extend((lo[0],hi[0]))
        if not invalid:
            stats['newborn_kept']+=bool(thin)
            stats['two_or_more_thin_dimensions']+=thin>=2
            boxes.append(tuple(coordinates))
    stats['kept_raw_boxes']=len(boxes)
    stats['distinct_boxes_after_exact_jet_dedup']=len(set(boxes))
    stats['unresolved_true_endpoint_order_jet_groups']=sum(
        len(keys)>1 for mapping in germ_groups for keys in mapping.values())
    stats['tie_interpretation']='higher-order endpoint order unresolved; derivative uses finite-union Lipschitz lemma'
    return boxes,stats


def directions():
    for j in range(6):
        for sign in (1,-1):
            h=[0]*6
            h[j]=sign
            yield f'spacing_{j+1}_{"plus" if sign>0 else "minus"}',h,0,0
    for sign in (1,-1):
        yield f'mass_{"plus" if sign>0 else "minus"}',[0]*6,sign,0
    for sign in (1,-1):
        yield f'gap_{"plus" if sign>0 else "minus"}',[0]*6,0,sign


def main():
    inputs=[MAIN/'rounds/007/frozen_lower_before_candidate.json',
            MAIN/'coordination/phase_20261025/frozen_thresholds.json',
            MAIN/'lower/round003/second_spacing_breakpoint_certificate_results.json',
            ARITHMETIC,HERE.parents[1]/'round006/astra_construction/report.md',
            HERE/'registration.md',HERE/'theory.md',Path(__file__)]
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    ledger=json.loads(inputs[0].read_text())
    certificate=json.loads(inputs[2].read_text())
    known=[Q(0)]*36
    for ij,c in certificate['coefficients'].items():
        i,j=map(int,ij.split(','))
        known[6*i+j]=Q(c)
    exact_baseline=F(a=known)
    out={'status':'running','pid':os.getpid(),
         'scope':'exact finite-family directional derivatives, no new C6 bound',
         'current_complete_upper':ledger['lower']['exact_interval'][1],
         'input_sha256':hashes,'limits':{'directions':16,'global_union_states':NODE_CAP,
                                        'global_wall_seconds':WALL_CAP},
         'results':[]}
    def save():
        out['elapsed_seconds']=time.monotonic()-STARTED
        out['union_states_global']=NODES
        (HERE/'jet_results.json').write_text(json.dumps(out,indent=2)+'\n')
    save()
    try:
        out['self_checks']=self_checks()
        raw,order,b,family=baseline_family(spacings())
        out['family']=family
        save()
        for name,h,D,E in directions():
            tick=time.monotonic()
            oldnodes=NODES
            boxes,stats=direction_boxes(raw,order,b,h,D,E)
            value=volume(boxes)*J(16)
            assert value.c==exact_baseline,'exact baseline constant coefficient mismatch'
            lo,hi=A.bounds(value.d.a)
            sign=0 if value.d==ZERO else (1 if lo>0 else (-1 if hi<0 else None))
            if sign is None:
                raise ArithmeticError('derivative sign not separated by field intervals')
            result={'name':name,'h':h,'D':D,'E':E,'active_axis_one_based':3,
                    'exact_baseline_coefficient_match':True,
                    'derivative_coefficients':coefficients(value.d),
                    'derivative_interval':[str(lo),str(hi)],
                    'derivative_float':float((lo+hi)/2),'derivative_sign':sign,
                    'geometry':stats,'union_states':NODES-oldnodes,
                    'elapsed_seconds':time.monotonic()-tick}
            out['results'].append(result)
            save()
            print(json.dumps({k:result[k] for k in ('name','derivative_float','derivative_sign',
                                                   'geometry','union_states','elapsed_seconds')}),flush=True)
            if sign>0:
                out['status']='positive_direction_stop_requires_one_finite_step_certificate'
                out['positive_direction']=result
                break
            # Direction caches are disposable; bound memory across 16 cases.
            A.add.cache_clear()
            A.neg.cache_clear()
            A.mul.cache_clear()
            A.bounds.cache_clear()
        else:
            out['status']='completed_no_positive_direction_in_registered_set'
    except (RuntimeError,TimeoutError,ArithmeticError,AssertionError) as exc:
        out['status']='stopped'
        out['stop_reason']=repr(exc)
    out['actual_process_status']='finished; no child process, MATLAB or agent'
    save()
    print(json.dumps({'status':out['status'],'directions_completed':len(out['results']),
                      'union_states':NODES,'elapsed_seconds':out['elapsed_seconds'],
                      'stop_reason':out.get('stop_reason')}),flush=True)


if __name__=='__main__':
    main()
