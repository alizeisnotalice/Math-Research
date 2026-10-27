#!/usr/bin/env python3
"""Portable manuscript regression check, no search.
The fixed-axis union kernel is repackaged from the completed independent audit.
Its implementation is reused, not presented as a third independent algorithm.
Only Python standard library is required. Inputs are never rewritten.
"""
import argparse
from collections import defaultdict
from fractions import Fraction as Q
from itertools import product
from math import isqrt, prod, comb
from pathlib import Path
import hashlib, json, time
HERE=Path(__file__).resolve().parent
D=10**35
S=4*D
LIMIT=2**22

def root8(m):
    scaled = m * D**8
    q = isqrt(isqrt(isqrt(scaled)))
    assert q**8 <= scaled < (q + 1)**8
    return q

class LimitReached(RuntimeError):
    pass

class FixedAxisSweep:
    """Each state is the set of physical boxes on the remaining fixed axes.

    Sweep first remaining axis. Event multiplicities maintain the union of
    suffix boxes. Fubini reduces to the next fixed axis; no orbit compression.
    """
    def __init__(self, budget):
        self.memo = {}
        self.budget = budget
        self.max_state_boxes = 0

    def volume(self, state):
        if not state:
            return 0
        if state in self.memo:
            return self.memo[state]
        if self.budget[0] >= LIMIT:
            raise LimitReached("2^22 total memo states exhausted; incomplete")
        self.budget[0] += 1
        self.max_state_boxes = max(self.max_state_boxes, len(state))
        if len(state[0]) == 2:
            total = 0
            left, right = state[0]
            for lo, hi in state[1:]:
                if lo > right:
                    total += right - left
                    left, right = lo, hi
                else:
                    right = max(right, hi)
            total += right - left
        else:
            events = defaultdict(list)
            for box in state:
                events[box[0]].append((box[2:], 1))
                events[box[1]].append((box[2:], -1))
            active = {}
            previous = None
            total = 0
            for x in sorted(events):
                if previous is not None and active:
                    child = tuple(sorted(active))
                    total += (x - previous) * self.volume(child)
                for suffix, delta in events[x]:
                    count = active.get(suffix, 0) + delta
                    assert count >= 0
                    if count:
                        active[suffix] = count
                    else:
                        del active[suffix]
                previous = x
            assert not active
        self.memo[state] = total
        return total

def upper_certificate(claims):
    t=Q(3837029,524288); p=1+1/(2*t); h=(1+1/t)/2
    alpha=1-h**8-4*sum(Q(comb(7,j))*p**(7-j)*Q(-1,2)**j*(1-t**(-j))/j for j in range(1,8))
    beta=8-4*p**7
    z=(t-1)/(t+1); term=z; partial=Q(0); N=400
    for j in range(N):
        partial+=term/(2*j+1); term*=z*z
    lo=2*partial; hi=lo+2*term/((2*N+1)*(1-z*z))
    blo,bhi=sorted((alpha+beta*lo,alpha+beta*hi))
    Bplus=Q(claims['upper_compact_B']); assert 0<blo<=bhi<=Bplus
    base=(1+1/t)**4; U=Q(claims['strict_upper']); d=U-base**2-Bplus
    residual=d*d-4*base*base*Bplus
    assert d>0 and residual>Q(claims['upper_squared_residual_lower'])>0
    return {'status':'PASS','t':str(t),'log_terms':N,'strict_upper':str(U),'positive_squared_residual_lower':claims['upper_squared_residual_lower']}

def finite_formula(L,claims):
    f=claims['finite_L1_original']; theta=Q(f['theta']); K=f['K'];loss=Q(f['loss_upper'])
    value=theta*(L-loss)*Q(2*K+1,2*K+5)**8
    assert value>Q(claims['strict_lower'])
    return value

def lower_certificate(claims):
    masses={prod(e) for e in product((1,2),*(range(1,5) for _ in range(7)))}
    roots={m:root8(m) for m in masses|{1944}}
    spacings=[(2*(D+roots[1944]),2*(D+roots[1944]+1))]+[(2*(D+roots[2]),2*(D+roots[2]+1))]*7
    assert roots[1]==D
    boxes={'inner':[],'outer':[]}; checks=0; widths=[];labels=0
    for e in product((1,2),*(range(1,5) for _ in range(7))):
        labels+=1; m=prod(e);rlo=2*roots[m];rhi=2*(roots[m]+1)
        ins=[];outs=[]; ranges=[]
        for ell,(alo,ahi) in zip(e,spacings):
            q=-(ell//2) if ell%2 else 1-ell//2; p=q+ell-1
            ranges.append(range(q,p+1))
            ins.extend((max(0,p*ahi-rlo),min(alo//2,q*ahi+rlo)))
            outs.extend((max(0,p*alo-rhi),min(ahi//2,q*alo+rhi)))
        for mode,box in [('inner',ins),('outer',outs)]:
            if all(box[2*i]<box[2*i+1] for i in range(8)): boxes[mode].append(tuple(box))
        if any(ins[2*i]>=ins[2*i+1] for i in range(8)):continue
        for i,((alo,ahi),indices) in enumerate(zip(spacings,ranges)):
            left,right=ins[2*i:2*i+2];assert 0<=left<right<=alo//2; widths.append(right-left)
            for index in indices:
                pl,ph=sorted((index*alo,index*ahi))
                assert right-pl<=rlo and ph-left<=rlo;checks+=2
                assert -2<=index<=2 and -2<=-index<=2
    assert labels==32768 and all(len(boxes[k])==882 for k in boxes)
    budget=[0];volumes={};stats={}
    for mode in ('inner','outer'):
        state=tuple(sorted(set(boxes[mode])));assert len(state)==644
        engine=FixedAxisSweep(budget);volume=Q(256*engine.volume(state),S**8)
        volumes[mode]=volume;stats[mode]={'physical_boxes':len(state),'memo_states':len(engine.memo)}
    assert [volumes[k] for k in ('inner','outer')]==list(map(Q,claims['full_union_interval']))
    eps=Q(claims['finite_L1_original']['epsilon']);assert Q(min(widths),S)>2*eps and all(Q(w,S)<2 for w in widths)
    loss=256*882*16*eps*2**7;assert loss==Q(claims['finite_L1_original']['loss_upper'])
    finite=finite_formula(volumes['inner'],claims);assert finite==Q(claims['finite_L1_original']['lower'])
    L0,U0=map(Q,claims['compact_union_interval']);assert L0<=volumes['inner']<=volumes['outer']<=U0
    compact=finite_formula(L0,claims);assert compact==Q(claims['finite_L1_from_compact_lower'])
    return {'status':'PASS','labels':labels,'union_states':budget[0],'coordinate_containment_checks':checks,'statistics':stats,'strict_lower':claims['strict_lower']}

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--scalar-only',action='store_true');args=parser.parse_args()
    started=time.perf_counter();claims=json.loads((HERE/'paper_claims.json').read_text())
    upper=upper_certificate(claims)
    finite_formula(Q(claims['compact_union_interval'][0]),claims)
    lower=({'status':'saved arithmetic only; union not recomputed'} if args.scalar_only else lower_certificate(claims))
    print(json.dumps({'status':'PASS','scope':'manuscript-export verification; no new parameter search','lower':lower,'upper':upper,'wall_seconds':time.perf_counter()-started,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2))
if __name__=='__main__': main()
