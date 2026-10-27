#!/usr/bin/env python3
"""Independent exact audit. Does not import or execute candidate mathematics.

Unlike candidate orbit recursion, this retains every positive labelled raw box,
sweeps one fixed coordinate, and strips only that coordinate. Integer Newton
roots and Fraction interval arithmetic independently construct all endpoints.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
import json, hashlib, math, time, sys, platform

HERE=Path(__file__).resolve().parent
SRC=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube7-lower/cube-dim-7-bounds/lower/round001')
BASE=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds/lower/round075/post_s8_overlap/results.json')
D=10**35
S=4*D
LIMIT=2**20

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def root(n):
    x=1 << ((n.bit_length()+6)//7)
    while True:
        y=(6*x+n//x**6)//7
        if y>=x: break
        x=y
    while x**7>n: x-=1
    while (x+1)**7<=n: x+=1
    assert x**7<=n<(x+1)**7
    return x

def add(a,b): return (a[0]+b[0],a[1]+b[1])
def mul(a,c): return tuple(sorted((a[0]*c,a[1]*c)))
def maximum(a,b): return (max(a[0],b[0]),max(a[1],b[1]))
def minimum(a,b): return (min(a[0],b[0]),min(a[1],b[1]))
def units(x):
    y=x*S
    assert y.denominator==1
    return y.numerator

class RawUnion:
    def __init__(self,budget): self.memo={}; self.budget=budget
    def volume(self,boxes):
        if not boxes: return 0
        if boxes in self.memo: return self.memo[boxes]
        self.budget[0]+=1
        assert self.budget[0]<=LIMIT, 'registered state limit'
        dims=len(boxes[0])
        assert all(len(b)==dims for b in boxes)
        if dims==1:
            total=0; right=None
            for ((a,b),) in sorted(boxes):
                if right is None or a>right: total+=b-a; right=b
                elif b>right: total+=b-right; right=b
        else:
            events={}
            for box in boxes:
                a,b=box[0]
                events.setdefault(a,[set(),set()])[0].add(box)
                events.setdefault(b,[set(),set()])[1].add(box)
            active=set(); last=None; total=0
            for x,(starts,ends) in sorted(events.items()):
                if last is not None and active:
                    section=tuple(sorted({b[1:] for b in active}))
                    total+=(x-last)*self.volume(section)
                active.difference_update(ends)
                active.update(starts)
                last=x
            assert not active
        self.memo[boxes]=total
        return total

def construct(power,claimed):
    lengths=list(product((1,2),*[(1,2,3,4)]*6))
    masses={math.prod(e) for e in lengths}
    roots={m:root(m*D**7) for m in masses}
    assert {str(m):str(v) for m,v in roots.items()}==claimed['root_lower_numerators']
    uu=(Q(roots[2],D),Q(roots[2]+1,D))
    p=root(2**power*D**7)
    assert str(p)==claimed['first_root_lower_numerator']
    ap=(1+Q(p,2*D),1+Q(p+1,2*D))
    at=mul(add((1,1),uu),Q(1,2))
    inward=[]; outward=[]; keys=[]; checks=0; minwidth=None
    for ell in lengths:
        mass=math.prod(ell)
        radius=(Q(roots[mass],2*D),Q(roots[mass]+1,2*D))
        ib=[]; ob=[]
        for j,count in enumerate(ell):
            spacing=ap if j==0 else at
            # Independent description: nearest block with right endpoint nonnegative.
            start=-(count-1)//2 if count%2 else 1-count//2
            end=start+count-1
            left=maximum((Q(0),Q(0)),add(mul(spacing,end),mul(radius,-1)))
            right=minimum(mul(spacing,Q(1,2)),add(mul(spacing,start),radius))
            ib.append((units(left[1]),units(right[0])))
            ob.append((units(left[0]),units(right[1])))
        if all(a<b for a,b in ob):
            assert all(a<b for a,b in ib), 'positive-box ambiguity'
            for j,count in enumerate(ell):
                spacing=ap if j==0 else at
                start=-(count-1)//2 if count%2 else 1-count//2
                end=start+count-1
                lo,hi=map(lambda x:Q(x,S),ib[j])
                # For every source index, worst possible position/radius enclosure.
                for index in range(start,end+1):
                    atom=mul(spacing,index)
                    assert lo>=atom[1]-radius[0]
                    assert hi<=atom[0]+radius[0]
                    checks+=2
                w=hi-lo
                minwidth=w if minwidth is None else min(minwidth,w)
            inward.append(tuple(ib)); outward.append(tuple(ob)); keys.append(ell)
    assert len(keys)==claimed['positive_raw_boxes']
    orbits={ (e[0],tuple(sorted(e[1:]))) for e in keys }
    assert len(orbits)==claimed['symmetric_box_orbits']
    # Explicit closure of all raw boxes under the adjacent transpositions generating S6.
    ik=set(inward); ok=set(outward)
    for boxes,known in ((inward,ik),(outward,ok)):
        for b in boxes:
            for j in range(1,6):
                c=list(b); c[j],c[j+1]=c[j+1],c[j]
                assert tuple(c) in known
    return tuple(sorted(ik)),tuple(sorted(ok)),checks,minwidth,max(ap[1],at[1]),len(keys)

def main():
    lock=json.loads((HERE/'input_lock.json').read_text())
    for name,h in lock['candidate_files'].items(): assert sha(SRC/name)==h, f'input changed: {name}'
    assert sha(BASE)==lock['baseline_sha256']=='dec62274abcf8204cf698cd7665dd6c754db3fee87a55beb5bfa1ece676e2817'
    ref=json.loads((SRC/'frozen_internal_reference.json').read_text())
    base=json.loads(BASE.read_text())['candidate']['exact_interval']
    assert list(map(Q,ref['interval']))==list(map(Q,base))
    reference_upper=Q(base[1])
    result={'input_lock_sha256':sha(HERE/'input_lock.json'),'check_sha256':sha(__file__),
            'algorithm':'full labelled raw boxes; fixed coordinate event sweep; no orbit quotient',
            'baseline_interval_verified_against_full_source':True,'seeds':[]}
    for power in (4,5):
        t=time.time(); cpu=time.process_time()
        claimed=json.loads((SRC/f'seed{power}_results.json').read_text())
        inner,outer,checks,minwidth,maxspacing,raw_count=construct(power,claimed)
        budget=[0]; values=[]; memo_counts=[]
        for direction,boxes in (('inner',inner),('outer',outer)):
            engine=RawUnion(budget)
            volume=Q(128*engine.volume(boxes),S**7)
            assert volume==Q(claimed['union_lower' if direction=='inner' else 'union_upper']), direction
            values.append(volume); memo_counts.append(len(engine.memo))
            print(json.dumps({'power':power,'direction':direction,'value':float(volume),'states_total':budget[0],'seconds':time.time()-t}),flush=True)
        assert values[0]<=values[1]
        eps=Q(1,10**25); eta=eps/2; theta=1-Q(1,10**18); K=10**19
        assert minwidth>2*eps and maxspacing/2<2
        # Candidate deliberately pays every raw box, including duplicate boxes.
        loss=128*raw_count*14*eps*2**6
        finite=theta*(values[0]-loss)*Q(2*K+1,2*K+5)**7
        assert loss==Q(claimed['finite_L1']['union_shrink_loss_upper'])
        assert finite==Q(claimed['finite_L1']['finite_ratio_lower'])
        assert finite>reference_upper
        if power==5: assert finite>Q(2135677,1000000)
        result['seeds'].append({'power':power,'raw_positive_boxes':raw_count,
            'distinct_positive_boxes':len(inner),
            'all_source_atom_containment_checks':checks,'memo_states_each_direction':memo_counts,
            'memo_states_both_directions':budget[0],'union_lower':str(values[0]),'union_upper':str(values[1]),
            'finite_L1_ratio_lower':str(finite),'finite_L1_above_complete_R075_upper':True,
            'strict_margin':str(finite-reference_upper),'wall_seconds':time.time()-t,'cpu_seconds':time.process_time()-cpu})
    assert Q(result['seeds'][1]['union_lower'])>Q(result['seeds'][0]['union_upper'])
    result['seed5_finite_ratio_strictly_above_2135677_over_1000000']=True
    result['seed5_ideal_lower_strictly_above_seed4_ideal_upper']=True
    for name,h in lock['candidate_files'].items(): assert sha(SRC/name)==h, f'input changed during run: {name}'
    assert sha(BASE)==lock['baseline_sha256'], 'baseline changed during run'
    result['all_locked_inputs_unchanged_after_run']=True
    result['status']='passed exact independent finite checker'
    result['runtime']={'executable':sys.executable,'python':sys.version,'platform':platform.platform()}
    (HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__': main()
