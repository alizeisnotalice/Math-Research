#!/usr/bin/env python3
"""Round 11: finite parity-dependent position perturbations, not weight search.

Infinite atoms have unit mass and coordinates n_j*a_j in the first two axes;
axis j>=2 is n_j*a_j + a_j*delta[n0%2,n1%2,j-2].  The four deltas in each
shifted axis have zero mean.  The actual fundamental period is
(2*a0,2*a1,a2,...), with full mass 4.  Every rectangular index block uses ALL
its actual atom positions to compute the common centre cube intersection.
Numerical discovery only. No exact lower-bound or global optimality claim.
"""
from __future__ import annotations
import functools
import itertools
import json
import math
from pathlib import Path
import platform
import sys
import time
import numpy as np

ROOT=Path(__file__).resolve().parent
BUDGET=80
LIM=4
BASELINE={3:1.7939332277947704,4:1.9017426111009557}
SEEDS={d:2026102411+d for d in (3,4)}
IMPROVEMENT_TOL=1e-10  # only a numerical candidate trigger, not certification
ENDPOINT_TOTAL_CELL_LIMIT=2_000_000_000
ENDPOINT_SLAB_CELL_LIMIT=32_000_000  # at most ~32 MB bool coverage per slab

def base_spacings(d):
    u=2.**(1./d);v=3.**(1./d)
    return np.array([1.+u**(d-2)/2.,(2.+u**(d-1)*v)/4.]
                    +[(1.+u)/2.]*(d-2))

def split_periodic(lo,hi,period):
    length=hi-lo
    if length<=0:return []
    if length>=period:return [(0.,float(period))]
    lo%=period
    if lo+length<=period:return [(float(lo),float(lo+length))]
    return [(float(lo),float(period)),(0.,float(lo+length-period))]

def prune(boxes):
    unique=list(dict.fromkeys(boxes))
    if not unique:return []
    arr=np.asarray(unique);keep=[]
    for i,b in enumerate(arr):
        contains=(np.all(arr[:,::2]<=b[::2],axis=1)
                  & np.all(arr[:,1::2]>=b[1::2],axis=1))
        contains[i]=False
        if not np.any(contains):keep.append(tuple(float(x) for x in b))
    return keep

def union_volume(boxes):
    """Recursive sweep used by discovery only."""
    @functools.lru_cache(maxsize=None)
    def rec(data):
        if not data:return 0.
        if len(data[0])==2:
            intervals=sorted(data);left,right=intervals[0];total=0.
            for a,b in intervals[1:]:
                if a>right:total+=right-left;left,right=a,b
                else:right=max(right,b)
            return total+right-left
        cuts=sorted({x for b in data for x in b[:2]});total=0.
        for left,right in zip(cuts,cuts[1:]):
            mid=(left+right)/2.
            rest=tuple(sorted(set(b[2:] for b in data if b[0]<mid<b[1])))
            total+=(right-left)*rec(rest)
        return total
    return rec(tuple(sorted(boxes)))

@functools.lru_cache(maxsize=None)
def block_indices(d):
    """Reusable integer blocks; no spacing-based extrema are precomputed."""
    data=[]
    for starts01 in itertools.product(range(2),repeat=2):
        starts=starts01+(0,)*(d-2)
        for lengths in itertools.product(range(1,LIM+1),repeat=d):
            points=np.array(list(itertools.product(*(range(s,s+ell) for s,ell in zip(starts,lengths)))),dtype=int)
            parity=2*(points[:,0]%2)+(points[:,1]%2)
            mass=math.prod(lengths)
            data.append((points,parity,mass**(1./d)))
    return tuple(data)

def boxes_for(spacings,deltas):
    d=len(spacings);periods=np.r_[2*spacings[:2],spacings[2:]];raw=[]
    for indices,parity,side in block_indices(d):
        # This is the essential difference from the old unperturbed lattice:
        # compute every real point, then min/max, including parity shifts.
        points=indices*spacings
        points[:,2:]+=deltas[parity,:]*spacings[2:]
        minimum=points.min(axis=0);maximum=points.max(axis=0)
        low=maximum-side/2.;high=minimum+side/2.
        if np.any(high<=low):continue
        pieces=[split_periodic(float(a),float(b),float(p))
                for a,b,p in zip(low,high,periods)]
        for product in itertools.product(*pieces):
            raw.append(tuple(x for pair in product for x in pair))
    return prune(raw)

def score(spacings,deltas):
    boxes=boxes_for(spacings,deltas)
    return union_volume(boxes)/4.,boxes

def endpoint_cell_check(boxes):
    """Independent cell marking, slabbed in one axis to bound memory.

    No recursive union/sweep call is used. Every elementary endpoint cell is
    either marked in a boolean slab or not, then integrated by cell widths.
    """
    if not boxes:return {'status':'checked','volume':0.,'score':0.}
    d=len(boxes[0])//2
    edges=[sorted({b[2*j+q] for b in boxes for q in (0,1)}) for j in range(d)]
    shape=tuple(len(e)-1 for e in edges)
    total_cells=math.prod(shape)
    # Slab the longest axis, minimizing the coverage array footprint.
    slab_axis=int(np.argmax(shape));others=[j for j in range(d) if j!=slab_axis]
    slab_shape=tuple(shape[j] for j in others);slab_cells=math.prod(slab_shape)
    meta={'cell_shape':list(shape),'total_cells':total_cells,'slab_axis':slab_axis,
          'peak_boolean_cells':slab_cells,'total_cell_limit':ENDPOINT_TOTAL_CELL_LIMIT,
          'slab_cell_limit':ENDPOINT_SLAB_CELL_LIMIT}
    if total_cells>ENDPOINT_TOTAL_CELL_LIMIT or slab_cells>ENDPOINT_SLAB_CELL_LIMIT:
        return {'status':'not_checked_resource_limit',**meta}
    index=[{x:i for i,x in enumerate(e)} for e in edges]
    widths=[np.diff(np.asarray(e)) for e in edges]
    indexed=[tuple((index[j][b[2*j]],index[j][b[2*j+1]]) for j in range(d)) for b in boxes]
    volume=0.
    for cut,width in enumerate(widths[slab_axis]):
        coverage=np.zeros(slab_shape,dtype=bool)
        for ib in indexed:
            if ib[slab_axis][0]<=cut<ib[slab_axis][1]:
                coverage[tuple(slice(*ib[j]) for j in others)]=True
        weighted=coverage
        # Contract the final axis repeatedly; only the bool slab is 3D in d=4.
        for j in reversed(others):
            weighted=np.tensordot(weighted,widths[j],axes=([-1],[0]))
        volume+=float(width)*float(weighted)
    return {'status':'checked','volume':volume,'score':volume/4.,**meta}

def normalized(raw):
    return raw-np.mean(raw,axis=0,keepdims=True)

def schedule(d):
    out=[('zero_control',np.zeros((4,d-2)))]
    amplitudes=(-.1,-.04,-.01,.01,.04,.1)
    for q in range(4):
        for t in amplitudes:
            raw=np.zeros((4,d-2));raw[q,0]=t
            out.append((f'axis3_single_parity_{q}_{t:+.2f}',normalized(raw)))
    modes=[np.array([1.,1.,-1.,-1.]),np.array([1.,-1.,1.,-1.]),np.array([1.,-1.,-1.,1.])]
    names=['row','column','checkerboard']
    for name,mode in zip(names,modes):
        for t in amplitudes:
            raw=np.zeros((4,d-2));raw[:,0]=t*mode
            out.append((f'axis3_{name}_{t:+.2f}',normalized(raw)))
    if d==4:
        for k in range(3):
            for t in amplitudes:
                raw=np.column_stack([t*modes[k],t*modes[(k+1)%3]])
                out.append((f'axes3and4_coupled_modes_{k}_{t:+.2f}',normalized(raw)))
    rng=np.random.default_rng(SEEDS[d]);random_index=0
    while len(out)<BUDGET:
        scale=(.01,.04,.1)[random_index%3]
        raw=rng.uniform(-scale,scale,size=(4,d-2))
        out.append((f'joint_position_uniform_scale_{scale:.2f}',normalized(raw)))
        random_index+=1
    assert len(out)==BUDGET
    return out

def item_for(d,number,name,deltas):
    spacings=base_spacings(d);value,boxes=score(spacings,deltas)
    return {'evaluation':number,'source':name,'score':value,'score_minus_reference':value-BASELINE[d],
            'spacings':spacings.tolist(),'delta_in_spacing_units':deltas.tolist(),
            'physical_delta':(deltas*spacings[2:]).tolist(),
            'box_count':len(boxes)},boxes

def save(data):
    (ROOT/'search_result.json').write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n')

def verify_item(item):
    boxes=boxes_for(np.array(item['spacings']),np.array(item['delta_in_spacing_units']))
    check=endpoint_cell_check(boxes)
    if check['status']=='checked':
        check['difference_from_recursive_score']=check['score']-item['score']
        check['agreement_tolerance']=1e-9
        check['agrees']=abs(check['difference_from_recursive_score'])<=1e-9
    return check

def main():
    if (ROOT/'search_result.json').exists():
        raise SystemExit('Recorded experiment already exists; refusing duplicate budget use.')
    started=time.monotonic()
    result={'round':11,'status':'running_numerical_discovery_only',
            'mechanism':'nonproduct parity-dependent periodic positions; masses=1; spacings fixed baseline',
            'period':'2x2x1^(d-2)','full_period_mass':4,
            'witness_index_lengths_each_coordinate':[1,LIM],
            'start_parities':[[0,0],[0,1],[1,0],[1,1]],
            'common_translation_removed':'columnwise arithmetic mean of four displacements is zero',
            'evaluation_budget_per_dimension':BUDGET,'improvement_trigger_tolerance':IMPROVEMENT_TOL,
            'environment':{'python':sys.version,'platform':platform.platform(),'numpy':np.__version__},
            'dimensions':[]}
    schedules={d:schedule(d) for d in (3,4)}
    # First reproduce both controls; these count as evaluation 1 in each budget.
    for d in (3,4):
        item,boxes=item_for(d,1,*schedules[d][0])
        dimension={'dimension':d,'seed':SEEDS[d],'baseline_reference':BASELINE[d],
                   'baseline':item,'evaluations_completed':1,'candidates':[item]}
        result['dimensions'].append(dimension)
        print(f'd={d} control={item["score"]:.16f} difference={item["score_minus_reference"]:+.3g} boxes={len(boxes)}',flush=True)
        save(result)
        if abs(item['score_minus_reference'])>1e-10:
            result['status']='control_failed';save(result);raise SystemExit('Zero control failed baseline reproduction.')
        dimension['baseline_endpoint_check']=verify_item(item)
        print(f'd={d} baseline endpoint check {dimension["baseline_endpoint_check"]}',flush=True)
        save(result)
    stopped=False
    for dimension in result['dimensions']:
        d=dimension['dimension'];best=None
        for number,(name,deltas) in enumerate(schedules[d][1:],2):
            item,boxes=item_for(d,number,name,deltas)
            dimension['candidates'].append(item);dimension['evaluations_completed']=number
            if best is None or item['score']>best['score']:
                best=item;dimension['best_noncontrol']=best
                print(f'd={d} eval={number} best_noncontrol={item["score"]:.16f} source={name} boxes={len(boxes)}',flush=True)
            if number%10==0:save(result)
            if item['score']>BASELINE[d]+IMPROVEMENT_TOL:
                dimension['improvement_candidate']=item
                result['status']='new_numerical_candidate_stop_expansion'
                stopped=True;save(result)
                print('NEW_CANDIDATE: expansion stopped immediately.',flush=True)
                break
        if best is not None:
            dimension['best_noncontrol_endpoint_check']=verify_item(best)
            dimension['best_noncontrol_minus_baseline']=best['score']-dimension['baseline']['score']
            dimension['new_candidate_triggered']=best['score']>BASELINE[d]+IMPROVEMENT_TOL
        save(result)
        if stopped:break
    if not stopped:result['status']='finite_experiment_no_improvement_above_trigger'
    result['elapsed_seconds']=time.monotonic()-started
    save(result)
    summary={k:v for k,v in result.items() if k!='dimensions'}
    summary['dimensions']=[{k:v for k,v in dim.items() if k!='candidates'} for dim in result['dimensions']]
    print(json.dumps(summary,indent=2,ensure_ascii=False),flush=True)
if __name__=='__main__':main()
