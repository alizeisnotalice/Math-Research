#!/usr/bin/env python3
"""Registered rational 2D/3D toys only; no C geometry or production imports."""
import time
START = time.monotonic()
from fractions import Fraction as Q
from itertools import product, combinations
from pathlib import Path
import hashlib, json, os, sys

HERE = Path(__file__).resolve().parent
TMAX = Q(1,16)
VOLUME_CALLS = 0


def interval(lo, hi, vlo=0, vhi=0):
    return ((Q(lo),Q(vlo)),(Q(hi),Q(vhi)))


FIXTURES = [
    ('moving_rectangle', [[interval(1,2,-1,2),interval(1,2,-2,1)]]),
    ('shared_face_gap', [[interval(0,1,0,-1),interval(0,2,0,1)],
                         [interval(1,2,1,0),interval(0,2,0,-1)]]),
    ('duplicate_faces', [[interval(1,2,-1,0),interval(1,2,0,1)],
                         [interval(1,2,0,1),interval(1,2,-1,0)]]),
    ('single_zero_birth', [[interval(1,1,-1,2),interval(1,2,0,1)]]),
    ('double_zero_birth', [[interval(1,1,-1,1),interval(2,2,-1,2)]]),
    ('clipped_boundary_birth', [[interval(0,0,-2,1),interval(3,4,-1,2)]]),
    ('triple_zero_birth', [[interval(1,1,-1,1)]*3]),
    ('three_dimensional_shared_faces', [
        [interval(0,1,0,1),interval(0,2,0,-1),interval(0,2,0,1)],
        [interval(1,2,-1,0),interval(0,2,0,1),interval(0,2,0,-1)],
        [interval(1,1,-2,2),interval(2,3,-1,1),interval(1,2,0,1)]]),
    ('covered_birth', [[interval(0,3),interval(0,3)],
                        [interval(1,1,-2,2),interval(1,2,-1,1)]]),
    ('strictly_inverted_box', [[interval(2,1,-1,1),interval(1,1,-1,1),interval(1,1,-1,1)],
                               [interval(0,1,0,1)]*3]),
]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def evaluated(boxes, t, active=None):
    # active=None changes all coordinates; integer active changes only that axis.
    return [tuple((max(Q(0),lo+vlo*(t if active is None or active==k else 0)),
                   min(Q(4),hi+vhi*(t if active is None or active==k else 0)))
                  for k,((lo,vlo),(hi,vhi)) in enumerate(box)) for box in boxes]


def cells(cuts):
    for pairs in product(*[list(zip(xs,xs[1:])) for xs in cuts]):
        width=Q(1)
        for a,b in pairs: width*=b-a
        yield tuple((a+b)/2 for a,b in pairs),width


def volume(boxes):
    global VOLUME_CALLS
    VOLUME_CALLS+=1
    d=len(boxes[0])
    good=[b for b in boxes if all(a<c for a,c in b)]
    if not good:return Q(0)
    cuts=[sorted({Q(0),Q(4)}|{z for b in good for z in b[k]}) for k in range(d)]
    return sum((w for x,w in cells(cuts)
                if any(all(a<xk<c for xk,(a,c) in zip(x,b)) for b in good)),Q(0))


def poly(xs,ys):
    # Exact Lagrange coefficients, with a separate topology certificate below.
    ans=[Q(0)]*len(xs)
    for i,(xi,yi) in enumerate(zip(xs,ys)):
        p=[Q(1)];den=Q(1)
        for j,xj in enumerate(xs):
            if i==j:continue
            out=[Q(0)]*(len(p)+1)
            for k,v in enumerate(p):out[k]-=xj*v;out[k+1]+=v
            p=out;den*=xi-xj
        for k,v in enumerate(p):ans[k]+=yi*v/den
    return ans


def clipped_jets(box):
    return tuple((max((Q(0),Q(0)),lo),min((Q(4),Q(0)),hi)) for lo,hi in box)


def local_line_derivative(intervals):
    """Integral of the local velocity-line union minus baseline half-lines."""
    intervals=[(lo,hi) for lo,hi in intervals if lo[0]<=hi[0]]
    ans=Q(0)
    for c in sorted({p[0] for row in intervals for p in row}):
        if any(lo[0]<c<hi[0] for lo,hi in intervals):continue
        left=[hi[1] for lo,hi in intervals if lo[0]<c==hi[0]]
        right=[lo[1] for lo,hi in intervals if lo[0]==c<hi[0]]
        born=[(lo[1],hi[1]) for lo,hi in intervals if lo[0]==c==hi[0]]
        cuts=sorted({Q(0)}|set(left)|set(right)|{z for row in born for z in row})
        for a,b in zip(cuts,cuts[1:]):
            z=(a+b)/2
            now=(any(z<v for v in left) or any(z>v for v in right)
                 or any(v<z<w for v,w in born))
            old=(bool(left) and z<0) or (bool(right) and z>0)
            ans+=(b-a)*(int(now)-int(old))
    return ans


def static_face_integral(boxes,axis):
    jets=[clipped_jets(b) for b in boxes]
    other=[k for k in range(len(boxes[0])) if k!=axis]
    cuts=[sorted({Q(0),Q(4)}|{max(Q(0),min(Q(4),p[0]))
                            for b in jets for p in b[k]}) for k in other]
    ans=Q(0)
    for x,w in cells(cuts):
        active=[b[axis] for b in jets if all(b[k][0][0]<xj<b[k][1][0]
                                            for k,xj in zip(other,x))]
        ans+=w*local_line_derivative(active)
    return ans


def changeset_length(before,after,axis):
    cuts=sorted({Q(0),Q(4)}|{max(Q(0),min(Q(4),z))
                          for box in before+after for z in box[axis]})
    ans=Q(0)
    for a,b in zip(cuts,cuts[1:]):
        x=(a+b)/2
        if any((p[axis][0]<x<p[axis][1]) != (q[axis][0]<x<q[axis][1])
               for p,q in zip(before,after)):ans+=b-a
    return ans


def main():
    reg=json.loads((HERE/'registration.json').read_text())
    assert reg['source_sha256']==sha(__file__)
    assert not (HERE/'results.json').exists()
    out={'pid':os.getpid(),'python':sys.executable,'status':'running',
         'source_sha256':sha(__file__),'registration_sha256':sha(HERE/'registration.json')}
    print(json.dumps(out),flush=True)
    try:
        for item in reg['readonly_inputs']:assert sha(item['path'])==item['sha256']
        rows=[];mixed_checks=0;crossing_pairs=0
        for name,boxes in FIXTURES:
            d=len(boxes[0]);xs=[Q(j,64) for j in range(d+1)]
            # No positive endpoint-order or clipping event before TMAX.
            for k in range(d):
                endpoints=[(Q(0),Q(0)),(Q(4),Q(0))]+[p for b in boxes for p in b[k]]
                for (a,v),(b,w) in combinations(endpoints,2):
                    crossing_pairs+=1
                    if v!=w:assert not (0<(b-a)/(v-w)<=TMAX)
            vals={};coeffs={}
            for axis in [None]+list(range(d)):
                vals[axis]=[volume(evaluated(boxes,t,axis)) for t in xs+[TMAX]]
                coeffs[axis]=poly(xs,vals[axis][:-1])
                assert sum(c*TMAX**j for j,c in enumerate(coeffs[axis]))==vals[axis][-1]
            assert coeffs[None][1]==sum(coeffs[k][1] for k in range(d))
            faces=[static_face_integral(boxes,k) for k in range(d)]
            assert faces==[coeffs[k][1] for k in range(d)]
            for j,t in enumerate(xs+[TMAX]):
                if not t:continue
                before=evaluated(boxes,0);after=evaluated(boxes,t)
                sizes=[changeset_length(before,after,k) for k in range(d)]
                remainder=vals[None][j]-sum(vals[k][j] for k in range(d))+(d-1)*vals[None][0]
                bound=d*sum(sizes[a]*sizes[b]*4**(d-2) for a,b in combinations(range(d),2))
                assert abs(remainder)<=bound
                mixed_checks+=1
            rows.append({'name':name,'dimension':d,'boxes':len(boxes),
                         'volume_polynomial':[str(c) for c in coeffs[None]],
                         'full_first_derivative':str(coeffs[None][1]),
                         'coordinate_first_derivatives':[str(x) for x in faces]})
        assert len(rows)==10 and VOLUME_CALLS==144 and mixed_checks==33
        for item in reg['readonly_inputs']:assert sha(item['path'])==item['sha256']
        out.update(status='qualified_toy_identity_and_face_formula_only',fixtures=rows,
                   exact_toy_volume_calls=VOLUME_CALLS,mixed_bound_checks=mixed_checks,
                   primitive_crossing_pairs_checked=crossing_pairs,
                   max_dimension=3,max_boxes=3,new_C_union_calls=0,old_baseline_calls=0,
                   production_descriptor_evaluations=0,real_direction_values=0)
    except Exception as error:
        out.update(status='stopped_failure_no_rerun',error=repr(error))
        raise
    finally:
        out['seconds_before_final_write']=time.monotonic()-START
        (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ('status','pid','exact_toy_volume_calls','mixed_bound_checks','seconds_before_final_write')}),flush=True)


if __name__=='__main__':main()
