"""Small independent rational checks of R035's new compression examples.

No LP solver, dimension/parameter search, or earlier baseline replay.
Usage: python verify.py /path/to/upper/round035/joint_cover/proof.md
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent

def cells(intervals):
    cuts = sorted({e for v in intervals for e in v})
    for lo,hi in zip(cuts,cuts[1:]):
        x = (lo+hi)/2
        yield hi-lo,{i for i,(a,b) in enumerate(intervals) if a<x<b}

def load(intervals,q):
    return max([F(0)]+[sum((q[i] for i in active),F(0)) for _,active in cells(intervals)])

def adjacent(v,w):
    return max(v[0],w[0])<min(v[1],w[1])

def check(last,module_expected,expected_new):
    vs = [(F(-1),F(5)),(F(-1,4),F(5,4)),(F(1),F(6)),last]
    h = [b-a for a,b in vs]
    qs = [((a+b)/2-(b-a)/6,(a+b)/2+(b-a)/6) for a,b in vs]
    order = sorted(range(4),key=lambda i:-h[i])
    d = [F(0)]*4
    for volume,active in cells(qs):
        if active:
            d[next(i for i in order if i in active)] += volume
    q = [v/b for v,b in zip(d,h)]
    assert q == [F(1,3),F(1,3),F(4,15),F(1,3)]
    assert sum(d,F(0)) == F(14,3)
    c = {0,1}
    module = all(all(adjacent(vs[j],vs[i]) for i in c) or
                 not any(adjacent(vs[j],vs[i]) for i in c) for j in {2,3})
    assert module == module_expected
    new = q.copy()
    new[0] = (h[0]*q[0]+h[1]*q[1])/h[0]
    new[1] = F(0)
    assert sum((x*y for x,y in zip(new,h)),F(0)) == F(14,3)
    old_peak,new_peak = load(vs,q),load(vs,new)
    assert old_peak == F(14,15)
    assert new_peak == expected_new
    good = {j for j in {2,3} if all(adjacent(vs[j],vs[i]) for i in c)}
    boundary = {j for j in {2,3} if adjacent(vs[j],vs[0]) and j not in good}
    def restricted_peak(indices):
        return load(vs,[q[i] if i in indices else F(0) for i in range(4)])
    savings = q[0]+q[1]-new[0]
    delta = restricted_peak(boundary)
    rest = restricted_peak({2,3})
    common = restricted_peak(good)
    bound = max(rest,new[0]+common+delta)
    assert new_peak == bound == max(rest,old_peak-savings+delta)
    assert savings == F(1,4)
    return {'module':module,'U':str(sum(d,F(0))),
            'q':[str(x) for x in q],'old_peak':str(old_peak),'new_peak':str(new_peak),
            'savings':str(savings),'boundary_bill':str(delta),
            'payment_conserved':True,'boundary_bound_attained':True}

rows = [check((F(11,2),F(8)),True,F(41,60)),
        check((F(9,2),F(7)),False,F(61,60))]
touch = [(F(-1),F(0)),(F(0),F(1))]
assert not adjacent(*touch)
assert load(touch,[F(1),F(1)]) == 1
# An absent positivity assumption would invalidate the clique-weight formula:
# the inner positive box can never be active without its outer negative box.
signed = [(F(1),F(2)),(F(0),F(3))]
assert load(signed,[F(1),F(-1)]) == 0
assert max(F(1),F(-1),F(0)) == 1

receipt = {'status':'passed; R035 new rational examples and assumption guards only',
           'examples':rows,
           'boundary_touching_not_edge':True,'nonnegative_q_is_essential':True,
           'source_proof_sha256':hashlib.sha256(Path(sys.argv[1]).read_bytes()).hexdigest(),
           'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           'scope':'No prior-baseline rerun, LP solve, n7/root/parameter search, MATLAB, or Git mutation.'}
(HERE/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
