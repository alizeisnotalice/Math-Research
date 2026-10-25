"""Bounded exact checks; independent endpoint-partition oracle, no differentiation."""
from fractions import Fraction as F
from itertools import product, combinations
from pathlib import Path
import hashlib
import json
from normal_intervals import signed_length, pure_satellite_term

P = Path(__file__).resolve().parent


def oracle(left=(), right=(), slabs=(), crossing=False):
    vals = list(left) + list(right) + [v for p in slabs for v in p]
    m = max([F(1)] + [abs(F(v)) + 1 for v in vals])
    cuts = sorted(set([-m, F(0), m] + [F(v) for v in vals]))
    total = F(0)
    for a, b in zip(cuts, cuts[1:]):
        x = (a + b) / 2
        new = (crossing or any(x < v for v in left)
               or any(x > v for v in right)
               or any(a0 < x < b0 for a0, b0 in slabs))
        old = crossing or (bool(left) and x < 0) or (bool(right) and x > 0)
        total += (b - a) * (int(new) - int(old))
    return total


def union_length(intervals):
    # Distinct from the production merge: endpoint cells, midpoint membership.
    cuts = sorted(set(v for p in intervals for v in p))
    return sum((b-a for a,b in zip(cuts,cuts[1:])
                if any(lo < (a+b)/2 < hi for lo,hi in intervals)), F(0))


def union_volume(boxes):
    # Inclusion-exclusion on only three boxes, exact and independent of face logic.
    ans = F(0)
    for k in range(1, len(boxes)+1):
        for subset in combinations(boxes, k):
            v = F(1)
            for j in range(4):
                v *= max(F(0), min(b[1][j] for b in subset)
                         - max(b[0][j] for b in subset))
            ans += (-1)**(k+1) * v
    return ans


def main():
    sides = [(), (F(-1),F(0)), (F(0),F(0)), (F(1),F(3,2))]
    catalog = [(), ((F(-2),F(-1)),), ((F(1),F(2)),),
               ((F(-1),F(1)),), ((F(0),F(0)),), ((F(1),F(-1)),),
               ((F(-2),F(-1)),(F(1),F(2))),
               ((F(-2),F(0)),(F(0),F(2))),
               ((F(-2),F(1)),(F(-1),F(2))),
               ((F(-1),F(0)),(F(-1),F(0))),
               ((F(-3,2),F(-1,2)),(F(1,2),F(3,2))),
               ((F(-2),F(2)),(F(-1),F(1))),
               ((F(0),F(1,2)),), ((F(-1,2),F(0)),),
               ((F(-1,3),F(2,3)),), ((F(-2),F(-2)),(F(2),F(2)))]
    count = 0
    for left,right,slabs,crossing in product(sides,sides,catalog,[False,True]):
        got = signed_length(left,right,slabs,crossing=crossing)
        assert got == oracle(left,right,slabs,crossing), (left,right,slabs,crossing)
        count += 1
    assert count == 512
    # Internal gap (-1,2), two disjoint newborn pieces, total fill=1, loss=2.
    assert signed_length([F(-1)],[F(2)],[(F(-1,2),F(0)),(F(1),F(3,2))]) == -2
    assert signed_length([F(-1)],[F(2)],[(F(-1,2),F(1,2))]) == -2
    group_a = []
    group_b = []
    for t in [F(1,8),F(1,16),F(1,32),F(1,64)]:
        a = union_length([(F(-4),-t),(2*t,F(4)),
                          (-t/2,F(0)),(t,3*t/2)])
        assert a - 8 == -2*t
        group_a.append(dict(epsilon=str(t), volume=str(a), difference=str(a-8),
                            exact_linear_term=str(-2*t)))
        boxes = [((F(-2),F(0),F(0),F(0)),
                  (-t,1+t,1+2*t,1+3*t)),
                 ((2*t,F(0),F(0),F(0)),
                  (F(2),1-t,1+t,1-2*t)),
                 ((-t/2,F(0),F(0),F(0)),
                  (t/2,F(1),F(1),F(1)))]
        b = union_volume(boxes)
        # Symbolic product identity, not an inferred derivative:
        # (2-t)(1+t)(1+2t)(1+3t)
        # +(2-2t)(1-t)(1+t)(1-2t)+t.
        exact = 4+6*t+18*t*t+7*t**3-10*t**4
        assert b == exact
        group_b.append(dict(epsilon=str(t),volume=str(b),difference=str(b-4),
                            linear_term=str(6*t),remainder=str(b-4-6*t),
                            exact_polynomial="4+6*t+18*t^2+7*t^3-10*t^4"))
    # Parent requested a separate N0=0 applicability check after first batch.
    # These are five fixed cone identities, not additional interval-grid search.
    center = (0,0,0,0)
    cone_cases = [[], [((-1,)*4,(1,)*4)], [((0,)*4,(1,)*4)],
                  [((-1,)*4,(0,1,1,1))],
                  [((-1,)*4,(0,1,1,1)),((0,-1,-1,-1),(1,)*4)]]
    cone_values = [pure_satellite_term(center,1,boxes) for boxes in cone_cases]
    assert cone_values == [F(1),F(0),F(15,16),F(1,2),F(0)]
    out = dict(status="exact finite checks passed; local theorem proof is separate",
               rational_configurations=count,group_a=group_a,group_b=group_b,
               pure_satellite_fixed_checks=list(map(str,cone_values)),
               production_sha256=hashlib.sha256((P/'normal_intervals.py').read_bytes()).hexdigest(),
               checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (P/'round005_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__ == '__main__':
    main()
