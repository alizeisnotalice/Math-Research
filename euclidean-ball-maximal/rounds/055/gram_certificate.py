#!/usr/bin/env python3
"""Exact certificates for a local Gram-normalization obstruction, not a global weak-type counterexample."""
from fractions import Fraction as F
import json


def radius(j):
    return F(4, 2**j)


def trace(atoms, t, depth):
    return {j: sum((w for a, w in atoms if abs(t-a) < radius(j)), F(0))/(2*radius(j)) for j in range(1, depth+1)}


def maximal(atoms, t):
    ds = sorted({abs(t-a) for a, w in atoms})
    assert ds[0] > 0
    return max(sum((w for a, w in atoms if abs(t-a) <= d), F(0))/(2*d) for d in ds)


def first(tr, threshold):
    return next(j for j, value in tr.items() if value > threshold)


def case(k, batch):
    r = radius(k)
    mass = F(23, 100)
    p = F(3, 2)*r
    atoms = [(F(0), p), (F(7, 20), mass), (F(10), 1-mass-p)]
    x, z = F(1, 4)-r/10, -F(3, 5)*r
    eps, smooth_h = r/1000, r/100000
    assert sum(w for a, w in atoms) == 1 and all(w > 0 for a, w in atoms)
    tx, tz = trace(atoms, x, k+1), trace(atoms, z, k+1)
    assert first(tx, F(1, 8)) == first(tz, F(1, 8)) == 3
    assert first(tx, F(1, 2)) == 5 and first(tz, F(1, 2)) == k
    assert first(tx, F(2, 5)) == 4 and first(tz, F(2, 5)) == k
    assert maximal(atoms, x) == F(23, 20)/(1+r)
    assert maximal(atoms, z) == F(5, 4)
    W = min(F(1, 2), tx[4], tz[k])-max(F(1, 4), max(tx[j] for j in range(1, 4)), max(tz[j] for j in range(1, k)))
    assert W == F(17, 200)+3*r
    common = sum((w for a, w in atoms if abs(x-a) < radius(4) and abs(z-a) < r), F(0))
    qx = sum((w for a, w in atoms if abs(x-a) < radius(4)), F(0))
    qf = sum((w for a, w in atoms if abs(z-a) < r), F(0))
    lost = sum((w for a, w in atoms if abs(x-a) < radius(4) and abs(z-a) >= radius(4)), F(0))
    assert common == qf == p and qx == mass+p and lost == mass
    sigma, tau, chi2 = common/qf, common/qx, common**2/(qx*qf)
    theta = W*(2*radius(4))/lost
    assert sigma == 1 and tau == chi2 == p/(mass+p)
    assert F(17, 92) <= theta <= 1
    H = common/(2*radius(4)*2*r)
    point_fee = W/F(1, 4)*H
    assert H == F(3, 2) and point_fee == 6*W >= F(51, 100)
    assert (2*radius(4)*H)**2/(chi2*tx[4]*tz[k]) == 2**(k-4)
    for t in (x, z):
        e = eps+smooth_h
        margin = min(abs(abs(t-a)-radius(j)) for a, w in atoms for j in range(1, k+2))
        assert margin > e
        dmin = min(abs(t-a) for a, w in atoms)
        q = e/dmin
        assert 0 < q < 1
        assert maximal(atoms, t)/(1+q) > 1
        assert maximal(atoms, t)/(1-q) < 2
        for v in (t-eps, t+eps):
            assert trace(atoms, v, k+1) == trace(atoms, t, k+1)
    assert x-z-2*eps > radius(4)
    assert x-z+2*eps < radius(4)+r
    assert abs(x)+eps+smooth_h < radius(4)
    assert abs(z)+eps+smooth_h < r
    # Fixed observer interval (6/25,13/50) retains original eligible status.
    for t in (F(6,25), F(1,4), F(13,50)):
        tr = trace(atoms, t, k)
        assert first(tr, F(1,8)) == 3 and first(tr, F(1,2)) == 5
        dmin = min(abs(t-a) for a, w in atoms)
        q = smooth_h/dmin
        assert maximal(atoms, t)/(1+q) > 1 and maximal(atoms, t)/(1-q) < 2
    fee = point_fee*4*eps**2
    assert fee >= F(51,25)*eps**2 > 0
    return dict(batch=batch, k=k, radius=str(r), W=str(W), sigma='1', theta=str(theta), gram_overlap_squared=str(chi2), original_point_fee=str(point_fee), observer_halfwidth=str(eps), smoothing_support=str(smooth_h), box_fee=str(fee), missing_volume_factor_squared=str(2**(k-4)))


def main():
    batches = ((9,10,12),(16,24,32),(64,128,256))
    results = [case(k, batch) for batch, ks in enumerate(batches) for k in ks]
    return dict(status='passed', cases=results, exact_arithmetic='Fraction', scope='Actual complete-band and first-clock local cells; observer box area tends to zero. No global O/X or weak-type counterexample.')


if __name__ == '__main__':
    print(json.dumps(main(), indent=2))
