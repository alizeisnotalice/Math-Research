#!/usr/bin/env python3
"""Independent exact checker for one fixed histogram selector.

Public entry points: histogram_coefficients(weights), selected_numerator(weights,
selector), prior_weights(n,K,R,selector), cube_values(tables,x,r,selector).
weights has shape [R][n][K]; Fraction inputs give exact Fraction outputs.
tables[b][k] is a sequence of (left,right,height) segments inside cell k/K.
All segments and all query coordinates are interpreted periodically with period 1.
Run this file to reproduce AUDIT.json. No main-program code is imported.
"""
from fractions import Fraction as F
from itertools import product
from math import factorial, prod
from pathlib import Path
import json


def compositions(n, k):
    if k == 1:
        yield (n,)
    else:
        for j in range(n + 1):
            for tail in compositions(n - j, k - 1):
                yield (j,) + tail


def histogram_coefficients(weights):
    """Coefficients of product_i sum_k weights[i][k] z_k, sparse DP."""
    if not weights:
        raise ValueError("provide at least one coordinate")
    k = len(weights[0])
    d = {(0,) * k: 1}
    for row in weights:
        if len(row) != k or any(w < 0 for w in row):
            raise ValueError("nonnegative rectangular weights required")
        nd = {}
        for a, v in d.items():
            for j, w in enumerate(row):
                aa = a[:j] + (a[j] + 1,) + a[j+1:]
                nd[aa] = nd.get(aa, 0) + v * w
        d = nd
    return d


def _table(selector, a):
    return selector(a) if callable(selector) else selector[a]


def selected_numerator(weights, selector):
    """Unnormalized cube integral; the SAME selector is used for every table."""
    ans = 0
    R = len(weights)
    for b, wb in enumerate(weights):
        for a, v in histogram_coefficients(wb).items():
            sb = _table(selector, a)
            if not isinstance(sb, int) or not 0 <= sb < R:
                raise ValueError("invalid selector table")
            if sb == b:
                ans += v
    return ans


def direct_numerator(weights, selector):
    """Independent exponential expansion for small n, with no polynomial DP."""
    R, n, K = len(weights), len(weights[0]), len(weights[0][0])
    ans = 0
    for labels in product(range(K), repeat=n):
        a = tuple(labels.count(k) for k in range(K))
        b = _table(selector, a)
        if not 0 <= b < R:
            raise ValueError("invalid selector table")
        ans += prod(weights[b][i][labels[i]] for i in range(n))
    return ans


def prior_weights(n, K, R, selector):
    pi = [F(0) for _ in range(R)]
    for a in compositions(n, K):
        b = _table(selector, a)
        if not isinstance(b, int) or not 0 <= b < R:
            raise ValueError("invalid selector table")
        m = factorial(n) // prod(factorial(v) for v in a)
        pi[b] += F(m, K ** n)
    return pi


def prior_numerator(weights, selector):
    R, n, K = len(weights), len(weights[0]), len(weights[0][0])
    pi = prior_weights(n, K, R, selector)
    return sum(pi[b] * prod(sum(row) for row in weights[b]) for b in range(R))


def periodic_segment_integral(segments, lo, hi):
    """Exact integral on [lo,hi] of periodized constant-height segments."""
    if hi < lo:
        raise ValueError("invalid interval")
    ans = 0
    for left, right, height in segments:
        for shift in range((lo-right).__floor__(), (hi-left).__ceil__()+1):
            overlap = min(hi, right+shift) - max(lo, left+shift)
            if overlap > 0:
                ans += overlap * height
    return ans


def validate_tables(tables):
    K = len(tables[0])
    for table in tables:
        if len(table) != K:
            raise ValueError("inconsistent K")
        for k, segments in enumerate(table):
            old = F(k, K)
            mass = 0
            for left, right, height in sorted(segments):
                if not F(k,K) <= left < right <= F(k+1,K) or height < 0 or left < old:
                    raise ValueError("segments must be nonnegative, inside cell, disjoint")
                old = right
                mass += (right-left) * height
            if mass != F(1,K):
                raise ValueError("each cell must have EXACT mass 1/K")
    return True


def cube_weights(tables, x, r):
    if r <= 0:
        raise ValueError("positive radius required")
    return [[[periodic_segment_integral(segments, xi-r, xi+r)
              for segments in table] for xi in x] for table in tables]


def cube_values(tables, x, r, selector):
    """Return selected and matched-prior cube averages."""
    w = cube_weights(tables, x, r)
    den = (2*r) ** len(x)
    return selected_numerator(w, selector)/den, prior_numerator(w, selector)/den


def coordinate_cell_min(segments, a, c, r):
    """Minimum moving-window integral for center x in [a,c], exact for boxes."""
    points = {a,c}
    for left,right,_ in segments:
        for edge in (left,right):
            for shift in range((a-r-edge).__floor__()-1, (c+r-edge).__ceil__()+2):
                for rr in (-r,r):
                    x = edge+shift+rr
                    if a <= x <= c:
                        points.add(x)
    return min(periodic_segment_integral(segments,x-r,x+r) for x in points)


def query_cell_lower(tables, cell, r, selector):
    """Certified uniform lower average throughout product_i cell[i].

    Takes coordinatewise minima before the positive polynomial DP. This can be
    conservative, but is exact rational evidence when inputs are Fraction.
    """
    w = [[[coordinate_cell_min(segments,a,c,r) for segments in table]
           for a,c in cell] for table in tables]
    return selected_numerator(w,selector)/(2*r)**len(cell)


def prototype_tables(K,R):
    # Narrow positive intervals stay strictly inside each source cell.
    out=[]
    for b in range(R):
        table=[]
        for k in range(K):
            left=F(k,K)+F(2+b,16*K)
            width=F(1,4*K)
            table.append([(left,left+width,F(4))])
        out.append(table)
    validate_tables(out)
    return out


def run_audit():
    import numpy as np
    rng=np.random.default_rng(420042)
    exact_count=0
    float_max=0.0
    cases=[]
    # Arbitrary nonnegative weights: DP vs independent K**n expansion.
    for K in (2,3):
        for n in (1,2,3,4):
            R=3
            selector={a:int(rng.integers(R)) for a in compositions(n,K)}
            pi=prior_weights(n,K,R,selector)
            assert sum(pi)==1
            mass_weights=[[[F(1,K)]*K for _ in range(n)] for _ in range(R)]
            mass=selected_numerator(mass_weights,selector)
            assert mass==direct_numerator(mass_weights,selector)==1
            for trial in range(7):
                w=[[[F(int(rng.integers(0,11)),13) for _ in range(K)]
                    for _ in range(n)] for _ in range(R)]
                got=selected_numerator(w,selector)
                direct=direct_numerator(w,selector)
                assert got==direct
                wf=np.asarray(w,dtype=float).tolist()
                err=abs(float(selected_numerator(wf,selector))-float(direct))
                float_max=max(float_max,err)
                exact_count+=1
            tables=prototype_tables(K,R)
            # Include boundary crossing, long intervals, and integer query shifts.
            x=[F(2*i+1,11) for i in range(n)]
            for r in (F(1,13),F(5,17),F(9,7)):
                value,prior=cube_values(tables,x,r,selector)
                assert cube_values(tables,[xi+3 for xi in x],r,selector)==(value,prior)
                tau=F(2,7)
                shifted=[[[ (l+tau,u+tau,h) for l,u,h in segments]
                          for segments in table] for table in tables]
                # Cell LABELS move with the sources, so stored per-k segments move.
                assert cube_values(shifted,[xi+tau for xi in x],r,selector)==(value,prior)
                weights=cube_weights(tables,x,r)
                assert selected_numerator(weights,selector)==direct_numerator(weights,selector)
                # Interchange table-mixture sum with direct label expansion.
                direct_prior=sum(pi[b]*sum(prod(weights[b][i][z[i]] for i in range(n))
                    for z in product(range(K),repeat=n)) for b in range(R))
                assert direct_prior==prior_numerator(weights,selector)
                exact_count+=1
            full_values=cube_values(tables,x,F(1),selector)
            assert full_values==(F(1),F(1))
            cell=[(xi-F(1,100),xi+F(1,100)) for xi in x]
            lower=query_cell_lower(tables,cell,F(5,17),selector)
            for corner in product(*[(a,(a+c)/2,c) for a,c in cell]):
                assert cube_values(tables,corner,F(5,17),selector)[0]>=lower
            cases.append(dict(K=K,n=n,R=R,mass=str(mass),prior=[str(p) for p in pi],
                              query_cell_lower=str(lower)))
    all_selector_checks=0
    for K in (2,3):
        n,R=2,2
        hist=list(compositions(n,K))
        mass_weights=[[[F(1,K)]*K for _ in range(n)] for _ in range(R)]
        for assignments in product(range(R),repeat=len(hist)):
            s=dict(zip(hist,assignments))
            assert selected_numerator(mass_weights,s)==direct_numerator(mass_weights,s)==1
            assert sum(prior_weights(n,K,R,s))==1
            all_selector_checks+=1
    # Intended prototype size, one independent floating enumeration of K**n labels.
    K,n,R=4,8,4
    s={a:int(rng.integers(R)) for a in compositions(n,K)}
    w=rng.uniform(0,1,size=(R,n,K)).tolist()
    proto_dp=selected_numerator(w,s)
    proto_direct=direct_numerator(w,s)
    proto_error=abs(proto_dp-proto_direct)
    assert proto_error <= 1e-10*max(1,abs(proto_direct))
    # An exactly computable, nonzero finite-grid weak lower bound.
    K,n,R=2,2,3
    tables=prototype_tables(K,R)
    selector={a:(a[0]+2*a[1])%R for a in compositions(n,K)}
    radii=[F(1,8),F(1,4),F(3,8)]
    threshold=F(1,2)
    grid=8
    good=0
    for index in product(range(grid),repeat=n):
        cell=[(F(j,grid),F(j+1,grid)) for j in index]
        if max(query_cell_lower(tables,cell,r,selector) for r in radii)>threshold:
            good+=1
    p=F(good,grid**n)
    # Interior integer periods for an L1 cutoff: N=10, ceil(max radii)=1.
    N,q=10,1
    cutoff_factor=F(N-q,N)**n
    cutoff_weights=cube_weights(tables,[F(0)]*n,F(N))
    cutoff_mass=selected_numerator(cutoff_weights,selector)
    assert cutoff_mass==(2*N)**n
    return dict(status="PASS",seed=420042,numpy_version=np.__version__,
        exact_dp_expansion_checks=exact_count,float_dp_max_absolute_error=float_max,
        all_selector_mass_checks=all_selector_checks,
        intended_shape_float_check=dict(K=4,n=8,R=4,
            dp=proto_dp,direct=proto_direct,absolute_error=proto_error,
            histogram_states=len(s)),
        cases=cases,grid_certificate=dict(n=n,K=K,R=R,grid=grid,
            radii=[str(r) for r in radii],threshold=str(threshold),
            certified_cells=good,total_cells=grid**n,measure_lower=str(p),
            periodic_weak_lower=str(threshold*p),cutoff_N=N,
            cutoff_mass=str((2*N)**n),cutoff_interior_factor=str(cutoff_factor),
            cutoff_weak_lower=str(threshold*p*cutoff_factor)),
        scope="Exact rational small-case formulas and one grid certificate; no asymptotic or log gain claim")


if __name__=="__main__":
    result=run_audit()
    path=Path(__file__).with_name("AUDIT.json")
    path.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2))
