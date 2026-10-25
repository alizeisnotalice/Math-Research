#!/usr/bin/env python3
"""Exact geometry audit for an ordered-simplex face cover; no kernel calls."""
from __future__ import annotations
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import sympy as sp

DIM, LEAF_TARGET = 5, 120

def dist2(a, b):
    return sum((x-y)**2 for x, y in zip(a, b))

def midpoint(a, b):
    return tuple((x+y)/2 for x, y in zip(a, b))

def vol(vs):
    b = sp.Matrix(vs[0])
    return abs(sp.Matrix.hstack(*[sp.Matrix(v)-b for v in vs[1:]]).det())/sp.factorial(DIM)

def bisect(vs):
    _, i, j = max(((dist2(vs[i],vs[j]),i,j) for i,j in combinations(range(DIM+1),2)),
                  key=lambda item: (item[0],-item[1],-item[2]))
    m = midpoint(vs[i],vs[j]); a, b = list(vs), list(vs)
    a[i], b[j] = m, m
    return tuple(a), tuple(b)

def dmax(vs):
    """Exact max of sum lam|v|^2-|sum lam v|^2 by all KKT supports."""
    m = DIM+1
    mat = sp.Matrix([[sp.Rational(x.numerator,x.denominator) for x in v] for v in vs])
    Q, d = mat*mat.T, sp.Matrix([(mat*mat.T)[i,i] for i in range(m)])
    best = None
    for count in range(1,m+1):
        for tup in combinations(range(m),count):
            S=list(tup); qs=Q.extract(S,S); one=sp.ones(count,1)
            K=(2*qs).row_join(one).col_join(one.T.row_join(sp.zeros(1,1)))
            sol=K.inv()*d.extract(S,[0]).col_join(sp.Matrix([1]))
            lam, mu=list(sol[:count,0]),sol[count]
            if any(x<0 for x in lam): continue
            full=sp.zeros(m,1)
            for a,i in enumerate(S): full[i]=lam[a]
            grad=d-2*Q*full
            if any(grad[i]>mu for i in range(m) if i not in S): continue
            value=(d.T*full)[0]-(full.T*Q*full)[0]
            best=value if best is None or value>best else best
    assert best is not None
    return sp.Rational(best)

def main():
    initial=tuple(tuple(Fraction(1) if j<i else Fraction(0) for j in range(DIM))
                  for i in range(DIM+1))
    leaves=[initial]
    while len(leaves)<LEAF_TARGET:
        loc=max(range(len(leaves)),key=lambda k:(vol(leaves[k]),-k))
        parent=leaves.pop(loc); a,b=bisect(parent)
        assert vol(a)+vol(b)==vol(parent)
        leaves.extend([a,b])
    dvs=[dmax(x) for x in leaves]; vols=[vol(x) for x in leaves]
    data={
      "status":"PASS_exact_ordered_simplex_cover",
      "dimension":DIM,
      "fundamental_vertices":[[str(x) for x in v] for v in initial],
      "initial_volume":str(vol(initial)),
      "coordinate_permutation_copies":int(sp.factorial(DIM)),
      "boundary_overlap":"coordinate-tie hyperplanes only (Lebesgue-null)",
      "leaves_in_fundamental_simplex":len(leaves),
      "bisections":len(leaves)-1,
      "physical_simplex_copies":len(leaves)*int(sp.factorial(DIM)),
      "sum_leaf_volumes":str(sum(vols)),
      "positive_face_volume":str(sp.factorial(DIM)*sum(vols)),
      "initial_exact_variance_max":str(dmax(initial)),
      "leaf_exact_variance_min":str(min(dvs)),
      "leaf_exact_variance_max":str(max(dvs)),
      "variance_method":"exact rational KKT enumeration over all nonempty vertex supports",
      "kernel_evaluations":0,
      "source_sha256":sha256(Path(__file__).read_bytes()).hexdigest()
    }
    Path(__file__).with_suffix(".json").write_text(json.dumps(data,indent=2,sort_keys=True)+"\n")
    print(json.dumps(data,indent=2,sort_keys=True))

if __name__=="__main__":
    main()
