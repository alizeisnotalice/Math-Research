"""Independent G30A replay for n=4,7 only, with exact leaf-partition checks.

Usage: python verify.py /path/to/upper/round030/compute/g30a_arb.py
Uses the reviewed Arb evaluator, a separate stack traversal, and exact Fractions.
Never executes the target module's __main__ (in particular, never scans n=64).
"""
import hashlib
import importlib.util
import json
from fractions import Fraction as Q
from pathlib import Path
import sys

import flint
from flint import arb, ctx
import sympy as sp

HERE = Path(__file__).resolve().parent
source = Path(sys.argv[1]).resolve()
spec = importlib.util.spec_from_file_location('g30a_reviewed', source)
g = importlib.util.module_from_spec(spec)
spec.loader.exec_module(g)
ctx.prec = 512

def replay(n):
    d = Q(1, 4*n*2**n)
    points = [1+d]
    x = d
    while x < 1:
        x = min(2*x, Q(1))
        points.append(1+x)
    while points[-1] < 4*n:
        points.append(min(2*points[-1], Q(4*n)))
    pending = list(zip(points, points[1:]))
    leaves = []
    count = 0
    while pending:
        assert count < 4096
        a, b = pending.pop()
        count += 1
        H, Hp = g.calc(n,a,b)
        s = g.sign(H)
        if s:
            leaves.append((a,b,'strict',s,s))
            continue
        derivative = g.sign(Hp)
        if derivative:
            left, right = g.sign(g.calc(n,a,a)[0]), g.sign(g.calc(n,b,b)[0])
            if left and right:
                if left != right:
                    assert (left,right,derivative) in ((1,-1,-1),(-1,1,1))
                    leaves.append((a,b,'root',left,right))
                else:
                    leaves.append((a,b,'monotone',left,right))
                continue
        mid = (a+b)/2
        pending.extend(((a,mid),(mid,b)))
    leaves.sort()
    assert leaves[0][0] == 1+d and leaves[-1][1] == 4*n
    assert all(a[1] == b[0] and a[4] == b[3] for a,b in zip(leaves,leaves[1:]))
    assert sum((b-a for a,b,*_ in leaves), Q(0)) == 4*n-(1+d)
    roots = [r for r in leaves if r[2] == 'root']
    assert len(roots) == 2
    assert [(r[3],r[4]) for r in roots] == [(1,-1),(-1,1)]
    assert leaves[0][3] == leaves[-1][4] == 1
    return {'n':n, 'boxes':count, 'leaf_count':len(leaves),
            'strict_leaves':sum(r[2]=='strict' for r in leaves),
            'same_sign_monotone_leaves':sum(r[2]=='monotone' for r in leaves),
            'roots':[[str(r[0]),str(r[1]),r[3],r[4]] for r in roots],
            'rational_coverage':'gap-free, no overlapping interiors, full initial interval',
            'unresolved':0}

rows = [replay(n) for n in (4,7)]

# Independent exact differentiation of the finite integral expression.
t = sp.symbols('t',positive=True)
symbolic = []
for n in (4,7):
    q = 1+1/(2*t)
    S = sum(sp.binomial(n-1,j)*q**(n-1-j)*(-sp.Rational(1,2))**j*(1-t**(-j))/j for j in range(1,n))
    B = 1-((1+1/t)/2)**n+n*sp.log(t)-n*(S+q**(n-1)*sp.log(t))/2
    A = (1+1/t)**n
    D, E = sp.diff(B,t), sp.diff(B,t,2)
    K = t*(t+1)*D/n
    Hp = 2*K*((2*t+1)*D+t*(t+1)*E)/n-sp.diff(A,t)*B-A*D
    assert sp.simplify(sp.diff(K*K-A*B,t)-Hp) == 0
    if n == 4:
        assert sp.simplify(B.subs(t,2)-(1841+72*sp.log(2))/768) == 0
    # Verify the Arb series evaluator against an independently differentiated
    # expression at fixed rational points; this is an algebra/implementation
    # check, not the source of the global root certificate.
    for point in (Q(3,2),Q(2),Q(7,2),Q(6)):
        Hball, Hpball = g.calc(n,point,point)
        for expr, actual in ((K*K-A*B,Hball),(Hp,Hpball)):
            poly = sp.Poly(sp.expand(expr),sp.log(t))
            expected = arb(0)
            for (degree,), coefficient in poly.terms():
                rational = coefficient.subs(t,sp.Rational(point.numerator,point.denominator))
                expected += g.ab(Q(int(sp.numer(rational)),int(sp.denom(rational))))*g.ab(point).log()**degree
            assert actual.overlaps(expected)
    symbolic.append({'n':n,'exact_Hprime_identity':True,'fixed_rational_points_overlap':4})

assert Q(1,4)*(1+Q(1,128))*Q(33,32)**2 < 1
assert Q(1913,768) < Q(49,16)
assert Q(1267,32)-18*Q(7,4) == Q(259,32) > 6

receipt = {'status':'ACCEPT n=4,7 parameter-family uniqueness only',
           'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
           'python_flint':flint.__version__,'sympy':sp.__version__,'precision_bits':ctx.prec,
           'runs':rows,'symbolic':symbolic,
           'n64':'not rerun; original 4041 residual boxes prevent root completeness',
           'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(HERE/'verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
