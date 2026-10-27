#!/usr/bin/env python3
"""One new exact near-boundary audit case, not the author's test batch."""
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
SOURCE = Path('/Users/zhengzhihao/Desktop/T/cube-r040-star-lp/cube-dim-7-bounds/upper/round040/star_lp/PROOF.md')
EXPECTED = '6aece373926c7dfc65e7373b7245d117ba3fe4ae8f17fcd09492d1acc55a2a95'
assert hashlib.sha256(SOURCE.read_bytes()).hexdigest() == EXPECTED
n, d = 2, 3
a, beta = F(267, 100), F(199, 100)
t = 1 / (beta - 1)
assert a > F(8, 3) and a > 2*t and t > 1
nodes = [p for k in range(d+1) for p in product((-1, 1), repeat=k)]
radius = {p: F(1, 8**len(p)) for p in nodes}
center = {(): F(0)}
for p in nodes[1:]:
    center[p] = center[p[:-1]] + p[-1] * radius[p[:-1]] / 2
is_ancestor = lambda p, q: len(p) < len(q) and q[:len(p)] == p
w = a**(-n)
tau = d*w
small_gaps, unrelated_beta_gaps = [], []
for p, q in combinations(nodes, 2):
    delta = abs(center[p]-center[q])
    gap = delta - (radius[p]+radius[q])/a
    assert gap > 0
    small_gaps.append(gap)
    if not is_ancestor(p,q) and not is_ancestor(q,p):
        gap = delta - beta*(radius[p]+radius[q])
        assert gap > 0
        unrelated_beta_gaps.append(gap)
for p in nodes[1:]:
    h = p[:-1]
    assert abs(center[p]-center[h]) + beta*radius[p] < beta*radius[h]
loads = []
for k,p in enumerate(nodes):
    hits = [q for q in nodes[:k] if abs(center[p]-center[q]) < beta*radius[q]]
    assert set(hits) == {p[:l] for l in range(len(p))}
    load = len(hits)*w
    assert load <= tau
    loads.append(load)

dual = {p: F(1, 2**(len(p)+1)) if len(p)<d else F(1,2**d) for p in nodes}
chain_count = 0
for deepest in nodes:
    ancestors = [deepest[:l] for l in range(len(deepest))]
    for mask in range(2**len(ancestors)):
        chain = [deepest] + [p for k,p in enumerate(ancestors) if mask & (1<<k)]
        rep = min(chain, key=len)
        assert sum((dual[p] for p in chain),F(0)) <= F(1,2**len(rep))
        chain_count += 1
assert sum(dual.values(),F(0)) == 1 + F(d,2)

columns = []
for v in nodes:
    for leaf in nodes:
        if len(leaf)==d and leaf[:len(v)]==v:
            chain = tuple(leaf[:k] for k in range(len(v),d+1))
            qv = F(1) if not v else F(1,2)
            columns.append((v,chain,qv/F(2**(d-len(v)))))
for p in nodes:
    assert sum((z for v,J,z in columns if p in J),F(0)) == 1
for leaf in (p for p in nodes if len(p)==d):
    assert sum((z for v,J,z in columns if leaf[:len(v)]==v),F(0)) == 1+F(d,2)

A = sum((F(2,8**n)**l for l in range(d+1)), F(0))
root_cost = a**(-n)*A
assert A < F(4,3)
assert root_cost < 1
result = dict(
    source_sha256=EXPECTED, arithmetic='fractions.Fraction exact rational',
    note='One new near-a=8/3, beta-near-2, equality-eligibility case; author batch not rerun.',
    n=n,d=d,a=str(a),beta=str(beta),t=str(t),tau=str(tau),
    nodes=len(nodes),all_chain_columns_checked=chain_count,
    primal_path_columns_checked=len(columns),K_star=str(1+F(d,2)),
    min_small_cube_gap=str(min(small_gaps)),
    min_unrelated_inflated_star_gap=str(min(unrelated_beta_gaps)),
    maximum_actual_eligibility=str(max(loads)),
    K_contain=root_cost.__str__(),all_assertions_passed=True)
(HERE/'boundary_result.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
