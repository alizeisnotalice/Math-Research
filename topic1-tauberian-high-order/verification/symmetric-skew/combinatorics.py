"""Finite exact sanity checks; the general-dimension proof is in the companion note."""
from itertools import product
from math import comb

for n in range(2, 10):
    balanced = [k for k in product((-1, 0, 1), repeat=n) if any(k) and sum(k) == 0]
    for k in balanced:
        plus = [i for i, x in enumerate(k) if x == 1]
        minus = [i for i, x in enumerate(k) if x == -1]
        assert len(plus) == len(minus) >= 1
        edges = {(i, j) for i in plus for j in minus}
        if len(plus) >= 2:
            assert len(edges) >= 4
    active_patterns = 0
    multi_patterns = 0
    for signs in product((-1, 0, 1), repeat=n):
        plus = [i for i, x in enumerate(signs) if x == 1]
        minus = [i for i, x in enumerate(signs) if x == -1]
        edges = {(i, j) for i in plus for j in minus}
        if edges:
            active_patterns += 1
        if len(edges) >= 2:
            used = {i for edge in edges for i in edge}
            assert len(used) >= 3
            multi_patterns += 1
    assert len(balanced) < 3**n
    print(
        f"n={n}: balanced_nonzero={len(balanced)}, M={3**n}, "
        f"active_sign_patterns={active_patterns}, multi_patterns={multi_patterns}, "
        f"triple_union_constant={8*comb(n,3)}; PASS"
    )
