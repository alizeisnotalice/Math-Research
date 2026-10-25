#!/usr/bin/env python3
"""Boolean endpoint-grid cover for the five cap-five boxes.

The old boxes live exactly in Q[u,v]/(u^6-2,v^6-3).  A new endpoint
involving 5^(1/6) is represented only by its D=10^60 rational enclosure.
Crucially, endpoints whose branch is proved to be a shared Q[u,v] endpoint
are retained as exact field elements, not demoted to an interval.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import prod
from pathlib import Path
import hashlib, json, time

D = 10**60
NODE_CAP = 10**6
WALL = 20 * 60
N = 6
TARGETS = ((3,5,5,5,5,5), (3,5,4,5,5,5), (3,5,5,4,5,5),
           (3,5,5,5,4,5), (3,5,5,5,5,4))

class Unresolved(Exception): pass
class Budget(Exception): pass

# Minimal exact Q[u,v]/(u^6-2,v^6-3) implementation, kept local so this
# verifier remains self-contained in the independent round-005 clone.
def floor_root_k(x, k):
    lo, hi = 0, 1
    while hi**k <= x: hi *= 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if mid**k <= x: lo = mid
        else: hi = mid
    return lo

EXACT_D = 10**45
RU, RV = floor_root_k(2 * EXACT_D**N, N), floor_root_k(3 * EXACT_D**N, N)
UL, VL, UH, VH = Q(RU, EXACT_D), Q(RV, EXACT_D), Q(RU+1, EXACT_D), Q(RV+1, EXACT_D)
SZ = N * N

@lru_cache(None)
def fadd(a, b): return tuple(x+y for x, y in zip(a, b))
@lru_cache(None)
def fneg(a): return tuple(-x for x in a)
@lru_cache(None)
def fmul(a, b):
    z = [Q(0)] * SZ
    for i, x in enumerate(a):
        if not x: continue
        ii, ij = divmod(i, N)
        for j, y in enumerate(b):
            if not y: continue
            ji, jj = divmod(j, N)
            z[((ii+ji) % N)*N + ((ij+jj) % N)] += x*y * 2**((ii+ji)//N) * 3**((ij+jj)//N)
    return tuple(z)
def exact_bounds(a):
    lo = hi = Q(0)
    for k, c in enumerate(a):
        if not c: continue
        i, j = divmod(k, N); x, y = UL**i * VL**j, UH**i * VH**j
        if c > 0: lo, hi = lo + c*x, hi + c*y
        else: lo, hi = lo + c*y, hi + c*x
    return lo, hi
def positive(a):
    if not any(a): return False
    lo, hi = exact_bounds(a)
    if lo > 0: return True
    if hi < 0: return False
    raise ArithmeticError('Q[u,v] isolating intervals insufficient')

class F:
    __slots__ = ('a',)
    def __init__(self, n=0, a=None): self.a = tuple(a) if a is not None else (Q(n),) + (Q(0),)*(SZ-1)
    def __hash__(self): return hash(self.a)
    def __eq__(self, b): return self.a == (b.a if isinstance(b, F) else F(b).a)
    def __add__(self, b): return F(a=fadd(self.a, (b if isinstance(b, F) else F(b)).a))
    __radd__ = __add__
    def __neg__(self): return F(a=fneg(self.a))
    def __sub__(self, b): return self + -b
    def __rsub__(self, b): return F(b) + -self
    def __mul__(self, b): return F(a=fmul(self.a, (b if isinstance(b, F) else F(b)).a))
    __rmul__ = __mul__
    def __truediv__(self, k): return F(a=tuple(x/Q(k) for x in self.a))
    def __pow__(self, k):
        z = F(1)
        for _ in range(k): z = z * self
        return z
    def __lt__(self, b): return positive(((b if isinstance(b, F) else F(b))-self).a)
    def __le__(self, b): return self == b or self < b

u = F(a=(Q(0),)*N + (Q(1),) + (Q(0),)*(SZ-N-1))
v = F(a=(Q(0), Q(1)) + (Q(0),)*(SZ-2))
ZERO = F(0)
SP = (F(1)+u**4/2, (u**3+u**4*v**3)/4, *((F(1)+u)/2 for _ in range(4)))

def old_root(k):
    a = b = 0
    while k % 2 == 0: a, k = a+1, k//2
    while k % 3 == 0: b, k = b+1, k//3
    assert k == 1
    return u**a * v**b

def floor_root(x): return floor_root_k(x, 6)
def rootI(n):
    z = floor_root(n * D**6)
    return Q(z, D), Q(z+1, D)
def mul(I, J):
    q = (I[0]*J[0], I[0]*J[1], I[1]*J[0], I[1]*J[1])
    return min(q), max(q)
def sub(I, J): return I[0]-J[1], I[1]-J[0]
def scale(I, c): return (I[0]*c, I[1]*c) if c >= 0 else (I[1]*c, I[0]*c)
def power(I, k): return I[0]**k, I[1]**k

UV, W = (rootI(2), rootI(3)), rootI(5)
def field_bounds(f):
    lo = hi = Q(0)
    for k, c in enumerate(f.a):
        if not c: continue
        i, j = divmod(k, 6); I = mul(power(UV[0], i), power(UV[1], j))
        if c >= 0: lo, hi = lo+c*I[0], hi+c*I[1]
        else: lo, hi = lo+c*I[1], hi+c*I[0]
    return lo, hi
SPI = tuple(field_bounds(a) for a in SP)

def mrI(ell):
    n = prod(ell); exponents = []
    for p in (2, 3, 5):
        e = 0
        while n % p == 0: e, n = e+1, n//p
        exponents.append(e)
    assert n == 1
    ans = (Q(1), Q(1))
    for I, e in zip((UV[0], UV[1], W), exponents): ans = mul(ans, power(I, e))
    return ans

@dataclass(frozen=True)
class Endpoint:
    exact: object
    lo: Q
    hi: Q
    label: str
    @classmethod
    def exact_f(cls, x, label):
        lo, hi = field_bounds(x)
        return cls(x, lo, hi, label)
    @classmethod
    def outer(cls, I, label): return cls(None, I[0], I[1], label)
    @property
    def kind(self): return 'exact_Q[u,v]' if self.exact is not None else 'rational_outer_interval'

def cmp(x, y):
    """Strictly certify order; exact equality is allowed only for exact F."""
    if x.exact is not None and y.exact is not None:
        if x.exact == y.exact: return 0
        return -1 if x.exact < y.exact else 1
    if x.hi < y.lo: return -1
    if x.lo > y.hi: return 1
    raise Unresolved('endpoint interval not separated: %s  versus  %s' % (x.label, y.label))

def old_box(ell):
    r = old_root(prod(ell)) / 2; box = []
    for j, (a, k) in enumerate(zip(SP, ell)):
        lo, hi = (ZERO, min(a/2, r-(k-1)*a/2)) if k % 2 else (max(ZERO, k*a/2-r), a/2)
        if hi <= lo: return None
        box += [Endpoint.exact_f(lo, f'old ell={ell} coordinate={j+1} lower'),
                Endpoint.exact_f(hi, f'old ell={ell} coordinate={j+1} upper')]
    return tuple(box)

def new_box(ell):
    """Outer new box, retaining every proved shared exact endpoint."""
    rr = scale(mrI(ell), Q(1, 2)); box = []
    for j, (a, aI, k) in enumerate(zip(SP, SPI, ell)):
        tag = f'new ell={ell} coordinate={j+1}'
        half, halfI = a/2, scale(aI, Q(1, 2))
        if k % 2:
            lo = Endpoint.exact_f(ZERO, tag + ' odd lower=0 (exact)')
            candidate = sub(rr, scale(aI, Q(k-1, 2)))
            if halfI[1] < candidate[0]:
                hi = Endpoint.exact_f(half, tag + ' odd upper=a/2 (exact proved min branch)')
            elif candidate[1] < halfI[0]:
                hi = Endpoint.outer(candidate, tag + ' odd upper=r-(ell-1)a/2 (outer)')
            else:
                raise Unresolved(tag + ' odd upper min branch interval not separated')
        else:
            # The even upper endpoint is physically a/2, independent of r.
            hi = Endpoint.exact_f(half, tag + ' even upper=a/2 (exact)')
            candidate = sub(scale(aI, Q(k, 2)), rr)
            if candidate[1] < 0:
                lo = Endpoint.exact_f(ZERO, tag + ' even lower=0 (exact proved max branch)')
            elif candidate[0] > 0:
                lo = Endpoint.outer(candidate, tag + ' even lower=ell*a/2-r (outer)')
            else:
                raise Unresolved(tag + ' even lower max branch interval not separated')
        if cmp(lo, hi) >= 0: raise Unresolved(tag + ' new interval nonpositive')
        box += [lo, hi]
    return tuple(box)

def sortcuts(cuts):
    ans = []
    for x in cuts:
        pos = 0
        while pos < len(ans) and cmp(ans[pos], x) < 0: pos += 1
        if pos < len(ans) and cmp(ans[pos], x) == 0: continue
        ans.insert(pos, x)
    return ans
def leq(x, y): return cmp(x, y) <= 0

def cover_one(new, old, started):
    nodes = 0
    @lru_cache(None)
    def rec(active, dim):
        nonlocal nodes
        nodes += 1
        if nodes > NODE_CAP: raise Budget('node cap')
        if time.monotonic() - started > WALL: raise Budget('wall time')
        if dim == 6: return bool(active)
        lo, hi = new[2*dim], new[2*dim+1]
        cuts = [lo, hi]
        for i in active:
            a, b = old[i][2*dim], old[i][2*dim+1]
            if cmp(lo, a) < 0 and cmp(a, hi) < 0: cuts.append(a)
            if cmp(lo, b) < 0 and cmp(b, hi) < 0: cuts.append(b)
        cuts = sortcuts(tuple(cuts))
        for left, right in zip(cuts, cuts[1:]):
            if cmp(left, right) >= 0: raise Unresolved('nonpositive sorted slab: %s / %s' % (left.label, right.label))
            nxt = tuple(i for i in active if leq(old[i][2*dim], left) and leq(right, old[i][2*dim+1]))
            if not nxt or not rec(nxt, dim+1): return False
        return True
    active = tuple(i for i, b in enumerate(old) if all(
        cmp(b[2*j+1], new[2*j]) >= 0 and cmp(b[2*j], new[2*j+1]) <= 0 for j in range(6)))
    ok = rec(active, 0)
    return ok, nodes, len(active)

def summary(box):
    return [{'coordinate': j+1, 'lower_kind': box[2*j].kind, 'upper_kind': box[2*j+1].kind,
             'lower_label': box[2*j].label, 'upper_label': box[2*j+1].label} for j in range(6)]

def main():
    started = time.monotonic()
    old = tuple(b for ell in product(range(1, 5), repeat=6) if (b := old_box(ell)) is not None)
    out = {'scope': 'symbolic endpoint-grid Boolean coverage; no volume arithmetic; no Q(u,v,w)',
           'D': str(D), 'node_cap': NODE_CAP, 'wall_seconds': WALL, 'old_box_count': len(old), 'targets': []}
    for ell in TARGETS:
        item = {'ell': list(ell)}
        try:
            new = new_box(ell); item['endpoint_tags'] = summary(new)
            ok, nodes, active = cover_one(new, old, started)
            item.update(status='covered' if ok else 'not_covered', nodes=nodes, initial_active_boxes=active)
        except (Unresolved, Budget) as exc:
            item.update(status='unresolved' if isinstance(exc, Unresolved) else 'budget_stop', reason=str(exc))
        out['targets'].append(item)
    statuses = {x['status'] for x in out['targets']}
    out['status'] = 'pass' if statuses == {'covered'} else ('budget_stop' if 'budget_stop' in statuses else ('unresolved' if 'unresolved' in statuses else 'not_covered'))
    out['elapsed_seconds'] = time.monotonic() - started
    here = Path(__file__).resolve(); out['code_sha256'] = hashlib.sha256(here.read_bytes()).hexdigest()
    here.with_name('endpoint_grid_cover_fixed_results.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(out, indent=2))

if __name__ == '__main__': main()
