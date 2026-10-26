"""Independent integer/rational and bounded-memory JSON helpers for audit130.

No fractions, SymPy, author checker, numerical kernel, logs, or H evaluation.
Geometry reads each axis's path subsequence directly rather than bisecting.
"""
from __future__ import annotations

import hashlib
import json
from math import gcd
import mmap
from pathlib import Path
import re


def rat(n, d=1):
    n, d = int(n), int(d)
    assert d != 0
    if d < 0:
        n, d = -n, -d
    k = gcd(n, d)
    return n // k, d // k


def parse(s):
    if isinstance(s, int):
        return s, 1
    fields = s.split('/')
    assert len(fields) in (1, 2)
    return rat(fields[0], fields[1] if len(fields) == 2 else 1)


def text(q):
    return str(q[0]) if q[1] == 1 else f'{q[0]}/{q[1]}'


def add(x, y):
    g = gcd(x[1], y[1])
    return rat(x[0]*(y[1]//g) + y[0]*(x[1]//g), x[1]*(y[1]//g))


def neg(x):
    return -x[0], x[1]


def sub(x, y):
    return add(x, neg(y))


def mul(x, y):
    g1, g2 = gcd(abs(x[0]), y[1]), gcd(abs(y[0]), x[1])
    return rat((x[0]//g1)*(y[0]//g2), (x[1]//g2)*(y[1]//g1))


def div(x, y):
    assert y[0]
    return mul(x, rat(y[1], y[0]))


def le(x, y):
    return x[0]*y[1] <= y[0]*x[1]


def minimum(xs):
    best = xs[0]
    for x in xs[1:]:
        if not le(best, x):
            best = x
    return best


def maximum(xs):
    best = xs[0]
    for x in xs[1:]:
        if not le(x, best):
            best = x
    return best


def sum_cleared(terms):
    """One integer numerator after clearing every positive denominator."""
    denominator = 1
    for _, d in terms:
        denominator = denominator // gcd(denominator, d) * d
    numerator = sum(n*(denominator//d) for n, d in terms)
    return rat(numerator, denominator)


def roots():
    return [(a,b,c,d) for a in range(32) for b in range(a+1)
            for c in range(b+1) for d in range(c+1)]


def node_key(node):
    rid, path = node
    assert type(rid) is int and 0 <= rid < 52360
    assert isinstance(path, str) and all(b in '01' for b in path)
    return rid, path


def box(node, root_list):
    rid, path = node_key(node)
    bounds = []
    for j, a in enumerate(root_list[rid]):
        subsequence = path[j::4]
        scale = 1 << len(subsequence)
        index = int(subsequence, 2) if subsequence else 0
        denominator = 32*scale
        low = a*scale + index
        bounds.append((rat(low, denominator), rat(low+1, denominator)))
    return bounds


def contains(bounds, point):
    assert len(point) == 5 and parse(point[0]) == (1, 1)
    return all(le(lo, parse(x)) and le(parse(x), hi)
               for (lo, hi), x in zip(bounds, point[1:]))


def split_partition(node, children, axis, root_list):
    parent = node_key(node)
    assert axis == len(parent[1]) % 4
    assert list(map(node_key, children)) == [(parent[0],parent[1]+'0'),(parent[0],parent[1]+'1')]
    outer = box(parent, root_list)
    left, right = (box(child, root_list) for child in children)
    for j in range(4):
        if j == axis:
            assert left[j][0] == outer[j][0] and right[j][1] == outer[j][1]
            assert left[j][1] == right[j][0]
            assert add(left[j][1], left[j][1]) == add(*outer[j])
        else:
            assert left[j] == right[j] == outer[j]


def frontier(leaves, splits):
    """Per-root dyadic interval tiling + exact prefix closure of internal nodes."""
    leafset = set(map(node_key, leaves))
    splitset = set(map(node_key, splits))
    assert len(leafset) == len(leaves) and len(splitset) == len(splits)
    assert not leafset.intersection(splitset)
    implied = set()
    by_root = [[] for _ in range(52360)]
    for rid, path in leafset:
        for i in range(len(path)):
            implied.add((rid,path[:i]))
        by_root[rid].append(path)
    assert implied == splitset
    maxdepth = 0
    for paths in by_root:
        assert paths
        depth = max(map(len,paths))
        maxdepth = max(maxdepth, depth)
        intervals = []
        for path in paths:
            width = 1 << (depth-len(path))
            low = (int(path,2) if path else 0)*width
            intervals.append((low,low+width))
        cursor = 0
        for low, high in sorted(intervals):
            assert low == cursor, (low,cursor)
            cursor = high
        assert cursor == 1 << depth
    assert len(leaves) == 52360+len(splits)
    return {'roots':52360,'leaves':len(leaves),'splits':len(splits),
            'prefix_free':True,'exact_unit_interval_tiling_and_Kraft':True,'maximum_path_depth':maxdepth}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        while chunk := f.read(8*1024*1024):
            h.update(chunk)
    return h.hexdigest()


def digest(obj):
    return hashlib.sha256(json.dumps(obj,sort_keys=True,separators=(',',':')).encode()).hexdigest()


_TOKENS = re.compile(rb'"(?:[^"\\]|\\.)*"|[{}\[\]]')
_STRING = re.compile(rb'"(?:[^"\\]|\\.)*"')
_SPACE = b' \n\r\t'


class JsonView:
    """Top-level selective JSON access; skip the 761MB seed cache without parsing it."""
    def __init__(self, path):
        self.file = Path(path).open('rb')
        self.data = mmap.mmap(self.file.fileno(),0,access=mmap.ACCESS_READ)
        self.fields = {}
        i = self._skip(0)
        assert self.data[i] == 123
        i = self._skip(i+1)
        while self.data[i] != 125:
            match = _STRING.match(self.data,i)
            assert match
            key = json.loads(self.data[i:match.end()])
            i = self._skip(match.end())
            assert self.data[i] == 58
            start = self._skip(i+1)
            end = self._end(start)
            assert key not in self.fields
            self.fields[key] = start,end
            i = self._skip(end)
            if self.data[i] == 44:
                i = self._skip(i+1)
            else:
                assert self.data[i] == 125
        assert self._skip(i+1) == len(self.data)

    def _skip(self,i):
        while i < len(self.data) and self.data[i] in _SPACE:
            i += 1
        return i

    def _end(self,i):
        first = self.data[i]
        if first == 34:
            match = _STRING.match(self.data,i)
            assert match
            return match.end()
        if first in (123,91):
            depth = 0
            for match in _TOKENS.finditer(self.data,i):
                ch = self.data[match.start()]
                if ch in (123,91): depth += 1
                elif ch in (125,93):
                    depth -= 1
                    if depth == 0: return match.end()
            raise AssertionError('Incomplete JSON container')
        end = i
        while self.data[end] not in b',}] \n\r\t': end += 1
        return end

    def get(self,key):
        lo,hi = self.fields[key]
        return json.loads(self.data[lo:hi])

    def array(self,key):
        lo,hi = self.fields[key]
        assert self.data[lo] == 91
        i = self._skip(lo+1)
        while self.data[i] != 93:
            end = self._end(i)
            yield json.loads(self.data[i:end])
            i = self._skip(end)
            if self.data[i] == 44: i = self._skip(i+1)
            else: assert self.data[i] == 93
        assert i+1 == hi

    def close(self):
        self.data.close()
        self.file.close()

    def __enter__(self): return self
    def __exit__(self,*exc): self.close()
