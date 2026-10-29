"""Implicit exact-local-phase oracle for the actual periodic source average.

Fractions enter before any global-coordinate floating point conversion. GF labels
are evaluated from g**j by binary exponentiation; no period-sized table is built.
The NPZ adapter exists only for independent historical small-instance validation.
"""
from fractions import Fraction as F
from pathlib import Path
import json
import os
import sys

# A normal installed package is preferred. Archived copies can point to any
# existing vendor root; the historical local project is only a convenience fallback.
LOCAL_NUMERICS = next((parent/'output/hl_cube_paper_20261013/numerics'
                       for parent in Path(__file__).resolve().parents
                       if (parent/'output/hl_cube_paper_20261013/numerics').is_dir()), None)
try:
    import mpmath as mp
except ModuleNotFoundError as exc:
    if exc.name != 'mpmath':
        raise
    vendor = os.environ.get('CUBE_MPMATH_VENDOR')
    if vendor is None and LOCAL_NUMERICS is not None:
        vendor = str(LOCAL_NUMERICS/'vendor')
    if vendor is None:
        raise ModuleNotFoundError('mpmath is required. Install it or set CUBE_MPMATH_VENDOR to an existing vendor root containing mpmath/.') from exc
    sys.path.insert(0, vendor)
    import mpmath as mp
mp.mp.dps = 70


def exact(x):
    if isinstance(x, float):
        raise TypeError('Pass Fraction/integer/decimal string; global float can lose peak phase.')
    return F(x)


def gf_mul(a, b, modulus, d):
    out = 0
    while b:
        if b & 1:
            out ^= a
        b >>= 1
        a <<= 1
        if a >> d:
            a ^= modulus
    return out


def gf_pow(a, k, modulus, d):
    out = 1
    while k:
        if k & 1:
            out = gf_mul(out, a, modulus, d)
        a = gf_mul(a, a, modulus, d)
        k >>= 1
    return out


def gf_trace(x, modulus, d):
    acc = 0
    for _ in range(d):
        acc ^= x
        x = gf_mul(x, x, modulus, d)
    if acc not in (0, 1):
        raise ValueError('Trace outside F2: invalid field representation')
    return acc


def xor_rank(rows):
    pivots = {}
    for x in rows:
        while x:
            k = x.bit_length() - 1
            if k in pivots:
                x ^= pivots[k]
            else:
                pivots[k] = x
                break
    return len(pivots)


class GFLabels:
    def __init__(self, d, b, modulus, generator, basis):
        self.d, self.b = int(d), int(b)
        self.modulus, self.generator = int(modulus), int(generator)
        self.basis = tuple(map(int, basis))
        self.K, self.B = (1 << self.d) - 1, 1 << self.b
        if len(self.basis) != self.b or xor_rank(self.basis) != self.b:
            raise ValueError('basis dimension mismatch')
        if self.modulus.bit_length() != self.d + 1:
            raise ValueError('modulus must include highest bit')
        # The trace pairing is compiled to b parity masks (O(b*d) storage).
        self.masks = tuple(sum(gf_trace(gf_mul(beta, 1 << k, self.modulus, self.d),
                                              self.modulus, self.d) << k
                               for k in range(self.d)) for beta in self.basis)
        self.calls = 0

    @classmethod
    def from_manifest(cls, path):
        obj = json.loads(Path(path).read_text())
        obj = obj.get('gf', obj)
        def integer(x):
            return int(x, 0) if isinstance(x, str) else int(x)
        mod = obj.get('modulus', obj.get('modulus_int'))
        gen = obj.get('generator', obj.get('primitive_element'))
        return cls(obj['d'], obj['b'], integer(mod), integer(gen),
                   [integer(x) for x in obj['basis']])

    def state(self, j):
        return gf_pow(self.generator, int(j) % self.K, self.modulus, self.d)

    def label(self, j):
        self.calls += 1
        state = self.state(j)
        return sum(((state & mask).bit_count() & 1) << i for i, mask in enumerate(self.masks))

    def trace_rows(self, shifts):
        return [gf_mul(self.state(s), beta, self.modulus, self.d)
                for s in shifts for beta in self.basis]

    def endpoint_rank(self, shifts):
        return xor_rank(self.trace_rows(shifts))


class NPZLabels:
    """Historical full table is loaded only in the bounded validation adapter."""
    def __init__(self, path, b):
        import numpy as np
        with np.load(path) as z:
            self.labels = z['labels'].copy()
        self.K, self.B = len(self.labels), 1 << b
        self.calls = 0
        if int(self.labels.max()) >= self.B:
            raise ValueError('label exceeds B')

    def label(self, j):
        self.calls += 1
        return int(self.labels[int(j) % self.K])


def mpf(x):
    x = exact(x)
    return mp.mpf(x.numerator) / x.denominator


class SourceOracle:
    def __init__(self, labels, delta, p=F(1, 4), M=32):
        self.labels, self.K, self.B = labels, labels.K, labels.B
        self.delta, self.p, self.M = exact(delta), exact(p), int(M)
        self.h = self.K * self.delta / 2
        if not (0 < self.h < F(1, 2 * self.B)):
            raise ValueError('Require box confined in its own cell')
        self.coeff = []
        for m in range(1, self.M + 1):
            bm = 4 * mp.mpf(729) / 29 * ((mp.sin(mp.pi*m*10/27)/(mp.pi*m))**4
                            + 3*(mp.sin(mp.pi*m/9)/(mp.pi*m))**4)
            self.coeff.append(mp.sqrt(2*bm/3))

    def point(self, coarse, phase):
        """Exact physical coordinate from unrestricted integer grid + local phase."""
        return (int(coarse) + exact(phase)) / self.K

    def primitive(self, x):
        """P_delta(x) by integrating P over at most two adjacent coarse cells."""
        z = exact(x) * self.K
        left, right = z - self.h, z + self.h
        start, stop = left.numerator // left.denominator, right.numerator // right.denominator
        area = F(0)
        for cell in range(start, stop + 1):
            a, b = max(left-cell, F(0)), min(right-cell, F(1))
            if b <= a:
                continue
            u = F(2*self.labels.label(cell)+1, 2*self.B)
            area += max(b-u, F(0)) - max(a-u, F(0)) - (b*b-a*a)/2
        return area/(2*self.h)

    def primitive_one_cell(self, x):
        """Equivalent single-label formula when every smoothed peak stays in-cell.

        The multi-cell integral above remains the implementation used by query,
        providing an independent check of the useful clipped-ramp expression.
        """
        z = exact(x)*self.K
        cell = z.numerator//z.denominator
        r = z-cell
        u = F(2*self.labels.label(cell)+1,2*self.B)
        return min(F(1),max(F(0),(r-u+self.h)/(2*self.h)))-r

    def fourier(self, x, L):
        x, L = exact(x), exact(L)
        return sum((c * mp.sin(mp.pi*mpf((m*L) % 2)) *
                    mp.cos(2*mp.pi*mpf((m*x) % 1))
                    for m, c in enumerate(self.coeff, 1)), mp.mpf(0))

    def query(self, x, lengths):
        """Every scale shares x, labels and the full Fourier component."""
        x = exact(x)
        result = []
        for L in map(exact, lengths):
            if L <= 0:
                raise ValueError('L must be positive')
            h = self.primitive(x+L/2)-self.primitive(x-L/2)
            b = self.fourier(x, L)
            d = mpf(self.p/(self.K*L))*(mpf(h)+b)
            result.append({'L': L, 'H': h, 'B': b, 'd': d, 'Y': 1+d})
        return result

    def direct_small(self, x, L):
        """Independent exact full-period peak-window overlaps; O(K) validation only."""
        x, L = exact(x), exact(L)
        lo, hi = x-L/2, x+L/2
        periods = L.numerator // L.denominator
        rem_hi = lo + (L-periods)
        mass = F(periods)
        # The remainder is shorter than one period. Include every shifted peak
        # which can meet it; work is independent of the number of full periods.
        first, last = lo.numerator//lo.denominator-1, rem_hi.numerator//rem_hi.denominator+1
        for j in range(self.K):
            center = F(j, self.K) + F(2*self.labels.label(j)+1, 2*self.B*self.K)
            for shift in range(first, last+1):
                overlap = max(F(0), min(rem_hi, center+shift+self.delta/2)
                                      -max(lo, center+shift-self.delta/2))
                mass += overlap/(self.K*self.delta)
        fourier_integral = mp.mpf(0)
        for m, c in enumerate(self.coeff, 1):
            sin_diff = mp.sin(2*mp.pi*mpf((m*hi) % 1))-mp.sin(2*mp.pi*mpf((m*lo) % 1))
            fourier_integral += c*sin_diff/2
        Y = mpf(1-self.p)+mpf(self.p*mass/L)+mpf(self.p/(self.K*L))*fourier_integral
        return {'mass': mass, 'Y': Y, 'd': Y-1}
