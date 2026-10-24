#!/usr/bin/env python3
"""One approved 20-box/60-second exact comparison; never a global proof."""
from time import monotonic
MODULE_STARTED = monotonic()
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
import math
from pathlib import Path
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
R = Q(15, 8)
E = {2: Q(69, 100), 6: Q(19, 50), 7: Q(11, 10), 8: Q(3, 4)}
BOX_LIMIT = 20
SECONDS = 60.
MODULE_IMPORT_SECONDS = monotonic()-MODULE_STARTED


class TimeLimit(Exception):
    pass


class Guard:
    def __init__(self):
        self.started = monotonic()
    def check(self):
        if monotonic()-self.started >= SECONDS:
            raise TimeLimit('60-second comparison budget reached')


class Poly:
    """Exact sparse bivariate polynomial in r,q, with rational coefficients."""
    def __init__(self, data=0):
        if isinstance(data, Poly):
            self.d = dict(data.d)
        elif isinstance(data, dict):
            self.d = {k: Q(v) for k, v in data.items() if v}
        else:
            self.d = {(0, 0): Q(data)} if data else {}
    def __add__(self, other):
        out = dict(self.d)
        for k, v in Poly(other).d.items():
            out[k] = out.get(k, Q(0))+v
        return Poly(out)
    __radd__ = __add__
    def __neg__(self):
        return Poly({k: -v for k, v in self.d.items()})
    def __sub__(self, other):
        return self+-Poly(other)
    def __rsub__(self, other):
        return Poly(other)+-self
    def __mul__(self, other):
        out = {}
        for (i, j), a in self.d.items():
            for (k, l), b in Poly(other).d.items():
                key = (i+k, j+l)
                out[key] = out.get(key, Q(0))+a*b
        return Poly(out)
    __rmul__ = __mul__
    def __truediv__(self, scalar):
        assert not isinstance(scalar, Poly)
        return Poly({k: v/Q(scalar) for k, v in self.d.items()})
    def __pow__(self, n):
        assert isinstance(n, int) and n >= 0
        out, base = Poly(1), self
        while n:
            if n & 1:
                out = out*base
            base = base*base
            n //= 2
        return out
    def diff(self, axis):
        out = {}
        for (i, j), a in self.d.items():
            exponent = i if axis == 0 else j
            if exponent:
                key = (i-1, j) if axis == 0 else (i, j-1)
                out[key] = a*exponent
        return Poly(out)
    def evaluate(self, r, q):
        return sum((a*r**i*q**j for (i, j), a in self.d.items()), Q(0))
    def shift(self, r0, q0, guard):
        out = {}
        for (i, j), a in self.d.items():
            guard.check()
            for k in range(i+1):
                ar = a*math.comb(i, k)*r0**(i-k)
                for l in range(j+1):
                    key = (k, l)
                    out[key] = out.get(key, Q(0))+ar*math.comb(j, l)*q0**(j-l)
        return Poly(out)


class AD:
    """Order-two interval AD; values and every derivative are enclosures."""
    I = None
    log_i = None
    def __init__(self, value, gradient=None, hessian=None):
        self.v = value
        zero = self.I.point(0)
        self.g = gradient if gradient is not None else [zero, zero]
        self.h = hessian if hessian is not None else [[zero, zero], [zero, zero]]
    @classmethod
    def constant(cls, value):
        return cls(cls.I.point(value))
    @classmethod
    def coerce(cls, value):
        return value if isinstance(value, cls) else cls.constant(value)
    @classmethod
    def variable(cls, interval, axis):
        g = [cls.I.point(0), cls.I.point(0)]
        g[axis] = cls.I.point(1)
        return cls(interval, g)
    def __add__(self, other):
        other = self.coerce(other)
        return AD(self.v+other.v, [self.g[i]+other.g[i] for i in range(2)],
                  [[self.h[i][j]+other.h[i][j] for j in range(2)] for i in range(2)])
    __radd__ = __add__
    def __neg__(self):
        return AD(-self.v, [-x for x in self.g], [[-x for x in row] for row in self.h])
    def __sub__(self, other):
        return self+-self.coerce(other)
    def __rsub__(self, other):
        return self.coerce(other)+-self
    def __mul__(self, other):
        other = self.coerce(other)
        return AD(self.v*other.v, [self.g[i]*other.v+self.v*other.g[i] for i in range(2)],
                  [[self.h[i][j]*other.v+self.g[i]*other.g[j]+self.g[j]*other.g[i]+self.v*other.h[i][j]
                    for j in range(2)] for i in range(2)])
    __rmul__ = __mul__
    def __truediv__(self, scalar):
        return self*(Q(1)/Q(scalar))
    def __pow__(self, n):
        assert isinstance(n, int) and n >= 0
        if n == 0:
            return AD.constant(1)
        if n == 1:
            return self
        # Used only at nonnegative u/z+ in this implementation.
        assert self.v.lo >= 0
        value = ipow(self.v, n)
        first = n*ipow(self.v, n-1)
        second = n*(n-1)*ipow(self.v, n-2)
        return AD(value, [first*x for x in self.g],
                  [[second*self.g[i]*self.g[j]+first*self.h[i][j] for j in range(2)] for i in range(2)])
    def log(self):
        assert self.v.lo > 0
        value = AD.log_i(self.v)
        return AD(value, [x/self.v for x in self.g],
                  [[self.h[i][j]/self.v-self.g[i]*self.g[j]/ipow(self.v, 2)
                    for j in range(2)] for i in range(2)])


def ipow(interval, n):
    if n == 0:
        return AD.I.point(1)
    lo, hi = interval.lo**n, interval.hi**n
    if n % 2:
        return AD.I(lo, hi)
    return AD.I(Q(0) if interval.lo <= 0 <= interval.hi else min(lo, hi), max(lo, hi))


def import_source(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def prepare():
    archive = ROOT.parent.parent/'Math-Research-lowdim-agent/cube-lowdim-bounds/experiments'
    core_path = next(p for p in [ROOT.parent/'round008/certify.py', archive/'round008/certify.py'] if p.exists())
    old_path = ROOT.parent/'round016/certify.py'
    c = import_source('round017_interval_core', core_path)
    old = import_source('round017_old_bound', old_path)
    assert old.R == R and old.E == E
    c.R = R
    AD.I, AD.log_i = c.I, staticmethod(c.log_i)
    r, q = Poly({(1, 0): 1}), Poly({(0, 1): 1})
    u, v, s = r*(1+q)/2, r*(1-q)/2, 1-r
    a, b, zp, zm = 1-v, 1-u, 1+r*q, 1-r*q
    polys = {}
    for p in E:
        n = p-1
        P = sum((a**(n-k)*(zp**k+s**k)+b**(n-k)*(zm**k+s**k))/Q(n-k) for k in range(n))
        polys[p] = (P, zp**n+s**n)
    T = 2/(R-1)*(2*sum(r**j/Q((j+2)*R**(j+1)) for j in range(29))
                   +2*sum((r*q)**(2*k)/Q((2*k+2)*R**(2*k+1)) for k in range(15)))
    components = [T, T.diff(0), T.diff(1), T.diff(0).diff(0), T.diff(0).diff(1), T.diff(1).diff(1)]
    assert all(all(a >= 0 for a in poly.d.values()) for poly in components)
    return c, old.make_lower(c), polys, components, core_path, old_path


def box_list():
    out = []
    for j in [4, 6, 8, 10]:
        out.append((f'corner_2^-{j}', (1-Q(1, 2**j), Q(1), 1-Q(1, 2**j), Q(1))))
    for i, (lo, hi) in enumerate([(Q(1, 5), Q(21, 100)), (Q(1, 2), Q(51, 100)),
                                   (Q(4, 5), Q(81, 100)), (Q(9, 10), Q(91, 100))]):
        out.append((f'axis_strip_{i+1}', (lo, hi, 1-Q(1, 256), Q(1))))
    for i, (lo, hi) in enumerate([(Q(0), Q(1, 16)), (Q(1, 4), Q(5, 16)),
                                   (Q(1, 2), Q(9, 16)), (Q(3, 4), Q(13, 16))]):
        out.append((f'support_boundary_{i+1}', (1-Q(1, 256), Q(1), lo, hi)))
    for i, (rc, qc) in enumerate([(Q(109, 250), Q(0)), (Q(43, 50), Q(4, 5)),
                                   (Q(108, 125), Q(4, 5)), (Q(443, 500), Q(13, 15))]):
        out.append((f'active_region_{i+1}', (rc-Q(1, 1000), rc+Q(1, 1000),
                                           Q(0) if i == 0 else qc-Q(1, 1000), qc+Q(1, 1000))))
    out.extend([('small_r', (Q(1, 100), Q(1, 50), Q(0), Q(1, 4))),
                ('bulk_1', (Q(1, 4), Q(13, 50), Q(1, 4), Q(13, 50))),
                ('bulk_2', (Q(3, 5), Q(61, 100), Q(2, 5), Q(41, 100))),
                ('bulk_3', (Q(3, 4), Q(19, 25), Q(3, 4), Q(19, 25)))])
    assert len(out) == BOX_LIMIT and len({box for _, box in out}) == BOX_LIMIT
    return out


def interval_poly(poly, hr, hq, guard):
    value = AD.I.point(0)
    for (i, j), coefficient in poly.d.items():
        guard.check()
        value += coefficient*ipow(hr, i)*ipow(hq, j)
    return value


def shifted_bundles(poly, box, guard):
    rl, rh, ql, qh = box
    rc, qc, rr, rq = (rl+rh)/2, (ql+qh)/2, (rh-rl)/2, (qh-ql)/2
    shifted = poly.shift(rc, qc, guard)
    components = [shifted, shifted.diff(0), shifted.diff(1), shifted.diff(0).diff(0),
                  shifted.diff(0).diff(1), shifted.diff(1).diff(1)]
    intervals = [interval_poly(p, AD.I(-rr, rr), AD.I(-rq, rq), guard) for p in components]
    points = [AD.I.point(p.d.get((0, 0), Q(0))) for p in components]
    def bundle(values):
        return AD(values[0], values[1:3], [[values[3], values[4]], [values[4], values[5]]])
    return bundle(points), bundle(intervals)


def tail_bundles(components, box, guard):
    rl, rh, ql, qh = box
    rc, qc = (rl+rh)/2, (ql+qh)/2
    intervals, points = [], []
    for p in components:
        guard.check()
        intervals.append(AD.I(p.evaluate(rl, ql), p.evaluate(rh, qh)))
        points.append(AD.I.point(p.evaluate(rc, qc)))
    def bundle(values):
        return AD(values[0], values[1:3], [[values[3], values[4]], [values[4], values[5]]])
    return bundle(points), bundle(intervals)


def build_A(c, box, polys, tail_components, guard):
    rl, rh, ql, qh = box
    rc, qc = (rl+rh)/2, (ql+qh)/2
    rp, qp = AD.variable(c.I.point(rc), 0), AD.variable(c.I.point(qc), 1)
    rb, qb = AD.variable(c.I(rl, rh), 0), AD.variable(c.I(ql, qh), 1)
    up, zp = rp*(1+qp)/2, 1+rp*qp
    ub, zb = rb*(1+qb)/2, 1+rb*qb
    point, uniform = tail_bundles(tail_components, box, guard)
    for p, (P, Cu) in polys.items():
        guard.check()
        Pp, Pb = shifted_bundles(P, box, guard)
        Cp, Cb = shifted_bundles(Cu, box, guard)
        n, alpha = p-1, 2*p*E[p]
        point -= alpha*(Pp+Cp*up.log()-zp**n*zp.log())
        uniform -= alpha*(Pb+Cb*ub.log()-zb**n*zb.log())
    return point, uniform


def singular_lower(c, box):
    rl, rh, ql, qh = box
    sl, sh, zml, zmh = 1-rh, 1-rl, 1-rh*qh, 1-rl*ql
    vh = rh*(1-ql)/2
    def w(lo, hi, n):
        if hi == 0:
            return Q(0)
        if lo == 0:
            return hi**n*c.log_point(hi).lo if hi <= Q(1, 4) else -Q(1, 2*n)
        ll, lh = c.log_point(lo).lo, c.log_point(hi).hi
        return min(lo**n*ll, lo**n*lh, hi**n*ll, hi**n*lh)
    C_lo, ws, wz = Q(0), {}, {}
    for p, e in E.items():
        n, alpha = p-1, 2*p*e
        C_lo += alpha*(zml**n+sl**n)
        ws[str(p)], wz[str(p)] = w(sl, sh, n), w(zml, zmh, n)
    favorable = -C_lo*min(Q(0), c.log_point(vh).hi) if C_lo and vh else Q(0)
    total = favorable+sum((2*p*e*(wz[str(p)]+2*ws[str(p)]) for p, e in E.items()), Q(0))
    return total, {'C_lower': str(C_lo), 'v_upper': str(vh), 'favorable_negative_log_lower': str(favorable),
                   'W_s_lower_by_power': {p: str(x) for p, x in ws.items()},
                   'W_zminus_lower_by_power': {p: str(x) for p, x in wz.items()}, 'total_lower': str(total)}


def maxabs(interval):
    return max(abs(interval.lo), abs(interval.hi))


def main():
    path = ROOT/'comparison.json'
    if path.exists():
        raise SystemExit('Comparison output exists; no automatic repeat or completion permitted.')
    prepare_start = monotonic()
    c, old_lower, polys, tail_components, core_path, old_path = prepare()
    construction_seconds = monotonic()-prepare_start
    boxes = box_list()
    counts = {'distinct_boxes_attempted': 0, 'completed_boxes': 0, 'old_bound_calls': 0,
              'new_bundle_calls': 0, 'center_bundles': 0, 'uniform_derivative_bundles': 0,
              'LPs': 0, 'global_certification_attempts': 0}
    rows = []
    guard = Guard()
    stop = 'all_20_boxes_completed'
    for index, (label, box) in enumerate(boxes):
        if counts['distinct_boxes_attempted'] >= BOX_LIMIT:
            stop = 'box_limit'
            break
        try:
            guard.check()
        except TimeLimit:
            stop = 'time_budget_exhausted'
            break
        counts['distinct_boxes_attempted'] += 1
        row = {'box_index': index, 'label': label, 'r': list(map(str, box[:2])), 'q': list(map(str, box[2:])),
               'center': [str((box[0]+box[1])/2), str((box[2]+box[3])/2)],
               'halfwidth': [str((box[1]-box[0])/2), str((box[3]-box[2])/2)], 'status': 'started'}
        box_start = monotonic()
        try:
            counts['old_bound_calls'] += 1
            old, tail, upper = old_lower(box)
            row['old_natural'] = {'lower': str(old), 'tail_lower': str(tail),
                                  'bump_upper_by_power': {str(p): str(x) for p, x in upper.items()}}
            guard.check()
            counts['new_bundle_calls'] += 1
            counts['center_bundles'] += 1
            counts['uniform_derivative_bundles'] += 1
            point, uniform = build_A(c, box, polys, tail_components, guard)
            guard.check()
            sing, sing_data = singular_lower(c, box)
            guard.check()
            rr, rq = (box[1]-box[0])/2, (box[3]-box[2])/2
            mean_A = point.v.lo-maxabs(uniform.g[0])*rr-maxabs(uniform.g[1])*rq
            second_A = (point.v.lo-maxabs(point.g[0])*rr-maxabs(point.g[1])*rq
                        -(maxabs(uniform.h[0][0])*rr*rr+2*maxabs(uniform.h[0][1])*rr*rq
                          +maxabs(uniform.h[1][1])*rq*rq)/2)
            mean, second = mean_A+sing, second_A+sing
            best = max(mean, second)
            row.update({'status': 'completed', 'smooth_center_value': point.v.pair(),
                        'smooth_center_gradient': [x.pair() for x in point.g],
                        'smooth_uniform_gradient': [x.pair() for x in uniform.g],
                        'smooth_uniform_hessian': [[x.pair() for x in h] for h in uniform.h],
                        'singular_remainder': sing_data, 'smooth_mean_lower': str(mean_A),
                        'smooth_second_order_lower': str(second_A), 'new_mean_lower': str(mean),
                        'new_second_order_lower': str(second), 'new_best_lower': str(best),
                        'new_minus_old': str(best-old), 'display': {'old': float(old), 'new_mean': float(mean),
                                                                  'new_second': float(second), 'new_best': float(best)}})
            counts['completed_boxes'] += 1
        except TimeLimit:
            row['status'] = 'incomplete_time_budget'
            stop = 'time_budget_exhausted'
        row['box_elapsed_seconds'] = monotonic()-box_start
        rows.append(row)
        print(json.dumps({'box': index+1, 'label': label, 'status': row['status'],
                          'display': row.get('display'), 'elapsed_seconds': monotonic()-guard.started}), flush=True)
        if row['status'] != 'completed':
            break
    elapsed = monotonic()-guard.started
    completed = [x for x in rows if x['status'] == 'completed']
    result = {'round': 17, 'status': 'finite_20_box_comparison_only' if len(completed) == BOX_LIMIT else 'incomplete_finite_comparison',
              'parameters': {'R': str(R), 'e_by_power': {str(p): str(e) for p, e in E.items()}},
              'limits': {'distinct_boxes': BOX_LIMIT, 'comparison_seconds': SECONDS},
              'timing': {'module_import_seconds': MODULE_IMPORT_SECONDS,
                         'one_time_import_and_exact_polynomial_construction_seconds': construction_seconds,
                         'comparison_seconds_including_all_box_computations': elapsed},
              'counts': counts, 'stop_reason': stop, 'boxes': rows,
              'summary': {'old_nonnegative_boxes': sum(Q(x['old_natural']['lower']) >= 0 for x in completed),
                          'new_nonnegative_boxes': sum(Q(x['new_best_lower']) >= 0 for x in completed),
                          'new_strictly_tighter_boxes': sum(Q(x['new_minus_old']) > 0 for x in completed)},
              'core_dependency': 'round008/certify.py', 'core_sha256': sha256(core_path.read_bytes()).hexdigest(),
              'old_bound_dependency': 'round016/certify.py', 'old_bound_sha256': sha256(old_path.read_bytes()).hexdigest(),
              'comparison_source_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
              'formula_sha256': sha256((ROOT/'centered_bounds.md').read_bytes()).hexdigest(),
              'scope': 'Strict lower inclusions on at most 20 prescribed boxes; no new grid, LP, adaptive certification, or new rigorous upper bound. No automatic retry after timeout.'}
    path.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'boxes'}, indent=2), flush=True)


if __name__ == '__main__':
    main()
