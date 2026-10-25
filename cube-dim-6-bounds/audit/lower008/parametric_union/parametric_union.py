#!/usr/bin/env python3
"""Certified affine-endpoint union polynomial on one supplied order interval.

Does not find events, optimize parameters, call MATLAB, or import D union code.
"""
from argparse import ArgumentParser
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
import hashlib
import importlib.util
import json
import os
import sys
import time

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[2]
NODE_CAP, WALL_CAP = 10**6, 1200
BOX_CAP, ENDPOINT_CAP = 10000, 100000


class InvalidCertificate(Exception):
    pass


def require(condition, message):
    if not condition:
        raise InvalidCertificate(message)


def rational(value):
    require(type(value) in (str, int) or isinstance(value, Q),
            'exact rational string/integer required; floats and booleans rejected')
    return Q(value)


class Engine:
    def __init__(self, arithmetic):
        spec = importlib.util.spec_from_file_location('audited_parametric_field', arithmetic)
        self.A = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.A)
        a = self.A
        a.D = 10**80
        a.RU, a.RV = a.floor_root(2*a.D**6, 6), a.floor_root(3*a.D**6, 6)
        a.UL, a.UH = Q(a.RU, a.D), Q(a.RU+1, a.D)
        a.VL, a.VH = Q(a.RV, a.D), Q(a.RV+1, a.D)
        require(a.UL**6 <= 2 < a.UH**6 and a.VL**6 <= 3 < a.VH**6,
                'root isolation failed')
        a.bounds.cache_clear()
        self.F, self.zero = a.F, a.F(0)
        self.started, self.nodes = time.monotonic(), 0

    def budget(self, node=False):
        self.nodes += int(node)
        require(self.nodes <= NODE_CAP, 'global million-node cap')
        require(time.monotonic()-self.started <= WALL_CAP, 'global 1200-second cap')

    def field(self, source):
        if isinstance(source, dict):
            coefficients = [Q(0)]*36
            for ij, value in source.items():
                i, j = map(int, ij.split(','))
                require(0 <= i < 6 and 0 <= j < 6, 'field index outside 0..5')
                coefficients[6*i+j] = rational(value)
        elif isinstance(source, list):
            require(len(source) == 36, 'dense field coefficient count must be 36')
            coefficients = list(map(rational, source))
        else:
            raise InvalidCertificate('field coefficient must be dictionary or 36-list')
        return self.F(a=coefficients)

    def sign(self, value):
        if value == self.zero:
            return 0
        lo, hi = self.A.bounds(value.a)
        if lo > 0:
            return 1
        if hi < 0:
            return -1
        raise InvalidCertificate('field sign unresolved; no tolerance merge permitted')

    def trim(self, p):
        p = tuple(p)
        while len(p) > 1 and p[-1] == self.zero:
            p = p[:-1]
        return p

    def add(self, p, q):
        return self.trim(tuple((p[k] if k < len(p) else self.zero)
                               +(q[k] if k < len(q) else self.zero)
                               for k in range(max(len(p), len(q)))))

    def subtract(self, p, q):
        return self.add(p, tuple(-x for x in q))

    def multiply(self, p, q):
        r = [self.zero]*(len(p)+len(q)-1)
        for i, a in enumerate(p):
            for j, b in enumerate(q):
                if a != self.zero and b != self.zero:
                    r[i+j] = r[i+j]+a*b
        return self.trim(r)

    def evaluate(self, p, t):
        answer = self.zero
        for coefficient in reversed(p):
            answer = answer*t+coefficient
        return answer

    def coefficients(self, value):
        return {f'{i},{j}': str(value.a[6*i+j]) for i in range(6) for j in range(6)
                if value.a[6*i+j]}

    def prepare(self, data):
        self.budget()
        dimension = data['dimension']
        require(type(dimension) is int and 1 <= dimension <= 6, 'dimension outside 1..6')
        left, right = map(rational, data['interval'])
        require(0 <= left < right, 'require finite rational 0<=L<U')
        check_t = rational(data.get('static_check_parameter', str(left+(right-left)/100)))
        require(left < check_t < right, 'static check must be strictly inside interval')
        endpoint_source = data['endpoints']
        require(len(endpoint_source) <= ENDPOINT_CAP, 'endpoint cap')
        endpoints = {}
        for key, source in endpoint_source.items():
            require(isinstance(key, str), 'endpoint ID must be string')
            require(set(source) <= {'c0', 'c1'}, 'only affine c0/c1 accepted')
            endpoints[key] = self.trim((self.field(source.get('c0', {})),
                                        self.field(source.get('c1', {}))))
        boxes = []
        require(len(data['boxes']) <= BOX_CAP, 'box cap')
        used = [set() for _ in range(dimension)]
        for source in data['boxes']:
            if isinstance(source, dict):
                source = source['intervals']
            if len(source) == 2*dimension and all(isinstance(x, str) for x in source):
                source = [source[2*h:2*h+2] for h in range(dimension)]
            require(len(source) == dimension, 'wrong box dimension')
            flat = []
            for h, pair in enumerate(source):
                require(len(pair) == 2, 'coordinate must have low/high IDs')
                for key in pair:
                    require(key in endpoints, 'unknown endpoint ID')
                    used[h].add(key)
                    flat.append(key)
            boxes.append(tuple(flat))
        certificate = data['order_certificate']
        require(len(certificate) == dimension, 'one order certificate per coordinate required')
        ranks, cuts, checked_gaps = [], [], 0
        boundary_ties = []
        for h, groups in enumerate(certificate):
            self.budget()
            ids = [key for group in groups for key in group]
            require(all(group for group in groups), 'empty tie group')
            require(set(ids) == used[h], 'certificate must partition all used axis IDs exactly')
            require(all(count == 1 for count in Counter(ids).values()), 'repeated endpoint ID in axis order')
            rank, axis_cuts = {}, []
            for k, group in enumerate(groups):
                representative = endpoints[group[0]]
                require(all(endpoints[key] == representative for key in group),
                        'tie group contains distinct affine functions')
                for key in group:
                    rank[key] = k
                axis_cuts.append(representative)
            for k, (low, high) in enumerate(zip(axis_cuts, axis_cuts[1:])):
                difference = self.subtract(high, low)
                sl = self.sign(self.evaluate(difference, left))
                sr = self.sign(self.evaluate(difference, right))
                require(sl >= 0 and sr >= 0, 'order certificate crosses or reverses inside interval')
                require(sl > 0 or sr > 0, 'identical affine groups must be exactly merged')
                if sl == 0 or sr == 0:
                    boundary_ties.append({'axis': h, 'adjacent_groups': [k, k+1],
                                          'left_equal': sl == 0, 'right_equal': sr == 0})
                checked_gaps += 1
            ranks.append(rank)
            cuts.append(tuple(axis_cuts))
        nonempty, empty, zero_width = [], 0, 0
        for box in boxes:
            ranked = tuple(ranks[h][box[2*h+k]] for h in range(dimension) for k in (0, 1))
            if any(ranked[2*h] > ranked[2*h+1] for h in range(dimension)):
                empty += 1
            elif any(ranked[2*h] == ranked[2*h+1] for h in range(dimension)):
                zero_width += 1
            else:
                nonempty.append(ranked)
        return {'dimension': dimension, 'left': left, 'right': right, 'check_t': check_t,
                'endpoints': endpoints, 'boxes': tuple(boxes), 'cuts': cuts,
                'ranked': tuple(sorted(set(nonempty))), 'empty_boxes': empty,
                'identically_zero_width_boxes': zero_width, 'checked_adjacent_gaps': checked_gaps,
                'boundary_ties': boundary_ties}

    def polynomial_union(self, prepared):
        dimension, cuts = prepared['dimension'], prepared['cuts']
        @lru_cache(None)
        def rec(boxes):
            self.budget(node=True)
            if not boxes:
                return (self.zero,)
            if not boxes[0]:
                return (self.F(1),)
            h = dimension-len(boxes[0])//2
            ranks = sorted({x for b in boxes for x in b[:2]})
            total = (self.zero,)
            for lo, hi in zip(ranks, ranks[1:]):
                child = tuple(sorted({b[2:] for b in boxes if b[0] <= lo and hi <= b[1]}))
                width = self.subtract(cuts[h][hi], cuts[h][lo])
                total = self.add(total, self.multiply(width, rec(child)))
            require(len(total)-1 <= dimension-h, 'polynomial degree exceeds dimension')
            return total
        return rec(prepared['ranked'])

    def static_union(self, prepared):
        """Independent direct F-coordinate union at one rational interior t."""
        t, endpoints = prepared['check_t'], prepared['endpoints']
        values = {key: self.evaluate(value, t) for key, value in endpoints.items()}
        boxes = []
        for source in prepared['boxes']:
            box = tuple(values[key] for key in source)
            if all(self.sign(box[2*h+1]-box[2*h]) > 0 for h in range(prepared['dimension'])):
                boxes.append(box)
        @lru_cache(None)
        def rec(active):
            self.budget(node=True)
            if not active:
                return self.zero
            if len(active[0]) == 2:
                ordered = sorted(active)
                lo, hi = ordered[0]
                total = self.zero
                for a, b in ordered[1:]:
                    if hi < a:
                        total, lo, hi = total+hi-lo, a, b
                    else:
                        hi = max(hi, b)
                return total+hi-lo
            cuts = sorted({x for b in active for x in b[:2]})
            total = self.zero
            for lo, hi in zip(cuts, cuts[1:]):
                child = tuple(sorted({b[2:] for b in active if b[0] <= lo and hi <= b[1]}))
                total = total+(hi-lo)*rec(child)
            return total
        return rec(tuple(sorted(set(boxes))))

    def certify(self, data, comparison_target=None):
        before = self.nodes
        prepared = self.prepare(data)
        polynomial = self.polynomial_union(prepared)
        static = self.static_union(prepared)
        at_check = self.evaluate(polynomial, prepared['check_t'])
        require(static == at_check, 'independent static union disagrees with polynomial')
        factor = rational(data.get('normalization_factor', '1'))
        require(factor > 0, 'normalization factor must be positive rational')
        candidate_t = rational(data.get('candidate_parameter', str(prepared['right'])))
        require(prepared['left'] <= candidate_t <= prepared['right'], 'candidate outside interval')
        candidate_value = self.evaluate(polynomial, candidate_t)*factor
        lo, hi = self.A.bounds(candidate_value.a)
        result = {'status': 'certified_interval_union_polynomial',
            'dimension': prepared['dimension'], 'interval': [str(prepared['left']), str(prepared['right'])],
            'degree': len(polynomial)-1,
            'volume_polynomial_coefficients': [self.coefficients(x) for x in polynomial],
            'normalization_factor': str(factor),
            'ratio_polynomial_coefficients': [self.coefficients(x*factor) for x in polynomial],
            'checked_adjacent_order_gaps': prepared['checked_adjacent_gaps'],
            'boundary_ties': prepared['boundary_ties'], 'input_boxes': len(prepared['boxes']),
            'empty_boxes_in_interval_interior': prepared['empty_boxes'],
            'identically_zero_width_boxes': prepared['identically_zero_width_boxes'],
            'retained_distinct_boxes': len(prepared['ranked']),
            'static_check_parameter': str(prepared['check_t']), 'static_check_exact_match': True,
            'static_volume_coefficients': self.coefficients(static),
            'candidate_parameter': str(candidate_t), 'candidate_ratio_coefficients': self.coefficients(candidate_value),
            'candidate_ratio_interval': [str(lo), str(hi)],
            'union_nodes_used': self.nodes-before,
            'order_event_search_performed': False, 'optimization_performed': False,
            'C6_status': 'pending source-geometry, finite-L1 and independent review'}
        if comparison_target is not None:
            result.update({'previous_complete_Rideal_upper': str(comparison_target),
                'candidate_lower_strictly_above_previous_complete_upper': lo > comparison_target,
                'candidate_lower_minus_previous_complete_upper': str(lo-comparison_target)})
        return result


def toy_input(boxes, interval=('0', '1/4')):
    """Finite known toys supply their own order certificate; verification is separate."""
    dimension = len(boxes[0])//2
    endpoints, encoded = {}, []
    for k, box in enumerate(boxes):
        names = []
        for j, (c0, c1) in enumerate(box):
            key = f'b{k}e{j}'
            endpoints[key] = {'c0': {'0,0': str(c0)} if c0 else {},
                              'c1': {'0,0': str(c1)} if c1 else {}}
            names.append(key)
        encoded.append([names[2*h:2*h+2] for h in range(dimension)])
    midpoint = sum(map(Q, interval), Q(0))/2
    certificates = []
    for h in range(dimension):
        groups = {}
        for k, box in enumerate(boxes):
            for j in (2*h, 2*h+1):
                groups.setdefault(tuple(map(Q, box[j])), []).append(f'b{k}e{j}')
        order = sorted(groups, key=lambda p: (p[0]+midpoint*p[1], p))
        certificates.append([groups[p] for p in order])
    return {'dimension': dimension, 'interval': list(interval), 'endpoints': endpoints,
            'boxes': encoded, 'order_certificate': certificates}


def selfcheck(engine):
    c = lambda n: (Q(n), Q(0))
    a = lambda n, m: (Q(n), Q(m))
    toys = [
        ('intersecting intervals', [(c(0), a(1, 1)), (a(1, -1), c(2))], ('0', '1/4'), [2]),
        ('coincident 2D boxes', [(c(0), a(1, 1), c(0), c(1))]*2, ('0', '1/4'), [1, 1]),
        ('inverted and identically zero-width empty boxes',
         [(c(0), c(1)), (a(1, 1), c(0)), (c(2), c(2))], ('0', '1/4'), [1]),
        ('one-width birth', [(c(0), c(1)), (c(2), a(2, 1))], ('0', '1/4'), [1, 1]),
        ('two-width birth', [(c(0), a(0, 1), c(0), a(0, 2))], ('0', '1/4'), [0, 0, 2]),
        ('six-dimensional degree-six cube', [tuple([c(0), a(1, 1)]*6)], ('0', '1/4'),
         [1, 6, 15, 20, 15, 6, 1]),
        ('order event on right boundary', [(c(0), a(1, 1)), (c(2), c(3))], ('0', '1'), [2, 1]),
        ('interior crossing must be rejected', [(c(0), a(1, 1)), (c(2), c(3))], ('0', '2'), None),
    ]
    output = []
    for name, boxes, interval, expected in toys:
        data = toy_input(boxes, interval)
        try:
            got = engine.certify(data)
            if expected is None:
                raise AssertionError('interior crossing certificate was incorrectly accepted')
            actual = got['volume_polynomial_coefficients']
            want = [engine.coefficients(engine.F(n)) for n in expected]
            require(actual == want, 'toy analytic polynomial disagreement')
            output.append({'name': name, 'status': 'pass', 'result': got})
        except InvalidCertificate as exc:
            if expected is not None:
                raise
            output.append({'name': name, 'status': 'rejected_as_expected', 'reason': str(exc)})
    return {'status': 'pass', 'scope': 'eight rational analytic interface toys only',
            'toy_count': len(toys), 'real_candidate_count': 0, 'new_C6_bound': False, 'results': output}


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--arithmetic', type=Path, default=PROJECT/'lower/round001_rectangular_certificate.py')
    parser.add_argument('--selfcheck', action='store_true')
    parser.add_argument('--input', type=Path)
    parser.add_argument('--comparison-results', type=Path)
    args = parser.parse_args()
    require(args.selfcheck != bool(args.input), 'choose exactly one mode')
    if args.input:
        require(args.comparison_results is not None, 'real candidate requires complete E007 comparison results')
    engine = Engine(args.arithmetic)
    inputs = [args.arithmetic]+([args.input, args.comparison_results] if args.input else [])
    hashes = {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs}
    try:
        if args.selfcheck:
            result = selfcheck(engine)
        else:
            comparison = json.loads(args.comparison_results.read_text())
            target = Q(comparison['exact_interval'][1])
            result = engine.certify(json.loads(args.input.read_text()), target)
            result['real_candidate_count'] = 1
            result['new_C6_bound'] = False
    except (InvalidCertificate, ArithmeticError, KeyError, ValueError, TypeError, AssertionError) as exc:
        result = {'status': 'stopped', 'reason': repr(exc), 'new_C6_bound': False}
    result.update({'pid': os.getpid(), 'input_sha256': hashes,
                   'root_denominator': str(engine.A.D),
                   'limits': {'union_nodes': NODE_CAP, 'wall_seconds': WALL_CAP,
                              'boxes': BOX_CAP, 'endpoints': ENDPOINT_CAP},
                   'global_union_nodes': engine.nodes,
                   'elapsed_seconds': time.monotonic()-engine.started,
                   'code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   'actual_process_status': 'finished', 'no_MATLAB': True})
    output = HERE/('toy_results.json' if args.selfcheck else 'candidate_results.json')
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ('results', 'volume_polynomial_coefficients', 'ratio_polynomial_coefficients',
                                   'candidate_ratio_interval', 'candidate_lower_minus_previous_complete_upper',
                                   'previous_complete_Rideal_upper')}, indent=2))


if __name__ == '__main__':
    main()
