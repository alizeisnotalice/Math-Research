#!/usr/bin/env python3
"""One independent saved-certificate audit; no geometric sampling or union.

Algebra and sixth-root enclosures below are independent of author field.py.
No source program is imported. All source files are read-only and hash-pinned.
"""
from pathlib import Path
from fractions import Fraction as Q
from math import prod, comb, factorial
import hashlib
import json
import resource
import signal
import sys
import time

sys.dont_write_bytecode = True
H = Path(__file__).resolve().parent
WALL, CPU = time.monotonic(), time.process_time()
resource.setrlimit(resource.RLIMIT_CPU, (15, 15))
signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(TimeoutError('wall limit')))
signal.alarm(165)
checks = 0
counts = {'source_hashes': 0, 'swap_targets': 0, 'formal_endpoint_catalogs': 0,
          'saved_gap_equalities_and_signs': 0, 'saved_width_identities': 0,
          'inherited_derivatives': 0, 'new_exact_parameter_configurations': 0,
          'new_union_cells': 0, 'new_integration_states': 0,
          'author_scientific_program_runs': 0, 'old_science_replays': 0}
out = {'status': 'running', 'scope': 'saved certificate only; conditional on E089 acceptance'}

def ck(ok, message):
    global checks
    checks += 1
    if not ok:
        raise AssertionError(message)
    if time.process_time() - CPU >= 13.5:
        raise TimeoutError('CPU prelimit')

Z = (Q(0),) * 36
def scalar(q):
    return (Q(q),) + Z[1:]
def poly(d):
    ans = list(Z)
    for ij, c in d.items():
        i, j = map(int, ij.split(','))
        ck(0 <= i < 6 and 0 <= j < 6, 'coefficient index')
        ans[6*i+j] = Q(c)
    return tuple(ans)
def add(a, b):
    return tuple(x+y for x, y in zip(a, b))
def neg(a):
    return tuple(-x for x in a)
def sub(a, b):
    return add(a, neg(b))
def scale(a, q):
    return tuple(x*q for x in a)
def radius(label):
    n = prod(label)
    exps = []
    for prime in (2, 3):
        p = 0
        while n % prime == 0:
            n //= prime
            p += 1
        exps.append(p)
    ck(n == 1, 'radius prime factors')
    i, j = exps
    ans = list(Z)
    ans[6*(i % 6)+(j % 6)] = Q(2**(i//6)*3**(j//6), 2)
    return tuple(ans)

def sixth_root_interval(n):
    # Different denominator from the author: exact integer bisection.
    D = 10**80
    lo, hi, target = D, 2*D, n*D**6
    while hi-lo > 1:
        m = (lo+hi)//2
        if m**6 <= target:
            lo = m
        else:
            hi = m
    ck(lo**6 <= target < hi**6 and hi == lo+1, 'sixth-root enclosure')
    return Q(lo, D), Q(hi, D)

UL, UH = sixth_root_interval(2)
VL, VH = sixth_root_interval(3)
PB = [(UL**i*VL**j, UH**i*VH**j) for i in range(6) for j in range(6)]
def bounds(a):
    lo = hi = Q(0)
    for c, (l, u) in zip(a, PB):
        lo += c*(l if c >= 0 else u)
        hi += c*(u if c >= 0 else l)
    return lo, hi
def nonnegative(a, message):
    ck(a == Z or bounds(a)[0] > 0, message)
def positive(a, message):
    ck(bounds(a)[0] > 0, message)

def load(path):
    return json.loads(Path(path).read_text())

try:
    gate = load(H/'source_gate.json')
    for path, expected in gate['source_pins'].items():
        ck(hashlib.sha256(Path(path).read_bytes()).hexdigest() == expected, 'source pin '+path)
        counts['source_hashes'] += 1
    author = Path(next(p for p in gate['source_pins'] if p.endswith('/round097/full_spacing_local_max/manifest.json'))).parent
    manifest, r, reg = [load(author/x) for x in ('manifest.json', 'results.json', 'registration.json')]
    ck(manifest['status'] == 'sealed_complete_pending_parent_acceptance_and_E089_dependency', 'author seal status')
    for ent in manifest['artifacts']:
        p = author/ent['path']
        ck(p.stat().st_size == ent['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest() == ent['sha256'], 'artifact '+ent['path'])
    for ent in manifest['frozen_dependencies']:
        ck(gate['source_pins'][ent['path']] == ent['sha256'], 'manifest dependency')
    paths = list(reg['source_pins'])
    labels = [tuple(x) for x in load(paths[0])['box_lengths']]
    ck(len(labels) == len(set(labels)) == 127, 'FULL127 unique labels')
    label_index = {x: k for k, x in enumerate(labels)}
    swaps = r['full127_S5_closure']['adjacent_swap_targets']
    ck(len(swaps) == 508 and {(k, j) for k, j, _ in swaps} == {(k, j) for k in range(127) for j in range(1, 5)}, 'swap coverage')
    for k, j, target in swaps:
        x = list(labels[k]); x[j], x[j+1] = x[j+1], x[j]
        ck(label_index.get(tuple(x)) == target and prod(labels[k]) == prod(labels[target]), 'swap target/radius')
        counts['swap_targets'] += 1
    covered = []
    for orbit in r['full127_S5_closure']['orbits']:
        tail = tuple(orbit['sorted_tail_lengths'])
        size = factorial(5)//prod(factorial(tail.count(x)) for x in set(tail))
        ids = orbit['raw_indices']; covered.extend(ids)
        ck(len(ids) == orbit['orbit_size'] == size, 'orbit cardinality')
        ck(all(labels[k][0] == orbit['first_length'] and tuple(sorted(labels[k][1:])) == tail for k in ids), 'orbit class')
    ck(sorted(covered) == list(range(127)), 'orbit partition')
    e79, r89 = load(paths[1]), load(paths[3])
    base = list(map(poly, r['base_spacing_coefficients']))
    ck(base[0] == poly(r89['a1']['coefficients']), 'a1 handoff')
    ck(all(x == poly(e79['s9']['coefficients']) == poly(r89['tail_spacing']['coefficients']) for x in base[1:]), 'tail handoff')
    ck(poly(e79['inherited_E075_complete_ideal']['coefficients']) == poly(r89['baseline']['coefficients']) == poly(r['inherited_ideal']['coefficients']), 'ideal handoff')
    inherited = {'a1_left': r89['sides']['left']['exact_derivative'],
                 'a1_right': r89['sides']['right']['exact_derivative'],
                 'tail_left': e79['left_slope'], 'tail_right': e79['right_slope']}
    slope_bounds = {}
    for name, dep in inherited.items():
        p = poly(r['inherited_derivatives'][name]['coefficients'])
        ck(p == poly(dep['coefficients']), 'derivative coefficient handoff '+name)
        lo, hi = bounds(p)
        ck(lo > Q(1, 25) if name.endswith('left') else hi < -Q(1, 25), 'strict slope margin '+name)
        # Independent rational enclosure; no decimal values enter a decision.
        slope_bounds[name] = {'lower': str(lo), 'upper': str(hi)}
        counts['inherited_derivatives'] += 1
    delta = Q(r['explicit_neighborhood']['delta'])
    ck(delta == Q(1, 10**17), 'delta')
    catalogs = []
    for axis, chamber in enumerate(r['endpoint_chambers']):
        fs = [(poly(x['c']), Q(x['m'])) for x in chamber['functions']]
        ck(len(fs) == len(set(fs)), 'distinct formal endpoint functions')
        expected = {(Z, Q(0)), (Z, Q(1, 2))}
        # Structural catalogue only: no evaluated geometric configuration.
        for ls in labels:
            rd = radius(ls); ell = ls[axis]
            expected.add((rd, -Q(ell-1, 2)) if ell % 2 else (neg(rd), Q(ell, 2)))
        ck(set(fs) == expected, 'formal endpoint catalog '+str(axis+1))
        ck(all(-1 <= m <= 2 for _, m in fs), 'endpoint velocity bound')
        catalogs.append(fs)
        counts['formal_endpoint_catalogs'] += 1
        orders = []
        for side, key in ((-1, 'minus'), (1, 'plus')):
            cert = chamber['sides'][key]
            order = cert['complete_order']; orders.append(order)
            ck(sorted(order) == list(range(len(fs))), 'complete endpoint permutation')
            ck([g['ids'] for g in cert['adjacent_gap_certificate']] == [list(p) for p in zip(order, order[1:])], 'adjacent gap coverage')
            for gap in cert['adjacent_gap_certificate']:
                i, j = gap['ids']
                c = sub(fs[j][0], fs[i][0]); dm = fs[j][1]-fs[i][1]
                g0 = add(c, scale(base[axis], dm))
                gd = add(g0, scalar(side*delta*dm))
                ck(g0 == poly(gap['at_center']) and gd == poly(gap['at_delta_endpoint']), 'saved affine gap identities')
                nonnegative(g0, 'center weak gap')
                positive(gd, 'single-side strict gap endpoint')
                counts['saved_gap_equalities_and_signs'] += 1
        # The two one-sided orders both give the same weak center order.
        # Recover selected widths by ranks, checking identities already saved.
        rank = {k: pos for pos, k in enumerate(orders[0])}
        zero, cap = fs.index((Z, Q(0))), fs.index((Z, Q(1, 2)))
        ax = r['configurations'][0]['per_axis'][axis]
        widths = list(map(poly, ax['all127_width_coefficients']))
        ck(len(widths) == 127, 'saved width coverage')
        for k, ls in enumerate(labels):
            rd, ell = radius(ls), ls[axis]
            f = (rd, -Q(ell-1, 2)) if ell % 2 else (neg(rd), Q(ell, 2))
            mov = fs.index(f)
            low = zero if ell % 2 else max((zero, mov), key=rank.get)
            high = min((cap, mov), key=rank.get) if ell % 2 else cap
            w = add(sub(fs[high][0], fs[low][0]), scale(base[axis], fs[high][1]-fs[low][1]))
            ck(w == widths[k], 'saved width identity')
            counts['saved_width_identities'] += 1
        minimum = poly(ax['minimum_width']['coefficients'])
        ck(minimum in widths, 'saved width minimum member')
        for w in widths:
            nonnegative(sub(w, minimum), 'saved minimum comparison')
        positive(sub(minimum, scalar(3*delta)), 'uniform positive width margin')
        positive(base[axis], 'positive base spacing')
        positive(sub(scalar(1), scale(base[axis], Q(1, 2))), 'half-cell width below one')
    ck([len(x) for x in catalogs] == [24, 37, 37, 37, 37, 37], 'catalog sizes')
    ck(all(set(x) == set(catalogs[1]) for x in catalogs[2:]), 'tail endpoint catalog identity')
    ck(counts['saved_gap_equalities_and_signs'] == 406, 'all saved gaps verified')
    N = prod(len(fs)-1 for fs in catalogs)
    S = sum(comb(6, k)*3**k for k in range(2, 7))
    B = 64*N*S
    e = r['explicit_neighborhood']
    ck(N == 23*36**5 == e['six_dimensional_elementary_cell_count_upper'], 'cell bound')
    ck(S == 4077 and B == 362878322540544 == e['uniform_remainder_B'], 'remainder constant')
    ck(B*delta < Q(1, 50) and Q(1, 25)-B*delta > Q(1, 50), 'explicit strict decrease')
    ck(Q(e['certified_decrease_coefficient']) == Q(1, 25)-B*delta, 'saved decrease coefficient')
    for path, expected in gate['source_pins'].items():
        ck(hashlib.sha256(Path(path).read_bytes()).hexdigest() == expected, 'post-run source immutability')
    out.update(status='PASS_saved_certificate_conditional_E089', independent_root_denominator='10^80',
               source_manifest_sha256=gate['author_manifest_sha256'], slope_bounds=slope_bounds,
               delta=str(delta), N=N, B=B, B_delta=str(B*delta),
               certified_decrease=str(Q(1,25)-B*delta), simple_decrease='1/50',
               mathematical_verdict='Full six-spacing theorem follows from the verified saved structure and accepted 079/089 derivatives; E089 acceptance remains a parent dependency.')
except Exception as exc:
    out.update(status='FAIL_first_failure', error=repr(exc))
finally:
    out.update(checks=checks, counts=counts, CPU_seconds_before_save=time.process_time()-CPU,
               wall_seconds_before_save=time.monotonic()-WALL,
               self_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               source_gate_sha256=hashlib.sha256((H/'source_gate.json').read_bytes()).hexdigest())
    with (H/'review_results.json').open('x') as f:
        json.dump(out, f, ensure_ascii=False, indent=2); f.write('\n')
    print(json.dumps({k: out[k] for k in ('status', 'checks', 'counts', 'CPU_seconds_before_save', 'wall_seconds_before_save')}))
    signal.alarm(0)
if out['status'] != 'PASS_saved_certificate_conditional_E089':
    sys.exit(1)
