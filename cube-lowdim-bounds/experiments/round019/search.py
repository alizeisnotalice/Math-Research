#!/usr/bin/env python3
"""One bounded, numerical four-atom experiment. Importing performs no evaluation.

The main entry point is deliberately separate: source review is not permission
to run it. A result file prevents accidental repetition of the finite run.
"""
from __future__ import annotations

from time import monotonic
_IMPORT_START = monotonic()
from hashlib import sha256
import itertools
import json
import math
import os
from pathlib import Path
import platform
import sys
import tempfile

import numpy as np
import scipy
from scipy.optimize import differential_evolution
_IMPORT_SECONDS = monotonic() - _IMPORT_START

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
SEED = 2026102419
BOUNDS = [(-.2, .2), (-.2, .2), (-.15, .15), (-.15, .15), (-.3, .3)]
ZERO = (0., 0., 0., 0., 0.)
NONZERO_CONTROL = (.1, -.1, .1, -.1, .2)
HARD_TOTAL = 1000
DE_LIMIT = 640
SEARCH_SECONDS = 45.
TOTAL_SECONDS = 60.
GAIN_TOL = 1e-9
CHECK_TOL = 1e-8


class StopBudget(Exception):
    pass


class ControlFailure(Exception):
    pass


class Clock:
    def __init__(self):
        self.start = monotonic()

    def elapsed(self):
        return monotonic() - self.start

    def check(self, stage):
        limit = SEARCH_SECONDS if stage == 'search' else TOTAL_SECONDS
        if self.elapsed() >= limit:
            raise StopBudget(f'{stage}_time_limit')


def lin_parameters():
    u = (2. + math.sqrt(22.)) / 3.
    h = (1. + u) / 2.
    return h, h + 1., u*u - 1.


def atom_patch(z, base, clock, stage):
    """Aggregate equal positions; deduplicating coordinates never drops mass."""
    clock.check(stage)
    xi0, xi1, eta0, eta1, beta = map(float, z)
    h, v, w = base
    px, py = 2*h, 2*v
    seeds = [(h/2, v/2, 1.),
             (h*(1.5+xi0), v*(.5+eta0), w*(1+beta)),
             (h/2, 1.5*v, 1.),
             (h*(1.5+xi1), v*(1.5+eta1), w*(1-beta))]
    seed_mass = math.fsum(t[2] for t in seeds)
    cell_mass = 2*(1+w)
    if abs(seed_mass-cell_mass) > 1e-12*cell_mass:
        raise ControlFailure('Complete period-mass identity did not reproduce.')
    merged = {}
    for kx, ky in itertools.product((-1, 0, 1), repeat=2):
        clock.check(stage)
        for x, y, mass in seeds:
            merged.setdefault((x+kx*px, y+ky*py), []).append(mass)
    points = [(x, y, math.fsum(masses)) for (x, y), masses in sorted(merged.items())]
    return points, px, py, cell_mass


def wrap_interval(lo, hi, period):
    width = hi-lo
    if width <= 0:
        return []
    if width >= period:
        return [(0., period)]
    left = lo % period
    right = left+width
    if right <= period:
        return [(left, right)]
    return [(left, period), (0., right-period)]


def discovery_boxes(z, base, clock, stage):
    points, px, py, mass = atom_patch(z, base, clock, stage)
    xs = sorted({p[0] for p in points})
    ys = sorted({p[1] for p in points})
    # Membership is based on all actual points in the endpoint rectangle.
    # Bit masks only avoid repeatedly constructing the same closure.
    x_masks = []
    y_masks = []
    for lo, hi in itertools.combinations_with_replacement(xs, 2):
        clock.check(stage)
        x_masks.append(sum(1 << i for i, p in enumerate(points) if lo <= p[0] <= hi))
    for lo, hi in itertools.combinations_with_replacement(ys, 2):
        clock.check(stage)
        y_masks.append(sum(1 << i for i, p in enumerate(points) if lo <= p[1] <= hi))
    seen, boxes = set(), set()
    nonempty_closures = positive_boxes = 0
    for xm in x_masks:
        clock.check(stage)
        for ym in y_masks:
            clock.check(stage)
            selected_mask = xm & ym
            if not selected_mask or selected_mask in seen:
                continue
            seen.add(selected_mask)
            chosen = []
            bits = selected_mask
            while bits:
                bit = bits & -bits
                chosen.append(points[bit.bit_length()-1])
                bits -= bit
            nonempty_closures += 1
            radius = math.sqrt(math.fsum(p[2] for p in chosen))/2
            xlo = max(p[0] for p in chosen)-radius
            xhi = min(p[0] for p in chosen)+radius
            ylo = max(p[1] for p in chosen)-radius
            yhi = min(p[1] for p in chosen)+radius
            if xhi <= xlo or yhi <= ylo:
                continue
            positive_boxes += 1
            for xp, yp in itertools.product(wrap_interval(xlo, xhi, px),
                                            wrap_interval(ylo, yhi, py)):
                boxes.add((xp[0], xp[1], yp[0], yp[1]))
    clock.check(stage)
    return sorted(boxes), mass, {
        'patch_atoms_before_aggregation': 36, 'patch_positions_after_aggregation': len(points),
        'unique_x_coordinates': len(xs), 'unique_y_coordinates': len(ys),
        'endpoint_rectangles': len(x_masks)*len(y_masks),
        'distinct_nonempty_closures': nonempty_closures, 'positive_unwrapped_boxes': positive_boxes,
        'unique_wrapped_boxes': len(boxes), 'periods': [px, py]}


def sweep_union(boxes, clock, stage):
    """True rectangle union by x events and a y-coverage segment tree."""
    clock.check(stage)
    if not boxes:
        return 0.
    ys = sorted({y for b in boxes for y in b[2:]})
    lookup = {y: i for i, y in enumerate(ys)}
    events = []
    for xlo, xhi, ylo, yhi in boxes:
        clock.check(stage)
        a, b = lookup[ylo], lookup[yhi]
        events.extend(((xlo, 1, a, b), (xhi, -1, a, b)))
    events.sort()
    n = len(ys)-1
    count, length = [0]*(4*n+8), [0.]*(4*n+8)

    def update(node, left, right, a, b, delta):
        if a <= left and right <= b:
            count[node] += delta
        else:
            mid = (left+right)//2
            if a < mid:
                update(2*node, left, mid, a, b, delta)
            if b > mid:
                update(2*node+1, mid, right, a, b, delta)
        if count[node] > 0:
            length[node] = ys[right]-ys[left]
        elif right-left == 1:
            length[node] = 0.
        else:
            length[node] = length[2*node]+length[2*node+1]

    contributions = []
    previous_x = events[0][0]
    for x, delta, a, b in events:
        clock.check(stage)
        contributions.append((x-previous_x)*length[1])
        update(1, 0, n, a, b, delta)
        previous_x = x
    area = math.fsum(contributions)
    clock.check(stage)
    return area


def independent_boxes(z, base, clock):
    """Independent raw-atom selection and integer-translate periodic clipping.

    Does not use discovery atom_patch, bit masks, or wrap_interval.
    """
    stage = 'audit'
    h, v, w = base
    px, py = 2*h, 2*v
    xi0, xi1, eta0, eta1, beta = map(float, z)
    locations = [(h*.5, v*.5, 1.),
                 (h*(1.5+xi0), v*(.5+eta0), w*(1.+beta)),
                 (h*.5, v*1.5, 1.),
                 (h*(1.5+xi1), v*(1.5+eta1), w*(1.-beta))]
    raw_points = []
    for j in (-1, 0, 1):
        for i in (-1, 0, 1):
            clock.check(stage)
            for x, y, weight in locations:
                raw_points.append((x+i*px, y+j*py, weight))
    # Keep coincident point records: direct fsum below includes every mass.
    xs = sorted({p[0] for p in raw_points})
    ys = sorted({p[1] for p in raw_points})
    boxes = set()

    def clip_translations(lo, hi, period):
        out = []
        for shift in range(math.floor(-hi/period)-1, math.ceil((period-lo)/period)+2):
            clock.check(stage)
            a, b = max(0., lo+shift*period), min(period, hi+shift*period)
            if b > a:
                out.append((a, b))
        return out

    for left, right in itertools.combinations_with_replacement(xs, 2):
        clock.check(stage)
        for bottom, top in itertools.combinations_with_replacement(ys, 2):
            clock.check(stage)
            chosen = [p for p in raw_points if left <= p[0] <= right and bottom <= p[1] <= top]
            if not chosen:
                continue
            side = math.sqrt(math.fsum(p[2] for p in chosen))
            x0, x1 = max(p[0] for p in chosen)-side*.5, min(p[0] for p in chosen)+side*.5
            y0, y1 = max(p[1] for p in chosen)-side*.5, min(p[1] for p in chosen)+side*.5
            if x1 <= x0 or y1 <= y0:
                continue
            for xp in clip_translations(x0, x1, px):
                for yp in clip_translations(y0, y1, py):
                    clock.check(stage)
                    boxes.add((*xp, *yp))
    clock.check(stage)
    return sorted(boxes), math.fsum(p[2] for p in locations)


def endpoint_cell_union(boxes, clock):
    """Independent endpoint-cell coverage; no segment tree or discovery boxes."""
    stage = 'audit'
    clock.check(stage)
    if not boxes:
        return 0.
    ys = sorted({y for b in boxes for y in b[2:]})
    index = {y: i for i, y in enumerate(ys)}
    changes = {}
    for x0, x1, y0, y1 in boxes:
        clock.check(stage)
        a, b = index[y0], index[y1]
        changes.setdefault(x0, []).append((a, b, 1))
        changes.setdefault(x1, []).append((a, b, -1))
    ywidth = np.diff(ys)
    coverage = np.zeros(len(ywidth), dtype=np.int64)
    pieces = []
    cuts = sorted(changes)
    for i, x in enumerate(cuts):
        clock.check(stage)
        for a, b, delta in changes[x]:
            clock.check(stage)
            coverage[a:b] += delta
        if i+1 < len(cuts):
            height = math.fsum(ywidth[coverage > 0])
            pieces.append((cuts[i+1]-x)*height)
    clock.check(stage)
    return math.fsum(pieces)


def atomic_json(path, value):
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    temporary.replace(path)


def main():
    result_path = ROOT/'result.json'
    if result_path.exists():
        raise SystemExit('result.json exists: refusing another finite run.')
    clock = Clock()  # Includes preparation, reference construction, controls, and search.
    counts = {k: {'started': 0, 'completed': 0, 'interrupted': 0, 'failed': 0}
              for k in ('controls', 'DE', 'audits')}
    trace_handle = None
    trace_path = None
    best = best_nonzero = best_de = None
    controls, checks = [], []
    generation_count = 0
    opt = None
    search_stop = None
    fatal = None
    base = None
    phi_record = None
    prep_seconds = search_end = None
    result = {'round': 19, 'status': 'started', 'seed': SEED,
              'scope': 'Finite numerical discovery only; no rigorous new lower bound.'}
    atomic_json(result_path, result)

    def total_started():
        return sum(c['started'] for c in counts.values())

    def begin(kind, stage):
        clock.check(stage)
        if total_started() >= HARD_TOTAL:
            raise StopBudget('hard_total_started_limit')
        limit = {'controls': 2, 'DE': DE_LIMIT, 'audits': 2}[kind]
        if counts[kind]['started'] >= limit:
            raise StopBudget(f'{kind}_started_limit')
        counts[kind]['started'] += 1
        return total_started()

    def emit(value):
        trace_handle.write(json.dumps(value, ensure_ascii=False, separators=(',', ':'), allow_nan=False)+'\n')
        trace_handle.flush()

    def evaluate(z, label, kind):
        nonlocal best, best_nonzero, best_de
        z = [float(x) for x in z]
        if len(z) != 5 or any(not a <= x <= b for x, (a, b) in zip(z, BOUNDS)):
            raise ControlFailure('Out-of-range candidate.')
        number = begin(kind, 'search')
        row = {'configuration': number, 'kind': kind, 'label': label, 'z': z,
               'started_at_seconds': clock.elapsed(), 'status': 'started'}
        emit(row)
        try:
            boxes, mass, stats = discovery_boxes(z, base, clock, 'search')
            area = sweep_union(boxes, clock, 'search')
            clock.check('search')
            if not math.isfinite(area) or area < 0 or area > stats['periods'][0]*stats['periods'][1]+1e-9:
                raise ControlFailure('Invalid union area.')
            row.update(status='completed', score=area/mass, union_area=area,
                       period_mass=mass, geometry=stats, completed_at_seconds=clock.elapsed())
            counts[kind]['completed'] += 1
            if best is None or row['score'] > best['score']:
                best = row
            if any(x != 0 for x in z) and (best_nonzero is None or row['score'] > best_nonzero['score']):
                best_nonzero = row
            if kind == 'DE' and (best_de is None or row['score'] > best_de['score']):
                best_de = row
            emit(row)
            if sum(counts[k]['completed'] for k in ('controls', 'DE')) % 20 == 0:
                print(json.dumps({'elapsed': clock.elapsed(), 'counts': counts,
                                  'best_score': best['score']}), flush=True)
            return row
        except StopBudget as exc:
            counts[kind]['interrupted'] += 1
            emit(dict(row, status='interrupted', reason=str(exc), ended_at_seconds=clock.elapsed()))
            raise
        except Exception as exc:
            counts[kind]['failed'] += 1
            emit(dict(row, status='failed', reason=f'{type(exc).__name__}: {exc}',
                      ended_at_seconds=clock.elapsed()))
            raise

    try:
        descriptor, trace_name = tempfile.mkstemp(prefix='cube-lowdim-round019-', suffix='.jsonl')
        trace_path = Path(trace_name)
        trace_handle = os.fdopen(descriptor, 'w')
        from phi_reference import reference
        phi = reference()
        phi_record = {'lower': str(phi.lo), 'upper': str(phi.hi),
                      'float_display': float(phi.hi), 'source': 'phi_reference.py',
                      'source_sha256': sha256((ROOT/'phi_reference.py').read_bytes()).hexdigest()}
        base = lin_parameters()
        clock.check('search')
        prep_seconds = clock.elapsed()
        controls.append(evaluate(ZERO, 'Lin_zero_supercell_control', 'controls'))
        if controls[0]['score'] < float(phi.lo)-GAIN_TOL:
            raise ControlFailure('Zero Lin control did not recover public Phi coverage.')
        controls.append(evaluate(NONZERO_CONTROL, 'fixed_four_level_nonrank1_control', 'controls'))

        def objective(z):
            return -evaluate(z, 'five_variable_DE', 'DE')['score']

        def callback(xk, convergence):
            nonlocal generation_count
            generation_count += 1
            clock.check('search')
            return False

        opt = differential_evolution(objective, BOUNDS, popsize=8, maxiter=15,
                                     polish=False, seed=SEED, workers=1, updating='immediate',
                                     tol=0., atol=0., x0=ZERO, callback=callback)
        search_stop = 'optimizer_returned'
    except StopBudget as exc:
        search_stop = str(exc)
    except Exception as exc:
        fatal = f'{type(exc).__name__}: {exc}'
        search_stop = 'failure'
    finally:
        search_end = clock.elapsed()

    # Fixed reconstruction only: never evaluates a new candidate or resumes DE.
    if fatal is None and base is not None:
        targets = []
        if controls:
            targets.append(('Lin_zero_supercell_control', controls[0]))
        if best_nonzero is not None:
            targets.append(('best_nonzero_completed', best_nonzero))
        for label, row in targets:
            try:
                number = begin('audits', 'audit')
            except StopBudget as exc:
                checks.append({'label': label, 'status': 'not_started', 'reason': str(exc)})
                break
            check = {'configuration': number, 'label': label, 'source_configuration': row['configuration'],
                     'z': row['z'], 'status': 'started', 'started_at_seconds': clock.elapsed()}
            emit(check)
            try:
                boxes, mass = independent_boxes(row['z'], base, clock)
                area = endpoint_cell_union(boxes, clock)
                clock.check('audit')
                difference = abs(area/mass-row['score'])
                mass_difference = abs(mass-row['period_mass'])
                mass_tolerance = 1e-12*row['period_mass']
                check.update(status='completed', independent_box_count=len(boxes),
                             independently_rebuilt_mass=mass, endpoint_cell_area=area,
                             endpoint_cell_score=area/mass, recorded_score=row['score'],
                             absolute_score_difference=difference, tolerance=CHECK_TOL,
                             mass_difference=mass_difference, mass_tolerance=mass_tolerance,
                             passed=difference <= CHECK_TOL and mass_difference <= mass_tolerance,
                             completed_at_seconds=clock.elapsed())
                counts['audits']['completed'] += 1
                checks.append(check)
                emit(check)
                if not check['passed']:
                    fatal = 'Independent geometry comparison failed.'
                    break
            except StopBudget as exc:
                counts['audits']['interrupted'] += 1
                check.update(status='interrupted', reason=str(exc), ended_at_seconds=clock.elapsed())
                checks.append(check)
                emit(check)
                break
            except Exception as exc:
                counts['audits']['failed'] += 1
                check.update(status='failed', reason=f'{type(exc).__name__}: {exc}',
                             ended_at_seconds=clock.elapsed())
                checks.append(check)
                emit(check)
                fatal = check['reason']
                break

    numerical_seconds = clock.elapsed()
    if trace_handle is not None:
        trace_handle.close()
    completed = sum(c['completed'] for c in counts.values())
    interrupted = sum(c['interrupted'] for c in counts.values())
    failed = sum(c['failed'] for c in counts.values())
    assert total_started() == completed+interrupted+failed <= 644 <= HARD_TOTAL
    assert counts['DE']['started'] <= DE_LIMIT
    phi_upper_float = float(phi_record['float_display']) if phi_record else None
    gain = best['score']-phi_upper_float if best is not None and phi_record else None
    nonzero_gain = best_nonzero['score']-phi_upper_float if best_nonzero is not None and phi_record else None
    result.update({
        'status': 'failed_saved_partial_evidence' if fatal else 'finite_numerical_run_closed',
        'failure': fatal, 'search_stop_reason': search_stop,
        'family': 'four_atom_2h_by_2V_independent_heavy_displacements_and_antisymmetric_heavy_masses',
        'variable_order': ['xi0', 'xi1', 'eta0', 'eta1', 'beta'], 'bounds': BOUNDS,
        'Lin_parameters': dict(zip(('h', 'V', 'w'), base)) if base else None,
        'complete_cell_mass': 2*(1+base[2]) if base else None,
        'patch': {'period_cells': [3, 3], 'atoms_before_aggregation': 36,
                  'rectangle_enumeration': 'all coordinate-endpoint closed rectangles; actual selected extrema',
                  'witness': 'common sqrt(selected_mass)/2 radius in both coordinates',
                  'union': 'actual periodic projection and geometric union, divided by complete cell mass'},
        'phi_reference': phi_record, 'numeric_gain_tolerance': GAIN_TOL,
        'controls': controls, 'best_completed': best, 'best_nonzero_completed': best_nonzero,
        'best_DE_completed': best_de, 'best_numerical_minus_Phi_upper_display': gain,
        'best_nonzero_numerical_minus_Phi_upper_display': nonzero_gain,
        'numerical_improvement_candidate': bool(gain is not None and gain > GAIN_TOL),
        'zero_control_possible_larger_witness_gain': bool(controls and phi_record and
                                                        controls[0]['score'] > phi_upper_float+GAIN_TOL),
        'optimizer': {'method': 'scipy.optimize.differential_evolution', 'seed': SEED,
                      'dimension': 5, 'popsize': 8, 'maxiter': 15, 'polish': False,
                      'workers': 1, 'updating': 'immediate', 'tol': 0., 'atol': 0., 'x0': ZERO,
                      'x0_in_population_and_counted': True, 'max_DE_started': DE_LIMIT,
                      'completed_generation_callbacks': generation_count,
                      'scipy_result_available': opt is not None,
                      'scipy_nfev': int(opt.nfev) if opt is not None else None,
                      'scipy_nit': int(opt.nit) if opt is not None else None,
                      'scipy_success': bool(opt.success) if opt is not None else None,
                      'scipy_message': str(opt.message) if opt is not None else None,
                      'interpretation': 'A generation cap or timeout is not convergence or a global optimum.'},
        'counts': counts, 'total_configurations_started': total_started(),
        'maximum_scheduled_configurations': 644, 'hard_total_started_limit': HARD_TOTAL,
        'independent_geometry_checks': checks,
        'independence_scope': 'Rebuilds atoms/closures/wrap; endpoint cells vs discovery segment tree. Shared mathematical definitions and float arithmetic; not a strict certificate.',
        'timing': {'prepare_seconds': prep_seconds, 'search_end_seconds': search_end,
                   'numerical_stage_seconds': numerical_seconds, 'prepare_controls_search_limit': SEARCH_SECONDS,
                   'total_numerical_limit': TOTAL_SECONDS, 'module_import_seconds_separate': _IMPORT_SECONDS,
                   'guards': 'Inside rectangle enumeration, wrap translation, and union event loops; interrupted partial areas are discarded.'},
        'search_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
        'full_trajectory': {'storage': 'local temporary JSONL; excluded from archive',
                            'basename': trace_path.name if trace_path else None,
                            'sha256': sha256(trace_path.read_bytes()).hexdigest() if trace_path else None,
                            'events': 'started and terminal records; every completed configuration retained'},
        'environment': {'python': sys.version, 'platform': platform.platform(),
                        'numpy': np.__version__, 'scipy': scipy.__version__},
        'scope': 'Numerical finite witnesses only. Any improvement requires a separately approved rigorous rational/interval and L1-transfer certificate.'})
    atomic_json(result_path, result)
    print(json.dumps({'status': result['status'], 'counts': counts, 'search_stop': search_stop,
                      'numerical_seconds': numerical_seconds, 'best_minus_Phi': gain,
                      'independent_checks': checks}, ensure_ascii=False), flush=True)


if __name__ == '__main__':
    main()
