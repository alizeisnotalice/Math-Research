#!/usr/bin/env python3
"""Finite six-variable periodic-atom search with complete mass accounting."""
from __future__ import annotations
from hashlib import sha256
import importlib.util
import json
import math
from pathlib import Path
import platform
import sys
import tempfile
from time import monotonic

import numpy as np
import scipy
from scipy.optimize import differential_evolution

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
SEED = 2026102415
BUDGET = 1000
PHI = 1.685509993355525
BOUNDS = [(.85, 1.15), (.85, 1.15), (.7, 1.3), (0., 1.), (0., 1.), (math.log(1e-4), 0.)]


def source_path():
    for path in [ROOT.parent/'round009/search_satellite.py',
                 ROOT.parent.parent/'Math-Research-lowdim-agent/cube-lowdim-bounds/experiments/round009/search_satellite.py']:
        if path.exists():
            return path
    raise FileNotFoundError('round009/search_satellite.py geometry core not found')


def load_geometry():
    path = source_path()
    spec = importlib.util.spec_from_file_location('round015_round009_geometry', path)
    core = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = core
    spec.loader.exec_module(core)
    return core, path


def parameters(core, ratios, x, y, mass):
    h0, v0, w0, _, _ = core.lin_parameters()
    h, v, w = h0*float(ratios[0]), v0*float(ratios[1]), w0*float(ratios[2])
    return {'h_ratio': float(ratios[0]), 'V_ratio': float(ratios[1]), 'w_ratio': float(ratios[2]),
            'h': h, 'V': v, 'w': w, 'x_fraction': float(x), 'y_fraction': float(y),
            'satellite_mass': float(mass), 'log_mass': math.log(mass) if mass > 0 else None,
            'periods': [2*h, v], 'period_mass': 1+w+mass}


def geometry(core, p):
    """Inject parameters only; reuse the unmodified full rectangle geometry."""
    h, v, w = p['h'], p['V'], p['w']
    px, py = 2*h, v
    def patch_points(x_fraction, y_fraction, satellite_mass):
        seeds = [(0., 0., 1.), (h, 0., w)]
        if satellite_mass > 0:
            seeds.append((px*x_fraction, py*y_fraction, satellite_mass))
        pts, weights = [], []
        for kx in (-1, 0, 1):
            for ky in (-1, 0, 1):
                for x, y, m in seeds:
                    pts.append((x+kx*px, y+ky*py))
                    weights.append(m)
        return np.asarray(pts), np.asarray(weights), px, py, 1+w+satellite_mass
    core.patch_points = patch_points
    boxes, mass = core.witness_boxes(p['x_fraction'], p['y_fraction'], p['satellite_mass'])
    assert mass == p['period_mass']
    area = core.union_area(boxes)
    return area/mass, boxes, mass, area


def main():
    if (ROOT/'result.json').exists():
        raise SystemExit('Result exists; refusing a duplicate finite-budget search.')
    start = monotonic()
    core, source = load_geometry()
    calls = 0
    best = None
    best_de = None
    trace = Path(tempfile.gettempdir())/'cube-lowdim-round015-full-trajectory.jsonl'
    trace_file = trace.open('x')
    def evaluate(p, label):
        nonlocal calls, best, best_de
        if calls >= BUDGET:
            raise RuntimeError('Hard objective-call budget reached.')
        calls += 1
        score, boxes, mass, area = geometry(core, p)
        row = {'evaluation': calls, 'source': label, 'parameters': p,
               'score': score, 'union_area': area, 'box_count': len(boxes), 'period_mass': mass}
        trace_file.write(json.dumps(row, separators=(',', ':'))+'\n')
        if p['satellite_mass'] > 0 and (best is None or score > best['score']):
            best = row
        if label == 'six_dimensional_DE' and (best_de is None or score > best_de['score']):
            best_de = row
        if calls % 100 == 0:
            trace_file.flush()
            print(json.dumps({'calls': calls, 'best_nonzero': best['score'],
                              'best_minus_phi': best['score']-PHI}), flush=True)
        return row
    zero = evaluate(parameters(core, [1., 1., 1.], 0., 0., 0.), 'zero_mass_Lin_baseline')
    assert abs(zero['score']-PHI) < 1e-10
    control = evaluate(parameters(core, [1., 1., 1.], .5, 0., 1e-4), 'round009_best_nonzero_control')
    assert abs(control['score']-1.6854857935889827) < 1e-10
    def objective(z):
        p = parameters(core, z[:3], z[3], z[4], math.exp(float(z[5])))
        return -evaluate(p, 'six_dimensional_DE')['score']
    opt = differential_evolution(objective, BOUNDS, popsize=10, maxiter=14, polish=False,
                                 seed=SEED, workers=1, updating='immediate', tol=0., atol=0.,
                                 x0=[1., 1., 1., .5, 0., math.log(1e-4)])
    trace_file.close()
    assert calls == int(opt.nfev)+2 <= BUDGET
    checks = []
    for label, row in [('zero_mass_Lin_baseline', zero), ('best_nonzero', best)]:
        value, boxes, mass, area = geometry(core, row['parameters'])
        edged = core.independent_edge_cells(boxes)
        check = {'label': label, 'box_count': len(boxes), 'period_mass': mass,
                 'sweep_area': area, 'endpoint_cell_area': edged,
                 'absolute_area_difference': abs(area-edged),
                 'endpoint_cell_score': edged/mass, 'recorded_score': row['score'],
                 'absolute_score_difference': abs(edged/mass-row['score'])}
        assert check['absolute_score_difference'] < 1e-10
        checks.append(check)
    h0, v0, w0, _, _ = core.lin_parameters()
    result = {'round': 15, 'status': 'finite_numerical_experiment_complete_no_certificate',
              'family': 'joint_h_V_w_and_free_positive_satellite_in_2h_by_V_period',
              'atoms_per_period': ['(0,0;1)', '(h,0;w)', '(2hx,Vy;m)'],
              'normalization': 'complete period mass 1+w+m', 'patch_periods': [3, 3],
              'enumeration': 'all coordinate-rectangular closed atom subsets in the true finite patch',
              'witness': 'common radius sqrt(M)/2 in both coordinates; true periodic split and geometric union',
              'Lin_parameters': {'h0': h0, 'V0': v0, 'w0': w0},
              'variable_order': ['h/h0', 'V/V0', 'w/w0', 'x', 'y', 'log(m)'], 'bounds': BOUNDS,
              'optimizer': {'method': 'scipy.optimize.differential_evolution', 'seed': SEED,
                            'popsize': 10, 'maxiter': 14, 'polish': False, 'tol': 0., 'atol': 0.,
                            'workers': 1, 'updating': 'immediate', 'actual_nfev': int(opt.nfev),
                            'actual_nit': int(opt.nit), 'success': bool(opt.success), 'message': str(opt.message)},
              'actual_counts': {'objective_calls': calls, 'objective_call_hard_limit': BUDGET,
                                'control_calls': 2, 'DE_calls': int(opt.nfev),
                                'independent_geometry_rebuilds': 2, 'total_geometry_builds': calls+2},
              'controls': {'zero_mass': zero, 'previous_best_nonzero': control, 'Phi': PHI,
                           'baseline_minus_Phi': zero['score']-PHI},
              'best_nonzero': best, 'best_DE': best_de, 'best_nonzero_minus_Phi': best['score']-PHI,
              'numerical_improvement_candidate': best['score'] > PHI+1e-10,
              'independent_endpoint_cell_checks': checks,
              'geometry_dependency': 'round009/search_satellite.py', 'geometry_sha256': sha256(source.read_bytes()).hexdigest(),
              'search_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
              'full_trajectory': 'local temporary JSONL only; excluded from archived result',
              'elapsed_seconds': monotonic()-start,
              'environment': {'python': sys.version, 'platform': platform.platform(),
                              'numpy': np.__version__, 'scipy': scipy.__version__},
              'scope': 'Numerical finite witnesses only. No rigorous lower-bound certificate, boundary expansion, or extra optimization.'}
    (ROOT/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
