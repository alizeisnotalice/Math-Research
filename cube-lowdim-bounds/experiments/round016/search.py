#!/usr/bin/env python3
"""Finite p=2..8 shape selection, derived integer slices and five LPs."""
from __future__ import annotations
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
import math
from pathlib import Path
import platform
import sys
from time import monotonic

import mpmath as mp
import numpy as np
import scipy
from scipy.optimize import linprog
import sympy as sp

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
POWERS = list(range(2, 9))
R_VALUES = [Q(15, 8), Q(151, 80), Q(19, 10), Q(153, 80)]
INCUMBENT = Q(411881, 106500)
LP_LIMIT = 5
POINT_LIMIT = 4000
COUNTS = {'distinct_interior_points': 0, 'high_precision_bump_values': 0,
          'audited_reference_bump_values': 0, 'tail_values': 0, 'LPs': 0,
          'direct_integral_points': 0, 'direct_quadrature_calls': 0}


def reference_source():
    for path in [ROOT.parent/'round012/search_quartic.py',
                 ROOT.parent.parent/'Math-Research-lowdim-agent/cube-lowdim-bounds/experiments/round012/search_quartic.py']:
        if path.exists():
            return path
    raise FileNotFoundError('round012/search_quartic.py reference missing')


def load_reference():
    source = reference_source()
    spec = importlib.util.spec_from_file_location('round016_audited_p234', source)
    c = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = c
    spec.loader.exec_module(c)
    return c, source


def hp(p, t):
    """Derived PV derivative integral; evaluated at current mpmath precision."""
    if not 0 < t < 1:
        raise ValueError('This implementation is the interior branch only.')
    n = p-1
    ap, am = 1+t, 1-t
    poly = mp.fsum((ap**k+am**k)/(n-k) for k in range(n))
    return 2*p*(poly+ap**n*(mp.log(t)-mp.log1p(t))
                +am**n*(mp.log(t)-mp.log1p(-t)))


def bump(p, u, v):
    u, v = mp.mpf(u), mp.mpf(v)
    a, b = 1-v, 1-u
    return a**(p-1)*hp(p, u/a)+b**(p-1)*hp(p, v/b)


def symbolic_checks():
    t = sp.symbols('t', positive=True)
    lt, lp, lm = sp.symbols('log_t log_1plus_t log_1minus_t')
    known = {2: 4*(2-(1-t)*lm-(1+t)*lp+2*lt),
             3: 6*(3+(1-t)**2*(lt-lm)+(1+t)**2*(lt-lp)),
             4: 8*(2*t*t+sp.Rational(11, 3)+(1+t)**3*(lt-lp)+(1-t)**3*(lt-lm))}
    checks = {}
    for p in [2, 3, 4]:
        n = p-1
        poly = sum(((1+t)**k+(1-t)**k)/sp.Integer(n-k) for k in range(n))
        generic = 2*p*(poly+(1+t)**n*(lt-lp)+(1-t)**n*(lt-lm))
        residual = sp.expand(generic-known[p])
        assert residual == 0
        checks[str(p)] = str(residual)
    return checks


def direct_checks():
    # Exactly ten predefined points, integrating the regularized PV derivative
    # rather than evaluating the polynomial-division closed form as integrand.
    points = [(p, s) for p in range(5, 9) for s in ['.2', '.8']]+[(5, '.0001'), (8, '.9999')]
    out = []
    with mp.workdps(80):
        for p, text in points:
            if COUNTS['direct_integral_points'] >= 10:
                raise RuntimeError('Direct-integral diagnostic budget reached')
            t = mp.mpf(text)
            n = p-1
            left = mp.quad(lambda y: (1+y)**n/(y-t), [-1, 0])
            def regular(y):
                if y == t:
                    return -n*(1-t)**(n-1)
                if y == 1:
                    return -(1-t)**n/(y-t)
                return (1-t)**n*mp.expm1(n*mp.log1p((t-y)/(1-t)))/(y-t)
            right = mp.quad(regular, [0, t, 1])+(1-t)**n*mp.log((1-t)/t)
            integral = 2*p*(left-right)
            closed = hp(p, t)
            error = abs(integral-closed)
            assert error < mp.mpf('1e-65')
            out.append({'p': p, 't': text, 'direct_integral': mp.nstr(integral, 60),
                        'closed_form': mp.nstr(closed, 60), 'absolute_error': mp.nstr(error, 8)})
            COUNTS['direct_integral_points'] += 1
            COUNTS['direct_quadrature_calls'] += 2
    return out


def endpoint_constants():
    out, values = {}, []
    for p in POWERS:
        a = Q(2, p-1)
        for k in range(p+1):
            if k != 1:
                a += 2*math.comb(p, k)*Q(2)**(p-k)*(-1)**k*(Q(2)**(k-1)-1)/Q(k-1)
        b = -p*2**p
        val = mp.mpf(a.numerator)/a.denominator+b*mp.log(2)
        assert val > 0
        values.append(float(val))
        out[str(p)] = {'rational_part': str(a), 'log2_coefficient': b,
                       'high_precision_value': mp.nstr(val, 55), 'LP_float': float(val)}
    return out, np.asarray(values)


def fixed_points():
    radii = np.sort(np.r_[np.linspace(.0025, .9975, 180),
                          [1e-6, 1e-5, 1e-4, .0005, .001, .8, .82, .84, .85, .86,
                           .88, .9, .998, .999, .9995, .9999, .99999, .999999, 1-1e-8, 1-1e-10]])
    angles = np.r_[np.arange(15)/15., .97, .99, .999, .99999, 1-1e-8]
    points = [(float(r), float(q)) for r in radii for q in angles]
    assert len(points) == len(set(points)) == POINT_LIMIT
    COUNTS['distinct_interior_points'] = len(points)
    return points


def tail_f(R, c):
    if abs(c/R) < .25:
        power, total = 1., 0.
        for n in range(40):
            total += power/((n+2)*R)
            power *= c/R
        return total
    return -(R/c**2)*math.log1p(-c/R)-1/c


def tail_rhs(Rq, points):
    R = float(Rq)
    COUNTS['tail_values'] += len(points)
    return np.asarray([2*(2*tail_f(R, r)+tail_f(R, r*q)+tail_f(R, -r*q))/(R-1) for r, q in points])


def solve(Rq, rows, rhs, points, H, control=False):
    if COUNTS['LPs'] >= LP_LIMIT:
        raise RuntimeError('Five-LP hard limit reached')
    COUNTS['LPs'] += 1
    R = float(Rq)
    columns = [0, 2] if control else list(range(7))
    selected_p = [POWERS[j] for j in columns]
    cost_terms = np.asarray([2/((p+1)*(p+2)) for p in selected_p])
    D1 = 2*(2*tail_f(R, 1)+1/R)/(R-1)
    all_rows = np.vstack([rows[:, columns], H[columns], np.full(len(columns), (R-1)/4)])
    all_rhs = np.r_[rhs, D1, R]
    opt = linprog(-cost_terms, A_ub=all_rows, b_ub=all_rhs, bounds=[(0, None)]*len(columns),
                  method='highs', options={'primal_feasibility_tolerance': 1e-9, 'dual_feasibility_tolerance': 1e-9})
    result = {'LP': COUNTS['LPs'], 'R_exact': str(Rq), 'R': R,
              'family': 'same_grid_old_p2_p4_control' if control else 'all_p2_through_p8',
              'powers_allowed': selected_p, 'success': bool(opt.success), 'solver_status': int(opt.status),
              'solver_message': opt.message, 'iterations': int(opt.nit),
              'fixed_points': len(points), 'constraints': len(all_rhs), 'D_R1': D1}
    if not opt.success:
        return result
    coeff = np.zeros(7)
    coeff[columns] = opt.x
    slack = all_rhs-all_rows@opt.x
    dual = -np.asarray(opt.ineqlin.marginals)
    chosen = set(map(int, np.argsort(slack)[:8])) | set(map(int, np.flatnonzero(dual > 1e-10)))
    active = []
    for i in sorted(chosen, key=lambda i: slack[i]):
        desc = {'kind': 'fixed_interior_generator', 'r': points[i][0], 'd_over_r': points[i][1]} if i < len(points) else {
            'kind': 'proved_exterior_sufficient' if i == len(points) else 'kernel_dominance_sufficient'}
        active.append({**desc, 'slack': float(slack[i]), 'benefit_dual': float(dual[i]),
                       'coefficients_in_allowed_power_order': all_rows[i].tolist(), 'rhs': float(all_rhs[i])})
    cost = R*R/(R-1)-float(cost_terms@opt.x)
    result.update({'e_by_power': {str(p): float(coeff[j]) for j, p in enumerate(POWERS)},
                   'selected_positive_powers': [p for j, p in enumerate(POWERS) if coeff[j] > 1e-10],
                   'cost_candidate': cost, 'cost_minus_current_strict_upper': cost-float(INCUMBENT),
                   'minimum_grid_slack': float(np.min(slack[:len(points)])),
                   'exterior_slack': float(slack[-2]), 'kernel_dominance_slack': float(slack[-1]),
                   'near_active_constraints': active, 'benefit_primal': float(-opt.fun),
                   'benefit_dual': float(all_rhs@dual), 'dual_gap': float(all_rhs@dual+opt.fun),
                   'lower_bound_reduced_costs': np.asarray(opt.lower.marginals).tolist()})
    return result


def main():
    if (ROOT/'result.json').exists():
        raise SystemExit('Result exists; refusing to repeat the finite search budget.')
    started = monotonic()
    mp.mp.dps = 60
    symbolic = symbolic_checks()
    direct = direct_checks()
    endpoint, H = endpoint_constants()
    points = fixed_points()
    reference, source = load_reference()
    rows = []
    discrepancies = {str(p): 0. for p in [2, 3, 4]}
    min_arg, max_arg = 1., 0.
    for i, (r, q) in enumerate(points):
        u, v = r*(1+q)/2, r*(1-q)/2
        args = [mp.mpf(u)/(1-mp.mpf(v)), mp.mpf(v)/(1-mp.mpf(u))]
        assert all(0 < t < 1 for t in args)
        min_arg, max_arg = min(min_arg, float(min(args))), max(max_arg, float(max(args)))
        row = [float(bump(p, u, v)) for p in POWERS]
        COUNTS['high_precision_bump_values'] += 7
        for j, p in enumerate([2, 3, 4]):
            old = reference.bump_generator(p, u, v)
            discrepancies[str(p)] = max(discrepancies[str(p)], abs(row[j]-old))
            COUNTS['audited_reference_bump_values'] += 1
        rows.append(row)
        if (i+1) % 1000 == 0:
            print(f'precomputed fixed points {i+1}/{len(points)}', flush=True)
    rows = np.asarray(rows)
    assert np.isfinite(rows).all()
    assert max(discrepancies.values()) < 1e-9
    cache, results = {}, []
    for Rq in R_VALUES:
        rhs = tail_rhs(Rq, points)
        cache[str(Rq)] = rhs
        entry = solve(Rq, rows, rhs, points, H)
        results.append(entry)
        print(json.dumps({k: entry.get(k) for k in ['R_exact', 'cost_candidate', 'selected_positive_powers', 'e_by_power']}), flush=True)
    feasible = [x for x in results if x['success']]
    if not feasible:
        raise RuntimeError('No feasible LP; no unplanned fifth solve will be run.')
    best = min(feasible, key=lambda x: x['cost_candidate'])
    Rq = Q(best['R_exact'])
    control = solve(Rq, rows, cache[str(Rq)], points, H, control=True)
    results.append(control)
    data = {'round': 16, 'status': 'finite_numerical_shape_selection_only_no_certificate',
            'current_strict_upper': str(INCUMBENT), 'powers': POWERS,
            'R_values_exact': list(map(str, R_VALUES)), 'LP_budget': LP_LIMIT,
            'distinct_point_budget': POINT_LIMIT, 'actual_counts': dict(COUNTS),
            'grid': {'definition': '200 prescribed radii times 20 prescribed d/r values; see PLAN.md and fixed_points()',
                     'point_digest_sha256': sha256(json.dumps(points).encode()).hexdigest(),
                     'minimum_slice_argument': min_arg, 'maximum_slice_argument': max_arg},
            'formula': {'branch': '0<t<1 only', 'n': 'p-1',
                        'J_n': 'sum(k=0..n-1) a^k/(n-k) + a^n log(abs((1-a)/a))',
                        'h_p': '2p*(J_n(1+t)+J_n(1-t))',
                        'precision': '60 decimal digits mpmath, then float for LP; diagnostic only',
                        'sympy_p234_exact_residuals': symbolic,
                        'same_grid_p234_max_absolute_difference_from_audited_float': discrepancies,
                        'direct_integral_precision': 80, 'direct_integral_checks': direct},
            'endpoint_Hp': endpoint, 'LPs': results, 'best_all_basis': best,
            'same_grid_old_basis_control': control,
            'control_cost_minus_all_basis_cost': control.get('cost_candidate', math.nan)-best['cost_candidate'],
            'geometry': 'true two-slice generator in every fixed interior point; nonnegative e_p',
            'exterior': 'sum e_p H_p <= D_R(1), using the proved round013 sufficient condition',
            'half_mass_discount': '2e_p/((p+1)(p+2))',
            'reference_dependency': 'round012/search_quartic.py', 'reference_sha256': sha256(source.read_bytes()).hexdigest(),
            'search_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
            'environment': {'python': sys.version, 'platform': platform.platform(), 'numpy': np.__version__,
                            'scipy': scipy.__version__, 'sympy': sp.__version__, 'mpmath': mp.__version__},
            'elapsed_seconds': monotonic()-started,
            'scope': 'No continuous angular certificate, new rigorous upper bound, extra R values, grid enrichment, or extra LPs.'}
    (ROOT/'result.json').write_text(json.dumps(data, indent=2)+'\n')
    print(json.dumps({'actual_counts': COUNTS, 'best_all_basis': best,
                      'same_grid_old_basis_control': control,
                      'control_cost_minus_all_basis_cost': data['control_cost_minus_all_basis_cost'],
                      'elapsed_seconds': data['elapsed_seconds']}, indent=2), flush=True)


if __name__ == '__main__':
    main()
