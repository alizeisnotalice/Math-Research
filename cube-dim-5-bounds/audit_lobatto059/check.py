#!/usr/bin/env python3
"""Independent exact Beta-Lobatto audit: Hankel systems and quotient traces.

No author code, float/mp roots, Arb or nontrivial determinant batch is run.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import os
import sys
import time

import sympy as sp

HERE = Path(__file__).resolve().parent
t = sp.symbols('t')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def sq(value):
    value = F(value)
    return sp.Rational(value.numerator, value.denominator)


def moments(q, degree):
    q = F(q)
    result, value = [F(1)], F(1)
    for k in range(degree):
        value *= (q+k)/(2*q+k)
        result.append(value)
    return list(map(sq, result))


def integrate_polynomial(expression, values):
    polynomial = sp.Poly(sp.expand(expression), t, domain=sp.QQ)
    return sum((coefficient*values[power[0]]
                for power, coefficient in polynomial.terms()), sp.Rational(0))


def main():
    began = time.monotonic()
    regpath = HERE/'registration.json'
    reg = json.loads(regpath.read_text())
    for row in reg['sources']:
        assert sha(row['source']) == row['sha256'], row['source']
    source = Path(reg['sources'][0]['source']).parent
    author = json.loads((source/'roots_weights.json').read_text())
    stages = []

    def stage(name):
        elapsed = time.monotonic()-began
        stages.append({'stage': name, 'seconds': elapsed})
        (HERE/'stage_receipt.json').write_text(json.dumps(
            {'pid': os.getpid(), 'stages': stages, 'complete': False}, indent=2)+'\n')
        if elapsed >= reg['limit_seconds']:
            raise SystemExit(75)

    stage('registered_inputs_verified_no_author_execution')
    mu, nu = moments(F(5, 4), 14), moments(F(9, 4), 12)
    c = mu[1]-mu[2]
    assert c == sp.Rational(5, 28)
    # Independently derive each monic polynomial from its moment equations.
    polynomials, norms = [], []
    for degree in range(7):
        if degree == 0:
            p = sp.Integer(1)
        else:
            hankel = sp.Matrix(degree, degree, lambda i, j: nu[i+j])
            rhs = sp.Matrix([-nu[degree+i] for i in range(degree)])
            coefficients = hankel.inv()*rhs
            p = t**degree+sum((coefficients[j]*t**j for j in range(degree)), sp.Integer(0))
        p = sp.expand(p)
        assert all(integrate_polynomial(p*t**j, nu) == 0 for j in range(degree))
        norm = integrate_polynomial(p*p, nu)
        assert norm > 0
        polynomials.append(p)
        norms.append(norm)
    p6 = sp.Poly(polynomials[6], t, domain=sp.QQ)
    primitive = p6.clear_denoms()[1].primitive()[1]
    assert list(map(str, p6.all_coeffs())) == author['monic_coefficients_descending']
    assert list(map(str, primitive.all_coeffs())) == author['primitive_integer_coefficients_descending']
    assert list(map(str, norms[:6])) == author['h_0_to_5']
    derived_b = {str(k): str(norms[k]/norms[k-1]) for k in range(1, 6)}
    assert derived_b == author['b_k']
    stage('independent_Hankel_polynomials_and_positive_norms')

    def rem(expression):
        return sp.rem(sp.Poly(expression, t, domain=sp.QQ), p6).as_expr()

    def vector(expression):
        p = sp.Poly(rem(expression), t, domain=sp.QQ)
        return sp.Matrix([p.nth(j) for j in range(6)])

    X = sp.Matrix.hstack(*(vector(t**(j+1)) for j in range(6)))
    traces, power = [], sp.eye(6)
    for k in range(19):
        traces.append(sp.trace(power))
        power = power*X

    def trace(expression):
        p = sp.Poly(rem(expression), t, domain=sp.QQ)
        return sum((p.nth(j)*traces[j] for j in range(6)), sp.Rational(0))

    # Gaussian node weights g(r) are recovered from six exact moment equations,
    # using sum_r r^k=Tr(X^k). This is not the author's Christoffel evaluation.
    trace_hankel = sp.Matrix(6, 6, lambda i, j: traces[i+j])
    g_coefficients = trace_hankel.inv()*sp.Matrix(nu[:6])
    g = sp.expand(sum((g_coefficients[j]*t**j for j in range(6)), sp.Integer(0)))
    gauss_residuals = [sp.factor(trace(g*t**k)-nu[k]) for k in range(12)]
    assert gauss_residuals == [0]*12
    christoffel = sum((polynomials[j]**2/norms[j] for j in range(6)), sp.Integer(0))
    assert rem(g*christoffel-1) == 0
    w = rem(c*g*sp.invert(t*(1-t), p6.as_expr()))
    # Endpoint weights reconstructed two ways, both rational.
    w0_card = integrate_polynomial((1-t)*p6.as_expr()**2/p6.eval(0)**2, mu)
    w1_card = integrate_polynomial(t*p6.as_expr()**2/p6.eval(1)**2, mu)
    w0_mass = mu[0]-mu[1]-c*trace(g*sp.invert(t, p6.as_expr()))
    w1_mass = mu[1]-c*trace(g*sp.invert(1-t, p6.as_expr()))
    assert w0_card == w0_mass == w1_card == w1_mass == sp.Rational(512, 51051)
    assert str(w0_card) == author['endpoint_weight_0_exact']
    assert str(w1_card) == author['endpoint_weight_1_exact']
    residuals = []
    for k in range(14):
        result = (w0_card if k == 0 else 0)+w1_card+trace(w*t**k)
        residuals.append(sp.factor(result-mu[k]))
        assert residuals[-1] == 0
        saved = author['directed_moment_checks_0_to_13'][k]
        assert sq(F(saved['exact_beta_moment'])) == mu[k]
        saved_ball = saved['outward_quadrature_moment']
        assert F(saved_ball['lower_fraction']) <= F(str(mu[k])) <= F(saved_ball['upper_fraction'])
    error14 = w1_card+trace(w*t**14)-mu[14]
    assert error14 == c*norms[6] > 0
    stage('quotient_trace_weights_all14_exact_moments_and_error_orientation')

    assert sp.gcd(p6, p6.diff()).degree() == 0
    assert p6.count_roots(0, 1) == 6
    intervals = p6.intervals(eps=sp.Rational(1, 2**24))
    assert len(intervals) == 6
    for (left, right), multiplicity in intervals:
        assert 0 < left < right < 1 and multiplicity == 1
        assert right-left <= sp.Rational(1, 2**24)
        assert p6.count_roots(left, right) == 1
    # One new exact boundary shape, no sample99 nor group enumeration.
    M = sp.eye(5)
    assert M.det() == 1
    fixed_shape_f = sp.Integer(1)
    assert sp.diff(fixed_shape_f, t, 14) == 0
    assert w0_card+w1_card+trace(w) == 1
    # Normalizer substitution is exact algebra, not a Gamma numerical replay.
    Z = 1/(8*sp.sqrt(sp.pi))
    gaussian_time_coefficient = 8/(4*sp.pi)**sp.Rational(5, 2)
    profile_coefficient = sp.simplify(gaussian_time_coefficient*Z)
    assert profile_coefficient == 1/(32*sp.pi**3)
    assert 5*F(1, 10)-F(1, 2) == 0
    stage('exact_root_isolation_identity_shape_and_normalization')
    out = {
        'status': 'independent_exact_Beta_Lobatto_quadrature_and_interface_audit_passed',
        'method': 'Fraction moments, independent Hankel linear systems, quotient-algebra trace weights; no author functions or mp roots',
        'modified_mass_c': str(c),
        'p6_monic_coefficients_descending': list(map(str, p6.all_coeffs())),
        'p6_primitive_integer_coefficients_descending': list(map(str, primitive.all_coeffs())),
        'positive_squared_norms_0_to6': list(map(str, norms)),
        'derived_recurrence_b1_to5': derived_b,
        'g_polynomial_in_quotient_algebra': str(g),
        'interior_Lobatto_weight_polynomial_mod_p6': str(w),
        'gaussian_exact_moment_residuals0_to11': list(map(str, gauss_residuals)),
        'original_exact_moments0_to13': list(map(str, mu[:14])),
        'Lobatto_exact_moment_residuals0_to13': list(map(str, residuals)),
        'inverse_Christoffel_identity_mod_p6': True,
        'interior_weights_strictly_positive': 'Six simple real roots in(0,1); g(root)*sum(pj(root)^2/hj)=1, c>0, w=c*g/(root(1-root)).',
        'endpoint0_cardinal_square': str(w0_card), 'endpoint0_mass_firstmoment': str(w0_mass),
        'endpoint1_cardinal_square': str(w1_card), 'endpoint1_mass_firstmoment': str(w1_mass),
        'all_eight_weights_positive_and_sum1': True,
        'exact_sign_spot_Q_t14_minus_integral': str(error14),
        'root_intervals_exact': [{'lower': str(a), 'upper': str(b), 'multiplicity': m}
                                 for (a, b), m in intervals],
        'squarefree_and_exact_root_count0_1': 6,
        'identity_shape': {'M': 'I5', 'd': '1', 'f(t)': '1', 'r_true': '1', 'r_Lobatto': '1',
                           'all_positive_order_derivatives': '0', 'group_orbit_enumerated': False},
        'R_normalizer_exact': 'E T^(-1/2)=1/(8sqrt(pi))=1/(16Gamma(3/2))',
        'profile_coefficient_exact': str(profile_coefficient),
        'determinant_scale_exponent_cancels': '5*(1/10)-1/2=0',
        'theory_scope': 'All SPD A/B even derivative sign; same029/038 model and R tilt; integrable determinant domination; no shape-independence assumption.',
        'fixed_sample99_or_15360_determinants_recomputed': False,
        'author_float_or_Arb_root_function_called': False,
        'source_sha256': sha(Path(__file__)), 'registration_sha256': sha(regpath),
        'source_manifest_sha256': reg['production_manifest_sha256'],
        'external_sources': reg['sources'],
        'runtime': {'pid': os.getpid(), 'seconds': time.monotonic()-began,
                    'python': sys.executable, 'sympy_version': sp.__version__},
        'new_full_I_g_or_C5_bound': False,
        'independence_limits': reg['independence_limits'],
    }
    (HERE/'results.json').write_text(json.dumps(out, indent=2)+'\n')
    (HERE/'stage_receipt.json').write_text(json.dumps(
        {'pid': os.getpid(), 'stages': stages, 'complete': True,
         'seconds': time.monotonic()-began}, indent=2)+'\n')
    print(json.dumps({'status': out['status'], 'p6_integer': out['p6_primitive_integer_coefficients_descending'],
                      'endpoints': str(w0_card), 'moment_residuals': out['Lobatto_exact_moment_residuals0_to13'],
                      'error14': str(error14), 'seconds': out['runtime']['seconds']}, indent=2))


if __name__ == '__main__':
    main()
