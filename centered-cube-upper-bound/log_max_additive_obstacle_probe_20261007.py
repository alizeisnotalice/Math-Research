#!/usr/bin/env python3
"""Certified rational bad-component constants; saved-cell additive-density guard.
Unknown near-optimal source remains symbolic. No old oracle or MC calls.
"""
from decimal import Decimal, localcontext
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import gzip
import hashlib
import json
import math
import time

BASE = Path(__file__).resolve().parent
PREFIX = 'log_max_additive_obstacle_probe_20261007'
BITS = 256
UNIT = 1 << BITS


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def lo(x):
    return F((x.numerator * UNIT) // x.denominator, UNIT)


def hi(x):
    return F(-((-x.numerator * UNIT) // x.denominator), UNIT)


def interval(a, b=None):
    b = a if b is None else b
    assert a <= b
    return lo(a), hi(b)


def plus(a, b):
    return interval(a[0] + b[0], a[1] + b[1])


def minus(a, b):
    return interval(a[0] - b[1], a[1] - b[0])


def scale(a, c):
    assert c >= 0
    return interval(c * a[0], c * a[1])


def times(a, b):
    vals = [x * y for x in a for y in b]
    return interval(min(vals), max(vals))


def pack(a):
    return dict(lower=str(a[0]), upper=str(a[1]), width=str(a[1] - a[0]),
                lower_float=float(a[0]), upper_float=float(a[1]),
                outward_rational_certified=True)


def ln2_bounds(N):
    x = F(1, 3)
    partial = 2 * sum((x**(2 * k + 1) / (2 * k + 1) for k in range(N)), F(0))
    remainder = 2 * x**(2 * N + 1) / ((2 * N + 1) * (1 - x*x))
    return interval(partial, partial + remainder)


def log1p_bounds(s, N):
    assert 0 <= s <= F(1, 10) and N % 2 == 0
    power, partial = F(1), F(0)
    for k in range(1, N + 1):
        power *= s
        partial += (1 if k % 2 else -1) * power / k
    return interval(partial, partial + power * s / (N + 1))


def e_bounds(N):
    partial = sum((F(1, math.factorial(k)) for k in range(N + 1)), F(0))
    tail = F(1, math.factorial(N + 1)) / (1 - F(1, N + 2))
    return interval(partial, partial + tail)


def bad_round(plan, reg, ln2, einv):
    n, eta, tau, H = plan['n'], F(plan['eta']), F(plan['tau']), plan['finite_H']
    B = F(2**n)
    Z = plus(interval(F(1)), scale(ln2, n))
    Kupper = times(Z, einv)  # Enclosure of Z/e, not an enclosure/value of unknown K.
    logs = log1p_bounds(eta / B, reg['interval_arithmetic']['log1p_terms'])
    logeta = log1p_bounds(eta, reg['interval_arithmetic']['log1p_terms'])
    S0 = minus(plus(plus(interval(1 / (1 + eta)), scale(ln2, n)), logs), logeta)
    gap = minus(S0, Kupper)
    nearZ = minus(S0, minus(Z, interval(2 * eta)))
    assert gap[0] > F(3, 10) and nearZ[0] >= 0
    # Closed form of the actual continuous bad Phi; background stays in Wbad.
    phi = minus(plus(scale(logs, B + eta), scale(ln2, eta * n)), scale(logeta, eta))
    phi = scale(phi, tau)
    Wbad = tau * (F(4**n) + eta)
    assert phi[0] > 0
    badratio = scale(phi, 1 / Wbad)
    shells = []
    finiteS = interval(1 / (1 + eta))
    finiteZ = interval(F(1))
    prev = F(1)
    for l in range(1, H + 1):
        R = 1 + F(l, H)
        delta = R**n - prev**n
        term, envelope = delta / (R**n + eta), delta / R**n
        finiteS = plus(finiteS, interval(term))
        finiteZ = plus(finiteZ, interval(envelope))
        shells.append(dict(l=l, R=str(R), previous_R=str(prev), shell_volume=str(delta),
                           exact_source_score_summand=str(term), exact_envelope_summand=str(envelope)))
        prev = R
    finitegap = minus(finiteS, Kupper)
    assert finitegap[0] > 0 and finiteZ[1] < Z[0]
    ell = F(1, 32 * n)
    factor = (1 - 2 * ell)**n
    assert factor >= F(15, 16) and 1 / (4 * ell) - n > 0
    l1lower = scale(S0, factor)
    l1gap = minus(l1lower, Kupper)
    assert l1gap[0] > F(1, 4)
    dilution = []
    for epsilon in map(F, reg['dilution_parameters']):
        ratio_threshold = 2 * Z[1] / epsilon  # N WA/Wbad sufficient threshold, WA unknown.
        bad_fraction = 1 / (1 + ratio_threshold)
        deficit_upper = epsilon / 2 + Kupper[1] * bad_fraction
        assert deficit_upper < epsilon
        dilution.append(dict(epsilon=str(epsilon), nearopt_WA=None, nearopt_Kg=None,
                             unknown_K=None, actual_integer_copy_count_N=None,
                             sufficient_N_times_WA_over_Wbad=str(ratio_threshold),
                             sufficient_N_times_WA_over_Wbad_float=float(ratio_threshold),
                             sufficient_N_ge='ceil(2*Z_upper*Wbad/(epsilon*WA))',
                             bad_full_component_mass_fraction_upper=str(bad_fraction),
                             bad_full_component_mass_fraction_upper_float=float(bad_fraction),
                             conditional_entropy_deficit_upper=str(deficit_upper),
                             conditional_deficit_upper_float=float(deficit_upper),
                             conditional_deficit_less_than_epsilon_certified=True,
                             atomic_hot_mass_fraction_upper=str(eta * tau * bad_fraction / Wbad),
                             unchanged_atomic_hot_profile=pack(S0),
                             dN_expression='N*WA/tau + 4^n + eta',
                             dN_lower_bound_from_sufficient_copy_ratio=str((1 + ratio_threshold) * Wbad / tau),
                             only_symbolic_conditional_algebra=True, real_nearmax_input_evaluated=False))
    return dict(n=n, tau=str(tau), eta=str(eta), full_Wbad=str(Wbad), Z=pack(Z),
                continuous_Z_over_e=pack(Kupper), unknown_K_value=None,
                continuous_atomic_S0=pack(S0), S0_minus_Z_over_e=pack(gap),
                S0_minus_Z_minus_2eta=pack(nearZ),
                actual_continuous_atomic_Phi_bad=pack(phi), bad_only_Phi_over_W=pack(badratio),
                finite_family=dict(H=H, scale_count=H + 1, scales_formula='1+l/H',
                    actual_finite_S0=pack(finiteS), finite_envelope_ZJ=pack(finiteZ),
                    actual_finite_S0_minus_continuous_Z_over_e=pack(finitegap),
                    full_true_finite_winner_formula='smallest common R_l>=max(1,2||x||inf)',
                    shells=shells, not_continuous_approximation_certificate=True),
                explicit_L1_analytic_factor_guard=dict(packet_halfwidth=str(ell), factor=str(factor),
                    factor_ge_15_over16=True, partial_logderivative_lower_minus_n_over_a=str(1 / (4 * ell) - n),
                    uniform_packet_profile_lower=pack(l1lower),
                    uniform_packet_profile_lower_minus_Z_over_e=pack(l1gap),
                    strict_gap_greater_than_1_over4=True, actual_positivewidth_max_simulated=False,
                    scope='Rational factor certificate using cited complete-packet truewinner proof; lower bound, not exact S.'),
                dilution=dilution, all_rational_interval_checks_passed=True)


def dec(x):
    return Decimal(x.numerator) / Decimal(x.denominator)


@lru_cache(maxsize=None)
def gainlog(s):
    return (1 + dec(s)).ln()


def nu_round(plan, reg):
    p = Path(plan['input_path'])
    assert sha(p) == plan['input_sha256']
    inp = json.loads(gzip.decompress(p.read_bytes()))
    n, tau = plan['n'], F(inp['input']['tau'])
    cells = inp['cells']
    b = max(map(F, inp['input']['J']))
    bounds = [(min(F(c['cell_bounds'][i][0]) for c in cells),
               max(F(c['cell_bounds'][i][1]) for c in cells)) for i in range(n)]
    G = [(low - b / 2, high + b / 2) for low, high in bounds]
    V = tau * math.prod(high - low for low, high in G)
    W = sum(map(F, inp['input']['weights']))
    D, Q = F(0), F(0)
    gains = [Decimal(0) for _ in reg['nu_t_values']]
    rows = []
    tvalues = list(map(F, reg['nu_t_values']))
    for index, c in enumerate(cells):
        if not c['E']:
            continue
        R, M, volume = F(c['winner_R']), F(c['M']), F(c['volume'])
        assert M > tau and R <= b
        for (low, high), (glo, ghi) in zip(c['cell_bounds'], G):
            assert F(low) - R / 2 >= glo and F(high) + R / 2 <= ghi
        ratio = tau / M
        assert 0 < ratio < 1
        D += tau * volume * ratio
        Q += tau * volume * ratio**2
        local = []
        for j, t in enumerate(tvalues):
            val = dec(tau * volume) * gainlog(t * ratio)
            gains[j] += val
            local.append(dict(t=str(t), frozen_E_loggain=str(val)))
        rows.append(dict(original_cell_index=index, original_cell_bounds=c['cell_bounds'],
                         original_volume=str(volume), original_winner_R=str(R), original_M=str(M),
                         A_R_nu=str(tau), p_nu=str(ratio), selected_query_inside_G_exact=True,
                         gains=local))
    assert 0 < Q <= D
    checks = []
    for t, gain in zip(tvalues, gains):
        rhs = t * D - t*t*Q/2
        residual = gain - dec(rhs)
        assert residual >= 0
        checks.append(dict(t=str(t), new_full_source_mass=str(W + t * V),
                           frozen_original_E_loggain=str(gain), rational_RHS=str(rhs),
                           decimal_RHS=str(dec(rhs)), numerical_residual=str(residual),
                           numerical_lower_bound_observed=True, logarithm_interval_certified=False,
                           no_new_full_maximum_evaluated=True))
    profile = BASE / (PREFIX + f'_n{n}_nu_frozen_profiles.json.gz')
    with gzip.open(profile, 'xt', encoding='utf8') as f:
        json.dump(dict(n=n, G_bounds=[[str(a), str(b)] for a, b in G], nu_density=str(tau),
                       full_nu_mass=str(V), rows=rows), f, separators=(',', ':'))
    return dict(n=n, original_input_path=str(p), original_input_sha256=plan['input_sha256'],
                original_saved_cell_count=len(cells), original_E_cell_count=len(rows),
                G_bounds=[[str(a), str(b)] for a, b in G], full_nu_density=str(tau), full_nu_mass=str(V),
                original_W=str(W), D_nu=str(D), Q_nu=str(Q), D_nu_minus_Q_nu=str(D - Q),
                Q_nu_le_D_nu_exact=True, full_selected_query_containment_exact=True,
                frozen_gain_checks=checks, profiles_path=str(profile), profiles_sha256=sha(profile),
                hot_set_A_evaluated=False, new_full_maximum_evaluated=False, old_oracle_called=False)


def main():
    start = time.perf_counter()
    regpath = BASE / (PREFIX + '_registration.json')
    amendment = BASE / (PREFIX + '_registration_amendment.json')
    reg, amend = json.loads(regpath.read_text()), json.loads(amendment.read_text())
    assert sha(regpath) == amend['original_registration_sha256']
    assert sha(Path(amend['analytic_proof_path'])) == amend['analytic_proof_sha256']
    ln2 = ln2_bounds(reg['interval_arithmetic']['ln2_terms'])
    e = e_bounds(reg['interval_arithmetic']['e_max_index'])
    assert e[0] > 2
    einv = interval(1 / e[1], 1 / e[0])
    bads = []
    for plan in reg['bad_component_rounds']:
        x = bad_round(plan, reg, ln2, einv)
        bads.append(x)
        print(json.dumps(dict(n=plan['n'], kind='bad_component',
                              S0_lower=x['continuous_atomic_S0']['lower_float'],
                              gap_lower=x['S0_minus_Z_over_e']['lower_float'],
                              finite_gap_lower=x['finite_family']['actual_finite_S0_minus_continuous_Z_over_e']['lower_float'],
                              L1_gap_lower=x['explicit_L1_analytic_factor_guard']['uniform_packet_profile_lower_minus_Z_over_e']['lower_float'])), flush=True)
    nus = []
    with localcontext() as ctx:
        ctx.prec = 160
        for plan in reg['nu_guard_rounds']:
            x = nu_round(plan, reg)
            nus.append(x)
            print(json.dumps(dict(n=plan['n'], kind='bounded_density_nu', D=x['D_nu'], Q=x['Q_nu'],
                                  E_cells=x['original_E_cell_count'], mass=x['full_nu_mass'])), flush=True)
    out = dict(status='passed', registration_sha256=sha(regpath), registration_amendment_sha256=sha(amendment),
               script_sha256=sha(Path(__file__)), analytic_proof_sha256=amend['analytic_proof_sha256'],
               ln2_interval=pack(ln2), e_interval=pack(e), inverse_e_interval=pack(einv),
               bad_component_rounds=bads, saved_cell_nu_rounds=nus, elapsed_seconds=time.perf_counter() - start,
               true_nearoptimal_space_input_created_or_evaluated=False,
               old_MC_or_oracle_rerun=False, no_parameter_or_sample_adaptation=True,
               scope='Explicit continuous atomic bad source and independent finite-family hot-point interval certificates; positive-width analytic lower-factor guard; symbolic conditional dilution; saved-cell bounded-density frozen-loggain diagnostics only.')
    target = BASE / (PREFIX + '_results.json')
    with target.open('x') as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
        f.write('\n')
    print(json.dumps(dict(status='passed', elapsed_seconds=out['elapsed_seconds'], results_sha256=sha(target))), flush=True)


if __name__ == '__main__':
    main()
