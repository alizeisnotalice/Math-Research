#!/usr/bin/env python3
"""Exact Lebesgue-cell finite-scale algebra, numerical Decimal log variation.
No imports/calls of any previous MC or arrival oracle.
"""
from decimal import Decimal, localcontext
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from pathlib import Path
import gzip
import hashlib
import json
import math
import time

BASE = Path(__file__).resolve().parent
PREFIX = 'log_max_variation_probe_20261007'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def dec(x):
    x = F(x)
    return Decimal(x.numerator) / Decimal(x.denominator)


@lru_cache(maxsize=None)
def logplus(ratio):
    return dec(ratio).ln() if ratio > 1 else Decimal(0)


def dump_new(p, data):
    with p.open('x') as f:
        json.dump(data, f, ensure_ascii=False, indent=2)
        f.write('\n')


def round_run(plan, tvalues, precision, q):
    n = plan['n']
    atoms = [list(map(F, row)) for row in plan['atoms']]
    a, v = list(map(F, plan['weights'])), list(map(F, plan['v']))
    J, tau = sorted(map(F, plan['J'])), F(plan['tau'])
    assert sum(a) == 1 and sum(x * y for x, y in zip(a, v)) == 0
    assert all(abs(x) <= 1 for x in v) and all(x > 0 for x in a)
    weights_t = [[x * (1 + t * y) for x, y in zip(a, v)] for t in tvalues]
    for pert in weights_t:
        assert all(x > 0 for x in pert) and sum(pert) == 1
    cuts = [sorted({atom[i] + sign * R / 2 for atom in atoms for R in J for sign in (-1, 1)}) for i in range(n)]
    intervals = [[(lo, hi) for lo, hi in zip(axis, axis[1:])] for axis in cuts]
    scale_power = [R**n for R in J]
    S, Q, Dcell, Evolume = [F(0)] * len(a), F(0), F(0), F(0)
    support_volume, bounding_volume, plateau_volume, qvolume = F(0), F(0), F(0), F(0)
    phi0, phit = Decimal(0), [Decimal(0) for _ in tvalues]
    records = []
    counts = dict(cells=0, positive_response=0, zero_response=0, E=0, threshold_plateau=0,
                  multiwinner=0, E_multiwinner=0, threshold_plateau_multiwinner=0,
                  strict_q_threshold=0, q_threshold_plateau=0)
    perturb_stats = [dict(t=str(t), E_cell_count=0, new_E_cell_count=0, lost_E_cell_count=0,
                         winner_changes=0, exact_frozen_winner_dominance_checks=0) for t in tvalues]
    for cell in product(*intervals):
        mid = [(lo + hi) / 2 for lo, hi in cell]
        volume = math.prod(hi - lo for lo, hi in cell)
        bounding_volume += volume
        captures = [[j for j, atom in enumerate(atoms) if all(abs(x - y) <= R / 2 for x, y in zip(mid, atom))] for R in J]
        mass = [sum((a[j] for j in cap), F(0)) for cap in captures]
        signed = [sum((a[j] * v[j] for j in cap), F(0)) for cap in captures]
        response = [m / power for m, power in zip(mass, scale_power)]
        M = max(response)
        winner = response.index(M)  # J sorted: minimum R among exact ties.
        E = M > tau
        p = signed[winner] / mass[winner] if M else F(0)
        assert abs(p) <= 1
        counts['cells'] += 1
        counts['positive_response'] += M > 0
        counts['zero_response'] += M == 0
        counts['E'] += E
        counts['threshold_plateau'] += M == tau
        nties = sum(r == M for r in response) if M else 0
        counts['multiwinner'] += nties > 1
        counts['E_multiwinner'] += E and nties > 1
        counts['threshold_plateau_multiwinner'] += M == tau and nties > 1
        counts['strict_q_threshold'] += M > q * tau
        counts['q_threshold_plateau'] += M == q * tau
        if M:
            support_volume += volume
        if E:
            Evolume += volume
            for j in captures[winner]:
                S[j] += tau * volume / mass[winner]
            Dcell += tau * volume * p
            Q += tau * volume * p * p
        if M == tau:
            plateau_volume += volume
        if M > q * tau:
            qvolume += volume
        phi0 += dec(tau * volume) * logplus(M / tau)
        perturb = []
        for ti, t in enumerate(tvalues):
            new_mass = [m + t * mv for m, mv in zip(mass, signed)]
            assert all(m >= 0 for m in new_mass)
            new_response = [m / power for m, power in zip(new_mass, scale_power)]
            Mt = max(new_response)
            wi = new_response.index(Mt)
            assert new_response[winner] == M * (1 + t * p)
            assert Mt >= M * (1 + t * p)
            perturb_stats[ti]['exact_frozen_winner_dominance_checks'] += 1
            Et = Mt > tau
            perturb_stats[ti]['E_cell_count'] += Et
            perturb_stats[ti]['new_E_cell_count'] += Et and not E
            perturb_stats[ti]['lost_E_cell_count'] += E and not Et
            perturb_stats[ti]['winner_changes'] += wi != winner
            phit[ti] += dec(tau * volume) * logplus(Mt / tau)
            perturb.append(dict(t=str(t), full_mass_by_scale=[str(m) for m in new_mass],
                                M=str(Mt), winner_R=str(J[wi]), E=Et))
        records.append(dict(cell_bounds=[[str(lo), str(hi)] for lo, hi in cell],
                            volume=str(volume), midpoint=[str(x) for x in mid],
                            full_capture_labels_by_scale=captures, full_mass_by_scale=[str(m) for m in mass],
                            full_response_by_scale=[str(r) for r in response], M=str(M), winner_R=str(J[winner]),
                            maximum_tie_count=nties, E=E, threshold_plateau=M == tau,
                            frozen_winner_p=str(p), perturbations=perturb))
    first = sum(x * s for x, s in zip(a, S))
    D = sum(x * y * s for x, y, s in zip(a, v, S))
    assert first == tau * Evolume and D == Dcell
    assert bounding_volume == math.prod(axis[-1] - axis[0] for axis in cuts)
    coincident_pairs = [(i, j) for i in range(len(a)) for j in range(i + 1, len(a)) if atoms[i] == atoms[j]]
    assert all(S[i] == S[j] for i, j in coincident_pairs)
    checks = []
    for t, phi, stats in zip(tvalues, phit, perturb_stats):
        rhs = t * D - t * t * Q / (2 * (1 - abs(t))**2)
        delta = phi - phi0
        residual = delta - dec(rhs)
        checks.append(dict(**stats, Phi_t=str(phi), Phi_base=str(phi0), Phi_delta=str(delta),
                           rational_RHS=str(rhs), RHS_decimal=str(dec(rhs)),
                           variation_lower_bound_residual=str(residual),
                           numerical_inequality_observed=residual >= 0, interval_certified=False))
        assert residual >= 0, '160-digit numerical lower-bound failure; retain raw data'
    optional = dec(tau * qvolume) * dec(q).ln()
    assert phi0 >= optional
    payload = dict(input=plan, finite_scale_only=True, cells=records)
    profile = BASE / (PREFIX + f'_n{n}_cell_profiles.json.gz')
    with gzip.open(profile, 'xt', encoding='utf8') as f:
        json.dump(payload, f, separators=(',', ':'))
    return dict(n=n, atom_labels=len(a), physical_atom_count=len({tuple(x) for x in atoms}),
                finite_scales=[str(R) for R in J], tau=str(tau), cut_count_per_axis=[len(axis) for axis in cuts],
                counts=counts, bounding_volume=str(bounding_volume), response_support_volume=str(support_volume),
                original_E_volume=str(Evolume), threshold_plateau_volume=str(plateau_volume),
                source_profile_S=[str(s) for s in S], coincident_label_pairs=coincident_pairs,
                sum_aS=str(first), tau_E_volume=str(tau * Evolume), D_source=str(D), D_receiver=str(Dcell),
                Q=str(Q), exact_algebra_checks_passed=True, Phi_base=str(phi0),
                variation_checks=checks, q_optional=dict(q=str(q), strict_level_volume=str(qvolume),
                    lower_bound=str(optional), residual=str(phi0 - optional),
                    numerical_inequality_observed=True, interval_certified=False),
                Decimal_precision=precision, cell_profiles_path=str(profile), cell_profiles_sha256=sha(profile))


def main():
    registration = BASE / (PREFIX + '_registration.json')
    reg = json.loads(registration.read_text())
    started = time.perf_counter()
    with localcontext() as ctx:
        ctx.prec = reg['Decimal_precision']
        tvalues = list(map(F, reg['t_values']))
        rounds = []
        for plan in reg['plans']:
            result = round_run(plan, tvalues, ctx.prec, F(reg['q_optional']))
            rounds.append(result)
            print(json.dumps(dict(n=plan['n'], cells=result['counts']['cells'], E_cells=result['counts']['E'],
                                  plateau_cells=result['counts']['threshold_plateau'], tie_cells=result['counts']['multiwinner'],
                                  minimum_variation_residual=min(Decimal(c['variation_lower_bound_residual']) for c in result['variation_checks']).__str__())), flush=True)
    out = dict(status='passed', registration_sha256=sha(registration), script_sha256=sha(Path(__file__)),
               rounds=rounds, elapsed_seconds=time.perf_counter() - started, old_MC_or_oracle_called=False,
               exact_Lebesgue_cell_algebra=True, logarithms_interval_certified=False,
               general_proof_or_dimension_order_claim=False,
               scope='Fixed full finite positive atoms and common finite scales only; exact Fraction algebra and non-interval Decimal variation diagnostics.')
    target = BASE / (PREFIX + '_results.json')
    dump_new(target, out)
    print(json.dumps(dict(status='passed', elapsed_seconds=out['elapsed_seconds'], results_sha256=sha(target))), flush=True)


if __name__ == '__main__':
    main()
