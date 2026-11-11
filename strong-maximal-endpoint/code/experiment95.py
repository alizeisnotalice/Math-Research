"""Exact finite checks for the r095 overlap and periodic-density contracts.

The overlap and m-by-k-copy models are finite independent Bernoulli point
models. They check the algebra used in core95; they are not simulations or a
high-dimensional physical source construction.
"""
import argparse
import hashlib
import json
import multiprocessing as mp
import os
import time
from fractions import Fraction as Q
from pathlib import Path

P_VALUES = (Q(1, 31), Q(1, 127))
M_VALUES = tuple(range(1, 17))
K_VALUES = tuple(range(1, 9))
CONTRACT_TEXT = 'theta/H_N(p) <= d_N(p,theta) <= -log(1-theta)/H_N(p)'


def qstr(value):
    return str(Q(value))


def _sha_rows(rows):
    text = '\n'.join('|'.join(map(str, row)) for row in rows)
    return hashlib.sha256(text.encode('ascii')).hexdigest()


def overlap_check():
    p = Q(1, 31)
    q = 1-p
    rows = []
    raw_mass = Q(0)
    surviving_mass = Q(0)
    lost_mass = Q(0)
    pair_bound_mass = Q(0)
    for mask in range(64):
        bits = tuple((mask >> shift) & 1 for shift in range(5, -1, -1))
        probability = Q(1)
        for bit in bits:
            probability *= p if bit else q
        raw = sum(bits)
        colors = tuple(bits[2*i]+bits[2*i+1] > 0 for i in range(3))
        surviving = sum(int(colors[i] and not any(colors[:i]+colors[i+1:]))
                        for i in range(3))
        lost = raw-surviving
        pair_count = raw*(raw-1)
        assert 0 <= lost <= pair_count
        raw_mass += probability*raw
        surviving_mass += probability*surviving
        lost_mass += probability*lost
        pair_bound_mass += probability*pair_count
        rows.append((''.join(map(str, bits)), qstr(probability), raw,
                     surviving, lost, pair_count))

    ordered_pair_bound = Q(30)*p*p
    assert len(rows) == 64
    assert raw_mass == Q(6, 31)
    assert surviving_mass == Q(148230000, 887503681)
    assert lost_mass == Q(23544906, 887503681)
    assert pair_bound_mass == ordered_pair_bound == Q(30, 961)
    assert lost_mass <= pair_bound_mass
    return {
        'model': '3 colors, 2 independent Bernoulli(p) copies per color at one point',
        'p': qstr(p),
        'patterns_enumerated': len(rows),
        'pattern_rows_sha256': _sha_rows(rows),
        'raw_mass_expectation': qstr(raw_mass),
        'surviving_mass_expectation': qstr(surviving_mass),
        'lost_mass_expectation': qstr(lost_mass),
        'expected_ordered_pair_bound': qstr(pair_bound_mass),
        'exact_loss_le_pair_bound': lost_mass <= pair_bound_mass,
        'pointwise_loss_le_l_times_l_minus_1': all(row[4] <= row[5] for row in rows),
    }


def simultaneous_margin_check():
    c, eta, delta = Q(1, 1024), Q(1, 128), Q(1, 1024)
    lam = Q(5, 4)*c*delta
    minimum_expected_success = Q(5, 8)*(1-Q(5, 4)*c)
    packing_success_lower = ((minimum_expected_success-Q(3, 10))/Q(7, 10)
                             -4*lam)
    guard_success_lower = 1-38*eta
    # T1=4P/eta, so the endpoint-fragment correction 2P/T1 is eta/2.
    long_edge_deleted_lower = 2*eta-eta/2
    assert minimum_expected_success == Q(20455, 32768)
    assert packing_success_lower == Q(485691, 1048576)
    assert guard_success_lower == Q(45, 64)
    assert long_edge_deleted_lower == Q(3, 256) and long_edge_deleted_lower > eta
    return {
        'c': qstr(c), 'eta': qstr(eta), 'delta': qstr(delta),
        'lambda_upper': qstr(lam),
        'minimum_expected_success_fraction': qstr(minimum_expected_success),
        'packing_simultaneous_success_probability_lower': qstr(packing_success_lower),
        'guard_simultaneous_success_probability_lower': qstr(guard_success_lower),
        'long_edge_deleted_fraction_lower': qstr(long_edge_deleted_lower),
        'strictly_exceeds_eta': long_edge_deleted_lower > eta,
        'derived_with_exact_fraction_arithmetic': True,
    }


def generalized_point_model():
    rows = []
    minimum_survival = None
    maximum_survival = None
    minimum_pair_margin = None
    maximum_loss_ratio = None
    for m in M_VALUES:
        for k in K_VALUES:
            for p in P_VALUES:
                raw = Q(m*k)*p
                survival = Q(m)*(1-(1-p)**k)*(1-p)**(k*(m-1))
                loss = raw-survival
                pair_bound = Q(m*k*(m*k-1))*p*p
                assert 0 <= loss <= pair_bound
                ratio = Q(0) if pair_bound == 0 else loss/pair_bound
                row = (m, k, qstr(p), qstr(raw), qstr(survival),
                       qstr(loss), qstr(pair_bound))
                rows.append(row)
                if minimum_survival is None or survival < minimum_survival['value']:
                    minimum_survival = {'value': survival, 'm': m, 'k': k, 'p': p}
                if maximum_survival is None or survival > maximum_survival['value']:
                    maximum_survival = {'value': survival, 'm': m, 'k': k, 'p': p}
                margin = pair_bound-loss
                if minimum_pair_margin is None or margin < minimum_pair_margin['value']:
                    minimum_pair_margin = {'value': margin, 'm': m, 'k': k, 'p': p}
                if maximum_loss_ratio is None or ratio > maximum_loss_ratio['value']:
                    maximum_loss_ratio = {'value': ratio, 'm': m, 'k': k, 'p': p}

    assert len(rows) == 16*8*2 == 256
    assert all(Q(row[6])-Q(row[5]) >= 0 for row in rows)
    def extremum(record):
        return {'value': qstr(record['value']), 'm': record['m'],
                'k': record['k'], 'p': qstr(record['p'])}
    return {
        'model': 'independent point coverage by m colors, k copies per color',
        'm_values': [min(M_VALUES), max(M_VALUES)],
        'k_values': [min(K_VALUES), max(K_VALUES)],
        'p_values': [qstr(p) for p in P_VALUES],
        'case_count': len(rows),
        'canonical_case_rows_sha256': _sha_rows(rows),
        'raw_mass_formula': 'm*k*p',
        'surviving_mass_formula': 'm*(1-(1-p)^k)*(1-p)^(k*(m-1))',
        'loss_upper_formula': 'm*k*(m*k-1)*p^2',
        'loss_le_pair_bound_all_cases': True,
        'minimum_surviving_mass': extremum(minimum_survival),
        'maximum_surviving_mass': extremum(maximum_survival),
        'minimum_pair_bound_minus_loss': extremum(minimum_pair_margin),
        'maximum_loss_over_pair_bound': extremum(maximum_loss_ratio),
    }


def log_contract():
    return {
        'statement': CONTRACT_TEXT,
        'domain': '0 < theta < 1; H_N(p) > 0',
        'exact_ordering_basis': '-log(1-theta)-theta >= 0, divided by positive H_N(p)',
        'proof': ('g(theta)=-log(1-theta)-theta has g(0)=0 and '
                  'g\'(theta)=theta/(1-theta)>0 on (0,1)'),
        'numeric_log_evaluation_used': False,
        'theorem_contract_only_not_a_finite_source_computation': True,
    }


def run_job(name):
    if name == 'overlap':
        return name, overlap_check(), os.getpid()
    if name == 'constants':
        return name, simultaneous_margin_check(), os.getpid()
    if name == 'generalization':
        return name, generalized_point_model(), os.getpid()
    raise ValueError(name)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('output_directory')
    ap.add_argument('--workers', type=int, default=2)
    args = ap.parse_args()
    assert 1 <= args.workers <= 2
    outdir = Path(args.output_directory)
    outdir.mkdir(parents=True, exist_ok=True)
    start = time.monotonic()
    with mp.get_context('spawn').Pool(args.workers) as pool:
        results = pool.map(run_job, ('overlap', 'constants', 'generalization'))
    sections = {name: value for name, value, _pid in results}
    pids = sorted({pid for _name, _value, pid in results})
    data = {
        'schema': 'r095-finite-exact-checks-v1',
        'pass_all': True,
        'configured_workers': args.workers,
        'actual_workers': len(pids),
        'worker_pids': pids,
        **sections,
        "log_contract": log_contract(),
        'seconds': time.monotonic()-start,
        'scope': ('Exact finite checks for the stated single-point independent overlap models, '
                  'the audit95 rational constants, and the 256-case algebraic loss bound. '
                  'These are not high-dimensional physical source constructions or onset estimates.'),
        'high_dimensional_source_construction_claimed': False,
        'finite_dimension_onset_claimed': False,
        'general_H_N_coefficient_claimed': False,
    }
    (outdir/'certificates95.json').write_text(json.dumps(data, indent=2)+'\n')
    print(json.dumps({key: value for key, value in data.items()
                      if key not in ('worker_pids',)}, indent=2))


if __name__ == '__main__':
    main()
