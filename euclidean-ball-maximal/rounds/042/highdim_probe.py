"""Exploratory axial-source Monte Carlo; not an interval or global-bound certificate.
Set alpha=1/v_n and proposal density q(x)=coverage(x)/(N v_n) on the
union of unit balls. Every strict superlevel witness has radius <1,
so this proposal covers the entire original band. Source integrals remain
weighted by the actual atomic masses. Canonical clocks use ratio4.
"""
import argparse
import json
import math
from pathlib import Path
import numpy as np


def evaluate(n, positions, masses, samples, seed):
    rng = np.random.default_rng(seed)
    positions = np.asarray(positions, dtype=float)
    masses = np.asarray(masses, dtype=float)
    assert abs(masses.sum()-1) < 1e-12 and masses.min() > 0
    count = len(masses)
    # Uniform-ball sampling uses independent Gaussian axis and chi-square
    # transverse norm; no independent-coordinate assumption on the input P.
    center = rng.integers(count, size=samples)
    axis = rng.normal(size=samples)
    transverse = rng.chisquare(n-1, size=samples)
    norm2 = axis*axis+transverse
    radius = np.exp(np.log(rng.uniform(size=samples))/n)
    z = positions[center]+radius*axis/np.sqrt(norm2)
    rho2 = radius*radius*transverse/norm2
    distance2 = (z[:, None]-positions[None, :])**2+rho2[:, None]
    log_r = 0.5*np.log(distance2)
    coverage = (distance2 < 1).sum(axis=1)
    assert coverage.min() >= 1
    importance = count/coverage
    order = np.argsort(distance2, axis=1)
    cumulative = np.cumsum(masses[order], axis=1)
    ordered_log_r = np.take_along_axis(log_r, order, axis=1)
    log_maximum = np.max(np.log(cumulative)-n*ordered_log_r, axis=1)
    band = (log_maximum > 0) & (log_maximum <= math.log(2))
    D = math.ceil(math.log2(8/masses.min()))+1
    J = np.zeros(samples, dtype=int)
    K = np.zeros(samples, dtype=int)
    radius_members = {}
    for j in range(1, D+1):
        volume = 8*2.0**(-j)
        member = log_r < math.log(volume)/n
        radius_members[j] = member
        mass = member @ masses
        J[(J == 0) & (mass > volume/8)] = j
        K[(K == 0) & (mass > volume/2)] = j
    assert np.all(J[band] > 0) and np.all(K[band] > 0)
    eligible = band & (J >= 2)
    h = np.zeros((count, D+1))
    X = np.zeros(D+1)
    observations = np.zeros(D+1, dtype=int)
    for k in range(2, D+1):
        mask = eligible & (K == k)
        observations[k] = int(mask.sum())
        weight = importance*mask
        X[k] = weight.mean()
        h[:, k] = np.sum(radius_members[k]*weight[:, None], axis=0)/samples/(8*2.0**(-k))
    gram = (h.T*masses)@h
    row = np.triu(gram, 1).sum(axis=1)
    # Ratios for poorly observed layers are explicitly excluded from ranking.
    accepted = np.flatnonzero(observations >= 100)
    best = int(max(accepted, key=lambda k: row[k]/X[k])) if len(accepted) else None
    return {'n': n, 'N': count, 'seed': seed, 'samples': samples, 'D': D,
            'alpha_normalization': 'alpha*v_n=1', 'a_over_alpha': 0.125, 'b_over_alpha': 0.5,
            'full_band_X_estimate': float(np.mean(importance*band)),
            'eligible_X_estimate': float(X.sum()),
            'Q_over_X_estimate': float(gram.sum()/X.sum()) if X.sum() else None,
            'best_observed_row': best,
            'best_row_ratio_estimate': float(row[best]/X[best]) if best is not None else None,
            'layer_rows': [{'j': int(j), 'observation_count': int(observations[j]),
                            'X_estimate': float(X[j]), 'row_estimate': float(row[j]),
                            'row_ratio_estimate': float(row[j]/X[j]) if X[j] else None,
                            'used_in_ranking': bool(observations[j] >= 100)}
                           for j in range(2, D+1) if observations[j]],
            'full_spatial_proposal_coverage_proved': True,
            'quadratic_estimator_has_finite_sample_product_bias': True,
            'confidence_or_interval_certificate': False}


def cases():
    for n in (2, 8, 32, 128, 512):
        yield 1, 'separated_dyadic_baseline', n, [0, 10, 20], [.25, .25, .5]
    for n in (2, 8, 32, 128, 512):
        for separation in (.5, 2.0):
            yield 2, f'two_near_c{separation}', n, [0, separation/math.sqrt(n), 10], [1/32, 13/32, 9/16]
    for n in (8, 32, 128, 512, 2048):
        for spacing in (.25, 1.0):
            near = [-spacing*i/math.sqrt(n) for i in range(8)]
            masses = [2.0**(-i-1) for i in range(8)]
            masses = [w*(7/16)/sum(masses) for w in masses]
            yield 3, f'correlated_cloud_c{spacing}', n, near+[10], masses+[9/16]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--samples', type=int, default=32768)
    args = parser.parse_args()
    rows = []
    for case_index, (batch, name, n, x, w) in enumerate(cases()):
        replicates = [evaluate(n, x, w, args.samples, 420711+100*case_index+rep) for rep in range(2)]
        record = {'batch': batch, 'family': name, 'positions': x, 'masses': w, 'replicates': replicates}
        if batch == 1:
            # Analytic baseline: X=1/2, Q/X=1/2, all off-diagonal pairs zero.
            assert all(abs(v['eligible_X_estimate']-.5) < .04 for v in replicates)
            assert all(abs(v['Q_over_X_estimate']-.5) < .08 for v in replicates)
            assert all(v['best_row_ratio_estimate'] == 0 for v in replicates)
            record['analytic_baseline'] = {'X': '1/2', 'Q_over_X': '1/2', 'all_rows': 0}
        rows.append(record)
        print(batch, name, n, [r['best_row_ratio_estimate'] for r in replicates], flush=True)
    result = {'status': 'passed', 'cases': rows, 'configuration_count': len(rows),
              'replicate_count': 2*len(rows), 'arithmetic': 'float64 exploratory Monte Carlo',
              'main_dimension_independent_bound_proved': False,
              'uniform_row_bound_proved': False,
              'warning': 'No confidence certificate; small observed layers excluded, unobserved layers not bounded. Product bias remains.'}
    args.output.write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
