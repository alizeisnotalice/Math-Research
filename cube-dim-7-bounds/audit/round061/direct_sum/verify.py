#!/usr/bin/env python3
"""New exact finite-cell measure examples for the direct-sum lemma only."""

from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import platform


def inspect(name, ambient, actual, layers, tau, expect_disjoint):
    atoms = tuple(ambient)
    assert set(actual) == set(atoms)
    assert all(0 <= actual[x] <= ambient[x] for x in atoms)
    nu = {x: ambient[x] - actual[x] for x in atoms}
    support_lists = [S for A, P in layers for S in (A, P)]
    disjoint = all(not (support_lists[i] & support_lists[j])
                   for i in range(len(support_lists))
                   for j in range(i + 1, len(support_lists)))
    assert disjoint == expect_disjoint
    rows = []
    lower_sum = F(0)

    for A, P in layers:
        K = {(w, z): F(w in A and z in P) for w in atoms for z in atoms}
        row_lambda = {w: sum(K[w, z] * ambient[z] for z in atoms) for w in A}
        row_mu = {w: sum(K[w, z] * actual[z] for z in atoms) for w in A}
        col_lambda = {z: sum(K[w, z] * ambient[w] for w in atoms) for z in P}
        r = min(row_lambda.values())
        D = max(tuple(row_lambda.values()) + tuple(col_lambda.values()))
        assert D > 0 and all(v <= tau for v in row_mu.values())

        def energy(left, right):
            return sum(K[w, z] * left[w] * right[z] for w in atoms for z in atoms)

        ell, emm = energy(ambient, ambient), energy(actual, actual)
        identity_rhs = energy(nu, ambient) + energy(actual, nu)
        local_deficit = sum(nu[x] for x in A | P)
        assert ell - emm == identity_rhs
        assert identity_rhs <= D * local_deficit
        lower = sum(ambient[w] for w in A) * (r - tau) / D
        assert lower <= local_deficit
        lower_sum += lower
        rows.append({'r': str(r), 'D': str(D), 'E_lambda_lambda': str(ell),
                     'E_mu_mu': str(emm), 'identity_rhs': str(identity_rhs),
                     'local_deficit': str(local_deficit), 'lower': str(lower)})

    total_deficit = sum(nu.values())
    if expect_disjoint:
        assert lower_sum <= total_deficit
    else:
        assert lower_sum > total_deficit
    return {'name': name, 'pairwise_point_supports_disjoint': disjoint,
            'tau': str(tau), 'layers': rows, 'lower_sum': str(lower_sum),
            'total_deficit': str(total_deficit)}


tau = F(1, 4)
layers = [({'a1'}, {'p1'}), ({'a2'}, {'p2'})]
ambient = {x: F(1) for x in ('a1', 'a2', 'p1', 'p2')}
actual = {'a1': F(1), 'a2': F(1), 'p1': tau, 'p2': tau}
base = inspect('disjoint; common items cross both bands', ambient, actual, layers, tau, True)
assert F(base['lower_sum']) == F(base['total_deficit']) == F(3, 2)

# One target item spans both A bands; one source item spans both P bands.
items = [{'a1': F(1), 'a2': F(1)}, {'p1': tau, 'p2': tau}]
for x in ambient:
    assert sum(item.get(x, F(0)) for item in items) == actual[x]
for item in items:
    restricted = sum(item.get(x, F(0)) for A, P in layers for x in A | P)
    assert restricted == sum(item.values())
assert all(sum(value > 0 for value in item.values()) == 2 for item in items)

c = F(3, 2)
scaled_ambient = {x: c for x in ambient}
scaled_actual = {'a1': c, 'a2': c, 'p1': tau, 'p2': tau}
scaled = inspect('disjoint finite-a scaling c=3/2', scaled_ambient, scaled_actual, layers, tau, True)
assert F(scaled['lower_sum']) == F(scaled['total_deficit']) == F(5, 2)
assert F(scaled['lower_sum']) >= F(base['lower_sum'])

overlap = inspect('source support overlap counterexample',
                  {x: F(1) for x in ('a1', 'a2', 'p')},
                  {'a1': F(1), 'a2': F(1), 'p': tau},
                  [({'a1'}, {'p'}), ({'a2'}, {'p'})], tau, False)
assert F(overlap['lower_sum']) == F(3, 2)
assert F(overlap['total_deficit']) == F(3, 4)

small_c = F(1, 2)
shrink = inspect('c<1 cannot preserve base coefficient',
                 {'a': small_c, 'p': small_c}, {'a': small_c, 'p': tau},
                 [({'a'}, {'p'})], tau, True)
assert F(shrink['lower_sum']) == F(shrink['total_deficit']) == F(1, 4)
assert F(1) * (1 - tau) > F(shrink['total_deficit'])

source = Path('/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round056/refined_inverse/PROOF.md')
source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
assert source_hash == '6986ed9ec48b46f3136a9a56722d0092db4f2942d52e35c9f4de37eb279f399a'
root = Path(__file__).resolve().parent
result = {'status': 'all new exact measure checks passed', 'python': platform.python_version(),
          'method': 'finite unit-cell/atom kernels; exact Fraction rows, columns and two-measure energies',
          'cases': [base, scaled, overlap, shrink],
          'cross_band_item_partition': 'checked exactly for two multi-band item measures',
          'source': {'path': str(source), 'sha256': source_hash, 'read_scope': 'Section 8 scaling identification only; not full R056 audit'},
          'artifact_sha256': {name: hashlib.sha256((root / name).read_bytes()).hexdigest()
                              for name in ('PROOF.md', 'verify.py')},
          'scope': 'abstract measure-lemma checks, not actual cube history realizations or root/tuple/operator certification'}
(root / 'verification.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({'status': result['status'], 'cases': len(result['cases']),
                  'disjoint_exact_deficit': base['total_deficit'],
                  'scaled_exact_deficit': scaled['total_deficit'],
                  'overlap_claimed_vs_actual': [overlap['lower_sum'], overlap['total_deficit']]}, indent=2))
