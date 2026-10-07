#!/usr/bin/env python3
"""Read-only reconstruction of saved source/receiver/capture and RN profiles.
Does not recompute an oracle maximum or sample any new point.
"""
import gzip
import hashlib
import json
import math
from fractions import Fraction as F
from pathlib import Path

BASE = Path(__file__).resolve().parent


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    reg = json.loads((BASE / 'registration.json').read_text())
    results = json.loads((BASE / 'results.json').read_text())
    rows_out = []
    for plan, summary in zip(reg['rounds'], results['rounds']):
        n = plan['n']
        inp = reg['source_inputs'][str(n)]['input']
        assert sha(Path(summary['source_profiles_path'])) == summary['source_profiles_sha256']
        Ns, Nr = plan['source_samples'], plan['receivers_per_replica']
        all1, all2, mean1, mean2 = [], [], [], []
        calls, captured, Ecount = 0, 0, 0
        coincidence = {}
        with gzip.open(summary['source_profiles_path'], 'rt') as stream:
            for expected, line in enumerate(stream):
                row = json.loads(line)
                assert row['source_index'] == expected
                source_k = row['source_k']
                codes = row['source_integer_codes']
                assert row['source_y'] == [str(F(c, source_k)) for c in codes]
                assert all(-inp['A'] * source_k <= c <= inp['A'] * source_k for c in codes)
                membership = [all(c * k % source_k == 0 for c in codes) for k in inp['ks']]
                assert membership == row['source_grid_membership']
                coincidence[str(sum(membership))] = coincidence.get(str(sum(membership)), 0) + 1
                repmean = []
                for replica in row['receiver_replicas']:
                    assert len(replica) == Nr
                    scores = []
                    for z in replica:
                        calls += 1
                        D, key = z['dyadic_denominator'], z['R_key']
                        assert D <= key <= 2 * D
                        dx = [F.from_float(float.fromhex(v)) for v in z['offsets_hex']]
                        assert (2 * max(abs(v) for v in dx) <= F(key, D)) == z['source_captured']
                        if z['core_or_shell'] == 'shell':
                            assert max(abs(v) for v in dx) == F.from_float(z['r']) / 2
                            assert dx[z['face']] == z['sign'] * F.from_float(z['r']) / 2
                        else:
                            assert max(abs(v) for v in dx) <= F(1, 2)
                        x = [F(c, source_k) + v for c, v in zip(codes, dx)]
                        R = F(key, D)
                        counts = []
                        for k in inp['ks']:
                            cr = []
                            for xi in x:
                                lo = k * (xi - R / 2)
                                hi = k * (xi + R / 2)
                                low = max(-inp['A'] * k, -((-lo.numerator) // lo.denominator))
                                high = min(inp['A'] * k, hi.numerator // hi.denominator)
                                cr.append(max(0, high - low + 1))
                            counts.append(cr)
                        assert counts == z['full_capture_counts']
                        m = sum(F(w) * F(math.prod(cr), q**n) for w, cr, q in zip(inp['weights'], counts, inp['qs']))
                        assert m == F(z['full_mass_numerator'], z['full_mass_denominator'])
                        M, tau = m / R**n, F(summary['tau'])
                        assert (M > tau) == z['E']
                        assert abs(float(tau / M) - z['g']) < 1e-14 if z['E'] else z['g'] == 0
                        free = z['free_coordinate_indices']
                        assert free == [i for i in range(n) if i != z['face']]
                        high = []
                        for i in free:
                            lo, hi = x[i] - F(3, 4), x[i] + F(3, 4)
                            low = max(-inp['A'], -((-lo.numerator) // lo.denominator))
                            up = min(inp['A'], hi.numerator // hi.denominator)
                            high.append(int(max(0, up - low + 1) == 2))
                        assert high == z['free_high_status']
                        T = float(F(plan['tilt_odds_T']))
                        logZ = sum(math.log(1 - p + p * T) for p in z['free_p_i'])
                        logL = sum(high) * math.log(T) - logZ
                        RN = 1 / (.25 + .75 * math.exp(logL))
                        assert math.isclose(logL, z['logL'], abs_tol=1e-12)
                        assert math.isclose(RN, z['RN'], rel_tol=1e-12)
                        geometric = (z['r'] / z['R_float'])**n if z['core_or_shell'] == 'shell' else z['R_float']**(-n)
                        score = z['g'] * summary['Cn'] * geometric * RN if z['E'] and z['source_captured'] else 0
                        assert math.isclose(score, z['score'], rel_tol=1e-12, abs_tol=1e-15)
                        assert 0 <= score <= 4 * summary['Cn'] * (1 + 1e-12)
                        scores.append(z['score'])
                        captured += z['source_captured']
                        Ecount += z['E']
                    repmean.append(sum(scores) / Nr)
                a, b = repmean
                assert a == row['S_replica_1'] and b == row['S_replica_2']
                assert a * b == row['I2_cross_score']
                assert (a + b) / 2 == row['S_estimate']
                all1.append(a)
                all2.append(b)
                mean1.append((a + b) / 2)
                mean2.append(a * b)
        assert len(mean1) == Ns and calls == 2 * Ns * Nr
        assert math.isclose(sum(mean1) / Ns, summary['I1_over_W']['mean'], rel_tol=1e-14)
        assert math.isclose(sum(mean2) / Ns, summary['I2_over_W']['mean'], rel_tol=1e-14)
        assert captured == summary['diagnostics']['source_captured']
        assert Ecount == summary['diagnostics']['E']
        rows_out.append(dict(n=n, complete_source_profiles=Ns, reconstructed_receivers=calls,
                             original_source_coincident_component_memberships=coincidence,
                             all_saved_capture_and_RN_checks_passed=True))
    out = dict(status='passed', results_sha256=sha(BASE / 'results.json'),
               review_script_sha256=sha(Path(__file__)), rounds=rows_out,
               oracle_reexecution=False, new_random_samples=False,
               scope='Exact saved receiver/capture/full-source/RN reconstruction, not a continuous floating-law interval certificate.')
    target = BASE / 'saved_profile_review.json'
    with target.open('x') as f:
        json.dump(out, f, indent=2)
        f.write('\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
