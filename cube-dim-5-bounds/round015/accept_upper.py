"""One registered exact acceptance/comparison; no H or old budget rerun.

Usage: python accept_upper.py REGISTRATION OUTPUT
"""
from fractions import Fraction as Q
import hashlib
import json
import os
from pathlib import Path
import sys
import time

P = Path(__file__).resolve().parent.parent


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    regpath, output = map(Path, sys.argv[1:])
    assert not output.exists()
    start = time.monotonic()
    reg = json.loads(regpath.read_bytes())
    assert reg['maximum_cases'] == 1
    for rel, digest in reg['source_sha256'].items():
        assert sha(P / rel) == digest, rel
    a = json.loads((P / 'round015/audit136_actual/results.json').read_bytes())
    gdoc = json.loads((P / 'round013/primary_globalg105_acceptance.json').read_bytes())
    old = json.loads((P / 'results/bounds_ledger.json').read_bytes())
    g = Q(gdoc['qualified_g'])
    assert g == Q(a['g']) == Q(21, 200000)
    I = Q(a['independent_complete_I_upper'])
    lo, hi = map(Q, a['independent_cost_interval'])
    fourth = 3125 * I**4 / g**5
    assert 0 < lo < hi < Q(799, 25) < 32
    assert lo**4 <= fourth <= hi**4
    assert old['upper']['exact'] == '32'
    assert I < Q(old['auxiliary_results_round014']['best_complete_same_alpha1_I_upper'])
    gate = 32**4 * g**5 - 3125 * I**4
    assert gate > 0 and gate == Q(a['independent_exact32_gate'])
    comparisons = {}
    for name in ('checked_public', 'opening_internal', 'continuation_reference'):
        comparisons[name] = {'old_upper': '32',
            'strict_fourth_power_margin': str(gate),
            'gain_interval': [str(32-hi), str(32-lo)],
            'simple_rational_upper': '799/25', 'simple_exact_gain': '1/25'}
    result = {'status': 'primary_accepts_complete_C5_upper',
        'dimension': 5, 'g': str(g), 'complete_I_upper': str(I),
        'best_upper': {'type': 'unique_positive_fourth_root',
            'fourth_power': str(fourth), 'interval': [str(lo), str(hi)],
            'display_only': a['cost_display_only'],
            'formula': '(3125*I_upper^4/(21/200000)^5)^(1/4)'},
        'simple_rational_upper': '799/25',
        'exact_I_improvement_vs_round014': str(Q(old['auxiliary_results_round014']['best_complete_same_alpha1_I_upper'])-I),
        'comparisons': comparisons,
        'same_model_general_L1_all_positive_common_radii': True,
        'primary_analytic_review': 'round015/primary_proof.md',
        'new_H_or_old_baseline_recomputations': 0,
        'lower_changed': False, 'lower_research': 'stopped_by_user',
        'dual_breakthrough_at': None, 'post_dual_rounds': 0,
        'source_sha256': reg['source_sha256'], 'registration_sha256': sha(regpath),
        'pid': os.getpid(), 'seconds': time.monotonic()-start}
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps({'status': result['status'], 'pid': result['pid'],
        'seconds': result['seconds'], 'upper_display_only': a['cost_display_only']}))


if __name__ == '__main__':
    main()
