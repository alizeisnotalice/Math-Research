#!/usr/bin/env python3
"""Pure B015 preparation checks. This file NEVER calls the native engine."""
import sys
sys.dont_write_bytecode = True
import ast
import json
import os
import subprocess
from fractions import Fraction as Q
from pathlib import Path
from time import perf_counter
import prepare as prep
import search

HERE = Path(__file__).resolve().parent

def main():
    tic = perf_counter()
    assert not (HERE/'release.json').exists(), 'This checker is only for the unreleased preparation state.'
    assert not (HERE/'checkpoint.json').exists()
    assert not (HERE/'native_outputs').exists()
    reg = json.loads((HERE/'preparation.json').read_text())
    for name in ('prepare.py', 'search.py', 'check_preparation.py'):
        ast.parse((HERE/name).read_text())
    for rec in reg['inputs']:
        assert prep.sha(rec['path']) == rec['sha256']
    cones, rays = prep.fan_data()
    forms = prep.all_zero_forms(rays)
    assert len(cones) == 21 and len(rays) == 18 and len(forms) == 9308
    for seed in reg['seeds']:
        pf = prep.preflight(seed['numerators'], seed['denominator'], forms, rays)
        assert pf == seed['preflight']
        assert pf['feasible'] and pf['int64_safe'] and pf['zero_form_count'] == 0
    s1 = reg['seeds'][0]
    assert all(abs(Q(x)) <= Q(1, 2*10**14) for x in s1['coordinate_rounding_errors'])
    assert any(Q(x) for x in s1['coordinate_rounding_errors'])
    assert prep.nearest(Q(1, 2)) == 0 and prep.nearest(Q(3, 2)) == 2
    assert prep.nearest(Q(-3, 2)) == -2 and prep.nearest(Q(-1, 2)) == 0
    assert search.canonical_key(reg['frame_sha256'], s1['numerators'], s1['denominator']) == \
        search.canonical_key(reg['frame_sha256'], [2*n for n in s1['numerators']], 2*s1['denominator'])
    shifted = s1['numerators'][:]
    shifted[-1] += 1
    assert search.canonical_key(reg['frame_sha256'], s1['numerators'], s1['denominator']) != \
        search.canonical_key(reg['frame_sha256'], shifted, s1['denominator'])
    known_zero = prep.preflight([prep.DEN]*6, prep.DEN, forms, rays)
    assert known_zero['zero_form_count'] > 0
    # Pure six-gap arithmetic: both end gaps on the lower bound move inward.
    boundary = [prep.DEN, prep.DEN-prep.GAP, 8*10**13, 5*10**13, 2*10**13, prep.GAP]
    repair_checks = []
    for rec in prep.gap_repairs(boundary, prep.DEN):
        assert sum(rec['gap_delta_numerators']) == 0
        assert sum(rec['candidate_gap_numerators']) == prep.DEN
        assert min(rec['candidate_gap_numerators']) >= prep.GAP
        assert rec['candidate_gap_numerators'][0] > prep.GAP
        assert rec['candidate_gap_numerators'][5] > prep.GAP
        pf = prep.preflight(rec['numerators'], prep.DEN, forms, rays)
        assert pf['feasible'] and pf['int64_safe']
        repair_checks.append({'k': rec['k'], 'largest_gap_index': rec['largest_gap_index'], 'preflight': pf})
    old_repair = next(prep.gap_repairs(reg['seeds'][1]['numerators'], reg['seeds'][1]['denominator']))
    assert not old_repair['base_is_same_exact_point']
    assert any(Q(v) for v in old_repair['base_rounding_errors'])
    # Decode an EXISTING complete output only; do not ask the engine to replay it.
    old = json.loads((prep.ENGINE/'native_result.json').read_text())
    seed = reg['seeds'][1]
    check = search.validate_native(old, seed['numerators'], seed['denominator'], Q(4, 125))
    assert check['PASS']
    # Runtime proof of the lock. The child is the Python gate, never the native.
    command = [sys.executable, str(HERE/'search.py'), '--run']
    result = subprocess.run(command, cwd=HERE, capture_output=True, text=True, timeout=10)
    assert result.returncode == 3 and result.stdout.strip() == 'NOT_RELEASED: no release.json; zero native calls'
    assert not (HERE/'checkpoint.json').exists() and not (HERE/'native_outputs').exists()
    assert not (HERE/'__pycache__').exists()
    output = {
        'status': 'PREPARATION_CHECKS_PASS_NOT_RELEASED',
        'pid': os.getpid(), 'native_evaluations': 0,
        'elapsed_seconds': perf_counter()-tic,
        'zero_form_count': len(forms), 'known_degenerate_point_zero_forms': known_zero['zero_form_count'],
        'seed_checks': [s['preflight'] for s in reg['seeds']],
        'six_gap_inward_repairs_pure_checks': repair_checks,
        'old_c_repair_explicitly_new_quantized_input': True,
        'existing_C014_output_decoder_and_exact_Euler_PASS': True,
        'existing_C014_output_sha256': prep.sha(prep.ENGINE/'native_result.json'),
        'gate_command': command, 'gate_exit_code': result.returncode,
        'gate_stdout': result.stdout, 'gate_stderr': result.stderr,
        'no_checkpoint_or_native_output_created': True,
        'sources': [{'path': str(HERE/name), 'sha256': prep.sha(HERE/name)}
                    for name in ('prepare.py', 'search.py', 'registration.md', 'check_preparation.py', 'preparation.json')]
    }
    prep.dump(HERE/'checked_preparation.json', output)
    print(json.dumps(output, indent=2))

if __name__ == '__main__':
    main()
