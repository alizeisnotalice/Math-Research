#!/usr/bin/env python3
"""Audit 057's saved bytes and rational endpoints; never import its producer."""
from fractions import Fraction
from pathlib import Path
import hashlib
import json
import os
import signal
import sys
import time

HERE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def timeout(_signum, _frame):
    raise TimeoutError('Registered 120-second saved-evidence audit limit')


def main():
    start = time.monotonic()
    signal.signal(signal.SIGALRM, timeout)
    signal.alarm(120)
    audit_reg = json.loads((HERE / 'recovery_registration.json').read_text())
    assert audit_reg['status'] == 'registered_before_audit'
    local = {}
    for name, key in [('results.json', 'input_results_sha256'),
                      ('registration.json', 'input_registration_sha256'),
                      ('check.py', 'input_check_sha256')]:
        local[name] = sha(HERE / name)
        assert local[name] == audit_reg[key]
    result = json.loads((HERE / 'results.json').read_text())
    registration = json.loads((HERE / 'registration.json').read_text())
    assert result['registration_sha256'] == local['registration.json']
    assert result['source_sha256'] == local['check.py']
    assert result['external_sources'] == registration['sources']
    source_checks = []
    for row in registration['sources']:
        digest = sha(row['path'])
        assert digest == row['sha256']
        source_checks.append({'path': row['path'], 'sha256': digest, 'matches': True})
    assert result['status'] == 'kernel_geometry_and_constants_certified_no_tube_integral_claim'
    assert result['new_Green_evaluations'] == 0
    assert result['new_I_or_g_or_C5_bound'] is False
    checks = {}
    constants = result['constants']['constants']
    for name, row in constants.items():
        lower, upper = Fraction(row['lower_dyadic']), Fraction(row['upper_dyadic'])
        assert 0 < lower < upper
        checks[name] = {'saved_endpoints_positive_ordered': True}
    for name, cap in registration['preselected_strict_rational_caps'].items():
        assert Fraction(constants[name]['upper_dyadic']) < Fraction(cap)
        checks[name]['saved_upper_strictly_below_preregistered_cap'] = cap
    checkpoint = Path(registration['checkpoint'])
    if checkpoint.exists():
        assert checkpoint.read_bytes() == (HERE / 'results.json').read_bytes()
        cp_state = 'present_and_byte_identical_to_results'
    else:
        cp_state = 'not_present_saved_results_remain_authoritative'
    receipt = {
        'id': audit_reg['id'], 'status': 'passed_saved_evidence_audit',
        'audit_pid': os.getpid(), 'python': sys.executable,
        'audit_registration_sha256': sha(HERE / 'recovery_registration.json'),
        'audit_source_sha256': sha(__file__), 'local_sha256': local,
        'external_sources_verified': source_checks, 'saved_constant_checks': checks,
        'checkpoint': cp_state, 'original_pid_recorded': result['pid'],
        'original_substage_seconds_recorded': {
            'geometry': result['geometry']['seconds'],
            'constants': result['constants']['seconds']},
        'original_shell_exit_code': None,
        'original_exit_code_status': 'not preserved in supplied evidence; not inferred',
        'producer_reruns': 0, 'new_Arb_calls': 0, 'new_Green_evaluations': 0,
        'new_tube_scalar_evaluations': 0,
        'analytic_proof_location': 'report.md; not claimed to be machine checked',
        'seconds': time.monotonic() - start
    }
    signal.alarm(0)
    (HERE / 'recovery_audit.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'status': receipt['status'], 'audit_pid': receipt['audit_pid'],
                      'seconds': receipt['seconds'], 'producer_reruns': 0}))


if __name__ == '__main__':
    main()
