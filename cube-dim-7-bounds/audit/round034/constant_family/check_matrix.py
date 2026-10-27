#!/usr/bin/env python3
"""Independently check only the 12 frozen C34A rational candidates.

No producer code is imported or run. Positive exponential partial sums give
an exact upper bound for exp(-r), which gives a lower bound for KP.
"""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json
import os
import platform
import resource
import sys
import time

started = time.perf_counter()
THEORY = Path('/Users/zhengzhihao/Desktop/T/cube-r030-parameters/cube-dim-7-bounds/upper/round034/constant_theory')
MATRIX = Path('/Users/zhengzhihao/Desktop/T/cube-r030-compute/cube-dim-7-bounds/upper/round034/compute/c34a_results.json')
HERE = Path(__file__).resolve().parent

def evaluate(Y, a, b, c, q, degree):
    Y, a, b, c, q = map(Q, (Y, a, b, c, q))
    g = b+c
    assert Y >= 8 and 0 < a < b < Y/2 and 0 < c < Y/2
    assert g*Y >= 1
    f1 = 1-4*g/Y
    f2 = 1-2*g*g/(1-b/Y)
    assert f1 > 0 and f2 > 0
    KA = c/(4*(1+c/Y))
    r = 9*a*Y/32
    term = total = Q(1)
    for k in range(1, degree+1):
        term *= r/k
        total += term
    exp_upper = 1/total
    KP_lower = (b-a)*(f1*f2-exp_upper/2)-1/Y-Q(3, 8)/Y**3
    assert 0 < q < KA and q < KP_lower
    return KA, KP_lower

data = json.loads(MATRIX.read_text())
assert len(data['results']) == 12
assert len({row['Y'] for row in data['results']}) == 12
assert len({row['pid'] for row in data['results']}) == 12
rows = []
for row in data['results']:
    cert = row['certificate']
    KA, KP = evaluate(row['Y'], cert['a'], cert['b'], cert['c'], cert['q'], 128)
    q = Q(cert['q'])
    assert KA == Q(cert['KA_lower'])
    assert KP == Q(cert['KP_lower'])
    assert KA-q == Q(cert['strict_KA_margin'])
    assert KP-q == Q(cert['strict_KP_margin'])
    assert row['status'] == 'certified'
    assert row['objective_calls'] <= row['max_calls']
    qstar = (q+min(KA,KP))/2
    assert q < qstar < KA and qstar < KP
    rows.append({
        'Y': row['Y'], 'L_min': max(64,row['Y']**2),
        'a': cert['a'], 'b': cert['b'], 'c': cert['c'], 'q': cert['q'],
        'domain_pass': True, 'positive_factors_pass': True,
        'q_lt_KA_pass': True, 'q_lt_independent_KP_lower_pass': True,
        'producer_lower_bounds_and_margins_match': True,
        'strict_qstar_exists_by_rational_midpoint': True,
        'producer_pid': row['pid'],
        'producer_objective_calls': row['objective_calls'],
        'producer_wall_seconds': row['wall_seconds'],
        'producer_cpu_seconds': row['cpu_seconds'],
        'producer_peak_rss_bytes': row['peak_rss_bytes'],
    })

# Independently check the hand-picked proof example with only six terms.
KA0, KP0 = evaluate(64, Q(1,10), Q(29,100), Q(21,100), Q(1,20), 6)
assert KA0 == Q(336,6421)
# The source uses an even coarser exp(-9/5)<1/6 bound.
coarse_KP0 = Q(19,100)*(Q(98301,203872)-Q(1,12))-Q(1,64)-Q(3,8*64**3)
assert coarse_KP0 == Q(12055437337,200414330880) > Q(1,20)
assert coarse_KP0-Q(1,20) == Q(2034720793,200414330880)

sources = {str(p): hashlib.sha256(p.read_bytes()).hexdigest()
           for p in [THEORY/'PROOF.md', THEORY/'verify.py',
                     THEORY/'verification.json', MATRIX]}
usage = resource.getrusage(resource.RUSAGE_SELF)
result = {
    'status': 'PASS', 'batch': data['batch'], 'rows_verified': len(rows),
    'arithmetic': 'stdlib fractions.Fraction; independent positive partial sums k=0..128',
    'no_producer_code_imported': True, 'no_search': True,
    'python': sys.version, 'platform': platform.platform(),
    'auditor_pid': os.getpid(), 'auditor_wall_seconds': time.perf_counter()-started,
    'auditor_cpu_seconds': usage.ru_utime+usage.ru_stime,
    'auditor_peak_rss_bytes': usage.ru_maxrss if sys.platform == 'darwin' else usage.ru_maxrss*1024,
    'source_sha256': sources,
    'producer_resources_status': 'reported by frozen input matrix; not independently observed timings',
    'producer_workers_max_reported': data['workers_max'],
    'producer_batch_wall_seconds_reported': data['wall_seconds'],
    'rows': rows,
    'scope': 'fixed candidate certificates only; no global-optimality certification',
}
(HERE/'matrix_receipt.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k != 'rows'}, indent=2))
for row in rows:
    print('PASS', row['Y'], row['q'], 'producer PID', row['producer_pid'])
