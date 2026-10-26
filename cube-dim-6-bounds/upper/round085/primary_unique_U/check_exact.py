"""One fixed composite exact proof-fact check; no pair/triple enumeration."""
import hashlib
import itertools
import json
import signal
import time
from fractions import Fraction as Q
from pathlib import Path

signal.alarm(170)
start = time.monotonic()
p = Path(__file__).resolve().parent
assert not (p/'results.json').exists()
sign_values = []
for a in (Q(1, 2), Q(-1, 2)):
    values = sorted({sum(s*c for s, c in zip(signs, (a, -a, Q(1))))
                     for signs in itertools.product((-1, 1), repeat=3)})
    assert values == list(map(Q, (-2, -1, 0, 1, 2)))
    assert Q(1, 2) not in values
    sign_values.append(list(map(str, values)))
assert Q(1, 864) > Q(1, 32)**2
assert Q(1, 144) == Q(1, 12)**2
assert 32**2*6 < 80**2
assert Q(1, 32)/16 == Q(1, 512)
assert Q(1, 512)/81 == Q(1, 41472)
r = dict(status='PASS_FIXED_SIGN_AND_RATIONAL_MARGIN_FACTS',
         sign_values=sign_values, rho='1/512', cutoff='1/41472',
         pair_distance_squared_lower='1/144',
         nonincident_triple_distance_squared_lower='1/864',
         configurations=1, body_pair_or_triple_enumerations=0,
         kernel_calls=0, old_baseline_evaluations=0,
         proof_sha256=hashlib.sha256((p/'proof.md').read_bytes()).hexdigest(),
         elapsed_seconds=time.monotonic()-start)
with (p/'results.json').open('x') as stream:
    stream.write(json.dumps(r, indent=2)+'\n')
print('PASS one fixed composite sign/margin verifier; unique-U proof remains analytic')
