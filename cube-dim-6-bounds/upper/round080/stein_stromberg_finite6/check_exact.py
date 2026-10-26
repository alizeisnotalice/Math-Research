"""Five fixed exact checks for the allocated covering formula; no grid/kernel."""
import hashlib
import json
import signal
import sys
import time
from fractions import Fraction as Q
from pathlib import Path

signal.alarm(170)
start = time.monotonic()
here = Path(__file__).resolve().parent
output = here/'results.json'
assert not output.exists(), 'Refuse to overwrite evidence'
registration = json.loads((here/'registration.json').read_text())
source = Path(registration['source_pdf'])
assert hashlib.sha256(source.read_bytes()).hexdigest() == registration['source_pdf_sha256']
# Exact finite prefixes; omitted terms are strictly positive analytically.
exp_prefix = sum((Q(1), Q(1), Q(1, 2), Q(1, 6)), Q(0))
assert exp_prefix == Q(8, 3)
log_two_first = 2*Q(1, 3)
assert log_two_first == Q(2, 3)
log_seven_prefix = 2*sum((Q(3, 4)**(2*k+1)/Q(2*k+1) for k in range(8)), Q(0))
assert log_seven_prefix > Q(97, 50)
middle = (1+Q(9, 7)**6)*Q(35, 3)
assert middle > Q(1602, 25) > 64
large = 2*(1+6*Q(8, 3)*Q(97, 50))
assert large == Q(1602, 25) and Q(65) > large
result = dict(status='PASS_FIVE_FIXED_EXACT_FACTS', configurations=5,
              exp_prefix=str(exp_prefix), log_two_prefix=str(log_two_first),
              log_seven_prefix=str(log_seven_prefix), middle_product_lower=str(middle),
              uniform_allocated_cost_strict_lower=str(large),
              new_C6_bound=False, true_covering_optimum_lower_claim=False,
              kernel_calls=0, parameter_grid_points=0, old_baseline_evaluations=0,
              python=sys.version, elapsed_seconds=time.monotonic()-start,
              source_sha256=registration['source_pdf_sha256'],
              code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              theorem_sha256=hashlib.sha256((here/'theorem.md').read_bytes()).hexdigest())
with output.open('x') as stream:
    stream.write(json.dumps(result, indent=2)+'\n')
print('PASS five exact facts: allocated formula >1602/25; no new C6 bound')
