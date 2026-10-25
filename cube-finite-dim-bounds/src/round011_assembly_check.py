"""Check exact conditional costs and the archived g3 certificate aggregation.

Does not certify I, rerun nodes, or update the bounds ledger.  The analytic
dependencies and remaining gates are in proofs/round011_candidate_assembly.md.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
from flint import arb, ctx

ROOT = Path(__file__).resolve().parents[1]
ctx.prec = 256
manifest = json.loads((ROOT / 'audit/033_manifest.json').read_text())
for name, expected in manifest['files'].items():
    assert hashlib.sha256((ROOT / 'audit' / name).read_bytes()).hexdigest() == expected, name
for name, expected in manifest['sources'].items():
    assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == expected, name
result = json.loads((ROOT / 'audit/033_results.json').read_text())
rows = result['patches']
assert len(rows) == 136
assert {tuple(r['patch']) for r in rows} == {(i,j) for i in range(16) for j in range(i+1)}
assert sum(r['fresh_nodes'] for r in rows) == 59976
assert min(F(r['H_lower_exact']) for r in rows) > F(131,25000)
cases = []
for n,g,I,U in [(3,F(131,25000),F(29,50000),F(8)),
                 (4,F(21,25000),F(19,100000),F(16))]:
    power = F(n)**n * I**(n-1) / g**n
    gap = U**(n-1) - power
    assert gap > 0
    cost = (arb(power.numerator)/power.denominator).root(n-1)
    cases.append(dict(n=n, g_target=str(g), I_target=str(I),
                      conditional_cost_power=str(power), exact_power_gap=str(gap),
                      conditional_cost_interval=str(cost)))
integral_path = ROOT / 'upper/continuous_integral_n4/ci4b_check_results.json'
assert hashlib.sha256(integral_path.read_bytes()).hexdigest() == 'ea99abf732e0ad6cb2981321f01e922f1a2251e5bdeb50271bad21a788ea0351'
actual_I = F(json.loads(integral_path.read_text())['I4_upper_rational'])
outward_I = F(92072351,500000000000)
assert actual_I <= outward_I
actual_power = 256*outward_I**3/F(21,25000)**4
simple_U = F(14753,1000)
assert actual_power < simple_U**3 < 16**3
print(json.dumps({'status':'PASS conditional arithmetic and g3 aggregation only',
                  'new_upper_bound_certified':False, 'cases':cases,
                  'n4_actual_integral_conditional_comparison': {
                      'outward_I':str(outward_I), 'candidate_U':str(simple_U),
                      'cost_cube':str(actual_power),
                      'strict_gap':str(simple_U**3-actual_power),
                      'remaining':'independent full representation/g/I audit'}}, indent=2))
