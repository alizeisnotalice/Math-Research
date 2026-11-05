"""Independent exact finite checks of the finite-depth bookkeeping.

The analytic reindexing proves arbitrary fixed N; these finite checks only
guard coefficients and indices and do not establish the geometric induction.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
from itertools import product

cD, b = F(1,2), F(3,7)
beta = b/cD
checks = []
for N in range(9,65):
    coefficients = {}
    # Entries are coefficients of f_k*v^(N-k); the initial f_0 term uses v^N.
    for h in range(1,N+1):
        coefficients[h] = coefficients.get(h,F(0))+cD*beta**(N-h)
        coefficients[h-1] = coefficients.get(h-1,F(0))-b*beta**(N-h)
    assert coefficients[N] == cD
    assert coefficients[0] == -b*beta**(N-1)
    assert all(coefficients[h] == 0 for h in range(1,N))
    checks.append(N)

out = {
    "rational_coefficients": {"cD": str(cD), "b": str(b), "beta": str(beta)},
    "telescoping_tested_depths": checks,
    "f0_A0_v_M_over_a_initial_power": "N+1",
    "sixth_little_o_weighted_power": "N",
    "examples": [],
    "source_geometry_certified": False,
}
reverse_cases = 0
for Av, a, Q, actual_ratio in product(
        (F(9,8),F(5,4),F(3,2)),
        (F(1,64),F(1,32),F(1,16)),
        (F(0),F(1,100),F(1,80)),
        (F(1),F(101,100),F(11,10))):
    V = F(7,3)
    M = V*(actual_ratio-1)
    total = V+M
    denominator = Av-a*(1+Q)
    assert denominator > 0
    benchmark = Av/denominator
    centered = total-Av*M/a+Q*total
    assert centered == denominator*V/a*(benchmark-actual_ratio)
    reverse_cases += 1
out["same_actual_M_reverse_identity_rational_cases"] = reverse_cases
for N in (9,10,11):
    # Existing seventh quantitative error O(v8); sixth leading remainder o(v6).
    # Existing local eighth/ninth-or-later errors O(v9) need not be improved.
    out["examples"].append({
        "N": N,
        "initial_O_power": N+1,
        "early_1_to_5_worst_O_power": N-5+6,
        "sixth_little_o_power": N-6+6,
        "seventh_quantitative_O_power": N-7+8,
        "eighth_local_O9_power": N-8+9,
        "penultimate_local_O9_power": 10,
        "terminal_local_O9_power": 9,
    })
out["all_assertions_passed"] = True
out["script_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
Path(__file__).with_suffix(".json").write_text(json.dumps(out,indent=2)+"\n")
print(json.dumps(out,indent=2))
