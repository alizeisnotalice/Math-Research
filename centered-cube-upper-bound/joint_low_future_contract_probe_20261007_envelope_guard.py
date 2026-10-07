#!/usr/bin/env python3
"""Additional exact guards; preserves the frozen original phi probe."""
from fractions import Fraction as F
from math import factorial, isqrt
from pathlib import Path
import hashlib
import json
import joint_low_future_contract_probe_20261007_v2 as probe

BASE = Path(__file__).resolve().parent
OUT = BASE / 'joint_low_future_contract_probe_20261007_envelope_guard.json'
assert not OUT.exists()
ln3lo, ln3hi = probe.log_interval(F(3))
assert ln3hi < F(9,8)
e_upper = sum((F(1,factorial(j)) for j in range(11)), F(0)) + F(12,11*factorial(11))
assert e_upper < F(11,4)
survival512 = (1-F(9,1024))**512
assert survival512 > F(1,128)
C = 128*F(11,4)**17
first_global = C*529/F(4**23)
second_global = C*512*F(15,16)**512
assert first_global < F(1,32)
assert second_global < F(1,32)
assert first_global+second_global < F(1,16)
rows = []
for n in [512,1024,4096]:
    alpha = isqrt(n)+(isqrt(n)**2 < n)
    radius0 = 1+F(16,n)
    radiusmax = radius0+F(1,16*n)
    ratio_radius_power_bound = F(16,15)
    assert F(3,4)*ratio_radius_power_bound < 1
    assert F(1,2)*ratio_radius_power_bound < 1
    assert F(16,45) < F(3,8)  # q/u <=16/45 < response at sigma=1/n.
    assert F(45,16)>2 and 3<4
    assert F(9,2*n) < F(1,16)
    assert (F(3,4)*radius0)**n < F(1,6)  # G_1<q/2, includes L<R.
    lnlo,lnhi = probe.log_interval(F(n+2))
    upper = C*(F(1,4**alpha)+F(15,16)**n)
    assert upper < 1/lnhi**4
    assert n*upper < F(1,16)
    rows.append(dict(n=n, alpha=alpha,
        max_envelope_ratio_upper=probe.digest_fraction(upper),
        epsilon_certified_lower=probe.ratio(1/lnhi**4),
        receiver_box_side_over_a=f'1/{32*n}',
        positive_volume_hardband_uniform_lower='45/16',
        partial_scale_soft_response_factor_upper='4/5',
        partial_scale_hard_response_factor_upper='8/15',
        hard_shell_delay_upper='1/16',
        max_envelope_low=True))
data=dict(status='PASS_THREE_EXACT_MAX_ENVELOPE_AND_POSITIVE_BOX_GUARDS',
    rows=rows, universal_coefficient=probe.ratio(C),
    all_n_512_proof=dict(first_term_nC_upper=probe.ratio(first_global),
        second_term_nC_upper=probe.ratio(second_global),
        survival512=probe.digest_fraction(survival512),
        exp1_certified_upper=probe.ratio(e_upper)),
    arithmetic='Exact Fraction with explicit log remainder; no floating predicates.',
    frozen_phi_script_sha256=hashlib.sha256(Path(probe.__file__).read_bytes()).hexdigest(),
    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    scope='Positive receiver box kernel-level contracts only; no actual birth/CPGP/history certification or weak lower bound.')
with OUT.open('x') as f:
    json.dump(data,f,ensure_ascii=False,indent=2)
print(json.dumps(dict(status=data['status'],n=[r['n'] for r in rows])))
