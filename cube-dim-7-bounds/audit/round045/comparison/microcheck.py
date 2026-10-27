#!/usr/bin/env python3
"""Independent compressed dyadic envelopes; does not run author checks."""

from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import platform


def two(exponent):
    return F(2 ** exponent) if exponent >= 0 else F(1, 2 ** (-exponent))


sources = {
    'comparison': '/Users/zhengzhihao/Desktop/T/cube-r030-parameters/cube-dim-7-bounds/upper/round045/comparison/PROOF.md',
    'geometry': '/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round045/soft_inverse_column/PROOF.md',
    'soft': '/Users/zhengzhihao/Desktop/T/cube-r030-soft/cube-dim-7-bounds/upper/round045/layered_source/PROOF.md',
    'cap': '/Users/zhengzhihao/Desktop/T/cube-r030-parameters/cube-dim-7-bounds/revision_v3/asymptotic/proof.md',
    'R042': '/Users/zhengzhihao/Desktop/T/cube-r030-parameters/cube-dim-7-bounds/upper/round042/explicit_threshold/THRESHOLD_PROOF.md',
    'old_comparison': '/Users/zhengzhihao/Desktop/T/cube-r030-parameters/cube-dim-7-bounds/upper/round042/comparison/PROOF.md',
}
hashes = {key: hashlib.sha256(Path(path).read_bytes()).hexdigest()
          for key, path in sources.items()}
assert hashes['comparison'] == 'ba7498d1ccdc2c9158b065ade728c60f7b83cc2706ce517580d866903451d2f9'
assert hashes['geometry'] == 'add435fdf2f55993eba9064bd50c493e4c42769c1c48463329c4015be2e66ab9'
assert hashes['soft'] == '22db93b275416b90e9f15a0118a88cd412fee0c2d0b20c4cce1ed32769f259ac'

# Use a single dominant dyadic term per E, not the author's exact sums.
L0 = two(240)
cases = {
    'I_display': (two(-55), two(-46), two(-246)),
    'S_display': (two(-106), two(-96), two(-144)),
    'I_endpoint': (two(-43), two(-40), two(-270)),
}
records = {}
for name, (sigma, delta, dominant) in cases.items():
    terms = (2 / L0 ** 2, 2 / (sigma ** 2 * L0 ** 2),
             16 / (sigma ** 2 * two(360)))
    assert all(term <= dominant for term in terms)
    margin = delta / 2 - 2 * sigma
    assert margin > 0
    assert sum(terms) < 3 * dominant < margin
    records[name] = {'dominant_E_term': str(dominant),
                     'coarse_E_upper': str(3 * dominant),
                     'remaining_margin': str(margin)}

# Strict complete-cost comparisons via coarser common gap certificates.
assert two(-238) + two(-111) < two(-110) < two(-97)
assert two(-238) + two(-115) < two(-114) < two(-97)
assert two(-238) + two(-96) < two(-95) < two(-47)
assert two(-238) + two(-46) < two(-45) < two(-41)

# Endpoint's actual budget vs positive scalar buffer, exactly at equality.
dmax, eps = two(-24), two(-20)
delta, sigma = dmax / 65536, 32 * dmax ** 2
assert delta == two(-40) and sigma == two(-43)
assert 2 * sigma == delta / 4
assert delta / 2 - 2 * sigma == delta / 4
assert 12 * dmax < F(1, 64)
assert dmax / 32 - 12 * dmax ** 2 > dmax / 64 > 0
assert 1280 * eps == F(5, 4096) < F(1, 8)

# A different, strictly larger probe L only checks simple algebraic guards;
# it is not asserted to be a full-operator activation threshold.
Lprobe = two(256)
guards = []
for d in (two(-32), two(-30), dmax):
    sig = 32 * d ** 2
    h = eps / (Lprobe * d)
    # n>=L^2/2 gives zeta*t/h <= this rational expression.
    zeta_t_over_h_upper = 4 * d / (eps * (1 - sig) * Lprobe ** 4)
    qmin = 1 - 3 * d / 4
    Rmin = qmin + d - h
    assert 0 < sig <= two(-43) < F(1, 2)
    assert h < d / 16 and zeta_t_over_h_upper < F(1, 2)
    assert Lprobe * d >= 1 and qmin > F(1, 2) and Rmin > 1
    assert d < Rmin / 2 and 1 + 3 * d / 4 < 3
    guards.append({'d': str(d), 'probe_L': str(Lprobe),
                   'h_over_d': str(h / d),
                   'zeta_t_over_h_upper': str(zeta_t_over_h_upper)})

out = {
    'status': 'all independent compressed exact checks passed',
    'python': platform.python_version(),
    'method': 'Fraction dyadic dominant-term envelopes and separate endpoint guards; no author verify/receipt or production program executed',
    'scalar_base_L': str(L0), 'coarse_certificates': records,
    'simple_guard_probe': guards, 'source_sha256': hashes,
    'limitations': 'No numerical exponentials, exact cap integration, new geometry certification, or full-operator activation threshold is claimed.',
}
path = Path(__file__).with_name('microcheck.json')
path.write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps({'status': out['status'], 'coarse_scalar_cases': len(records),
                  'endpoint_budget': str(delta), 'endpoint_display': str(delta / 2),
                  'simple_guard_probes': len(guards)}, indent=2))
