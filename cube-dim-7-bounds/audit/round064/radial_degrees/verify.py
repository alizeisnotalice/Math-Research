#!/usr/bin/env python3
"""Seven new bounded rational checks of radial degree bookkeeping only."""

from fractions import Fraction as F
from math import prod
from pathlib import Path
import hashlib
import json
import platform

source = Path('/Users/zhengzhihao/Desktop/T/cube-r040-cover-selector/cube-dim-7-bounds/upper/round064/radial_calibration/PROOF.md')
source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
assert source_hash == '17ce38ab5af46c0358046bdcc2aa57a05844c0799e9a97280d959887f59615d6'

s, q, t = F(3, 2), F(4, 5), F(20)
h, zeta = F(1, 2 ** 72), F(1, 2 ** 100)
epsilon, rho = zeta * t, t * q
kappa = 1 + 1 / t - zeta
R, a0 = s * q - h, (1 - epsilon) * q
g = R - a0
e, delta = F(1, 1024), F(1, 2 ** 30)
rlo, rhi = 1 + delta, F(2727, 2000) - delta
qmin, qmax = F(2, 3), F(909, 1000)
assert 0 < g <= R / 3 and R > 1
assert epsilon * qmax <= h / s
assert a0 == (kappa - 1) * rho


def direct_fraction(radius, coord):
    width = max(F(0), min(radius, coord + kappa * rho)
                - max(-radius, coord - kappa * rho))
    return width / (2 * radius)


records = []
for dim, offset, theta in ((1, F(1, 4), F(1, 4)),
                            (2, F(1, 2), F(7, 8)),
                            (3, F(1), F(1, 2)),
                            (4, F(7, 4), F(7, 8))):
    r = R + offset * h
    assert rlo < r < rhi and abs(r - s * q) < h and rho >= 4 * r
    boundary = kappa * rho - r
    if dim == 1:
        coords = [rho]
    elif dim == 2:
        coords = [rho, F(0)]
    elif dim == 3:
        coords = [rho, rho - F(1, 4), -boundary]
    else:
        coords = [rho, -rho, F(0), boundary]
    assert max(map(abs, coords)) == rho
    ell, ell_ref, pins = [], [], []
    for z in coords:
        ai = kappa * rho - abs(z)
        assert ai > 0
        frac = direct_fraction(r, z)
        ref = direct_fraction(R, z)
        assert frac == min(F(1), (r + ai) / (2 * r))
        assert ref == min(F(1), (R + ai) / (2 * R))
        Ti = sum(abs(z - sign * r) < kappa * rho for sign in (-1, 1))
        assert Ti / (2 * frac) <= 1
        if frac < 1:
            assert Ti == 1
        ell.append(frac)
        ell_ref.append(ref)
        pins.append(Ti)
    if dim >= 3:
        assert ell[-1] == 1 and pins[-1] == 1  # strict-capture tangency
    G, Gref = prod(ell), prod(ell_ref)
    face = sum(F(pins[i], 2) * prod(ell[j] for j in range(dim) if i != j)
               for i in range(dim))
    assert face <= dim * G <= dim * Gref
    bsoft = min(F(1), theta / Gref)
    matched_width = min(rhi, s * q + h) - max(rlo, s * q - h)
    assert matched_width == 2 * h
    column_upper = bsoft * dim * Gref * matched_width / rhi
    k = dim * h * theta / rhi
    assert column_upper <= 2 * k
    records.append({'check': f'complete-column-{dim}D', 'n': dim,
                    'pinned_target_sign_counts': pins,
                    'face_over_nG': str(face / (dim * G)),
                    'soft_column_upper_over_2k': str(column_upper / (2 * k))})

# 5. No floating logarithm: use exact atanh-series lower/upper envelopes.
r = F(6, 5)
qminus, qplus = (r - h) / s, (r + h) / s
assert qmin <= qminus < qplus <= qmax
assert s * qplus - h == r and s * qminus - h == r - 2 * h
assert qplus - qminus == 2 * h / s
u = h / r
width_lower = (r / rhi) * 2 * u
width_upper = width_lower / (1 - u ** 2)
assert width_lower == 2 * h / rhi
assert width_upper <= (1 + e) * 2 * h / rhi
records.append({'check': 'radial-log-width', 'lower_over_2h_rhi': '1',
                'upper_over_2h_rhi': str(1 / (1 - u ** 2))})

# 6. Constant-size count computation, not a 4096-coordinate integration.
n, free_low_count = 4096, 2730
Nsum = F(n + free_low_count)
endpoint = []
for qq in (qminus, qplus):
    RR, aa = s * qq - h, (1 - epsilon) * qq
    assert 0 < aa < r
    # Low coordinates are zero; other free coordinates equal 5r/6.
    assert F(5, 6) * r > aa
    mean = n + (n - 1) * aa / r
    assert abs(Nsum - mean) <= e * n / 4
    pin = (1 + aa / RR) / 2
    ratio = Nsum / (2 * n * pin)
    endpoint.append((pin, ratio))
assert 2 * h + F(1, n) <= 3 * e / 4
assert endpoint[0][0] > endpoint[1][0]
assert 1 - e <= endpoint[0][1] <= endpoint[1][1] <= 1 + e
records.append({'check': 'endpoint-pin-count-calibration', 'n': n,
                'Nsum': str(Nsum), 'qminus_ratio': str(endpoint[0][1]),
                'qplus_ratio': str(endpoint[1][1]),
                'scope': 'count condition only; no J integration/product event claimed'})

# 7. Radial k cancels exactly in the correctly scaled deficit, and cap/k
# keeps rhi/h. This is algebra, not a full-history/probability realization.
P0 = 1 - e - F(1, 2 ** 50)
Ulo = P0 * (1 - e) / (1 + e)
Uhi = (1 + e) / (P0 * (1 - e))
M = 2 * (1 + e) * Uhi
rawlo = 2 * (1 - e) ** 2 * Ulo
rowlo, bbudget = rawlo - e, rawlo - 2 * e
assert M > 2 and bbudget > 0
mass_A, sample_n, theta, c = F(1239, 2000), F(2), F(1, 4), F(3, 2)
k = sample_n * h * theta / rhi
tau, energy, degree = e * k, mass_A * rowlo * k, M * k
scaled = (c ** 2 * energy - tau * c * mass_A) / (c * degree)
base = mass_A * bbudget / M
assert scaled == mass_A * (c * rowlo - e) / M > base
assert (energy - tau * mass_A) / degree == base
Csrc, formal_tail = F(727, 2000), F(1, 11)
cap_mass = Csrc * sample_n * formal_tail
assert cap_mass / k == Csrc * rhi / h * formal_tail / theta
records.append({'check': 'finite-a-and-cap-normalization', 'M_gt_2': M > 2,
                'scaled_deficit_over_base': str(scaled / base),
                'cap_prefactor': 'Csrc*rhi/h',
                'scope': 'formal cap tail variable, no exponent or scalar gate evaluation'})

assert len(records) == 7
root = Path(__file__).resolve().parent
result = {'status': 'all seven new bounded exact checks passed',
          'python': platform.python_version(), 'checks': records,
          'source': {'path': str(source), 'sha256': source_hash},
          'coordinate_case_scope': 'local capture guards only; n,t,L not asserted to form the author fixed operator tuple',
          'excluded': 'relative-J proof, soft probability/retention, full tuple/scalar and C_n theorem',
          'artifact_sha256': {name: hashlib.sha256((root / name).read_bytes()).hexdigest()
                              for name in ('PROOF.md', 'verify.py')}}
(root / 'verification.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({'status': result['status'], 'checks': len(records),
                  'column_dimensions': [1, 2, 3, 4],
                  'count_calibration': 'passed by coordinate multiplicities',
                  'finite_a_normalization': 'passed'}, indent=2))
