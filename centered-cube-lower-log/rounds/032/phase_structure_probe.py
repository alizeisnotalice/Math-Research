"""Deterministic continuous-scale n=2 phase diagnostic, not matched-source MC.

rho(x)=1+eta*cos(2*pi*x), eta=.9.  On I=[10/9,11/9], its
window average is 1+b(L)cos(2*pi*x), with b monotone.  A quadratic
maximization in b gives the exact continuous-scale maximum at each
query grid point.  The query grid is quadrature, not probability sampling.
"""
from pathlib import Path
from fractions import Fraction
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
ETA = .9
GRID = 256
I = (10 / 9, 11 / 9)
b0, b1 = (ETA * math.sin(math.pi * L) / (math.pi * L) for L in I)
low, high = sorted((b0, b1))


def maximum_quadratic(s, product):
    left = 1 + low * s + low * low * product
    right = 1 + high * s + high * high * product
    result = np.maximum(left, right)
    concave = product < 0
    critical = np.divide(-s, 2 * product, out=np.zeros_like(s), where=concave)
    valid = concave & (critical >= low) & (critical <= high)
    return np.where(valid, np.maximum(result, 1 + critical * s
                                     + critical * critical * product), result)


def empirical_weak(values):
    ordered = np.sort(values.ravel())[::-1]
    return float(np.max(ordered * np.arange(1, len(ordered) + 1) / len(ordered)))


x = (np.arange(GRID) + .5) / GRID
x1, x2 = np.meshgrid(x, x, indexing="ij")
c1, c2 = np.cos(2 * np.pi * x1), np.cos(2 * np.pi * x2)
base = maximum_quadratic(c1 + c2, c1 * c2)
rows = []
for R in (1, 2, 3, 4, 8, 16, 64):
    envelope = np.zeros_like(base)
    for r in range(R):
        ca = np.cos(2 * np.pi * (x1 - r / R))
        cb = np.cos(2 * np.pi * (x2 - r / R))
        envelope += maximum_quadratic(ca + cb, ca * cb) / R
    if R == 1:
        shared = base.copy()
    elif R == 2:
        coefficient = c1 * c2
        shared = 1 + np.where(coefficient >= 0, low * low, high * high) * coefficient
    else:
        coefficient = .5 * np.cos(2 * np.pi * (x1 - x2))
        shared = 1 + np.where(coefficient >= 0, low * low, high * high) * coefficient
    assert np.min(envelope - shared) > -1e-12
    rows.append(dict(R=R, shared_weak_grid=empirical_weak(shared),
                     shared_mean_grid=float(shared.mean()),
                     envelope_weak_grid=empirical_weak(envelope),
                     envelope_mean_grid=float(envelope.mean()),
                     mean_common_scale_gap=float((envelope-shared).mean()),
                     max_common_scale_gap=float((envelope-shared).max())))

# For R>=3 the entire mixed source stabilizes by Fourier filtering.
minor = high * high / 2
major = low * low / 2
exact_shared_weak = 1 - minor  # minor<1/(1+pi), as verified below.
assert minor < 1 / (1 + math.pi)
exact_shared_mean = 1 + (major-minor) / math.pi

B = 512
N = B-1
T = N*(N+1)//2
D = N*(N-1)
H = N + T*(T-1) + 2*D*D + 4*D
output = dict(numpy_version=np.__version__, source="positive pure first harmonic",
              n=2, eta=ETA, query_grid=GRID, scales="exact continuous I",
              interval=I, b_interval=(low, high), rows=rows,
              R_ge_3_exact_weak=exact_shared_weak,
              R_ge_3_exact_mean=exact_shared_mean,
              label_B512=dict(nonzero_labels=N, distinct_products_required=T,
                              distinct_ratios_required=D, H=H,
                              smallest_d_with_2d_gt_H=H.bit_length()))
(HERE / "phase_structure_results.json").write_text(json.dumps(output, indent=2)+"\n")
print(json.dumps(output, indent=2))
