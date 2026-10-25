#!/usr/bin/env python3
"""One finite exact/Arb batch for complete spatial elimination of I.

No spatial values, quadrature, sampling, or peak estimates are computed.
"""
from fractions import Fraction as F
from itertools import product, permutations
from math import factorial
from pathlib import Path
import hashlib
import json
import time

import flint
from flint import arb, ctx

ctx.prec = 256
began = time.monotonic()
root = Path(__file__).parent
V = [(1,)+s for s in product((-1, 1), repeat=4)]
frame = [[sum(v[i]*v[j] for v in V) for j in range(5)] for i in range(5)]
assert frame == [[16*int(i == j) for j in range(5)] for i in range(5)]
for p in permutations(range(5)):
    for s in product((-1, 1), repeat=5):
        image = [tuple(s[j]*v[p[j]] for j in range(5)) for v in V]
        normalized = {w if w[0] == 1 else tuple(-x for x in w) for w in image}
        assert normalized == set(V)

beta, b = F(3, 4), F(3, 14)
a, gamma = 5*b/16, 5*b-beta
assert a == F(15, 224) and gamma == F(9, 28)
assert 0 < gamma < beta and 5*b*(1/beta-1) == F(5, 14)
rising = F(1)
for j in range(5):
    rising *= beta+j
Ed = beta**5*factorial(5)/rising


def ab(x):
    x = F(x)
    return arb(x.numerator)/x.denominator


def show(x):
    return x.str(65)


pi = arb.pi()
c = 3/(arb(2)**ab(F(15, 2))*pi**ab(F(5, 2)))
negative_s_moment = (1+ab(a/beta)).gamma()/(1+ab(a)).gamma()
positive_s_moment = (1-ab(gamma/beta)).gamma()/(1-ab(gamma)).gamma()
normalization = 1/(16*(1+ab(beta)).gamma())
Q = (arb(16)**ab(F(5, 14))*(1+ab(beta)).gamma()
     *(1-ab(gamma/beta)).gamma()/(1-ab(gamma)).gamma()
     *negative_s_moment**16)
Q_unsimplified = (arb(16)**(-5*ab(b))/normalization
                  *arb(16)**ab(gamma/beta)*positive_s_moment
                  *negative_s_moment**16)
assert (Q-Q_unsimplified).contains(0)
Q_lower = ab(Ed)**(-ab(b))
spatial_factor = pi**2/60
I_upper = c**ab(F(10, 7))*spatial_factor*Q
I_old_plane_upper = 5*I_upper
target = ab(F(11, 1000000))
Q_target_gate = target/(c**ab(F(10, 7))*spatial_factor)
checks = {"negative_moment_upper < 5/2": Q < ab(F(5, 2)),
          "I_upper < 1/50000": I_upper < ab(F(1, 50000)),
          "I_upper < 11/1000000": I_upper < target}
inputs = [root/'round011_shape_integral_registration.json',
          Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-upper/cube-dim-5-bounds/astra_upper001/report.md'),
          Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-upper/cube-dim-5-bounds/astra_upper001/verify.py'),
          Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube5-astra-upper/cube-dim-5-bounds/astra_upper001/verification.json')]
out = {
    "status": "exact_spatial_elimination_and_complete_I_upper_bound_target_comparison",
    "C5_improved": False,
    "model": "16-body alpha=3/2 independent stable sum, Lebesgue measure, Q=[-1,1]^5",
    "exact_geometry": {"frame": frame, "signed_permutation_maps_checked": 120*32,
                       "face_orthant_pieces": 10*16, "sphere_area": "8*pi^2/3",
                       "exact_spatial_factor": "pi^2/60"},
    "rational_parameters": {"beta": str(beta), "b": str(b), "a": str(a),
                            "gamma": str(gamma), "first_determinant_moment": str(Ed)},
    "arb_256bit": {"shape_prefactor_c": show(c),
                   "normalization_ET_minus_beta": show(normalization),
                   "one_variable_negative_moment": show(negative_s_moment),
                   "one_variable_positive_moment": show(positive_s_moment),
                   "negative_determinant_moment_upper_Q": show(Q),
                   "negative_determinant_moment_Jensen_lower": show(Q_lower),
                   "target_sufficient_negative_moment_upper_gate": show(Q_target_gate),
                   "complete_I_upper": show(I_upper),
                   "old_full_plane_bound_with_same_moment_upper": show(I_old_plane_upper),
                   "I_upper_over_target": show(I_upper/target)},
    "fixed_comparisons": checks,
    "proof_boundary": "Q is an upper bound for E_R det(M)^(-3/14); Q exceeding the sufficient target gate does not rule out a sharper moment bound, and Jensen may also lose information. No global g or C5 improvement is proved.",
    "cost": {"flint_version": flint.__version__, "precision_bits": ctx.prec,
             "elapsed_seconds": time.monotonic()-began,
             "spatial_points": 0, "quadrature_calls": 0, "sampling_calls": 0},
    "source_sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
    "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
}
raw = (json.dumps(out, indent=2)+"\n").encode()
destination = root/'round011_shape_integral_results.json'
destination.write_bytes(raw)
print(json.dumps({"sha256": hashlib.sha256(raw).hexdigest(),
                  "rational_parameters": out["rational_parameters"],
                  "arb_256bit": out["arb_256bit"],
                  "fixed_comparisons": checks, "cost": out["cost"]}, indent=2))
