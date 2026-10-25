"""Independent scalar/direction audit; it performs no kernel evaluations."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json

from flint import arb, ctx

ctx.prec = 384
MAIN = Path(__file__).resolve().parents[3]


def a(x):
    return arb(x.numerator) / x.denominator


mesh_path = MAIN / 'upper/round005/cover/simplex_cover_audit.json'
spectral_path = MAIN / 'upper/round004/primary/spectral_moment.json'
point_path = MAIN / 'upper/round003/direct_kernel/local_minimum/exact_local_reference_results.json'
mesh = json.loads(mesh_path.read_text())
spectral = json.loads(spectral_path.read_text())
point = json.loads(point_path.read_text())

# All 120 leaves have D_S >= Dmin; their volumes sum to 1/120.
Dmin = Q(mesh['leaf_exact_variance_min'])
assert Dmin == Q(25, 192)
assert Q(mesh['sum_leaf_volumes']) == Q(1, 120)
assert mesh['coordinate_permutation_copies'] == 120

# The previously certified uniform spectral lower bound is strictly above L.
L = a(Q(6, 10**8))
assert arb(spectral['uniform_g_lower']) > L

# The point is certified off S2 in its report; its continuity neighborhood makes
# this also an essential-infimum upper bound.  Use only the deliberately loose U.
U = a(Q(23, 10**6))
G_point = arb(point['G_arb_384'])
assert G_point < U

q = a(Q(6, 5))
split = a(Q(1, 100))
threshold = a(Q(16, 3)) * U
K = a(Q(6, 7)) * a(Q(1, 243)) ** a(Q(1, 6)) * arb.pi() ** (-a(Q(11, 6)))

# Formula (3) on every leaf with every vertex upper bound U_i >= H_eps(v_i)
# >= L has total right-hand L^q exp(q Dmin/(32 eps^2)); after 1/q it is this.
coarse_at_split = L * (a(Dmin) / (32 * split**2)).exp()
tail_at_split = K * split ** a(Q(7, 6))
assert coarse_at_split > threshold
assert tail_at_split > threshold

out = {
    'status': 'PASS_fixed_method_only',
    'no_kernel_evaluations': True,
    'Dmin': str(Dmin),
    'epsilon_split': '1/100',
    'L': str(L),
    'U': str(U),
    'point_G_upper_used_for_essential_g': str(G_point),
    'necessary_budget_ceiling_from_g_le_U': str(threshold),
    'coarse_upper_certificate_lower_envelope_at_split': str(coarse_at_split),
    'tail_upper_certificate_floor_at_split': str(tail_at_split),
    'coarse_to_ceiling_ratio': str(coarse_at_split / threshold),
    'tail_to_ceiling_ratio': str(tail_at_split / threshold),
    'direction_check': {
        'small_epsilon': 'The mesh formula is an upper certificate, not a lower bound on actual norm. Its RHS is forced >= L exp(Dmin/(32 eps^2)) and increases as eps decreases.',
        'large_epsilon': 'The whole-face tail is an upper certificate, not a lower bound on actual tail. Its allocated amount K eps^(7/6) increases as eps increases.',
        'cost_gate': 'For q=6/5, 6^q(B/g)^q < 64 iff B < 16g/3. Essential g <= U makes B < 16U/3 necessary for this fixed-budget proof route.',
    },
    'scope': 'Only excludes a proof using this fixed 120-leaf log-convex upper certificate, this spectral L insertion, this essential-g point cap, and K_new tail allocation. It says nothing about the actual H norm, actual tail, I, model cost, or C6.',
    'inputs_sha256': {
        'mesh': hashlib.sha256(mesh_path.read_bytes()).hexdigest(),
        'spectral': hashlib.sha256(spectral_path.read_bytes()).hexdigest(),
        'point': hashlib.sha256(point_path.read_bytes()).hexdigest(),
    },
    'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
}
Path(__file__).with_suffix('.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps(out, indent=2))
