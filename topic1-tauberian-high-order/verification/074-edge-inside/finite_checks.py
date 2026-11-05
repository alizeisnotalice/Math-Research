"""Finite rational checks accompanying the geometric R074 proof.

No floating-point experiments; these checks do not prove the covering lemma.
"""
from fractions import Fraction as F
from itertools import product, combinations
from pathlib import Path
import hashlib
import json

out = {}
K = 8
out['tube_large_radius_fraction_bound'] = str(F(12, K*K))
assert F(12, K*K) <= F(1, 4)
H = 10*(K+1)
out['main_edge_separation_margins'] = {
    'free_endpoint_lower': H-K,
    'transverse_upper': 2*K+1,
    'midpoint_gap_each_side': H-K,
}
assert H-K > 2*K+1

# Exact finite tests of the elementary common-orthant inequality behind the
# 2^3 overlap bound. The analytic proof permits all real coordinates/radii.
coords = list(product([F(0), F(1, 3), F(2, 3), F(1)], repeat=3))
checks = 0
for x,y in combinations(coords,2):
    rx,ry = max(x),max(y)
    assert max(abs(x[i]-y[i]) for i in range(3)) <= max(rx,ry)
    checks += 1
out['same_orthant_exact_checks'] = checks

# Canonical whole-pair layer: source cost / ell at v=0 has coefficient
# p*q*(p+q)/n, p=q=n/3, per undirected normal edge. Actual slot capacity is
# choose(n,2) at leading order; n=3 gives 2K_3=2.
out['canonical_volume_coefficients'] = []
for n in range(2,8):
    p = F(n,3)
    per_edge = p*p*(p+p)/n
    slots = n*(n-1)//2
    Kn = F(slots*n*n,27)
    assert per_edge*slots == 2*Kn
    out['canonical_volume_coefficients'].append({
        'n':n, 'per_edge':str(per_edge), 'slot_count':slots,
        'whole':str(per_edge*slots), '2Kn':str(2*Kn)})

# Ledger units: physical deletion is priced by v/a; inside bound contributes
# a times each missing physical volume. L=vR, h=v^3R.
n = 3
out['n3_exact_power_accounting'] = {
    'normal_physical':n+5,
    'normal_inside_deletion':2*n+5,
    'normal_inside_ledger':n+6,
    'collar_physical':9,
    'collar_inside_deletion':n+9,
    'collar_inside_ledger':10,
    'raw_physical_from_old_Rstar8':n+7,
    'extension_physical_from_F8':9,
    'raw_inside_ledger':n+8,
    'extension_inside_ledger':10,
    'bad_weight_from_threshold_av5':9-(n+5),
    'bad_normal_inside_ledger':n+7,
}
assert out['n3_exact_power_accounting']['normal_inside_ledger'] == 9
assert min(out['n3_exact_power_accounting'][x] for x in
           ['collar_inside_ledger','raw_inside_ledger',
            'extension_inside_ledger','bad_normal_inside_ledger']) == 10
out['all_assertions_passed'] = True

base=Path(__file__).resolve().parent
out['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
(base/'finite_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
