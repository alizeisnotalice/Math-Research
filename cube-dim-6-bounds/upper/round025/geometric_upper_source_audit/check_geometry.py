#!/usr/bin/env python3
"""Eight bounded exact Fraction checks; no kernel or spatial enumeration."""
from fractions import Fraction as F
from pathlib import Path
import json, hashlib, datetime
root=Path(__file__).resolve().parent
plan=json.loads((root/'geometry_plan.json').read_text())
e=F(plan['epsilon']); eta=F(plan['eta'])
results=[]
def check(name, predicate, **data):
    assert predicate, name
    results.append({'name':name, 'pass':True, **{k:str(v) for k,v in data.items()}})
check('01_epsilon_and_center_exclusion', 0<e<F(1,2) and 2*(1-e)>1,
      epsilon=e, center_distance=2*(1-e))
check('02_common_intersection', (2*e)**6/F(64)==e**6,
      common_volume=(2*e)**6, ratio_to_one_cube=e**6)
ratios=[e**k for k in range(1,7)]
check('03_pair_overlap', all(0<r<=e for r in ratios), ratios=ratios)
red=1-(1-e/2)**6
check('04_union_redundancy', (4-2*e)**6/F(4096)==1-red and 0<red<3*e,
      union_volume=(4-2*e)**6, redundancy_ratio=red)
t=(1-eta)**6
check('05_corner_mass_threshold', 21*t<1<22*t,
      single_corner_threshold=t, twenty_one_times=21*t, twenty_two_times=22*t)
delta=43*eta**6/F(64); gamma=64-delta
check('06_empty_corner_gap', delta==F(43,15625) and gamma==F(999957,15625),
      empty_corners=43, delta=delta, gamma=gamma)
check('07_large_mass_branch', F(64)/(22*t)<gamma,
      large_mass_bound=F(64)/(22*t), gamma=gamma)
theta=delta/256
check('08_conditional_leakage_arithmetic', gamma/(1-theta)+delta/4<=64-delta/2,
      theta=theta, leakage_budget=delta/4, conditional_global_constant=64-delta/2)
assert len(results)==8
out={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
     'arithmetic':'Python fractions.Fraction; exact integers and rationals',
     'named_checks':len(results), 'kernel_calls':0, 'matlab_calls':0,
     'new_spatial_grid_points':0, 'proof_status':'General statements proved analytically in theorem.md; finite checks are sanity checks.',
     'input_sha256':hashlib.sha256((root/'geometry_plan.json').read_bytes()).hexdigest(),
     'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), 'results':results}
(root/'geometry_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
