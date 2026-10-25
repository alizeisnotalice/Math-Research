"""Exact Q1 counterexample and source-hash record; stdlib arithmetic only.

Does not solve a continuous kernel problem or re-solve the main finite LP.
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import hashlib, json

OUT = Path(__file__).resolve().parent
MAIN = Path('/Users/zhengzhihao/Desktop/T/Math-Research-finite-dim/cube-finite-dim-bounds')
AUDIT = Path('/Users/zhengzhihao/Desktop/T/Math-Research-finite-audit/cube-finite-dim-bounds')

def source_record():
    names = ['rounds/005-expansion.md', 'results/round005_expansion_snapshot.json',
             'rounds/005-registration.md', 'STATE.md', 'results/bounds_ledger.json',
             'src/free_kernel_lp.py', 'src/check_free_kernel_dual.py',
             'proofs/finite_grid_kernel_diagnostic.md', 'results/round005_lp_dual.json']
    paths = [MAIN / name for name in names] + [AUDIT / 'audit/007b_results.json']
    return {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}

def value(k):
    if max(map(abs, k)) >= 2:
        return F(0)
    return [F(11,4), F(5,8), F(1,4), F(1,8)][sum(x != 0 for x in k)]

def d2(k, axis):
    a, b = list(k), list(k)
    a[axis] -= 1
    b[axis] += 1
    return value(a) + value(b) - 2 * value(k)

def run():
    # Complete support plus its one-step exterior collar, with exact arithmetic.
    nodes = list(product(range(-2, 3), repeat=3))
    lap = {k: sum((d2(k,j) for j in range(3)), F(0)) for k in nodes}
    assert lap[(0,0,0)] == -F(51,4)
    assert all(v >= 0 for k,v in lap.items() if k != (0,0,0))
    assert all(lap[k] == 0 for k in product(range(-1,2),repeat=3) if k != (0,0,0))
    assert sum(lap.values(), F(0)) == 0
    assert min(value(k) for k in product(range(-1,2),repeat=3)) == F(1,8)
    mass = sum((value(k) for k in nodes), F(0))
    assert mass == F(21,2)  # exact integral of tensor hats of grid width h=1
    # On x=0, 0<y,z<1, the jump density is bilinear in these four values.
    face = [d2((0,y,z),0) for y,z in product((0,1),repeat=2)]
    assert face == [-F(17,4), -F(3,4), -F(3,4), -F(1,4)]
    center_density = sum(face, F(0))/4
    assert center_density == -F(3,2)
    # Exact test: phi = tent(x/(1/4))*tent((y-1/2)/(1/4))*tent((z-1/2)/(1/4)).
    # Support meets only plane x=0 and avoids every other grid plane and origin.
    # Symmetry makes integration of a bilinear density equal its midpoint value
    # times (integral tent)^2 = (1/4)^2. phi(0,.,.) is the tangential product.
    test_value = center_density/F(16)
    assert test_value == -F(3,32)
    # Finite checker equivalent to the declared main model at q=1,R=3:
    model_nodes = list(product(range(-2,3),repeat=3))
    assert all(value(k) >= (F(1,8) if max(map(abs,k))<=1 else 0) for k in model_nodes)
    assert all(lap[k] >= 0 for k in model_nodes if k != (0,0,0))
    rec = {
      'status': 'PASS exact counterexample',
      'scope': 'Q1 bridge failure only; no new maximal weak-type bound or continuous-family obstruction',
      'arithmetic': 'Python stdlib Fraction; derivation and checker share formulas',
      'dimension': 3, 'h': '1', 'zero_at_max_abs_index': 2,
      'nonzero_orbit_values': {'0_nonzero': '11/4', '1_nonzero': '5/8', '2_nonzero': '1/4', '3_nonzero': '1/8'},
      'checked_nodes_including_collar': len(nodes),
      'finite_off_origin_min_laplacian': str(min(v for k,v in lap.items() if k!=(0,0,0))),
      'origin_laplacian': str(lap[(0,0,0)]),
      'finite_laplacian_total': str(sum(lap.values(),F(0))),
      'cube_obstacle_minimum': '1/8', 'exact_Q1_integral': str(mass),
      'negative_face': 'x=0, 0<y,z<1', 'negative_face_corner_jump_densities': list(map(str,face)),
      'negative_face_midpoint_jump_density': str(center_density),
      'nonnegative_compact_tent_test_distribution_pairing': str(test_value),
      'test_support': '[-1/4,1/4] x [1/4,3/4] x [1/4,3/4]; avoids origin',
      'continuous_laplacian_nonnegative_off_origin': False,
      'new_U3': False, 'new_U4': False}
    (OUT/'cb5_results.json').write_text(json.dumps(rec,indent=2)+'\n')
    (OUT/'cb5_sources.json').write_text(json.dumps(source_record(),indent=2)+'\n')
    print(json.dumps(rec,indent=2))

if __name__ == '__main__': run()
