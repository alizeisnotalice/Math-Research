"""Independent reconstruction of saved capacity-tree arithmetic, not spatial tests."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import hashlib
import json

BASE = Path(__file__).parent
source = BASE / 'capacity_pair_occupancy_exact_guard_results_20261007.json'
data = json.loads(source.read_text())
script = BASE / 'capacity_pair_occupancy_exact_guard_20261007.py'
assert hashlib.sha256(script.read_bytes()).hexdigest() == data['script_sha256']
counts = {'threshold_rows': 0, 'capacity_levels': 0, 'compressed_chains': 0}

def examine(rows, saved):
    union = sum((v*(p+r-p*r) for p,r,v,g in rows), F())
    failure = sum((v*g*(1-p)*(1-r) for p,r,v,g in rows), F())
    for entry in saved:
        delta = F(entry['delta'])
        hi = [(p,r,v,g) for p,r,v,g in rows if p+r-p*r >= delta]
        high = sum((v*g*(1-p)*(1-r) for p,r,v,g in hi), F())
        assert F(entry['union']) == union
        assert F(entry['failure']) == failure
        assert F(entry['high_failure']) == high
        assert F(entry['low_failure']) == failure-high
        assert entry['high_count'] == len(hi)
        assert entry['low_count'] == len(rows)-len(hi)
        assert entry['equality_count'] == sum(p+r-p*r == delta for p,r,_,_ in rows)
        assert high <= (1/delta-1)*union
        assert union+high <= union/delta
        counts['threshold_rows'] += 1

for run in data['rounds']:
    n, root = run['n'], run['sqrt_n']
    assert root**2 == n
    k = 0
    while 2**k < n+2:
        k += 1
    assert k == run['k_n'] and F(run['delta_n']) == F(1,k**4)
    for model in run['models'][:2]:
        b, depth = model['arity'], model['depth']
        paths = list(product(range(b), repeat=depth))
        skew = model['model'].startswith('skew')
        raw = [2**path[-1] if skew else 1 for path in paths]
        masses = [F(v,sum(raw)) for v in raw]
        born = [F(0) if skew and i%3 == 0 else v for i,v in enumerate(masses)]
        def mass(prefix, vector):
            return sum((v for path,v in zip(paths,vector) if path[:len(prefix)] == prefix),F())
        def coins(prefix, child):
            cm,ch = mass(child,masses),mass(child,born)
            sm,sh = mass(prefix,masses)-cm,mass(prefix,born)-ch
            return (min(F(1),cm/(root*sh)) if sh else F(1),
                    min(F(1),ch/(n*sm)) if sm else F(1),sm,sh)
        for saved in model['capacities']:
            length = depth-saved['level']
            forward = reverse = F()
            for prefix in product(range(b), repeat=length):
                for c in range(b):
                    p,r,sm,sh = coins(prefix,prefix+(c,))
                    forward += p*sh
                    reverse += r*sm
            assert forward == F(saved['forward']) <= F(1,root)
            assert reverse == F(saved['reverse']) <= sum(born)/n
            counts['capacity_levels'] += 1
        rows = []
        for i,y in enumerate(paths):
            for j,z in enumerate(paths):
                if i == j or not born[j]:
                    continue
                split = next(t for t in range(depth) if y[t] != z[t])
                parent = y[:split]
                p = coins(parent,y[:split+1])[0]
                r = coins(parent,z[:split+1])[1]
                rows.append((p,r,masses[i]*born[j]*F((i+3*j)%5+1,5),F((2*i+j)%7,6)))
        assert len(rows) == model['ordered_positive_distinct_pairs']
        assert sum(masses) == F(model['whole_mass']) and sum(born) == F(model['born_mass'])
        examine(rows,model['thresholds'])
    model = run['models'][2]
    b,D = model['occupied_children'],model['chain_depth']
    p,r = F(1,root*(b-1)),F(1,n*(b-1))
    ell,w = p+r-p*r,(1-p)*(1-r)
    for name,v in [('p',p),('r',r),('ell',ell),('w',w)]:
        assert v == F(model[name])
    assert ell < F(run['delta_n'])
    examine([(p,r,F(b*b//2-b,b*b),F(1,2)),
             (p,r,F(1,2),F(1,3))],model['star_thresholds'])
    rings = [F(b-1,b**j) for j in range(D,0,-1)]
    assert rings == list(map(F,model['ring_masses']))
    assert sum(rings)+F(model['inner_leaf_mass']) == 1
    assert sum(rings) == F(model['ring_total'])
    assert w*sum(rings) == F(model['failure_weighted_ring_total']) > F(99,100)
    counts['compressed_chains'] += 1

out = {'status':'PASS','counts':counts,
       'scope':'Independent prefix-path reconstruction of saved finite trees and exact fractions; no actual FIRST or spatial certification',
       'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
       'author_script_sha256':data['script_sha256']}
(BASE/'root_capacity_saved_review_20261007.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out))
