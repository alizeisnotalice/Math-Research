"""Read saved exact category masses; recompute Z and actual weights, no MC."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
P=ROOT/'theory'
r=json.loads((P/'imbalance_results.json').read_text())
src=ROOT/'compensation/source_jitters.json'
assert hashlib.sha256(src.read_bytes()).hexdigest()==r['source_jitters_sha256']
p=list(map(F,r['rho_class_mass']));R=len(p)
assert sum(p,F(0))==1
d=sum((abs(x-F(1,R))for x in p),F(0))
assert d==F(r['d_L1_bound'])<1
def conv(a,b):return [sum((a[j]*b[(q-j)%R]for j in range(R)),F(0))for q in range(R)]
rows=[]
for saved in r['rows']:
    n=saved['n'];m=[]
    for b in range(R):
        pb=[p[(a-b)%R]for a in range(R)]
        dist=[F(int(a==0))for a in range(R)]
        for i in range(n):dist=conv(dist,pb)
        m.append(dist[b])
    Z=sum(m,F(0));zeta=(R-1)*d**n;eps=(R-1)*((1+d)/2)**n
    assert zeta<1 and abs(Z-1)<=zeta and all(abs(R*x-1)<=zeta for x in m)
    assert float(Z-1)==saved['Z_minus_one_float']
    assert float(zeta)==saved['candidate_zeta_float']
    weights=[x/Z for x in m]
    assert sum(weights,F(0))==1 and min(weights)>0
    rows.append({'n':n,'Z_minus_one_exact':str(Z-1),'zeta_exact':str(zeta),
                 'actual_weights_exact':list(map(str,weights)),
                 'all_source_and_normalizer_bounds_passed':True,
                 'macro_lower_factor_display':float(max(F(0),1-eps)/(1+zeta)),
                 'macro_upper_factor_display':float((1+eps)/(1-zeta))})
out={'scope':'Exact saved-category-mass convolution arithmetic only; no cube/tail simulation.',
     'source_jitters_hash_passed':True,'d_L1_bound_exact':str(d),'rows':rows}
(HERE/'imbalance_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print('Six exact Z and actual-prior checks passed; source hash agrees.')
