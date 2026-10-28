"""Exact mass-normalization controls from the actually saved dyadic source."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib
OUT=Path(__file__).resolve().parent
src=OUT.parent/'compensation/source_jitters.json'
if not src.exists():src=OUT.parent/'round036_compensation/source_jitters.json'
sources=json.loads(src.read_text());one=next(s for s in sources if s['K']==63)
K=63;R=3;v=8;delta=F(1,4*K*v);u=list(map(F,one['U']));atoms=[(j+z)/K for j,z in enumerate(u)]
# Exact periodic class masses of the smoothed peak law, with genuine cross-cell transfers.
mu=[F(0)]*R
for a in atoms:
 for j in range(K):
  for sh in [-1,0,1]:
   ell=max(F(0),min(a+delta/2,F(j+1,K)+sh)-max(a-delta/2,F(j,K)+sh))
   mu[j%R]+=ell/(K*delta)
assert sum(mu)==1
p=[F(1,2*R)+z/2 for z in mu] # Low modes M3<K/R integrate to zero over each class.
d=sum(abs(z-F(1,R))for z in p)
def conv(a,b):return [sum((a[j]*b[(r-j)%R]for j in range(R)),F(0))for r in range(R)]
rows=[]
for n in [1,4,7,15,31,63]:
 z=F(0)
 for b in range(R):
  pb=[p[(r-b)%R]for r in range(R)]
  coeff=[F(1),F(0),F(0)]
  for i in range(n):coeff=conv(coeff,pb)
  z+=coeff[b]
 bound=(R-1)*d**n;assert abs(z-1)<=bound
 rows.append({'n':n,'Z_minus_one_float':float(z-1),'candidate_zeta_float':float(bound),
 'exact_inequality_passed':True,'normalizer_positive':z>0})
res={'source_jitters_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'K':K,'R':R,'v':v,
 'rho_class_mass':[str(z)for z in p],'d_L1_bound':str(d),'rows':rows,
 'scope':'Exact source mass controls for cyclically shifted tables; general unbalanced theorem separate.'}
(OUT/'imbalance_results.json').write_text(json.dumps(res,indent=2)+'\n')
for r in rows:print(r)
