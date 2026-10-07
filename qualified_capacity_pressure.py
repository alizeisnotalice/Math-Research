"""Three rounds of exact qualified-capacity cutoff diagnostics.
Geometric hard shells are encoded by integer log-volume j: 2||z||inf=2**(j/n).
Soft rows and gates are algebraic subprobability fixtures, NOT original FIRST.
No spatial quadrature or endpoint certification is performed.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json
rows=[]
alpha=F(1,1024) # alpha = 2 delta lambda; a=lambda=1
for round_id,s in enumerate((4,16,64),1):
 n=s*s
 js=(0,s//2,s,s+1,2*s,n)
 for pattern in ('uniform','skew','zeros'):
  p0=[1]*8 if pattern=='uniform' else ([1,2,4,8,16,32,64,128] if pattern=='skew' else [0,0,1,0,3,0,4,0])
  ps=[F(v,sum(p0)) for v in p0]
  # Unique LCA for 8 fixed binary leaf labels. Actual geometric soft FIRST unasserted.
  gs=[]
  for h in range(6):
   val=F(0)
   for i,p in enumerate(ps):
    if i==h:continue
    lca_height=(i^h).bit_length()
    gate=F(((i+2*h+lca_height)%5),4)
    val+=p*gate
   assert 0<=val<=1
   gs.append(val)
  for size in (F(1,4),F(1),F(4)):
   masses=[alpha*2**s*size*F(j+1,21) for j in range(6)]
   for k in (0,s//2,s,s+1,2*s,n):
    # R**n=2**k and hard membership j<=k, including closed boundaries.
    C=sum((m for j,m,g in zip(js,masses,gs) if j<=k and g>0),F(0))
    T=sum((m*g/F(2**k) for j,m,g in zip(js,masses,gs) if j<=k),F(0))
    small=C<=alpha*2**s
    near=k<=s
    envelope=sum((m/F(2**max(0,j)) for j,m in zip(js,masses) if j<=s),F(0))
    if near: assert T<=envelope
    if small and not near: assert T<=alpha
    rows.append(dict(round=round_id,n=n,pattern=pattern,mass_factor=str(size),radius_log2_volume=k,capacity=str(C),traffic=str(T),small_capacity=small,near_window=near,fixed_kernel_column_bound=str(envelope)))
result={'status':'PASS_EXACT_QUALIFIED_CAPACITY_CUTOFF','rows':rows,'row_count':len(rows),'rounds':3,'scope':'Exact cube membership and finite conditional LCA gate arithmetic; no original soft FIRST, source-high split or continuous early-exit certification; no integrated geom estimate certified.','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
out=Path(__file__).with_name('qualified_capacity_pressure_results.json')
if out.exists():raise SystemExit('Refuse overwrite')
out.write_text(json.dumps(result,indent=2))
print(json.dumps({k:v for k,v in result.items() if k!='rows'}))
