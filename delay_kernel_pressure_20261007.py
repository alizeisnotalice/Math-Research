"""Exact component pressure for actual-source long geometric delay envelope.
Integer dyadic radial volumes avoid floating gate comparisons. No FIRST claim.
"""
from fractions import Fraction as F
from pathlib import Path
import json,hashlib
out=Path(__file__).with_name('delay_kernel_pressure_20261007_results.json')
if out.exists():raise SystemExit('Refuse overwrite of frozen results')
def pow2(j):return F(2)**j
rows=[]; boundaries=0
for stage,q in enumerate((3,5,7),1):
 n=2**q
 for k in sorted(set((0,n//4,n//2,n-1,n))):
  # a=1,b=2, R^n=2^k, exp(delay cutoff)=n=2^q.
  js=sorted(set((-q-1,-q,0,k-q-1,k-q,k-q+1,n-q,n-q+1)))
  for pattern in ('uniform','skew','zero_transmission'):
   raw=[1 if pattern=='uniform' else 2**i for i in range(len(js))]
   ms=[F(x,sum(raw)) for x in raw]
   gs=[F((3*i+1)%5,4) if pattern=='zero_transmission' else F(1) for i in range(len(js))]
   full=long=short=majorant=F(0)
   for j,m,g in zip(js,ms,gs):
    hard=pow2(-k) if j<=k else F(0)
    is_long=k-j>=q
    part=hard if is_long else F(0)
    env=F(1) if j<=-q else (pow2(-j-q) if j<=n-q else F(0))
    assert part<=env
    if j==k-q:
     assert part==env and part>0
     boundaries+=1
    full+=m*g*hard;long+=m*g*part;short+=m*g*(hard-part);majorant+=m*env
   assert full==long+short and long<=majorant
   rows.append({'stage':stage,'n':n,'q':q,'R_log2_volume':k,'pattern':pattern,'radial_log2_volumes':js,'source_masses':list(map(str,ms)),'transmission':list(map(str,gs)),'full':str(full),'long_delay':str(long),'short_delay':str(short),'fixed_column_majorant':str(majorant)})
# ln 2 = 2 sum_{k>=0} (1/3)^(2k+1)/(2k+1); positive tail upper bound.
integrals=[]
for stage,(n,terms) in enumerate(zip((8,32,128),(8,16,32)),1):
 lo=2*sum((F(1,3)**(2*j+1)/F(2*j+1) for j in range(terms)),F(0))
 tail=2*F(1,3)**(2*terms+1)/(F(2*terms+1)*(1-F(1,9)))
 hi=lo+tail
 assert F(0)<lo<hi<F(1)
 # Exact envelope norm = (1+n ln2)/n; hence <1+1/n uniformly.
 integrals.append({'stage':stage,'n':n,'series_terms':terms,'norm_lower':str(F(1,n)+lo),'norm_upper':str(F(1,n)+hi),'uniform_upper':str(1+F(1,n))})
res={'status':'PASS_EXACT_DELAY_ENVELOPE_COMPONENT','rows':rows,'row_count':len(rows),'equality_boundary_checks':boundaries,'integral_intervals':integrals,'scope':'Finite positive-source cube kernel, frozen subunit transmission and exact delay split; no original FIRST or short-shell budget certified.','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
out.write_text(json.dumps(res,indent=2));print({k:v for k,v in res.items() if k not in ('rows','integral_intervals')})
