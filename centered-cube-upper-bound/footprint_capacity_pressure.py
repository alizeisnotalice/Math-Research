"""Exact pressure of universal full-footprint + small-mass sufficient condition.
Uniform unit cube, a=1/4, b=1/2; no original geom gate certification.
"""
from fractions import Fraction as F
from pathlib import Path
import json,hashlib
rows=[]
for stage,(n,s) in enumerate(((16,4),(256,16),(4096,64)),1):
 assert s*s==n
 # exp(D)=2**sqrt(n), hence D=log(2)*sqrt(n); alpha=1 is a generous bound.
 allowed=F(2**s);required=F(4,3)**n
 assert allowed<required
 rows.append(dict(stage=stage,n=n,D=f'{s}*log(2)',alpha='1',allowed=str(allowed),required=str(required),necessary_capacity_fails=True))
p=Path(__file__).with_name('footprint_capacity_pressure_results.json')
if p.exists():raise SystemExit('Refuse overwrite')
p.write_text(json.dumps({'status':'PASS_EXACT_OBSTRUCTION_CHECK','rows':rows,'scope':'No uniform raw-full-footprint construction under stated small-mass cap; not a cube maximal or actual geom counterexample.','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2))
print('3 exact capacity comparisons passed; no original traffic evaluated')
