#!/usr/bin/env python3
"""Read exact frozen q fields and check only the comparison's small gates."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import os

root=Path(__file__).resolve().parent
base=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-7/cube-dim-7-bounds')
matrix=base/'upper/round034/compute/c34a_results.json'
assert hashlib.sha256(matrix.read_bytes()).hexdigest()=='87ac7a68bf3979f6b654f7954044feac4d0d26031e14dc38152c8e48684d4a3d'
data=json.loads(matrix.read_text())
rows=[{'Y':r['Y'],'q':str(F(r['certificate']['q'])),'L_min':max(64,r['Y']**2)} for r in data['results']]
assert len(rows)==12 and all(0<F(r['q'])<F(1,8) for r in rows)
best=max(rows,key=lambda r:F(r['q']))
assert best['Y']==256 and F(best['q'])==F(74275164113,10**12)
assert F(39,320)<F(1,8)
assert 240*F(2,3)>128 and 240<16**2
gap=F(1,2**32)-F(1,2**35)-F(1,2**158)
assert gap>0
assert F(32,2**40)==F(1,2**35)
assert F(4,2**160)==F(1,2**158)
assert 2**240>=2**108 and all(2**240>=r['L_min'] for r in rows)
theory=Path('/Users/zhengzhihao/Desktop/T/cube-r030-parameters/cube-dim-7-bounds/upper/round042/comparison')
sources=[theory/'PROOF.md',theory/'verify.py',theory/'verification.json',matrix,
 Path('/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round040/energy_refinement/PROOF.md'),
 base/'upper/round030/parameter_theory/proof.md',base/'upper/round030/geometry_theory/REPORT.md',
 base/'upper/round034/constant_theory/PROOF.md',
 Path('/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round042/second_moment_cut/PROOF.md'),
 Path('/Users/zhengzhihao/Desktop/T/cube-r030-parameters/cube-dim-7-bounds/upper/round042/explicit_threshold/THRESHOLD_PROOF.md')]
hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
assert hashes[str(sources[0])]=='e9dc339e64164f4afb51e20ee81dd03d38f49132353469849de906ec84844c6a'
receipt={'status':'PASS','auditor_pid':os.getpid(),'arithmetic':'stdlib Fraction and exact integer powers',
 'exact_q_rows':rows,'largest_displayed_certified_q':best,'normalized_power_gap':str(gap),
 'source_sha256':hashes,
 'scope':'read q fields only; no parameter recertification, LP rerun, production or decimal reconstruction'}
(root/'exact_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:v for k,v in receipt.items() if k!='exact_q_rows'},indent=2))
