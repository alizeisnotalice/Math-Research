#!/usr/bin/env python3
"""Small exact algebra checks; no search, LP solve, dimension grid or MATLAB."""
from fractions import Fraction as Q
from math import comb
from pathlib import Path
import hashlib
import json
import platform
import sys

HERE = Path(__file__).resolve().parent
exp_partial = 1+Q(7,10)+Q(7,10)**2/2+Q(7,10)**3/6
assert exp_partial == Q(12013,6000) > 2
base = 40*(Q(236,5)/2048+Q(53386,25)/2048**3)
assert base < 1
for n,q in ((2,1),(3,2),(16,7),(64,32),(256,255)):
    assert comb(n-1,q)*q >= n-1
    assert 2**q >= 2*q
    assert 2**q*comb(n-1,q) >= 2*(n-1)
# A two-witness incidence LP primal and dual equality, without a solver.
# Left/common/right cells have unit volume; each witness has volume two.
m = (Q(0),Q(2),Q(0))
z = (Q(1,2),Q(1,2))
assert m[0]+m[1] >= 2 and m[1]+m[2] >= 2
assert z[0] <= 1 and sum(z) <= 1 and z[1] <= 1
assert sum(m) == 2*sum(z) == 2
paths = [
 Path('/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round032/multiscale_theory/REPORT.md'),
 Path('/Users/zhengzhihao/Desktop/T/cube-r030-soft/cube-dim-7-bounds/upper/round032/joint_cover_theory/proof.md'),
 Path('/Users/zhengzhihao/Desktop/T/cube-r030-parameters/cube-dim-7-bounds/upper/round032/anisotropic_cover/proof.md'),
 Path('/Users/zhengzhihao/Desktop/T/cube-r030-parameters/cube-dim-7-bounds/upper/round034/constant_theory/PROOF.md'),
]
result = {'status':'PASS','python':sys.version,'platform':platform.platform(),
          'arithmetic':'stdlib fractions.Fraction and exact integer combinatorics',
          'exp_7_over_10_partial_sum':str(exp_partial),
          'normalized_W_base_upper_bound_at_L_2pow22':str(base),
          'source_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
          'scope':'finite algebra only; analytical proofs are in the review files'}
(HERE/'small_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
