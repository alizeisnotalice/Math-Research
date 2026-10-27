#!/usr/bin/env python3
"""Small exact arithmetic check, not numerical validation of the theorem."""
import hashlib
import json
import platform
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent
d, eps, k, b, hr = F(1,2**64),F(1,2**10),2**20,2**22,2**24
c = F(1,256*hr)
checks = {
    'fixed_gap_small': d<=F(1,16),
    'scalar_buffer_positive': b*d*d<F(1,2),
    'buffer_pays_weight_exponent': b==4*k and F(b,2)-k==k,
    'row_product_lower_uses_exp_minus_a_ge_1_minus_a': 24*eps<=F(1,2),
    'pinned_face_sum_second_moment_constant': 16*16*12<4096,
    'four_energy_losses_leave_one_over_128': F(1,64)-F(3,4096)-F(8,k)>=F(1,128),
    'deficit_constant': c==F(1,2**32),
    'scalar_overhead_relative_to_deficit': 2*b*d*d/(c*d)==F(1,512),
    'final_linear_coefficient': c*d/2==F(1,2**97),
    'pinned_source_fraction_positive': 1-d/2>=F(1,2),
}
assert all(checks.values())
sources = [HERE/'PROOF.md',HERE/'verify_constants.py',
    Path('/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round042/second_moment_cut/PROOF.md'),
    Path('/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round043/fractional_moment_cut/PROOF.md')]
receipt = {'status':'passed','python':platform.python_version(),'arithmetic':'standard-library Fraction',
    'checks':checks,'constant_values':{'d':str(d),'epsilon':str(eps),'K':k,'B':b,'H_R':hr,
        'raw_energy_lower':'1/64','good_energy_lower':'1/128','deficit_coefficient':str(c),
        'final_linear_coefficient':str(c*d/2)},
    'sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
    'scope':'Exact rational arithmetic only. All-coordinate moment estimates, eventual guards, actual-history transfer and operator conclusion require the analytic proof and independent review. No production integration, parameter search, MATLAB or old baseline rerun.'}
(HERE/'arithmetic_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':receipt['status'],'checks':checks},indent=2))
