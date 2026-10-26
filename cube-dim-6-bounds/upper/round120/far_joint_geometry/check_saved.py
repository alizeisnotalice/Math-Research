#!/usr/bin/env python3
"""Minimal exact saved-data checker; no author generator import or rerun."""
import sys
sys.set_int_max_str_digits(0)
import hashlib
import json
from pathlib import Path
from fractions import Fraction as F
r=Path(__file__).resolve().parent
inp=json.loads((r/'upstream.json').read_text())
out=json.loads((r/'results.json').read_text())
assert out['status']=='COMPLETE_TWO_FIXED_CONFIGURATIONS'
g,b=out['configurations']
assert F(g['tail_coefficient_over_pi'])==F(3,2)
assert F(g['tube_coefficient_over_pi_squared'])==F(1519,1152)
P0=F(3665,10368)
assert F(g['P_far_lower_coarse'])==P0==1-F(1,2)-F(1519,10368)
W5=F(45,2**86)*P0**5
RL,RU=F(inp['far_allowance']['lower']),F(inp['far_allowance']['upper'])
assert F(b['J_far_lower_fifth_exact'])==W5
assert F(b['residual_vs_RL_fifth'])==W5-RL**5>0
assert F(b['residual_vs_RU_fifth'])==W5-RU**5>0
lo,hi=map(F,b['J_far_lower_interval'])
assert lo**5<=W5<=hi**5
ratioL,ratioU=map(F,b['intrinsic_lower_over_RU_interval'])
assert ratioL**5<=W5/RU**5<=ratioU**5
up=F(b['accepted108_complete_actual_far_I_upper_interval'][1])
assert up**5>=F(b['accepted108_complete_actual_far_upper_fifth_input'])
assert not b['actual_I_far_lower_bound'] and not out['full_far_competitive']
assert not out['new_C6_bound'] and not b['global_g0_proved']
print(json.dumps({'status':'PASS_SAVED_ONLY','results_sha256':hashlib.sha256((r/'results.json').read_bytes()).hexdigest(),'new_configurations':0,'native_calls':0}))
