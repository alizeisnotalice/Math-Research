#!/usr/bin/env python3
"""Exact saved-only checker, including original-failure/continuation lineage."""
import sys
sys.set_int_max_str_digits(0)
import hashlib
import json
from pathlib import Path
from fractions import Fraction as F
root = Path(__file__).resolve().parent
def read(name): return json.loads((root/name).read_text())
def sha(name): return hashlib.sha256((root/name).read_bytes()).hexdigest()
inp, out = read('upstream.json'), read('results.json')
assert out['status'] == 'COMPLETED_2_FIXED_CONFIGURATIONS_AFTER_OUTPUT_ONLY_V2'
assert out['cfg1_original_checkpoint_sha256'] == sha('checkpoint.json')
assert out['original_failure_receipt_sha256'] == sha('run_receipt.json')
assert read('run_receipt.json')['status'] == 'FAILED'
assert read('run_receipt.json')['completed'] == 1
assert out['cfg_window_calls_v2'] == 0 and out['cfg_tail_calls_v2'] == 1
assert read('observed_tool_receipt_v2.json')['scientific_child_exit_code'] == 0
w, t = out['configurations']
assert w == read('checkpoint.json')['completed_configurations'][0]
RL, RU = F(inp['far_allowance']['lower']), F(inp['far_allowance']['upper'])
assert 0 < RL <= RU
assert all(F(d) > F(1, 5) for d in w['D_at_b_certified_lower'])
assert F(w['source_probability_lower']) == F(23, 24)**31-F(32, 180**30) > F(15, 64)
assert F(w['source_mass_window_lower']) > F(7, 1760)
assert F(w['far_volume_lower']) == 1-150*F(1, 64)**4 > F(1, 2)
W5 = F(1, 64)**6 * F(7, 1760)**5 / 5
assert F(w['I_envelope_lower_to_fifth']) == W5
assert F(w['residual_to_RL_fifth']) == W5-RL**5 > 0
assert F(w['residual_to_RU_fifth']) == W5-RU**5 > 0
lo, hi = map(F, w['I_envelope_lower_interval'])
assert lo**5 <= W5 <= hi**5
tailL, tailU = map(F, t['tail_norm_interval'])
pL, pU = map(F, t['P6_interval'])
assert tailL == pL/(5*F(1, 8)**5) and tailU == pU/(5*F(1, 8)**5)
assert F(t['new_over_old_exact']) == F(135, 4194304)
assert F(t['tail_below_one_twentieth_remaining_norm_residual']) == RL**5/20**6-tailU**6 > 0
sl, su = map(F, t['short_time_norm_allocation_interval'])
assert 0 < sl and (sl+tailU)**6 <= RL**5 <= (su+tailL)**6
assert not w['actual_I_far_lower_bound'] and not out['full_far_paid']
assert not out['new_C6_bound'] and not t['global_g0_proved']
print(json.dumps({'status': 'PASS_SAVED_EXACT_ONLY_V2_LINEAGE',
                  'results_sha256': sha('results.json'),
                  'new_scientific_configurations': 0, 'kernel_peak_grid_runs': 0}))
