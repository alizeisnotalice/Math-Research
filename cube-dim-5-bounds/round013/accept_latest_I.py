#!/usr/bin/env python3
"""Primary exact root-budget acceptance, after analytic/source review.

This binds the full094 audit and complete092 budget; no H or old baseline is
recalculated. The cap comparisons are comparisons of certified upper budgets.
"""
from fractions import Fraction as F
from pathlib import Path
import json
import hashlib
import os
import time

HERE = Path(__file__).resolve().parent


def main():
    start = time.monotonic()
    regpath = HERE/'primary_I092_registration.json'
    reg = json.loads(regpath.read_text())
    for entry in reg['sources'].values():
        assert hashlib.sha256(Path(entry['path']).read_bytes()).hexdigest() == entry['sha256']
    audit = json.loads(Path(reg['sources']['094_whole_results']['path']).read_text())
    assert audit['status'] == 'pass' and audit['nodes'] == 56139
    assert audit['bands'] == 4098147 and audit['blocks'] == 14
    assert audit['all_new_bins_exact_match'] and audit['rational_sum_fallbacks'] == 0
    assert audit['physical_volume'] == '1' and audit['physical_frontier'] == 55819
    assert audit['new_H_calls'] == audit['Arb_log_exp_asinh_calls'] == 0
    summary = json.loads(Path(reg['sources']['092_summary']['path']).read_text())
    budget = summary['budget']
    U, B, E, I = [F(budget[key]) for key in
                  ('UH_truncated_integral', 'B_upper', 'fixed_E076_upper', 'I_upper')]
    cap = F(summary['I_upper_simple'])
    assert cap == F(audit['I_upper_simple']) == F(6010267, 10**11)
    assert B > E > 0 and (B-E)**5 >= U**4
    assert I**4 >= B**5 and cap >= I
    old_088, old_010 = F(70301737, 10**12), F(71555899, 10**12)
    assert cap < old_088 < old_010
    g = F(budget['valid_g'])
    assert g == F(1, 10000) and budget['no_tube_added']
    assert 3125*cap**4 > 32**4*g**5
    result = {
        'status': 'primary_accepts_complete_same_alpha1_I092',
        'I_upper': str(cap),
        'previous088_cap_minus_new_exact': str(old_088-cap),
        'previous010_cap_minus_new_exact': str(old_010-cap),
        'full_space_input': 'single positive face F, complete G=H512+shorttime, no tube added',
        'positive_bands_and_full_face': True,
        'shared_accepted_H_and_M2': True,
        'independent_new_ratio_range_mean_chord_hierarchy': '094 full56139 and4098147 bands',
        'primary_analytical_review': '075,076,085,091 and094 analysis/verify_whole read in full; same029 generalL1/allpositive scales retained',
        'qualified_g': str(g),
        'paired_cap_cost_fourth_power': str(3125*cap**4/g**5),
        'cost_still_above32': True,
        'formal_C5_upper': '32',
        'formal_C5_lower_changed': False,
        'formal_C5_upper_changed': False,
        'dual_breakthrough_time': None,
        'post_dual_rounds': 0,
        'new_H_or_special_function_evaluations': 0,
        'pid': os.getpid(),
        'seconds': time.monotonic()-start,
        'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'registration_sha256': hashlib.sha256(regpath.read_bytes()).hexdigest(),
        'frozen_inputs': reg['sources'],
    }
    out = HERE/'primary_I092_acceptance.json'
    assert not out.exists()
    out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('frozen_inputs','paired_cap_cost_fourth_power')},ensure_ascii=False))


if __name__ == '__main__':
    main()
