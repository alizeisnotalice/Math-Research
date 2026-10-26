#!/usr/bin/env python3
"""Exact comparisons of full budgets; these are not lower bounds for I."""
from pathlib import Path
from fractions import Fraction as F
import json
import hashlib
import time
import os

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent


def main():
    started = time.monotonic()
    summary = json.loads((PROJECT / 'upper_adapt082/summary.json').read_text())
    b = summary['budget']
    U, B, E, I = [F(b[k]) for k in
                  ('UH_truncated_integral', 'B_upper', 'fixed_E076_upper', 'I_upper')]
    old_U = F(1730073, 8000000000)
    old_080_I = F(220694961, 10**12)
    best_I_cap = F(71555899, 10**12)
    g = F(1, 10000)
    assert 0 < U < old_U and 0 < I < old_080_I
    assert B > E > 0 and (B-E)**5 >= U**4
    assert I**4 >= B**5
    assert I > best_I_cap
    assert 5**5 * B**5 > 32**4 * g**5
    assert 3125 * best_I_cap**4 > 32**4 * g**5
    result = {
        'status': 'exact_full_budget_comparisons_pass',
        '080_UH_cap_minus082_UH': str(old_U-U),
        '080_I_cap_minus082_I': str(old_080_I-I),
        '082_I_cap_minus_inherited_best_I_cap': str(I-best_I_cap),
        '082_B_norm_recipe_root_checks': True,
        '082_I_cap_fourth_power_covers_B_fifth_power': True,
        'old_best_I_with_qualifiedg100_cost_fourth_power': str(3125*best_I_cap**4/g**5),
        'both_budget_costs_strictly_above32': True,
        'warning_on_interpretation': 'Only exact differences of guarantee budgets; no lower bound on actual I or C5 is inferred',
        'formal_C5_upper_bound': '32',
        'no_new_C5_bounds': True,
        'new_H_or_special_function_calls': 0,
        'pid': os.getpid(),
        'seconds': time.monotonic()-started,
        'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        '082_summary_sha256': hashlib.sha256((PROJECT/'upper_adapt082/summary.json').read_bytes()).hexdigest(),
    }
    output = HERE/'budget_comparison.json'
    assert not output.exists()
    output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if 'minus' not in k and 'fourth_power' not in k},ensure_ascii=False))


if __name__ == '__main__':
    main()
