"""Four new root enclosures: add independently paid high region to C/E."""
import argparse
import hashlib
import importlib.util
import json
import time
from fractions import Fraction as Q
from pathlib import Path


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--old-account', type=Path, required=True)
    p.add_argument('--high-audit', type=Path, required=True)
    p.add_argument('--output-dir', type=Path, required=True)
    a = p.parse_args()
    a.output_dir.mkdir(parents=True, exist_ok=False)
    tic = time.monotonic()
    old = json.loads(a.old_account.read_text())
    high = json.loads(a.high_audit.read_text())
    assert high['all_parameter_events_paid'] and high['all_positive_time_paid']
    b = Q(high['strict_high_norm_upper'])
    assert b == Q(27779670841319, 10**18)
    core = Path(__file__).resolve().parents[2]/'round059/global_region_account/check_account.py'
    spec = importlib.util.spec_from_file_location('accepted_integer_root_enclosure', core)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    rows = []

    def root(x, n, label):
        lower, upper = module.enclosure(x, n)
        rows.append(dict(label=label, radicand=str(x), degree=n,
                         lower=str(lower), upper=str(upper)))
        return lower, upper

    h_i = root(b**6, 5, 'new_full_high_B_to_I')[1]
    paid_i = Q(old['paid_I_conservative_upper'])+h_i
    paid_b = root(paid_i**5, 6, 'disjoint_C_E_H_B')[1]
    target_i_lower = Q(next(row['lower'] for row in old['root_certificates']
                           if row['label'] == 'conditional_target_I'))
    allowance_i = target_i_lower-paid_i
    assert allowance_i > 0
    allowance_b = root(allowance_i**5, 6, 'remaining_D_F_sufficient_B')[0]
    # Independently pay the positive margin in the I slot too.
    paid_fraction = root((paid_b/Q(old['conditional_Bcrit']))**6, 5,
                         'conservative_fraction_of_conditional_I_budget')[1]
    result = dict(status='PASS_THREE_DISJOINT_FULL_TIME_REGIONS_ARITHMETIC',
                  new_configurations=4, root_certificates=rows,
                  new_high_B_upper=str(b), paid_C_E_H_I_upper=str(paid_i),
                  paid_C_E_H_B_upper=str(paid_b),
                  sufficient_remaining_D_F_I_allowance=str(allowance_i),
                  sufficient_remaining_D_F_B_allowance=str(allowance_b),
                  conservative_fraction_of_conditional_I_budget=str(paid_fraction),
                  unpaid=['middle full time: short paid, long still missing',
                          'far full Green integral', 'global g>=1/65000'],
                  global_I=False, global_g=False, new_C6_bound=False,
                  old_scalar_recalculations=0, kernel_calls=0,
                  elapsed_seconds=time.monotonic()-tic,
                  input_sha256={str(f):hashlib.sha256(f.read_bytes()).hexdigest()
                                for f in [a.old_account, a.high_audit, core]})
    (a.output_dir/'results.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(dict(status=result['status'], paid_B_display=float(paid_b),
                         sufficient_remaining_B_display=float(allowance_b),
                         no_global_I_or_g=True)))


if __name__ == '__main__':
    main()
