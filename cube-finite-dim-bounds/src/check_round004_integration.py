"""Check saved exact comparisons and give an explicit finite smooth-input margin.

This does not rerun discovery or geometric union algorithms. The separate audit
checker does that. All arithmetic below is Fraction/integer arithmetic.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import argparse

ROOT = Path(__file__).resolve().parents[1]


def read(rel):
    return json.loads((ROOT / rel).read_text())


def run():
    inp = read('lower/n3/round004_certificate_input.json')
    author = read('lower/n3/round004_certificate_result.json')
    grid = read('lower/n3/round004_grid_result.json')
    audit = read('audit/006-satellite_n3_results.json')
    originals = {x['dimension']: x for x in read('baseline/original_certificates.json')}
    oldhi = F(originals[3]['lower_bound_interval'][1])
    old4lo = F(originals[4]['lower_bound_interval'][0])
    intervals = [tuple(map(F, x['ratio_interval'])) for x in (author, grid, audit)]
    assert all(lo <= hi and lo > oldhi for lo, hi in intervals)
    lo, hi = max(x[0] for x in intervals), min(x[1] for x in intervals)
    assert lo < hi
    simple = F('1.7939726044040153')
    assert oldhi < simple < lo
    expected = '5f7c4e028a4b3d34d119410c24db9a569f96146427fdb60bd7eefed0e5515970'
    assert hashlib.sha256((ROOT/'lower/n3/round004_certificate_input.json').read_bytes()).hexdigest() == expected
    assert all(x['input_sha256'] == expected for x in (author, grid, audit))
    assert len(inp['witness_keys']) == 170 and author['b'] == audit['b'] == 2
    periods = list(map(F, inp['periods']))
    mass = sum(map(F, inp['weights']))
    assert mass == F(500280708211, 500000000000)
    N, tau = 10**6, F(1, 10**12)
    # Shrinking each of 170 clipped boxes by tau on all faces loses at most
    # 2*tau*sum(pairwise period products) per box, by product telescoping.
    loss = 2*len(inp['witness_keys'])*tau*sum(periods[j]*periods[k] for j in range(3) for k in range(j+1, 3))
    finite = (intervals[0][0] - loss/mass)*F(2*N+1, 2*N+5)**3
    assert finite > oldhi
    public3 = F(read('results/bounds_ledger.json')['tasks'][0]['public'])
    assert simple > public3
    n4 = read('lower/n4/round004_certificate.json')
    # Its field schema is recorded rather than inferring a whole-input bound.
    n4_upper = F('1902811183099396243470782409473350989894373508837/1001000000000000000000000000000000000000000000000')
    assert n4_upper < old4lo
    assert str(n4_upper) in json.dumps(n4)
    return {
        'date': '2026-10-25',
        'status': 'PASS: stored exact comparisons and finite smooth-input margin',
        'scope': 'No discovery, old-baseline rerun, or geometric union rerun in this script.',
        'input_sha256': expected,
        'combined_new_T3_interval': [str(lo), str(hi)],
        'simple_rational_lower': str(simple),
        'gain_above_original_upper_lower': str(lo-oldhi),
        'above_public': True,
        'above_original_internal': True,
        'n4_selected_finite_union_upper_below_original_lower': True,
        'finite_smooth_input': {'N': N, 'tau': str(tau), 'eta': str(tau/4),
            'cell_volume_loss_upper': str(loss), 'ratio_lower': str(finite),
            'decimal_display_only': float(finite), 'above_original_upper': True},
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--save', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.save:
        (ROOT/'results/round004_integration.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))
