"""Portable saved-data root enclosure and finite-L1 comparison; no union."""
import argparse
import hashlib
import json
import time
from fractions import Fraction as Q
from pathlib import Path


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--certificate', type=Path, required=True)
    p.add_argument('--old-certificate', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    assert not a.output.exists(), 'Refuse to overwrite a certificate'
    tic = time.monotonic()
    new = json.loads(a.certificate.read_text())
    old = json.loads(a.old_certificate.read_text())
    ul, uh = map(Q, old['root_enclosures']['u'])
    vl, vh = map(Q, old['root_enclosures']['v'])
    assert 0 < ul < uh and ul**6 < 2 < uh**6
    assert 0 < vl < vh and vl**6 < 3 < vh**6
    lower = upper = Q(0)
    for key, coefficient in new['new_ideal_ratio']['coefficients'].items():
        i, j = map(int, key.split(','))
        assert 0 <= i < 6 and 0 <= j < 6
        c = Q(coefficient)
        l, h = ul**i * vl**j, uh**i * vh**j
        lower += c * (l if c >= 0 else h)
        upper += c * (h if c >= 0 else l)
    assert [lower, upper] == list(map(Q, new['new_ideal_ratio']['exact_interval']))
    f = new['finite_L1']
    epsilon = Q(f['epsilon'])
    assert epsilon == Q(1, 10**22) and Q(f['bump_halfwidth']) == epsilon/2
    assert f['all_half_cell_widths_below_one'] and f['max_relative_source_index'] == 2
    assert min(Q(x) for x in new['finite_L1']['minimum_new_box_width']['exact_interval']) > 2*epsilon
    loss = 64*127*12*epsilon
    assert loss == Q(f['union_loss_upper'])
    inward = lower-loss
    assert inward == Q(f['inward_ratio_lower'])
    k = f['receiver_index_radius']
    assert k == 10**16 and f['source_index_radius'] == k+2
    theta = Q(f['threshold'])
    assert theta == 1-Q(1, 10**16)
    boundary = Q(2*k+1, 2*k+5)**6
    assert boundary == Q(f['boundary_factor'])
    finite = theta*inward*boundary
    assert finite == Q(f['weak_ratio_lower'])
    target = Q(2068316441022, 10**12)
    assert target == Q(f['strict_rational_target'])
    old_upper = Q(new['complete_current058_upper'])
    assert old_upper >= Q(old['new_ideal_ratio']['exact_interval'][1])
    assert finite > target > old_upper
    # Geometry is supplied by the new independent audit, not rerun here.
    assert new['new_endpoint_all_127_boxes_valid']
    assert len(new['complete_order_gaps']) == 36 and len(new['effective_order_gaps']) == 21
    assert not new['old_666_inventory_replayed'] and not new['old_22_crossing_checks_replayed']
    assert new['resource']['full_union_calls'] == 0
    result = dict(status='PASS_ROOT_AND_FINITE_L1_EXACT_COMPARISON',
                  ideal_interval=[str(lower), str(upper)],
                  finite_weak_ratio_lower=str(finite), strict_target=str(target),
                  target_strictly_above_complete_old_upper=True,
                  finite_margin=str(finite-target),
                  geometric_coverage_dependency='audit/round064/a2_one_more_event_review/report.md',
                  geometry_recomputed=False, union_calls=0, kernel_calls=0,
                  elapsed_seconds=time.monotonic()-tic,
                  input_sha256={str(x):hashlib.sha256(x.read_bytes()).hexdigest()
                                for x in [a.certificate, a.old_certificate]})
    a.output.write_text(json.dumps(result, indent=2)+'\n')
    print('PASS exact sixth-root, complete ideal value, finite L1, and full old-upper comparisons')


if __name__ == '__main__':
    main()
