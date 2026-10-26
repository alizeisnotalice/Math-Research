"""Read saved exact coefficients/roots and finite-L1 arithmetic; no geometry/union."""
import argparse
import hashlib
import json
import time
from fractions import Fraction as Q
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--round', type=int, choices=(73, 75), required=True)
    parser.add_argument('--certificate', type=Path, required=True)
    parser.add_argument('--old-certificate', type=Path, required=True)
    parser.add_argument('--root-certificate', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists(), 'Refuse to overwrite evidence'
    start = time.monotonic()
    new = json.loads(args.certificate.read_text())
    old = json.loads(args.old_certificate.read_text())
    roots = json.loads(args.root_certificate.read_text())
    assert Q(str(new['root_denominator'])) == Q(str(roots['root_denominator']))
    ul, uh = map(Q, roots['root_enclosures']['u'])
    vl, vh = map(Q, roots['root_enclosures']['v'])
    assert 0 < ul < uh and ul**6 < 2 < uh**6
    assert 0 < vl < vh and vl**6 < 3 < vh**6
    lower = upper = Q(0)
    for key, value in new['new_ideal_ratio']['coefficients'].items():
        i, j = map(int, key.split(','))
        assert 0 <= i < 6 and 0 <= j < 6
        coefficient = Q(value)
        lo, hi = ul**i * vl**j, uh**i * vh**j
        lower += coefficient * (lo if coefficient >= 0 else hi)
        upper += coefficient * (hi if coefficient >= 0 else lo)
    assert [lower, upper] == list(map(Q, new['new_ideal_ratio']['exact_interval']))
    finite = new['finite_L1']
    epsilon = Q(finite['epsilon'])
    assert epsilon == Q(1, 10**22)
    assert Q(finite['bump_halfwidth']) == epsilon / 2
    assert min(map(Q, finite['minimum_new_box_width']['exact_interval'])) > 2*epsilon
    assert finite['max_relative_source_index'] == 2
    width_key = 'all_box_widths_below_one'
    assert finite[width_key]
    loss = 64*127*12*epsilon
    assert loss == Q(finite['union_loss_upper'])
    inward = lower-loss
    assert inward == Q(finite['inward_ratio_lower'])
    k = finite['receiver_index_radius']
    assert k == 10**16 and finite['source_index_radius'] == k+2
    theta = Q(finite['threshold'])
    assert theta == 1-Q(1, 10**16)
    boundary = Q(2*k+1, 2*k+5)**6
    assert boundary == Q(finite['boundary_factor'])
    ratio = theta*inward*boundary
    assert ratio == Q(finite['weak_ratio_lower'])
    previous = 70 if args.round == 73 else 73
    target = Q(2072477059470 if args.round == 73 else 2073015624830, 10**12)
    assert target == Q(finite['strict_rational_target'])
    old_upper = Q(new[f'complete_current{previous:03d}_upper'])
    assert old_upper >= Q(old['new_ideal_ratio']['exact_interval'][1])
    assert ratio > target > old_upper
    assert new['new_endpoint_all_127_boxes_valid']
    assert new['all_127_minmax_choices_unchanged']
    assert len(new['six_disjoint_projection_pairs']) == 6
    assert len(new['new_complete_order_gaps']) == 36
    assert len(new['new_effective_gaps']) == 21
    proof = ('audit/round073/post_s7_mask_review/report.md' if args.round == 73
             else 'audit/round075/post_s8_overlap_review/report.md')
    assert new['resource']['full_union_calls'] == 0
    result = dict(status='PASS_EXACT_ROOT_VALUE_FINITE_L1_AND_FULL_PREVIOUS_UPPER',
                  round=args.round, ideal_interval=[str(lower), str(upper)],
                  finite_weak_ratio_lower=str(ratio), strict_target=str(target),
                  finite_margin=str(ratio-target), old_complete_upper=str(old_upper),
                  geometry_proof_dependency=proof, geometry_recomputed=False,
                  union_calls=0, kernel_calls=0, elapsed_seconds=time.monotonic()-start,
                  input_sha256={str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in (args.certificate, args.old_certificate, args.root_certificate)})
    with args.output.open('x') as output:
        output.write(json.dumps(result, indent=2)+'\n')
    print(f'PASS {args.round:03d}: exact roots, full coefficients, finite L1 and complete previous upper')


if __name__ == '__main__':
    main()
