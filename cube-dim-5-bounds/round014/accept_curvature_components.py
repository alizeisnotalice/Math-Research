"""Primary consumption of frozen 123/128/129/130 evidence; zero kernel calls.

Check archive identities, exact rational signs and source associations. This
does not rerun a successful numerical experiment or replace the analytic proof.
Pass a new output path to preserve the original acceptance receipt.
"""
import argparse
from fractions import Fraction as F
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import tarfile
import time

PROJECT = Path(__file__).resolve().parent.parent


def sha(data):
    return hashlib.sha256(data).hexdigest()


def read(relative):
    return json.loads((PROJECT / relative).read_bytes())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--registration', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists()
    started = time.monotonic()
    reg = json.loads(args.registration.read_bytes())
    for path, expected in reg['input_locks'].items():
        assert sha((PROJECT / path).read_bytes()) == expected, path
    one = read('audit_g123_sol/structure_batch0001.json')
    two = read('audit_g123_sol/structure_batch0002.json')
    for result, n, p, s in ((one, 4096, 307, 3789), (two, 2141, 265, 1876)):
        assert result['status'] == 'new_values_corners_journal_and_partial_tree_qualified'
        assert result['new_points'] == n and result['new_components'] == 5*n
        assert result['new_Delta_checks'] == n and result['integer_corners'] == 16*n
        assert result['passed_this_batch'] == p and result['geometric_splits_this_batch'] == s
        assert len(result['certificates']) == n and not result['global_g110_qualified']
    archive = read_gzip('audit_g123_sol/evidence/numeric_batches.json.gz')
    originals = 0
    for batch in archive['batches'].values():
        for item in batch['files'].values():
            assert 'input_locks' not in item['json']
            value = {**item['json'], 'input_locks': batch['input_locks']}
            raw = (json.dumps(value, sort_keys=True, indent=2) + '\n').encode()
            assert sha(raw) == item['original_sha256'] and len(raw) == item['original_bytes']
            originals += 1
    assert originals == 53
    final = read('audit_curvature130/final_qualification.json')
    assert final['status'] == 'finite11330_new_curvature_component_independently_qualified_partial'
    assert final['source123_final_dependency_closed'] is True
    assert not final['global_g110_qualified'] and final['pending'] == 5427
    packed = {}
    with tarfile.open(fileobj=io.BytesIO((PROJECT / 'audit_curvature130/evidence.tar.gz').read_bytes()), mode='r:gz') as tar:
        for member in tar:
            assert member.isfile() and member.name not in packed
            packed[member.name] = tar.extractfile(member).read()
    for item in final['lossless_archive']['decoded_originals_verified']:
        raw = packed[item['name']]
        assert len(raw) == item['original_bytes'] and sha(raw) == item['sha256']
    candidates = [json.loads(line) for line in packed['structure001/independent_candidates.jsonl'].splitlines()]
    results = [json.loads(line) for line in packed['integers001/integer_results.jsonl'].splitlines()]
    pending = json.loads(packed['integers001/pending_nodes.json'])
    assert len(candidates) == len(results) == 11330
    seen = set()
    passed = []
    failed = []
    for i, (candidate, result) in enumerate(zip(candidates, results)):
        assert candidate['candidate_ordinal'] == result['candidate_ordinal'] == i
        node = tuple(candidate['child_node'])
        assert node == tuple(result['child_node']) and node not in seen
        seen.add(node)
        assert candidate['record_sha256'] == result['source_record_sha256']
        assert candidate['saved_certificate_sha256'] == result['source_certificate_sha256']
        assert candidate['H64_lower'] == '21/200000'
        assert tuple(candidate['direct_parent_node']) == tuple(result['parent_node'])
        assert node[0] == result['parent_node'][0] and node[1][:-1] == result['parent_node'][1]
        q = F(int(result['q_integer_numerator']), int(result['q_integer_denominator']))
        margin = F(int(result['minimum_margin_numerator']), int(result['minimum_margin_denominator']))
        assert 256 <= q < 1024 and result['checked_points'] == 16
        assert result['passed'] == (margin > 0) and result['author_q_margin_bits_status_exactly_match']
        (passed if margin > 0 else failed).append((node, margin))
    assert len(passed) == 5903 and len(failed) == 5427
    assert [list(node) for node, _ in failed] == pending
    result = {'status': 'primary_accepts_123_values_and128_130_partial_curvature_component',
              'pid': os.getpid(), 'seconds': time.monotonic()-started,
              'registration_sha256': sha(args.registration.read_bytes()),
              'independent_new_jet_points': 6237, 'independent_components': 31185,
              'archive_original_jsons_byte_verified': originals,
              'curvature_candidates': 11330, 'closed_corners': 181280,
              'new_passed': len(passed), 'pending': len(failed),
              'minimum_strict_pass_margin': str(min(margin for _, margin in passed)),
              'passed_frontier': 104482, 'splits': 57549, 'roots': 52360,
              'primary_analytic_proof': 'round014/primary_anchor_curvature_proof.md',
              'qualified_global_g': '21/200000', 'global_g110_accepted': False,
              'formal_C5_upper': '32', 'formal_C5_changed': False,
              'new_H_gradient_log_or_old_baseline_calls': 0,
              'lower_research': 'stopped_by_user', 'Astra_dispatch': False,
              'independence_limits': 'shared017fan/081formula/075theory;123 directPF and130 int-gcd paths independently qualified; primary reviewed full math/code',
              'dual_breakthrough_time': None, 'post_dual_rounds': 0,
              'input_locks': reg['input_locks']}
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps({key: result[key] for key in ('status','pid','seconds','new_passed','pending','formal_C5_upper')}))


def read_gzip(relative):
    return json.loads(gzip.decompress((PROJECT / relative).read_bytes()))


if __name__ == '__main__':
    main()
