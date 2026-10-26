"""Prepared new q / 181280 closed-corner audit. NOT YET AUTHORIZED TO RUN.

Independent Python integer pairs + common-denominator sign tests. No Fraction,
author code imports, H, gradient, log, Arb, kernel evaluation, or new splitting.
Successful arithmetic remains provisional until source123 numerical receipt is
locked in a separate final qualification; this entry point cannot claim globalg.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import time

from integer_core import (add, box, digest, div, le, maximum, minimum, mul,
                          neg, node_key, parse, rat, roots, sha, sub,
                          sum_cleared, text)

HERE = Path(__file__).resolve().parent


def exclusive(path,obj):
    with Path(path).open('x') as f:
        json.dump(obj,f,ensure_ascii=False,sort_keys=True,separators=(',',':'))
        f.write('\n')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--reviewed-run-plan-sha',required=True)
    ap.add_argument('--run-plan',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args = ap.parse_args()
    started = time.monotonic()
    assert sha(args.run_plan) == args.reviewed_run_plan_sha
    plan = json.loads(args.run_plan.read_bytes())
    assert plan['new_q_and_corners_authorized'] is True
    assert plan['finite_children'] == 11330 and plan['finite_corners'] == 181280
    for name,expected in plan['code_locks'].items(): assert sha(HERE/name) == expected
    for item in plan['input_locks'].values(): assert sha(item['path']) == item['sha256']
    locks = plan['input_locks']
    receipt = json.loads(Path(locks['structure_receipt']['path']).read_bytes())
    assert receipt['status'] == 'structure_and_ancestry_checked_no_q_or_corner_test'
    assert receipt['children_once_each'] == 11330
    assert receipt['candidate_packet_sha256'] == locks['independent_candidates']['sha256']
    state = json.loads(Path(locks['state123']['path']).read_bytes())
    assert receipt['source_state123_sha256'] == locks['state123']['sha256']
    assert len(state['queue']) == 11330 and len(state['new_splits']) == 5665
    assert not state['faults'] and not state['inflight'] and not state['obstructions_H128']
    assert not args.output.exists()
    args.output.mkdir()
    exclusive(args.output/'run_registration.json',{'task':'audit_curvature130_integer_corners',
        'pid':os.getpid(),'registered_utc':datetime.now(timezone.utc).isoformat(),
        'reviewed_run_plan_sha256':args.reviewed_run_plan_sha,'workers':1,
        'finite_children':11330,'finite_corners':181280,'new_H_gradient_log_kernel_calls':0,
        'new_splits':0,'qualification':'new_curvature_arithmetic_only;source123_final_numeric_receipt_separate'})
    author = {}
    with Path(locks['author_results']['path']).open() as f:
        for line in f:
            row = json.loads(line)
            key = node_key(row['child_node'])
            assert key not in author
            author[key] = row
    queue = set(map(node_key,state['queue']))
    assert len(author) == 11330 and set(author) == queue
    seen = set()
    passed = 0
    pending_nodes = []
    qbounds = None
    corners = 0
    R = roots()
    cache = state['new_cache_by_epsilon']['1/128']
    with Path(locks['independent_candidates']['path']).open() as source, (args.output/'integer_results.jsonl').open('x') as output:
        for ordinal,line in enumerate(source):
            item = json.loads(line)
            node = node_key(item['child_node'])
            assert item['candidate_ordinal'] == ordinal and node in queue and node not in seen
            seen.add(node)
            split = state['new_splits'][item['parent_split_index']]
            assert split['node'] == item['direct_parent_node']
            assert split['children'][item['side']] == item['child_node']
            cert = split['failed_certificate']
            rec = cache[item['point_key']]
            assert digest(rec) == item['record_sha256'] and digest(cert) == item['saved_certificate_sha256']
            assert cert['point_key'] == rec['pair_cache_key'][1] == item['point_key']
            assert rec['point'] == item['point_key'].split(',')
            assert rec['model'] == 'n5 alpha1 body16 unnormalised inverse Fourier(2pi)^-5'
            assert rec['kernel_sha256'] == '128d98c94e55632336edafe024aa3a27c1486898670a52be89cb39cc22ed6cc0'
            assert rec['epsilon'] == cert['leaf_epsilon'] == '1/128'
            assert rec['qualification_target'] == cert['target'] == '11/100000'
            assert item['H64_lower'] == '21/200000'
            hl,hu = map(parse,rec['values'][0])
            assert hl[0] > 0 and le(hl,hu) and le((21,200000),hu)
            # q=(1024*200000*U_num - 768*21*U_den)/(200000*U_num).
            # This explicit denominator formula is independent of author Fraction arithmetic.
            q = rat(1024*200000*hu[0]-768*21*hu[1],200000*hu[0])
            assert le((256,1),q) and not le((1024,1),q)
            if qbounds is None: qbounds = [q,q]
            else: qbounds = [minimum([qbounds[0],q]),maximum([qbounds[1],q])]
            gradient = [tuple(map(parse,pair)) for pair in cert['loggradient_directed']]
            assert len(gradient) == 4
            for (lo,hi),raw in zip(gradient,rec['values'][1:]):
                dl,du = map(parse,raw)
                assert le(dl,du) and le(lo,hi)
                ratios = [div(dl,hl),div(dl,hu),div(du,hl),div(du,hu)]
                assert le(lo,minimum(ratios)) and le(maximum(ratios),hi)
            delta = parse(cert['Delta_lower'])
            anchor = list(map(parse,rec['point'][1:]))
            bounds = box(node,R)
            best = None
            worst = None
            for mask in range(16):
                bits = [(mask>>(3-j)) & 1 for j in range(4)]
                disp = [sub(bounds[j][bits[j]],anchor[j]) for j in range(4)]
                terms = [delta]
                for u,(gl,gu) in zip(disp,gradient):
                    terms.append(mul(gl if u[0] >= 0 else gu,u))
                    terms.append(neg(mul(q,mul(u,u))))
                margin = sum_cleared(terms)
                if best is None or not le(best,margin): best,worst = margin,bits
                corners += 1
            ok = best[0] > 0
            passed += ok
            if not ok: pending_nodes.append(item['child_node'])
            a = author[node]
            for field in ('candidate_ordinal','child_node','direct_parent_node','point_key',
                          'record_sha256','saved_certificate_sha256','H64_lower_source_node'):
                assert a[field] == item[field],(ordinal,field)
            assert a['epsilon'] == '1/128' and a['target'] == '11/100000'
            assert a['H64_lower'] == '21/200000' and parse(a['H128_upper']) == hu
            assert parse(a['penalty_coefficient']) == q
            assert parse(a['minimum_log_margin']) == best
            assert a['worst_corner_bits'] == worst and a['passed'] == ok
            assert a['checked_points'] == 16 and a['Delta_loggradient_recomputed'] is False
            result = {'candidate_ordinal':ordinal,'child_node':item['child_node'],
                'parent_node':item['direct_parent_node'],'source_record_sha256':item['record_sha256'],
                'source_certificate_sha256':item['saved_certificate_sha256'],
                'q_integer_numerator':str(q[0]),'q_integer_denominator':str(q[1]),
                'minimum_margin_numerator':str(best[0]),'minimum_margin_denominator':str(best[1]),
                'worst_corner_bits':worst,'checked_points':16,'passed':ok,
                'author_q_margin_bits_status_exactly_match':True}
            output.write(json.dumps(result,sort_keys=True,separators=(',',':'))+'\n')
            if (ordinal+1)%256 == 0:
                output.flush()
                os.fsync(output.fileno())
    assert seen == queue and len(seen) == 11330 and corners == 181280
    exclusive(args.output/'pending_nodes.json',pending_nodes)
    summary = {'task':'audit_curvature130_integer_corners',
        'status':'all_new_q_and181280_corners_compared_source123_final_qualification_still_required',
        'pid':os.getpid(),'tested_once_each':len(seen),'checked_integer_closed_corners':corners,
        'passed':passed,'pending':len(pending_nodes),'retained_old_and123_passed':98579,
        'prospective_passed_frontier':98579+passed,'frontier_splits':57549,'frontier_roots':52360,
        'frontier_geometry_unchanged':True,'prospective_complete_frontier':not pending_nodes,
        'penalty_interval':[text(x) for x in qbounds],
        'new_H_gradient_log_kernel_calls':0,'new_splits':0,
        'author_outputs_matched_all':True,'integer_results_sha256':sha(args.output/'integer_results.jsonl'),
        'source_numeric123_final_receipt_required':True,'global_g110_qualified':False,
        'Cn_C5_I_updated':False,'wall_seconds':time.monotonic()-started}
    exclusive(args.output/'integer_summary.json',summary)
    print(json.dumps(summary,ensure_ascii=False))


if __name__ == '__main__': main()
