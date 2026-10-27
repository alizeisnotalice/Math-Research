"""Prepared168 fixed actual158 audit; no q/log or old corner evaluation."""
import argparse
from datetime import datetime,timezone
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path
import queue
import time
import traceback
from integer168 import one_source,initialize,accept_test
from rational168 import Q
from sources168 import sha,digest,read,source_gates

HERE = Path(__file__).resolve().parent


def put(path,obj):
    with Path(path).open('x') as handle:
        json.dump(obj,handle,sort_keys=True,separators=(',',':'));handle.write('\n')


def worker(spec,binding,log_upper_text,destination,result_queue):
    dest = Path(destination);dest.mkdir();start = time.monotonic();evaluated = skipped = 0
    numeric = hashlib.sha256()
    put(dest/'registration.json',{'pid':os.getpid(),'shard_id':spec['shard_id'],
        'fixed_pairs':spec['pairs'],'receivers':spec['children'],'new_q_log_H_gradient_calls':0})
    try:
        data = read(spec['input_path']);sid = spec['shard_id']
        folder = Path(binding['output_directory'])/f'shard{sid:02d}'
        reg,cp = read(folder/'registration.json'),read(folder/'checkpoint.json')
        fp = reg['fingerprint'];assert cp['status'] == 'complete' and cp['fingerprint'] == fp
        assert fp['shard_id'] == sid and fp['input_sha256'] == spec['input_sha256']
        assert fp['reviewed_manifest_sha256'] == binding['prepared158_manifest_sha256']
        assert fp['log155_sha256'] == binding['existing155_log_sha256']
        assert fp['worker_code_sha256'] == binding['worker158_sha256']
        upper = Q(log_upper_text);best = {k:initialize(c) for k,c in data['children'].items()}
        for value in best.values():assert Q(value['best_margin115']) == Q(value['best_margin110'])-upper
        previous = '0'*64;count = 0;negative110 = 0
        with (folder/'pair_results.jsonl').open() as handle:
            for pair,line in zip(data['pairs'],handle,strict=True):
                assert line.endswith('\n');event = json.loads(line);h = event.pop('hash')
                assert event['seq'] == count+1 and event['previous'] == previous and digest(event) == h
                previous = h;count += 1;ci = str(pair['child_index'])
                assert event['pair_index'] == pair['pair_index'] and event['child_index'] == pair['child_index']
                assert event['canonical_id'] == pair['canonical_id']
                if event['kind'] == 'skip_after_first115PASS':
                    assert best[ci]['first_pass115'] is not None;skipped += 1
                    continue
                assert event['kind'] == 'evaluated_single_source' and best[ci]['first_pass115'] is None
                source = data['sources'][pair['canonical_id']];child = data['children'][ci]
                margin115,margin110,bits = one_source(source,child,upper)
                certificate = event['certificate']
                for field,value in [('child_index',pair['child_index']),('leaf_index',pair['leaf_index']),
                    ('node',pair['node']),('canonical_id',pair['canonical_id']),('new_core_sha256',pair['new_core_sha256']),
                    ('epsilon',source['core']['epsilon']),('point',source['core']['point']),('q_new',source['core']['q_new'])]:
                    assert certificate[field] == value
                assert Q(certificate['margin115']) == margin115 and Q(certificate['margin110']) == margin110
                assert certificate['worst_corner_bits'] == bits and certificate['checked_points'] == 16
                assert certificate['single_anchor_all16'] and certificate['new_q_log_H_gradient_calls'] == 0
                assert certificate['passed115'] == (margin115 > Q(0))
                # All tested cases, including negative110 or failed115, update max.
                accept_test(best[ci],certificate);evaluated += 1;negative110 += margin110 < Q(0)
                numeric.update(json.dumps({'pair_index':pair['pair_index'],'margin115':str(margin115),
                    'margin110':str(margin110),'bits':bits},sort_keys=True,separators=(',',':')).encode())
        assert count == spec['pairs'] == cp['next_pair'] and previous == cp['last_hash']
        assert evaluated == cp['evaluated'] and skipped == cp['skipped'] and digest(best) == cp['best_sha256']
        assert best == read(folder/'receiver_best.json')
        summary = read(folder/'summary.json')
        assert summary['fingerprint'] == fp and summary['status'] == 'complete'
        assert summary['evaluated'] == evaluated and summary['skipped'] == skipped
        assert summary['passed115_receivers'] == sum(b['first_pass115'] is not None for b in best.values())
        put(dest/'receiver_best.json',best)
        result = {'status':'complete','pid':os.getpid(),'shard_id':sid,'metadata_pairs':count,
            'new_16corner_cases':evaluated,'new_corners':16*evaluated,'skipped':skipped,
            'negative110_source_cases_allowed_and_matched':negative110,'receivers':len(best),
            'new115_passes':sum(b['first_pass115'] is not None for b in best.values()),
            'best_sha256':sha(dest/'receiver_best.json'),'numeric_digest':numeric.hexdigest(),
            'seconds':time.monotonic()-start,'new_q_log_or_old_corner_calls':0}
        put(dest/'result.json',result);result_queue.put(result)
    except BaseException as exc:
        result = {'status':'failed_preserved','pid':os.getpid(),'shard_id':spec['shard_id'],
                  'completed_new_cases':evaluated,'skipped':skipped,'numeric_digest':numeric.hexdigest(),
                  'exception':repr(exc),'traceback':traceback.format_exc(),'automatic_retry':False}
        put(dest/'fault.json',result);result_queue.put(result)
        raise


def main():
    ap = argparse.ArgumentParser();ap.add_argument('--reviewed-plan-sha',required=True)
    ap.add_argument('--binding',type=Path,required=True);ap.add_argument('--reviewed-binding-sha',required=True)
    ap.add_argument('--output',type=Path,required=True);args = ap.parse_args()
    assert sha(HERE/'plan.json') == args.reviewed_plan_sha and sha(args.binding) == args.reviewed_binding_sha
    plan,binding = read(HERE/'plan.json'),read(args.binding)
    assert binding['execution_ready'] is True and binding['plan_sha256'] == args.reviewed_plan_sha
    for name,h in plan['code_locks'].items():assert sha(HERE/name) == h
    for path,h in plan['input_locks'].items():assert sha(path) == h,path
    for path,h in binding['candidate_locks'].items():
        assert plan['input_locks'].get(path,h) == h and sha(path) == h,path
    folder = Path(binding['output_directory'])
    required = [folder/name for name in ['complete_domain_candidate.json','receiver_best_updates.json','target115_pending.json']]
    required.append(Path(binding['producer_summary_path']))
    for sid in range(15):
        required += [folder/f'shard{sid:02d}'/name for name in ['registration.json','checkpoint.json','pair_results.jsonl','receiver_best.json','summary.json']]
    assert all(str(path) in binding['candidate_locks'] for path in required)
    assert args.output.is_absolute() and args.output.is_relative_to(HERE) and not args.output.exists()
    args.output.mkdir();start = time.monotonic()
    put(args.output/'registration.json',{'pid':os.getpid(),'utc':datetime.now(timezone.utc).isoformat(),
        'plan_sha256':args.reviewed_plan_sha,'binding_sha256':args.reviewed_binding_sha,'workers':15,
        'metadata_pairs':21494,'new_actual_16corner_cases':binding['evaluated_pairs'],
        'new_corners':16*binding['evaluated_pairs'],'new_q_log_H_gradient_calls':0})
    fixed,children,old_rows,upper = source_gates(plan)
    assert binding['existing155_log_sha256'] == fixed['existing155_log']['sha256']
    put(args.output/'source_gate.json',{'status':'all4769_sources_copied_from_qualified155_and_fixed_pairs_closed_receivers_bound',
        'sources':4769,'receivers':1240,'pairs':21494,'old_complete_leaves':111201,'old_log_inherited_hash_only':True})
    ctx = mp.get_context('spawn');rq = ctx.Queue();processes = []
    for spec in fixed['shards']:
        p = ctx.Process(target=worker,args=(spec,binding,upper,str(args.output/f'shard{spec["shard_id"]:02d}'),rq))
        p.start();processes.append(p)
    results = {}
    while len(results) < 15:
        try:
            result = rq.get(timeout=.2);assert result['shard_id'] not in results;results[result['shard_id']] = result
        except queue.Empty:
            assert not any(p.exitcode not in (None,0) for p in processes),'worker fault retained;no automatic retry'
    for p in processes:p.join()
    rq.close();rq.join_thread()
    assert all(p.exitcode == 0 for p in processes) and all(r['status'] == 'complete' for r in results.values())
    evaluated = sum(r['new_16corner_cases'] for r in results.values());skipped = sum(r['skipped'] for r in results.values())
    assert evaluated == binding['evaluated_pairs'] and skipped == binding['skipped_pairs'] and evaluated+skipped == 21494
    receivers = {};by_leaf = {}
    for i in range(15):
        path = args.output/f'shard{i:02d}/receiver_best.json';assert sha(path) == results[i]['best_sha256']
        for text,best in read(path).items():
            ci = int(text);assert ci not in receivers;receivers[ci] = best;by_leaf[children[ci]['leaf_index']] = (children[ci],best)
    assert set(receivers) == set(range(1240)) and len(by_leaf) == 1240
    minimum = None;old_minimum = None;newpasses = 0;pending = [];best_updates = [];unchanged = 0
    improved_receivers = improved_without_pass115 = 0
    for index in range(111201):
        row = old_rows[index];old = Q(row['margin110_new']);margin = old
        old_minimum = old if old_minimum is None else min(old_minimum,old)
        witness = {'kind':'unchanged155','canonical_id':row['canonical_id']}
        if index in by_leaf:
            child,best = by_leaf[index]
            assert child['node'] == row['node'] and Q(best['best_margin110']) >= old
            assert Q(best['best_margin115']) == Q(best['best_margin110'])-Q(upper)
            margin = Q(best['best_margin110']);witness = {'kind':best['best_source_kind'],'canonical_id':best['best_canonical_id']}
            improved_receivers += margin > old
            improved_without_pass115 += margin > old and best['first_pass115'] is None
            best_updates.append({'leaf_index':index,'node':row['node'],'best':best})
            newpasses += best['first_pass115'] is not None
            if best['first_pass115'] is None:
                pending.append({'leaf_index':index,'node':row['node'],'best_margin110':best['best_margin110'],
                    'best_margin115':best['best_margin115'],'best_canonical_id':best['best_canonical_id']})
        else:assert row['passed115'];unchanged += 1
        key = (margin,index)
        if minimum is None or key < minimum[0]:minimum = (key,row['node'],witness)
    assert unchanged == 109961 and len(pending)+newpasses == 1240 and minimum[0][0] >= old_minimum > Q(0)
    assert best_updates == read(folder/'receiver_best_updates.json') and pending == read(folder/'target115_pending.json')
    candidate = read(folder/'complete_domain_candidate.json')
    assert candidate['full111201_closed_domain_retained'] and candidate['unchanged155_leaves'] == 109961
    assert candidate['updated_receivers'] == 1240 and Q(candidate['minimum_margin110_updated']) == minimum[0][0]
    assert Q(candidate['old155_minimum_margin110']) == old_minimum
    assert candidate['witness_leaf_index'] == minimum[0][1] and candidate['witness_node'] == minimum[1] and candidate['witness_source'] == minimum[2]
    assert Q(candidate['minimum_margin115']) == minimum[0][0]-Q(upper)
    assert Q(candidate['candidate_linear_floor']) == Q(11,100000)*(1+minimum[0][0])
    assert candidate['full115_passed'] == 109961+newpasses and candidate['target115_pending'] == len(pending)
    producer = read(binding['producer_summary_path'])
    assert producer['status'] == 'complete_fixed21494_bestmargin_candidate' and producer['fault'] is None
    assert producer['master_pid'] == binding['producer_master_pid']
    assert len(producer['worker_evidence']) == 15 and all(p['exitcode'] == 0 and not p['alive'] for p in producer['worker_evidence'])
    assert producer['evaluated_pairs'] == evaluated and producer['skipped_after_firstPASS'] == skipped
    assert producer['new115_passes'] == newpasses and producer['remaining115'] == len(pending)
    assert producer['new_q_log_H_gradient_or_split_calls'] == 0 and not producer['formal_g_C5_update']
    for path,h in {**plan['input_locks'],**binding['candidate_locks']}.items():assert sha(path) == h,path
    put(args.output/'result.json',{'status':'all_actual158_corners_firstpass_best_and_complete_domain_independently_qualified',
        'pid':os.getpid(),'worker_pids':[p.pid for p in processes],'evaluated_pairs':evaluated,'new_corners':16*evaluated,
        'skip_metadata_pairs':skipped,'receivers':1240,'bound_sources':4769,'unchanged155_leaves':109961,
        'new115_passes':newpasses,'pending115':len(pending),'full115_passed':109961+newpasses,
        'receiver_margin110_improvements':improved_receivers,'improved_without_PASS115':improved_without_pass115,
        'minimum_margin110_updated':str(minimum[0][0]),'witness_leaf_index':minimum[0][1],
        'witness_node':minimum[1],'witness_source':minimum[2],'complete_old111201_domain_preserved':True,
        'all_source_hashes_before_after_unchanged':True,'new_q_log_H_or_old_corner_calls':0,
        'P8_or_Cn_calls':0,'seconds':time.monotonic()-start})


if __name__ == '__main__':main()
