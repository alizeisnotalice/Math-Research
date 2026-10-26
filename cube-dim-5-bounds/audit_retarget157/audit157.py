"""Prepared finite15-shard independent155 audit. No producer imports."""
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
from integer157 import coefficient,corners,public_log,roots
from rational157 import Q
from sources157 import digest,load,sha,source_gates

HERE = Path(__file__).resolve().parent


def write(path,obj):
    with Path(path).open('x') as handle:
        json.dump(obj,handle,sort_keys=True,separators=(',',':'));handle.write('\n')


def journal(path):
    previous = '0'*64
    with Path(path).open() as handle:
        for seq,line in enumerate(handle,1):
            assert line.endswith('\n')
            event = json.loads(line);h = event.pop('hash')
            assert event['seq'] == seq and event['previous'] == previous and digest(event) == h
            previous = h
            yield event['body'],previous


def worker(spec,binding,log_upper_text,output,result_queue):
    output = Path(output);output.mkdir()
    write(output/'registration.json',{'pid':os.getpid(),'shard_id':spec['shard_id'],
          'objects':spec['objects'],'leaves':spec['leaves'],'new_H_log_gradient_calls':0})
    count_q = count_leaf = 0;numeric = hashlib.sha256();start = time.monotonic()
    try:
        data = load(spec['input_path']);sid = spec['shard_id'];folder = Path(binding['output_directory'])/f'shard{sid:02d}'
        cp,reg = load(folder/'checkpoint.json'),load(folder/'registration.json')
        assert cp['status'] == 'complete' and cp['next_q'] == spec['objects'] and cp['next_leaf'] == spec['leaves']
        assert cp['fingerprint'] == reg['fingerprint']
        assert cp['fingerprint']['shard_id'] == sid and cp['fingerprint']['input_sha256'] == spec['input_sha256']
        assert cp['fingerprint']['reviewed_manifest_sha256'] == binding['prepared155_manifest_sha256']
        assert cp['fingerprint']['worker_sha256'] == binding['worker155_sha256']
        assert cp['common_log_sha256'] == binding['common_log_sha256']
        upper = Q(log_upper_text);parameters = {};below256 = 0;last_q = '0'*64
        for source,(row,last_q) in zip(data['sources'],journal(folder/'q_updates.jsonl'),strict=True):
            c = source['old_core'];cid = source['canonical_id'];q,caps = coefficient(source)
            assert row['canonical_id'] == cid and cid not in parameters
            for field,value in [('epsilon',c['epsilon']),('point',c['point']),('target110','11/100000'),
                                ('target115','23/200000'),('Delta110_lower',c['Delta110_lower']),
                                ('loggradient_directed',c['loggradient_directed']),
                                ('mass_qualification',source['mass_qualification']),
                                ('raw_record_sha256',source['raw_record_sha256'])]:
                assert row[field] == value
            assert Q(row['q_old']) == Q(c['penalty_coefficient']) and Q(row['q_new']) == q
            assert Q(row['Delta115_lower']) == Q(c['Delta110_lower'])-upper
            assert row['q_new_rule'] == '153_time_shape_uniform_or_original_g105_mass'
            assert row['source_penalty_sha256'] == binding['penalty153_sha256']
            for key,value in caps.items():
                assert Q(row['caps153'][key]) == value
            parameters[cid] = (c,q);below256 += q < Q(256);count_q += 1
            numeric.update(json.dumps({'canonical_id':cid,'q':str(q)},sort_keys=True).encode())
        assert count_q == spec['objects'] and last_q == cp['q_hash']
        table = roots();minimum = None;passed = 0;pending = [];last_leaf = '0'*64
        for leaf,(row,last_leaf) in zip(data['leaves'],journal(folder/'leaf_results.jsonl'),strict=True):
            c,q = parameters[leaf['canonical_id']]
            m,shifted,bits = corners(c,q,leaf,upper,table)
            assert row['leaf_index'] == leaf['leaf_index'] and row['node'] == leaf['node'] and row['canonical_id'] == leaf['canonical_id']
            assert row['epsilon'] == leaf['epsilon'] == c['epsilon'] and row['checked_points'] == 16
            assert row['q_new_rule'] == '153_time_shape_uniform_or_original_g105_mass'
            assert Q(row['q_new']) == q and Q(row['q_old']) == Q(c['penalty_coefficient'])
            assert Q(row['margin110_new']) == m and Q(row['margin115']) == shifted
            assert row['saved_old_margin110'] == leaf['saved_g110_margin'] and row['worst_corner_bits'] == bits
            assert row['margin110_new_ge_old_positive'] and row['full110_leaf_preserved'] and row['one_anchor_all16corners']
            assert row['passed115'] == (shifted > Q(0))
            key = (m,leaf['leaf_index']);minimum = key if minimum is None else min(minimum,key)
            passed += shifted > Q(0);count_leaf += 1
            if not shifted > Q(0):
                pending.append({'leaf_index':leaf['leaf_index'],'node':leaf['node'],'canonical_id':leaf['canonical_id'],
                                'margin115':str(shifted),'margin110_new':str(m)})
            numeric.update(json.dumps({'leaf_index':leaf['leaf_index'],'margin110':str(m),
                                       'margin115':str(shifted),'bits':bits},sort_keys=True).encode())
        assert count_leaf == spec['leaves'] and last_leaf == cp['leaf_hash']
        assert cp['passed115'] == passed and cp['minimum110'] == [str(minimum[0]),minimum[1]]
        result = {'status':'complete','pid':os.getpid(),'shard_id':sid,'objects':count_q,'leaves':count_leaf,
                  'corners':16*count_leaf,'q_below256':below256,'passed115':passed,'pending115':pending,
                  'minimum110':[str(minimum[0]),minimum[1]],'numeric_digest':numeric.hexdigest(),
                  'seconds':time.monotonic()-start,'H_gradient_or_old_corner_calls':0}
        write(output/'result.json',result);result_queue.put(result)
    except BaseException as exc:
        result = {'status':'failed_preserved','pid':os.getpid(),'shard_id':spec['shard_id'],
                  'completed_q':count_q,'completed_leaves':count_leaf,'numeric_digest':numeric.hexdigest(),
                  'exception':repr(exc),'traceback':traceback.format_exc(),'automatic_retry':False}
        write(output/'fault.json',result);result_queue.put(result)
        raise


def main():
    ap = argparse.ArgumentParser();ap.add_argument('--reviewed-plan-sha',required=True)
    ap.add_argument('--binding',type=Path,required=True);ap.add_argument('--reviewed-binding-sha',required=True)
    ap.add_argument('--output',type=Path,required=True);args = ap.parse_args()
    assert sha(HERE/'plan.json') == args.reviewed_plan_sha and sha(args.binding) == args.reviewed_binding_sha
    plan,binding = load(HERE/'plan.json'),load(args.binding)
    assert binding['execution_ready'] is True and binding['plan_sha256'] == args.reviewed_plan_sha
    assert args.output.is_absolute() and args.output.is_relative_to(HERE) and not args.output.exists()
    for name,h in plan['code_locks'].items():
        assert sha(HERE/name) == h
    for path,h in plan['input_locks'].items():
        assert sha(path) == h,path
    for path,h in binding['candidate_locks'].items():
        assert plan['input_locks'].get(path,h) == h and sha(path) == h,path
    folder = Path(binding['output_directory'])
    required = [folder/'common_log_certificate.json',folder/'math_start.json',folder/'full_domain_candidate.json',folder/'target115_pending.json',Path(binding['producer_summary_path'])]
    for sid in range(15):
        required += [folder/f'shard{sid:02d}'/name for name in ['q_updates.jsonl','leaf_results.jsonl','checkpoint.json','registration.json','summary.json']]
    assert all(str(path) in binding['candidate_locks'] for path in required)
    args.output.mkdir();write(args.output/'registration.json',{
        'pid':os.getpid(),'utc':datetime.now(timezone.utc).isoformat(),'plan_sha256':args.reviewed_plan_sha,
        'binding_sha256':args.reviewed_binding_sha,'maximum_q':107530,'maximum_leaves':111201,
        'maximum_corners':1779216,'new_public_log_objects':1,'workers':15,'old_H_gradient_or_corner_calls':0})
    start = time.monotonic();specs,source_result = source_gates(plan)
    write(args.output/'source_gates.json',source_result)
    common = load(folder/'common_log_certificate.json');upper = public_log(common)
    assert sha(folder/'common_log_certificate.json') == binding['common_log_sha256']
    mathstart = load(folder/'math_start.json')
    assert mathstart['common_log_sha256'] == binding['common_log_sha256'] and mathstart['reviewed_manifest_sha256'] == binding['prepared155_manifest_sha256']
    write(args.output/'public_log_check.json',{'status':'positive_polynomial_integral_and_positive_tail_bound_passed',
          'upper':str(upper),'terms':64,'new_public_constants':1})
    ctx = mp.get_context('spawn');rq = ctx.Queue();processes = []
    for spec in specs:
        process = ctx.Process(target=worker,args=(spec,binding,str(upper),str(args.output/f'shard{spec["shard_id"]:02d}'),rq))
        process.start();processes.append(process)
    results = {}
    while len(results) < 15:
        try:
            result = rq.get(timeout=.2);sid = result['shard_id'];assert sid not in results;results[sid] = result
        except queue.Empty:
            assert not any(p.exitcode not in (None,0) for p in processes),'worker failure retained;no automatic retry'
    for process in processes:
        process.join()
    rq.close();rq.join_thread()
    assert all(p.exitcode == 0 for p in processes) and all(r['status'] == 'complete' for r in results.values())
    assert sum(r['objects'] for r in results.values()) == 107530 and sum(r['leaves'] for r in results.values()) == 111201
    minimum = min((Q(r['minimum110'][0]),r['minimum110'][1]) for r in results.values())
    passed = sum(r['passed115'] for r in results.values())
    pending = sorted((x for r in results.values() for x in r['pending115']),key=lambda x:x['leaf_index'])
    assert pending == load(folder/'target115_pending.json') and passed+len(pending) == 111201
    candidate = load(folder/'full_domain_candidate.json');floor = Q(11,100000)*(1+minimum[0])
    assert candidate['full110_leaf_count'] == 111201 and candidate['complete_old110_domain_preserved']
    assert Q(candidate['minimum_margin110_new']) == minimum[0] and candidate['minimum_witness_leaf_index'] == minimum[1]
    assert Q(candidate['minimum_margin115']) == minimum[0]-upper and Q(candidate['candidate_floor']) == floor
    assert candidate['target115_passed'] == passed and candidate['target115_pending'] == len(pending)
    producer = load(binding['producer_summary_path'])
    assert producer['status'] == 'complete_fixed153_jointretarget_candidate' and producer['fault'] is None
    assert producer['master_pid'] == binding['producer_master_pid']
    assert len(producer['worker_evidence']) == 15 and len({r['pid'] for r in producer['worker_evidence']}) == 15
    assert all(r['exitcode'] == 0 and not r['alive'] for r in producer['worker_evidence'])
    assert producer['q_objects'] == 107530 and producer['leaf_cases'] == 111201
    assert producer['target115_passed'] == passed and producer['target115_pending'] == len(pending)
    assert producer['new_H_gradient_or_split_calls'] == producer['new_exp_calls'] == producer['oldq_corner_baseline_runs'] == 0
    for path,h in {**plan['input_locks'],**binding['candidate_locks']}.items():
        assert sha(path) == h,path
    write(args.output/'result.json',{'status':'finite155_newq_log_jointcorners_and_complete110_cover_independently_qualified',
          'pid':os.getpid(),'worker_pids':[p.pid for p in processes],'objects':107530,'leaves':111201,
          'new_corners':1779216,'new_public_log_objects':1,'minimum110':[str(minimum[0]),minimum[1]],
          'passed115':passed,'pending115':len(pending),'old110_complete_domain_preserved':True,
          'producer_linear_floor':str(floor),'source_hashes_before_after_unchanged':True,
          'old_H_gradient_log_or_corner_baselines_replayed':False,'new_Cn_claim':False,
          'seconds':time.monotonic()-start})


if __name__ == '__main__':
    main()
