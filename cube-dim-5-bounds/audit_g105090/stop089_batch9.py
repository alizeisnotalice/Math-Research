#!/usr/bin/env python3
"""Independent reconstruction of finite batch9 with preserved counter leaves."""
import sys
sys.dont_write_bytecode = True
import json, gzip, os, time
from pathlib import Path
from fractions import Fraction as F
from interfaces import REG, HERE, common as C
from structural089 import load_base, choose, certify, write


def main():
    start = time.monotonic(); registration_path = HERE/'stop089_batch9_registration.json'
    registration = json.loads(registration_path.read_text())
    for path, digest in registration['input_locks'].items():
        C.require(C.file_sha(path) == digest, 'frozen input changed: '+path)
    for path, digest in REG['protected_old_locks'].items():
        C.require(C.file_sha(path) == digest, 'protected old file changed')
    root = Path(registration['batch_directory'])
    manifest = json.loads((root/'manifest.json').read_text())
    for path, digest in manifest['files'].items(): C.require(C.file_sha(path) == digest, 'manifest file changed')
    previous_result = json.loads((HERE/'structural_through0008_result.json').read_text())
    C.require(previous_result['status'] == 'all_frozen089_structure_values_qualified_partial_g105' and previous_result['through_batch'] == 8, 'through8 not accepted')
    previous = json.loads((root.parent/'batch0008/state.json').read_text())
    initial = json.loads((root/'initial_overlay.json').read_text())
    final = json.loads((root/'state.json').read_text())
    br = json.loads((root/'registration.json').read_text())
    C.require(br['parent_overlay_sha256'] == C.file_sha(root.parent/'batch0008/state.json'), 'wrong prior overlay')
    C.require(br['source_state_sha256'] == initial['source_sha256'] == REG['source086_sha256'], 'wrong fixed source')
    for key in ['new_cache', 'new_passed', 'new_splits', 'queue', 'hardness', 'counts', 'pre_split_count']:
        C.require(initial[key] == previous[key], 'initial does not retain previous prefix: '+key)
    C.require(initial['counterexample'] is None and initial['faults'] == [] and initial['inflight'] == initial['resume_jobs'] == {}, 'initial work is not clean')
    C.require(br['selected_nodes'] == initial['queue'] and len(initial['queue']) == 186 and C.digest(br['selected_nodes']) == br['selected_sha256'], 'registered selection is not exact186')
    C.require(br['workers'] == 12 and br['wall_deadline'] is None, 'resource registration differs')
    base = load_base(set())
    C.require(br['base_cache_canonical_sha256'] == base['cache_sha256'], 'source cache reference changed')
    fan = json.loads(gzip.decompress(Path(REG['fan']).read_bytes()))
    numeric = json.loads((HERE/'obstruction9_values_result.json').read_text())
    C.require(numeric['status'] == 'all12_H_plus4grad_verified_four_fixed_H64_obstructions', 'new12 not independently verified')
    numeric_registration = json.loads((HERE/'obstruction9_values_registration.json').read_text())
    C.require(numeric['registration_sha256'] == C.file_sha(HERE/'obstruction9_values_registration.json'), 'numeric registration changed')
    for path, digest in numeric_registration['input_locks'].items(): C.require(C.file_sha(path) == digest, 'numeric dependency changed')
    pm = json.loads((root/'point_manifest.json').read_text())
    entries = {entry['path']: entry for entry in pm['points']}
    proofs = {proof['point_manifest_index']: proof for proof in numeric['points']}
    expected = json.loads(json.dumps(initial)); queue = list(initial['queue']); hard = dict(initial['hardness'])
    inflight = {}; dispatch = {}; received = {}; merge_count = 0; counters = []; complete = []; certs = []
    available = base['keys'] | set(initial['new_cache']); seen_keys = set(); prefix = '0'*64; seq = 0; stop_seen = False
    for line in (root/'journal.jsonl').read_text().splitlines():
        event = json.loads(line); body = event['body']; kind = event['kind']
        C.require(event['seq'] == seq+1 and event['previous'] == prefix and C.digest({k:v for k,v in event.items() if k!='hash'}) == event['hash'], 'journal chain differs')
        seq += 1; prefix = event['hash']
        if kind == 'dispatch':
            ordinal = body['ordinal']
            C.require(not stop_seen and ordinal == len(dispatch) and ordinal < 12, 'dispatch after obstruction or beyond first pool')
            C.require(body['node'] == br['selected_nodes'][ordinal] == queue.pop(0), 'dispatch not exact selection prefix')
            point, candidate, tests = choose(C.node(body['node']), base['roots'], fan['vertices'])
            C.require(body['point'] == point and body['candidate'] == candidate and body['geometry_tests'] == tests, 'independent physical anchor differs')
            key = ','.join(point)
            C.require(body['cached'] is False and key not in available and key not in seen_keys, 'new anchor reuses accepted source point')
            seen_keys.add(key); dispatch[ordinal] = body; inflight[str(ordinal)] = dict(body)
        elif kind == 'received':
            ordinal = body['ordinal']; got = body['result']; task = dispatch[ordinal]
            C.require(ordinal not in received and got['kind'] in ('complete','counterexample'), 'fault or duplicate result')
            entry = entries[task['point_checkpoint']]
            C.require(C.file_sha(entry['path']) == entry['sha256'], 'raw SHA differs')
            raw = json.loads(Path(entry['path']).read_text()); record = got['record']
            C.require(raw['stage'] == entry['outcome'] == got['kind'] and raw['partial'] is None, 'raw outcome differs')
            C.require(raw['record'] == entry['record'] == record and raw['point'] == record['point'] == task['point'], 'raw/return/dispatch record differs')
            C.require(raw['task_sha256'] == C.digest(task) and raw['kernel_segment_calls'] == got['calls'] == 1, 'raw task/call association differs')
            C.require(raw['pid'] == got['pid'] == entry['pid'] and record['bits'] == entry['bits'] == 256 and record['attempts'] == raw['attempts'] == got['attempts'], 'precision/worker association differs')
            C.require(all(a['all_cones'] == 3504 and a['numeric_quality'] is True for a in record['attempts']), 'not complete cone evaluation')
            C.require(record['candidate'] == task['candidate'] and record['geometry_tests'] == task['geometry_tests'], 'generic geometry record differs')
            proof = proofs[ordinal]
            C.require(proof['point'] == record['point'] and proof['raw_checkpoint_sha256'] == entry['sha256'], 'independent numeric point differs')
            for value, ball in zip(proof['independent_values_512_dyadic'], record['values']):
                lo, hi = map(F,value); L,U = map(F,ball); C.require(L < lo <= hi < U, 'independent value outside original ball')
            if got['kind'] == 'counterexample':
                C.require(F(record['values'][0][1]) < F(21,200000) and F(proof['independent_values_512_dyadic'][0][1]) < F(21,200000), 'counter upper does not prove fixed-H obstruction')
                stop_seen = True
            received[ordinal] = got; inflight[str(ordinal)]['returned'] = got
        elif kind == 'merge':
            ordinal = body['ordinal']; C.require(ordinal == merge_count and ordinal in received, 'merge ordering differs')
            merge_count += 1; task = dispatch[ordinal]; got = received[ordinal]; record = got['record']; node = C.node(task['node'])
            if got['kind'] == 'counterexample':
                C.require(set(body) == {'ordinal'}, 'counter leaf received an invalid positive certificate')
                counters.append(ordinal)
                expected['counterexample'] = {'H_interval':record['values'][0], 'node':task['node'], 'point':task['point'], 'record':record,
                    'scope':'H_epsilon target only; no fullG or C5 counterexample claim', 'target':'21/200000'}
            else:
                cert = body['certificate']; margin = certify(node, cert, record, base['roots'])
                certs.append({'ordinal':ordinal, 'node':list(node), 'minimum_log_margin':margin, 'passed':cert['passed']})
                expected['new_cache'][','.join(record['point'])] = record
                hard.pop(str(node[0])+':'+node[1],None); inflight.pop(str(ordinal)); complete.append(ordinal)
                if cert['passed']: expected['new_passed'].append({'node':list(node),'certificate':cert})
                else:
                    children = [[node[0],node[1]+'0'], [node[0],node[1]+'1']]
                    expected['new_splits'].append({'node':list(node),'axis':len(node[1])%4,'children':children,'failed_certificate':cert})
                    queue.extend(children)
                    for child in children: hard[str(child[0])+':'+child[1]] = [record['values'][0][0], cert['minimum_log_margin']]
        elif kind == 'finish':
            pool = body['pool_evidence']
            C.require(body['status'] == 'counterexample' and len(dispatch) == len(received) == merge_count == 12 and len(counters) == 4 and len(complete) == 8, 'finite stop sets differ')
            C.require(pool['close_called'] is True and pool['join_returned'] is True and pool['terminated_on_actual_fault'] is False and pool['master_pid'] == br['master_pid'], 'pool was not cleanly joined')
            C.require(len(pool['worker_pids']) == 12 and all(x['exitcode'] == 0 and x['is_alive'] is False for x in pool['exits']), 'worker pool did not exit0')
            C.require({x['pid'] for x in received.values()} <= set(pool['worker_pids']), 'returned PID outside pool')
            expected['pool_evidence'] = pool; expected['status'] = body['status']
        else: raise AssertionError('unknown event kind')
    queue.sort(key=lambda n:(*map(F,hard[str(n[0])+':'+n[1]]),*n))
    expected.update({'queue':queue, 'hardness':hard, 'inflight':inflight, 'resume_jobs':{}, 'next_merge':merge_count, 'journal_seq':seq, 'journal_hash':prefix,
        'counts':{'cache_hits':initial['counts']['cache_hits'], 'kernel_segment_calls':initial['counts']['kernel_segment_calls']+12,
                  'new_anchor_starts':initial['counts']['new_anchor_starts']+12, 'new_complete':initial['counts']['new_complete']+8}})
    C.require(expected == final, 'exact independent counter-overlay reconstruction differs')
    C.require(len(entries) == len(seen_keys) == 12 and len(queue) == 190 and len(inflight) == 4 and final['faults'] == [], 'terminal finite counts differ')
    C.require(all(','.join(dispatch[i]['point']) not in final['new_cache'] for i in counters), 'counter was put in positive cache')
    for split in final['new_splits']:
        node = C.node(split['node'])
        C.require(split['axis'] == len(node[1])%4 and split['children'] == [[node[0],node[1]+'0'],[node[0],node[1]+'1']], 'closed child partition differs')
    leaves = base['passed']+[C.node(x['node']) for x in final['new_passed']]+list(map(C.node,queue))+[C.node(x['node']) for x in inflight.values()]
    tree = C.validate_tree(base['roots'],leaves,set(base['splits'])|{C.node(x['node']) for x in final['new_splits']},'089terminalIncluding4CounterLeaves')
    pids = [br['master_pid'],*final['pool_evidence']['worker_pids']]
    for pid in pids:
        try: os.kill(pid,0)
        except ProcessLookupError: continue
        raise AssertionError('production PID still present or reused')
    for path, digest in registration['input_locks'].items(): C.require(C.file_sha(path) == digest, 'input changed during audit')
    for path, digest in REG['protected_old_locks'].items(): C.require(C.file_sha(path) == digest, 'old lock changed during audit')
    result = {'task':'audit_g105090_terminal089', 'status':'fixed_H64_obstruction_confirmed_terminal_partial_tree_qualified',
        'pid':os.getpid(), 'wall_seconds':time.monotonic()-start, 'registration_sha256':C.file_sha(registration_path),
        'journal_events':seq, 'selected':186, 'dispatched_and_received':12, 'positive_complete':8, 'new_positive_passed':sum(x['passed'] for x in certs),
        'new_positive_split':sum(not x['passed'] for x in certs), 'fresh_fullQ_corners':16*len(certs), 'fresh_certificates':certs,
        'counter_ordinals':counters, 'all12_all5_values_strictly_verified':True, 'independent_overlay_exact':True,
        'full_pending_frontier':194, 'queue':190, 'retained_counter_inflight':4, 'passed':len(base['passed'])+len(final['new_passed']),
        'splits':len(base['splits'])+len(final['new_splits']), 'positive_cache':len(base['keys'])+len(final['new_cache']), 'tree':tree,
        'production_pids_os_absent':pids, 'no_dispatch_after_counter_received':True, 'old_locks_unchanged':True,
        'source180040_recomputed':False, 'prior4744_recomputed':False, 'numeric_values_total':4756, 'numeric_components_total':23780,
        'global_g105_qualified':False, 'full_G_counterexample':False, 'Cn_C5_I_updated':False}
    path = HERE/'stop089_batch9_result.json'; write(path,result)
    print(json.dumps({k:result[k] for k in ['pid','status','wall_seconds','journal_events','passed','splits','positive_cache','full_pending_frontier']}),flush=True)
    print(json.dumps({'result_sha256':C.file_sha(path)}),flush=True)


if __name__ == '__main__': main()
