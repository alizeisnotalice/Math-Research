#!/usr/bin/env python3
"""A090 saved-only audit. Never launch, import reconstruct, certify, or write source.

The one production-core recover uses a clearly labelled reviewer Context whose
only adaptation replaces oracle/Arb reconstruction by exact saved-byte checks.
"""
import sys
sys.dont_write_bytecode = True
import signal
signal.signal(signal.SIGALRM, signal.SIG_DFL)
signal.setitimer(signal.ITIMER_REAL, 179.5)
import resource
resource.setrlimit(resource.RLIMIT_CPU, (20, 20))
import hashlib
import json
import os
from pathlib import Path
import time
from fractions import Fraction as Q
from datetime import datetime, timezone

START = time.monotonic()
CPU = time.process_time()
OUT = Path(__file__).resolve().parent
ROOT = Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds')
C = ROOT/'upper/round068/cover_batch001'
R = ROOT/'upper/round082/runtime_revision'
A = ROOT/'upper/round090/cover_batch001_actual'
sys.path.insert(0, str(C))
import core_v2 as V
import math_contract as M
from batch_context import Context


def need(ok, text):
    if not ok:
        raise ValueError(text)


def read(p):
    return json.loads(Path(p).read_text())


READ_HASHES = {}


def sha(p):
    p = Path(p)
    value = hashlib.sha256(p.read_bytes()).hexdigest()
    old = READ_HASHES.setdefault(str(p), value)
    need(old == value, 'Source changed during audit: '+str(p))
    return value


def ref(item):
    need(sha(item['path']) == item['sha256'], 'Hash reference mismatch: '+item['path'])
    return read(item['path'])


def write_new(name, data):
    with (OUT/name).open('x') as stream:
        json.dump(data, stream, indent=2, sort_keys=True)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())


SUMMARIES = []


class SavedOnlyContext(Context):
    def verified_observation(self, previous, token, path, expected_sha):
        saved = self.bound_file(token, path, expected_sha, 'verification.json')
        folder = self.attempt_directory(token)
        attempt = previous['nodes'][previous['inflight']['node_id']]['attempts'][-1]
        expected = attempt['input']
        native, receipt = read(folder/'native.json'), read(folder/'native_receipt.json')
        intent, original = read(folder/'intent.json'), read(folder/'original_process.json')
        reserved, attached = read(folder/'reserved_input.json'), read(folder/'original_handle.json')
        ns, rs, ss, ins = (sha(folder/name) for name in
                           ('native.json','native_receipt.json','stderr.txt','intent.json'))
        need(reserved['input'] == expected and reserved['token'] == token and
             reserved['input_sha256'] == attempt['input_sha256'], 'Reserved input binding')
        need(intent['input'] == expected and intent['token'] == token and
             intent['input_sha256'] == attempt['input_sha256'], 'Original intent binding')
        need(receipt['intent_sha256'] == ins and receipt['release_sha256'] == previous['release']['release_sha256']
             and intent['release_sha256'] == previous['release']['release_sha256'], 'Release/intent binding')
        need(intent['command'] == receipt['command'] and intent['synthetic_fixture'] is False,
             'Intent native command/production binding')
        handle = attempt['execution_handle']
        need(handle == attached['handle'] and handle['native_pid'] == receipt['owned_native_pid'] == original['owned_native_pid']
             and handle['wrapper_pid'] == receipt['owned_wrapper_pid'] == original['owned_wrapper_pid'] == intent['owned_wrapper_pid'],
             'Original native/wrapper binding')
        need(original['token'] == token and original['input_sha256'] == attempt['input_sha256'] and
             original['owned_process_group_id'] == original['runtime_envelope_worker_pid'] == handle['wrapper_pid'] and
             receipt['native_inherited_worker_pgid'] == handle['wrapper_pid'], 'Original native PG binding')
        need(receipt['original_process_sha256'] == sha(folder/'original_process.json'), 'Original process hash')
        need(receipt['runtime_release_sha256'] == sha(R/'runtime_release.json'), 'Runtime release binding')
        need(receipt['automatic_retries'] == 0 and 0 <= receipt['wall_seconds'] <= 180, 'Native retry/time')
        need(not (folder/'failure.json').exists(), 'Failure artifact in successful attempt')
        checked = M.validate_native(self, expected, token, native, receipt, ns, ss)
        # The following arithmetic decodes original dyadics and validates Euler;
        # it does not evaluate any point, gradient, oracle, or exponential.
        proof = f'C068_ORIGINAL_NATIVE:{folder/"native.json"}#sha256={ns};receipt={rs}'
        jet = {**{key: expected[key] for key in ('model','model_definition','window','center')},
               'F_interval': list(map(str, checked['F'])),
               'face_gradient_intervals': [list(map(str,pair)) for pair in checked['physical_gradients'][1:]],
               'source_proof': proof,
               'provenance': dict(kind='C016_NATIVE_SAME_WINDOW_POSITIVE_BAND', precision_bits=512,
                    full_group_terms=967680, computed_sign_pairs=483840, synthetic_fixture=False,
                    native_sha256=ns, receipt_sha256=rs, gradient_gate_sha256=self.f061_manifest_sha256,
                    independent_full_gradient_crosschecked_for_this_new_attempt=False)}
        self.base.validate_jet(jet, expected, test_only=False)
        observation = dict(input_sha256=attempt['input_sha256'],receipt_exit_code=0,stderr_empty=True,jet=jet)
        classification = ('CENTER_STRICTLY_ABOVE' if checked['F'][0] > self.base.TARGET else
                          'CENTER_STRICTLY_BELOW' if checked['F'][1] < self.base.TARGET else
                          'CENTER_OVERLAP_OR_EQUALITY_STOP')
        box_pass = None
        box_lower = None
        if checked['F'][0] > self.base.TARGET:
            box = saved['observation']['box_certificate']
            box_pass = self.base.check_saved_box(jet, expected, box)
            need(box['precision_bits'] == 384 and box['kernel_calls'] == 0 and box['global_g_certificate'] is False,
                 'Saved box precision/scope')
            need(box['scope'] == 'one whole box on x1=1; no global coverage assertion', 'Box scope')
            need(Q(box['tangent_bound_upper']) == min(Q(x['certificate_upper']) for x in box['corners']),
                 'Saved upper minimum')
            for row in box['corners']:
                need(row['passes_target'] == (Q(row['certificate_lower']) >= self.base.TARGET), 'Corner pass flag')
            observation.update(box_certificate=box, box_certificate_sha256=self.base.digest(box),
                               box_checker_sha256=self.base.BOX_CHECKER_SHA256)
            box_lower = box['tangent_bound_lower']
        expected_saved = dict(schema='C068-reconstructed-native-and-box-v1', token=token,
            input_sha256=attempt['input_sha256'],native_artifact_sha256=ns,receipt_artifact_sha256=rs,
            stderr_artifact_sha256=ss,intent_artifact_sha256=ins,source_sha256=receipt['source_sha256'],
            binary_sha256=receipt['binary_sha256'],fan_sha256=receipt['fan_sha256'],
            complete_physical_gradient_intervals=[list(map(str,pair)) for pair in checked['physical_gradients']],
            endpoint_intervals={'p_a':list(map(str,checked['p_a'])),'p_T':list(map(str,checked['p_T']))},
            independent_Euler_interval=list(map(str,checked['independent_Euler'])),
            native_Euler_interval=list(map(str,checked['native_Euler'])),native_counts=checked['counts'],
            independent_zero_repeat_count_recomputation=False,center_classification=classification,
            observation=observation,exp_reconstructed_with_pinned_A057_checker=checked['F'][0] > self.base.TARGET,
            synthetic_fixture=False,universal_gradient_crosscheck_claim=False,
            global_g_certificate=False,new_C6_bound=False)
        need(saved == expected_saved, 'Saved native/box fields differ from exact original-byte check')
        SUMMARIES.append(dict(token=token,center=expected['center'],window=expected['window'],
            face_box=expected['face_box'],input_sha256=attempt['input_sha256'],
            handle=handle,center_classification=classification,F_interval=list(map(str,checked['F'])),
            whole_box_pass=box_pass,whole_box_lower=box_lower,native_counts=checked['counts'],
            kernel_seconds=native['kernel_seconds'],native_wrapper_wall_seconds=receipt['wall_seconds'],
            verification_sha256=expected_sha,receipt_sha256=rs))
        for name in ('reserved_input.json','intent.json','original_process.json','original_handle.json',
                     'native.json','stderr.txt','native_receipt.json','verification.json'):
            sha(folder/name)
        return observation


def main():
    need(len(sys.argv)==3, 'Need final seal and independently observed session table')
    seal_path, observed_path = map(Path,sys.argv[1:])
    need(seal_path.parent == A, 'Final seal outside actual 090 directory')
    seal = read(seal_path)
    need(sha(seal_path)=='5b9f004437d34fc0fd58be76ef762f1a7f975fa7cc4ce0e61dac2a89cc02f82d',
         'Main final scientific freeze identity')
    need(seal['schema']=='090-scientific-artifact-freeze-v1' and seal['native_calls']==8 and
         seal['new_C6_bound'] is False and seal['actual_exit_code']==0,'Main final freeze fields')
    for item in seal['files']:
        p=ROOT/item['path']
        need(p.is_relative_to(ROOT) and '..' not in Path(item['path']).parts,'Seal path scope')
        need(p.stat().st_size==item['bytes'] and sha(p)==item['sha256'],'Final seal artifact mismatch '+str(p))
    observed = read(observed_path)
    sha(observed_path)
    # Do not infer table shape or fields: parent supplies a normalized observation
    # table with provenance referring to original tool executions.
    need(observed['schema']=='A090-main-observed-original-sessions-v1', 'Observed table schema')
    need(len(observed['attempts'])==8, 'Exactly eight original observations')
    observed_by_token = {x['token']:x for x in observed['attempts']}
    need(len(observed_by_token)==8, 'Unique observed tokens')
    for name in ('registration.json','core_acceptance.json','runtime_acceptance.json','preflight_receipt.json',
                 'preflight_environment_failure.json'):
        sha(A/name)
    for folder in (C,R):
        for name in ('input_manifest.json','preparation_manifest.json'):
            sha(folder/name)
    context = SavedOnlyContext()
    # Fail closed if anything attempts transcendental/old oracle recomputation.
    def forbidden(*args,**kwargs):
        raise RuntimeError('Forbidden oracle/Arb recomputation in saved-only audit')
    M.reconstruct = forbidden
    context.checker.certify = forbidden
    latest = V.recover(C/'state', context)  # only call, original event chain
    need(len(SUMMARIES)==8, 'Exactly eight imports in one core recover')
    state=latest['state']
    need(state['batch']['native_reservations']==8 and len(state['batch']['started_nodes'])==8 and
         state['production_kernel_calls']==8 and state['next_attempt_serial']==9,'8/8 consumed, no refund')
    need(state['inflight'] is None and state['halted'] is None and latest['pending_files']==[], 'No unresolved inflight/halt/pending')
    need(len(state['queue'])==7 and state['passed_leaves']==['r000'] and not state['empty_leaves'],
         'Actual seven-node unresolved frontier and one passed leaf expected')
    need(all(x['center_classification']=='CENTER_STRICTLY_ABOVE' for x in SUMMARIES) and
         [x['token'] for x in SUMMARIES if x['whole_box_pass']]==['r000.a0.n000008'],
         'Eight center passes, only r000 whole-box pass')
    need(state['nodes']['r000']['box']==[['0','1/8']]*5 and
         sum(x['status']=='split' for x in state['nodes'].values())==7,'Passed local box and split count')
    need(state['global_g_certificate'] is False and state['new_C6_bound'] is False,'No global theorem flags')
    for item in SUMMARIES:
        obs=observed_by_token[item['token']]
        for key in ('tool_session_id','native_pid','wrapper_pid'):
            need(obs[key]==item['handle'][key], 'Observed tool/PID binding mismatch: '+item['token'])
        need(obs['native_receipt_sha256']==item['receipt_sha256'] and
             obs['original_handle_sha256']==sha(C/'attempts'/item['token']/'original_handle.json'),
             'Observed tool receipt/handle bytes binding')
    for path in sorted((C/'state').iterdir()):
        if path.is_file() and path.name!='.writer.lock':
            sha(path)
    sha(C/'state_initialized_once.json')
    sha(C/'main_release.json')
    sha(R/'runtime_release.json')
    # Record all pinned source/input evidence, not only the dynamic batch files.
    for path in context.pins['readonly_sha256']:
        sha(path)
    for name in context.pins['owned_source_sha256']:
        sha(C/name)
    for acceptance in (read(A/'core_acceptance.json'),read(A/'runtime_acceptance.json')):
        for item in acceptance['independent_review_artifacts']:
            ref(item)
    for name,digest in read(R/'preparation_manifest.json')['artifact_sha256'].items():
        need(sha(R/name)==digest,'Runtime seal artifact mismatch')
    # Imported state exactly retains ordered-root topology. Runtime and PID audit
    # is the separate second saved-only interface group.
    result=dict(status='PASS_SAVED_ONLY_CORE_AND_EXACT_CERTIFICATE_BINDINGS',
        final_seal_path=str(seal_path),final_seal_sha256=sha(seal_path),
        observed_original_sessions_path=str(observed_path),observed_original_sessions_sha256=sha(observed_path),
        source_root=str(ROOT),sequence=latest['sequence'],snapshot_sha256=latest['sha256'],
        queue=state['queue'],inflight=state['inflight'],halted=state['halted'],
        passed_leaves=state['passed_leaves'],empty_leaves=state['empty_leaves'],
        reservations=8,distinct_started_nodes=state['batch']['started_nodes'],
        original_native_calls_observed=state['production_kernel_calls'],attempts=SUMMARIES,
        pure_core_recover_calls=1,reviewer_context_adaptation='verified_observation saved-only exact check',
        saved_corner_rows_checked=sum(32 for x in SUMMARIES if x['whole_box_pass'] is not None),
        independent_exponential_recalculation=False,independent_full_gradient_recalculation=False,
        independent_zero_repeat_enumeration=False,native_calls=0,Arb_corner_calls=0,
        source_hashes=READ_HASHES,global_g_certificate=False,new_C6_bound=False)
    # Re-read every byte hash to detect concurrent or accidental source mutation.
    for p,digest in list(READ_HASHES.items()):
        need(hashlib.sha256(Path(p).read_bytes()).hexdigest()==digest,'Final source hash changed: '+p)
    write_new('core_results.json',result)
    write_new('core_run_receipt.json',dict(status='PASS',pid=os.getpid(),exit_code=0,
        utc=datetime.now(timezone.utc).isoformat(),python=sys.executable,
        wall_seconds=time.monotonic()-START,cpu_seconds=time.process_time()-CPU,
        cpu_limit_seconds=20,wall_limit_seconds=180,cleanup_and_receipt_included=True,
        native_calls=0,Arb_corner_calls=0,pure_core_recover_calls=1,automatic_retries=0,
        source_writes=0,result_sha256=sha(OUT/'core_results.json')))
    print(json.dumps({k:result[k] for k in ('status','sequence','snapshot_sha256','queue','reservations','saved_corner_rows_checked')}))


if __name__=='__main__':
    try:
        main()
    except BaseException as error:
        write_new('core_failure_receipt.json',dict(status='STOP_NO_RETRY',pid=os.getpid(),error=repr(error),
            wall_seconds=time.monotonic()-START,cpu_seconds=time.process_time()-CPU,native_calls=0,
            Arb_corner_calls=0,source_writes=0))
        raise
