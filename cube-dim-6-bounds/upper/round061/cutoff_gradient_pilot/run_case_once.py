"""C061 one-case runner. NO invocation is authorized by the preparation.

Usage after explicit main release only:
    python -B run_case_once.py CASE_ID MAIN_RELEASE.json

Each invocation is a different predeclared case; existing case directories are
never resumed or replayed. A 180 s alarm covers this entire program, including
the native child and box postprocessing. There is no compiler or kernel retry.
"""
import sys
sys.dont_write_bytecode = True
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import signal
import subprocess
import time
from datetime import datetime, timezone
from fractions import Fraction as Q
from jet_adapter import (MODEL, MODEL_DEFINITION, CASE_IDS, registered_case,
                         proposed_command, dyadic, build_jet, classification,
                         box_inputs)

P = Path(__file__).resolve().parent
LIMIT = 180
ABOVE = 'CENTER_STRICTLY_ABOVE_TARGET'
BELOW = 'CENTER_STRICTLY_BELOW_TARGET'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def utc():
    return datetime.now(timezone.utc).isoformat()


def write_new(path, data):
    with path.open('x') as stream:
        json.dump(data, stream, indent=2)
        stream.write('\n')


def require(test, message):
    if not test:
        raise ValueError(message)


def release_gate(release_path):
    release = read(release_path)
    require(release.get('schema') == 'C061-main-production-release-v1', 'Wrong release schema')
    require(release.get('production_authorized') is True, 'No explicit production release')
    require(release.get('released_by') == 'main', 'Main must own release')
    require(isinstance(release.get('authorization_message'), str) and release['authorization_message'],
            'Exact main authorization message is required')
    require(release.get('per_program_wall_timeout_seconds') == LIMIT, '180 s policy required')
    require(release.get('max_new_native_calls') in (2, 3), 'Only the registered 2/3-call batch')
    require(release.get('F057_final_PASS_accepted_by_main') is True, 'Main must accept final F057 PASS')
    required = ('registration.json', 'input_manifest.json', 'policy_amendment_180s.json',
                'preparation_results.json', 'preparation_execution_receipt.json',
                'jet_adapter.py', 'run_case_once.py')
    for name in required:
        require(sha(P / name) == release['C061_sha256'][name], f'Release pin mismatch: {name}')
    policy = read(P / 'policy_amendment_180s.json')
    require(policy['per_program_wall_timeout_seconds'] == LIMIT, 'Policy mismatch')
    require(sha(P / 'registration.json') == policy['original_registration_sha256'], 'Original registration changed')
    require(read(P / 'preparation_results.json')['status'] ==
            'PASS_TWELVE_PREPARATION_CONFIGURATIONS_ZERO_PRODUCTION_KERNELS', 'Preparation not passed')
    require(read(P / 'preparation_execution_receipt.json')['status'] ==
            'PASS_PREPARATION_ONLY_NO_PRODUCTION_RELEASE', 'Preparation receipt not passed')
    audit_path = Path(release['F057_final_PASS_path'])
    require(sha(audit_path) == release['F057_final_PASS_sha256'], 'F057 PASS pin changed')
    audit = read(audit_path)
    require(audit.get('status') == release['F057_expected_status'] and
            str(audit['status']).startswith('PASS'), 'F057 evidence is not the expected final PASS')
    a057_path = Path(policy['A057_sealed_manifest_path'])
    require(sha(a057_path) == policy['A057_sealed_manifest_sha256'], 'A057 manifest changed')
    require(release['A057_sealed_manifest_sha256'] == policy['A057_sealed_manifest_sha256'],
            'F057/main release must concern the pinned A057 seal')
    a057 = read(a057_path)
    for name, item in a057['artifacts'].items():
        require(sha(a057_path.parent / name) == item['sha256'], f'A057 artifact changed: {name}')
    pins = read(P / 'input_manifest.json')
    for path, digest in pins['source_sha256'].items():
        require(sha(Path(path)) == digest, f'Source/runtime/input changed: {path}')
    return release, a057_path.parent / 'box_certificate.py'


def branch_gate(case_id, release):
    root = P / 'production'
    expected = {'low-primary': [], 'low-fallback': ['low-primary'],
                'S3-primary': ['low-primary'], 'S3-fallback': ['low-primary', 'low-fallback']}[case_id]
    actual = sorted(path.name for path in root.iterdir()) if root.exists() else []
    require(actual == sorted(expected), 'Existing cases do not match the unique permitted branch')
    require(len(actual) < release['max_new_native_calls'], 'Released native budget exhausted')
    prior = {}
    for name in expected:
        folder = root / name
        receipt = read(folder / 'receipt.json')
        result = read(folder / 'result.json')
        require(receipt['status'] == 'COMPLETED_NO_RETRY' and
                receipt['result_sha256'] == sha(folder / 'result.json') and
                receipt['native_sha256'] == sha(folder / 'native.json'), 'Previous case failed or changed')
        require(receipt['release_sha256'] == sha(Path(release['_path'])), 'Release changed across cases')
        prior[name] = result
    if case_id == 'low-fallback':
        require(prior['low-primary']['center_classification'] == BELOW, 'Fallback only after strict primary failure')
    if case_id == 'S3-primary':
        require(prior['low-primary']['center_classification'] == ABOVE, 'S3 primary requires primary center pass')
    if case_id == 'S3-fallback':
        require(prior['low-primary']['center_classification'] == BELOW and
                prior['low-fallback']['center_classification'] == ABOVE, 'Wrong fallback branch')
    if case_id.startswith('S3'):
        low = prior['low-' + case_id.split('-')[1]]
        require(low['whole_box_status'] == 'PASS' or
                release.get('allow_S3_after_center_pass_with_unresolved_local_radius') is True,
                'Local radius unresolved; pause unless main explicitly released S3 after center pass')


def independent_euler(native, window, center):
    def scale(interval, factor):
        x, y = (factor * value for value in interval)
        return min(x, y), max(x, y)
    terms = [scale(dyadic(native['band']), Q(5)),
             scale(dyadic(native['p_T']), Q(window['T'])),
             scale(dyadic(native['p_a']), -Q(window['a']))]
    terms.extend(scale(dyadic(raw), Q(coordinate)) for raw, coordinate in
                 zip(native['gradient'], center))
    bounds = [sum(pair[index] for pair in terms) for index in (0, 1)]
    require(bounds[0] <= 0 <= bounds[1], 'Independent exact interval Euler misses zero')
    require(dyadic(native['p_a'])[0] > 0 and dyadic(native['p_T'])[0] > 0,
            'Endpoint density positivity unresolved')
    return list(map(str, bounds))


def main():
    if len(sys.argv) != 3 or sys.argv[1] not in CASE_IDS:
        raise SystemExit('Usage: run_case_once.py CASE_ID MAIN_RELEASE.json; explicit main release required')
    tic = time.monotonic()
    case_id, release_path = sys.argv[1], Path(sys.argv[2]).resolve()

    def timed_out(signum, frame):
        raise TimeoutError('180-second whole-program wall-clock limit')

    signal.signal(signal.SIGALRM, timed_out)
    signal.setitimer(signal.ITIMER_REAL, LIMIT)
    release, checker_path = release_gate(release_path)
    release['_path'] = str(release_path)
    branch_gate(case_id, release)
    registration, _, window, center = registered_case(case_id)
    root = P / 'production'
    root.mkdir(exist_ok=True)
    folder = root / case_id
    folder.mkdir()  # exclusive: an attempted case is never replayed
    receipt_path = folder / 'receipt.json'
    native_path, stderr_path = folder / 'native.json', folder / 'stderr.txt'
    command = proposed_command(case_id)
    receipt = dict(status='REGISTERED_NOT_YET_LAUNCHED', case_id=case_id, command=command,
                   model=MODEL, model_definition=MODEL_DEFINITION, window=window, center=center,
                   start_utc=utc(), timeout_seconds=LIMIT, release_path=str(release_path),
                   release_sha256=sha(release_path), source_sha256=sha(Path(registration['source_path'])),
                   binary_sha256=sha(Path(registration['binary_path'])),
                   fan_sha256=sha(Path(registration['fan_path'])),
                   production_attempt_count=1, new_native_invocations=0, automatic_retries=0,
                   new_compiles=0, old_seed_recalculations=0, C055_reruns=0,
                   global_g_certificate=False, new_C6_bound=False)
    write_new(receipt_path, receipt)

    def save_receipt():
        receipt['wall_seconds'] = time.monotonic() - tic
        temp = folder / 'receipt.json.tmp'
        temp.write_text(json.dumps(receipt, indent=2) + '\n')
        temp.replace(receipt_path)

    child = None
    try:
        with native_path.open('xb') as out, stderr_path.open('xb') as err:
            child = subprocess.Popen(command, cwd=P, stdout=out, stderr=err, start_new_session=True)
            receipt.update(status='NATIVE_RUNNING', owned_child_pid=child.pid,
                           owned_process_group_id=child.pid, new_native_invocations=1, launched_utc=utc())
            save_receipt()
            code = child.wait(timeout=max(0.001, LIMIT - (time.monotonic() - tic)))
        receipt.update(status='NATIVE_RETURNED', native_exit_code=code, native_end_utc=utc(),
                       native_sha256=sha(native_path), stderr_sha256=sha(stderr_path),
                       stderr_empty=stderr_path.stat().st_size == 0)
        save_receipt()
        require(code == 0 and receipt['stderr_empty'], 'Native failure: stop, no retry')
        native = read(native_path)
        source_proof = f'C061_SAME_WINDOW_NATIVE:{native_path}#sha256={receipt["native_sha256"]}'
        jet = build_jet(native, receipt, case_id, source_proof)
        euler = independent_euler(native, window, center)
        write_new(folder / 'jet.json', jet)
        center_class = classification(jet)
        result = dict(status='CASE_CHECKED', case_id=case_id, center_classification=center_class,
                      F_interval=jet['F_interval'], independent_Euler_interval=euler,
                      full_gradient_independently_crosschecked=False,
                      inherited_validation='C016 exact scalar DD5/6 and physical chain; Euler does not independently prove six components',
                      whole_box_status='NOT_APPLICABLE_CENTER_NOT_STRICTLY_ABOVE',
                      box_candidates=[], additional_kernel_calls_for_boxes=0,
                      native_node_counts={name: native[name] for name in
                                          ('full_group_terms', 'computed_sign_pairs', 'confluent_terms',
                                           'zero_node_terms', 'unique_node_values')},
                      native_kernel_seconds=native['kernel_seconds'],
                      global_g_certificate=False, new_C6_bound=False)
        if center_class == ABOVE:
            spec = importlib.util.spec_from_file_location('c061_pinned_a057_box', checker_path)
            checker = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(checker)
            result['whole_box_status'] = 'UNRESOLVED_AFTER_EIGHT_PREDECLARED_WIDTHS'
            for data in box_inputs(jet):
                checked = checker.certify(data, precision_bits=384)
                item = dict(halfwidth_index=data['halfwidth_index'], halfwidth_exact=data['halfwidth_exact'],
                            face_box=data['face_box'], certificate=checked)
                result['box_candidates'].append(item)
                write_new(folder / f'box_{data["halfwidth_index"]}.json', item)
                if checked['status'] == 'PASS_CONDITIONAL_ON_MATCHING_CENTER_JET':
                    result['whole_box_status'] = 'PASS'
                    result['passing_halfwidth_exact'] = data['halfwidth_exact']
                    break
        write_new(folder / 'result.json', result)
        receipt.update(status='COMPLETED_NO_RETRY', ended_utc=utc(),
                       result_sha256=sha(folder / 'result.json'), jet_sha256=sha(folder / 'jet.json'))
        save_receipt()
        print(json.dumps(dict(case_id=case_id, center_classification=center_class,
                              whole_box_status=result['whole_box_status'],
                              native_calls_this_invocation=1, global_g_certificate=False)))
    except BaseException as error:
        # Disable the expired timer only for owned-process reaping and durable receipt.
        signal.setitimer(signal.ITIMER_REAL, 0)
        if child is not None and child.poll() is None:
            try:
                os.killpg(child.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass  # The owned child exited between poll and kill.
            child.wait()
            receipt['killed_only_owned_process_group'] = True
        receipt.update(status='TIMEOUT_STOP_NO_RETRY' if isinstance(error, (TimeoutError, subprocess.TimeoutExpired))
                       else 'FAILED_STOP_NO_RETRY', error=repr(error), ended_utc=utc())
        for name in ('native.json', 'stderr.txt', 'jet.json'):
            path = folder / name
            if path.exists():
                receipt[name + '_retained_sha256'] = sha(path)
        save_receipt()
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


if __name__ == '__main__':
    main()
