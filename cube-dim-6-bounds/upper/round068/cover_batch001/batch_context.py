"""Pinned C068 context and authenticated artifact import; no process launches."""
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import sys
sys.dont_write_bytecode = True
import math_contract as M

P = Path(__file__).resolve().parent


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def import_pinned(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_new(path, data):
    raw = (json.dumps(data, sort_keys=True, separators=(',', ':'), ensure_ascii=False) + '\n').encode()
    with Path(path).open('xb') as stream:
        stream.write(raw)
        stream.flush()
        os.fsync(stream.fileno())
    return hashlib.sha256(raw).hexdigest()


class Context:
    def __init__(self, *, test_only=False, fixture_root=None):
        self.test_only = test_only
        require(test_only or fixture_root is None, 'Production context cannot redirect artifacts')
        self.root = Path(fixture_root) if test_only else P
        self.pins = read(P / 'input_manifest.json')
        self.execution_input_manifest_sha256 = sha(P / 'input_manifest.json')
        self.verify_pins()
        self.paths = self.pins['paths']
        self.fan_bounds = self.pins['fan_bounds']
        self.runtime_pins = self.pins['runtime_sha256']
        self.f061_manifest_sha256 = self.pins['F061_manifest_sha256']
        self.base = import_pinned('c068_frozen_c063', self.paths['C063_core'])
        self.checker = import_pinned('c068_frozen_a057', self.paths['A057_checker'])
        import flint
        require(flint.__FLINT_VERSION__ == '3.6.0', 'Box checker FLINT runtime version mismatch')
        self.source_snapshot_sha256 = sha(self.paths['C063_initial_snapshot'])
        require(self.source_snapshot_sha256 == self.pins['C063_initial_snapshot_sha256'], 'Frozen source genesis changed')
        original = read(self.paths['C063_initial_snapshot'])
        require(original['sequence'] == 0 and original['previous_sha256'] is None, 'Source must be original C063 genesis')
        self.source_state = original['state']
        require(self.source_state['environment'] == 'PREPARE_ONLY' and self.source_state['queue'] == ['r'] and
                self.source_state['inflight'] is None and self.source_state['batch']['native_reservations'] == 0,
                'Source genesis is not the unexecuted root')

    def external_sha(self, path):
        return self.pins['readonly_sha256'][str(path)]

    def verify_pins(self):
        for path, expected in self.pins['readonly_sha256'].items():
            require(sha(path) == expected, f'Immutable external input changed: {path}')
        for name, expected in self.pins['owned_source_sha256'].items():
            require(sha(P / name) == expected, f'Immutable C068 source changed: {name}')

    def attempt_directory(self, token):
        require(isinstance(token, str) and re.fullmatch(r'r[01]*\.a[0-4]\.n[0-9]{6}', token) is not None,
                'Malformed reservation token')
        return self.root / 'attempts' / token

    def bound_file(self, token, path, digest, name):
        path = Path(path)
        require(path.resolve() == (self.attempt_directory(token) / name).resolve(), 'Artifact escaped its unique attempt/name')
        require(sha(path) == digest, f'Original artifact changed: {name}')
        return read(path)

    def verified_release(self, path, expected_sha):
        require(not self.test_only, 'Fixture context cannot accept production release')
        path = Path(path)
        require(path.resolve() == (P / 'main_release.json').resolve() and sha(path) == expected_sha, 'Main release identity mismatch')
        release = read(path)
        require(release.get('schema') == 'C068-main-finite-batch-release-v1' and
                release.get('production_authorized') is True and release.get('owner') == 'main', 'No main production release')
        require(release.get('max_new_nodes') == 8 and release.get('max_native_reservations') == 8 and
                release.get('timeout_seconds') == 180, 'Finite batch/runtime policy mismatch')
        require(release.get('input_manifest_sha256') == self.execution_input_manifest_sha256 and
                release.get('preparation_manifest_sha256') == sha(P / 'preparation_manifest.json'), 'Release is not for this frozen revision')
        require(release.get('F061_manifest_sha256') == self.f061_manifest_sha256 and
                release.get('future_native_jets_accepted_under_engine_proof') is True and
                release.get('F061_is_not_a_universal_new_point_crosscheck') is True, 'Native proof scope not accepted')
        acceptance_ref = release['revision_acceptance']
        require(sha(acceptance_ref['path']) == acceptance_ref['sha256'], 'Revision acceptance changed')
        acceptance = read(acceptance_ref['path'])
        require(acceptance.get('status') == 'PASS_C068_REVISION_ACCEPTED_BY_MAIN' and
                acceptance.get('C068_preparation_manifest_sha256') == release['preparation_manifest_sha256'],
                'Independent revision review plus main acceptance missing')
        require(bool(acceptance.get('independent_review_artifacts')), 'No independent review evidence')
        for item in acceptance['independent_review_artifacts']:
            require(sha(item['path']) == item['sha256'], 'Independent revision review artifact changed')
        f061 = read(self.paths['F061_manifest'])
        require(f061['full_gradient_gap_closed_for_registered_cases'] is True and
                f061['universal_all_points_windows_gradient_crosscheck_claim'] is False, 'F061 scope mismatch')
        for name, digest in f061['artifact_sha256'].items():
            require(sha(Path(self.paths['F061_manifest']).parent / name) == digest, 'F061 sealed artifact changed')
        return dict(main_authorized=True, F061_gradient_gate_accepted=True,
                    F061_gradient_gate_sha256=self.f061_manifest_sha256,
                    input_manifest_sha256=self.source_state['input_manifest_sha256'],
                    timeout_seconds=180, release_sha256=expected_sha)

    def verified_handle(self, previous, token, path, expected_sha):
        require(previous['inflight'] is not None and previous['inflight']['token'] == token, 'No matching inflight handle')
        record = self.bound_file(token, path, expected_sha, 'original_handle.json')
        attempt = previous['nodes'][previous['inflight']['node_id']]['attempts'][-1]
        require(record['token'] == token and record['input_sha256'] == attempt['input_sha256'], 'Handle reservation mismatch')
        require(record['synthetic_fixture'] is self.test_only, 'Fixture handle in production')
        if not self.test_only:
            folder = self.attempt_directory(token)
            intent = read(folder / 'intent.json')
            require(record['intent_sha256'] == sha(folder / 'intent.json') and
                    record['handle']['wrapper_pid'] == intent['owned_wrapper_pid'], 'Original launcher mismatch')
            require(type(record['handle'].get('tool_session_id')) is int and record['handle']['tool_session_id'] > 0,
                    'Observed original tool session required')
            require(type(record['handle'].get('native_pid')) is int and record['handle']['native_pid'] > 0,
                    'Observed original native PID required')
            original = record['observed_original_process_record']
            require(original['token'] == token and original['owned_native_pid'] == record['handle']['native_pid'],
                    'Native PID was not recorded for this attempt')
            require(record['original_process_sha256'] == sha(folder/'original_process.json') and
                    original == read(folder/'original_process.json'), 'Original process record changed')
        return copy.deepcopy(record['handle'])

    def verified_observation(self, previous, token, path, expected_sha):
        saved = self.bound_file(token, path, expected_sha, 'verification.json')
        rebuilt = M.reconstruct(self, previous, token, self.attempt_directory(token))
        require(saved == rebuilt, 'Saved mathematical result differs from native/Arb reconstruction')
        return rebuilt['observation']

    def verified_failure(self, previous, token, path, expected_sha):
        require(previous['inflight'] is not None and previous['inflight']['token'] == token, 'No matching failed/orphan attempt')
        failure = self.bound_file(token, path, expected_sha, 'failure.json')
        attempt = previous['nodes'][previous['inflight']['node_id']]['attempts'][-1]
        require(failure['token'] == token and failure['input_sha256'] == attempt['input_sha256'], 'Failure identity mismatch')
        require(failure['synthetic_fixture'] is self.test_only, 'Fixture failure in production')
        require(failure['reason'] in ('NATIVE_ERROR_STOP', 'PROGRAM_TIMEOUT_STOP', 'UNKNOWN_ORPHAN_STOP',
                                      'VALIDATION_ERROR_STOP'), 'Unregistered stop reason')
        if not self.test_only:
            require(bool(failure.get('retained_artifacts')), 'Failure/orphan observations must be retained')
            for item in failure['retained_artifacts']:
                require(sha(item['path']) == item['sha256'], 'Failure evidence changed')
            if failure['reason'] == 'UNKNOWN_ORPHAN_STOP':
                require(failure['native_launch_observed'] is None, 'Unknown orphan cannot assert no launch')
        return failure
