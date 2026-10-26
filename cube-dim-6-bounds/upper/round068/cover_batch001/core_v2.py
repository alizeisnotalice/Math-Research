"""C068 event-replayed cover state. No kernel/process launcher.

The frozen C063 module supplies exact geometry and deterministic transitions;
this module replaces permissive snapshot acceptance by full typed-event replay.
The context authenticates artifacts and reconstructs mathematical observations.
"""
import copy
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re

SCHEMA = 'C068-cover-event-snapshot-v2'
ROOT_BOX = [['0', '1']] * 5


def require(ok, message):
    if not ok:
        raise ValueError(message)


def canonical(data):
    return (json.dumps(data, sort_keys=True, separators=(',', ':'), ensure_ascii=False) + '\n').encode()


def digest(data):
    return hashlib.sha256(canonical(data)).hexdigest()


def exact_keys(data, names):
    require(isinstance(data, dict) and set(data) == set(names), 'Unexpected/missing event fields')


def validate_domain(state, context):
    context.base.validate_state(state)
    root = state['nodes']['r']
    require(root['id'] == 'r' and root['parent'] is None and root['depth'] == 0,
            'Immutable root identity changed')
    require(root['raw_box'] == ROOT_BOX and root['box'] == ROOT_BOX,
            'Root must be the complete closed ordered-face bounding box')
    require(root['tightening'] == context.base.tighten(ROOT_BOX), 'Root order domain changed')
    require(state['input_manifest_sha256'] == context.source_state['input_manifest_sha256'], 'C063 input identity changed')
    require(state['model'] == context.base.MODEL and state['target'] == '1/65000', 'Model/target changed')
    require(state['batch']['id'] == 'first-batch-proposed' and state['batch']['max_new_nodes'] == 8 and
            state['batch']['max_native_reservations'] == 8, 'This finite batch has immutable 8/8 caps')
    blocked = [node['id'] for node in state['nodes'].values() if node['status'] == 'blocked']
    require((not blocked and state['halted'] is None) or
            (len(blocked) == 1 and state['halted'] is not None and state['halted']['node_id'] == blocked[0]),
            'Blocked frontier and halt identity disagree')
    if context.test_only:
        require(state['environment'] == 'TEST_ONLY' and state['release'] is None, 'Fixture state cannot become production')
    else:
        require(state['environment'] in ('PREPARE_ONLY', 'MAIN_RELEASED'), 'Fixture state forbidden in production history')
        if state['environment'] == 'PREPARE_ONLY':
            require(state == context.source_state, 'Unreleased preparation state cannot consume/change the root tree')
        else:
            require(isinstance(state['release'], dict) and isinstance(state['release'].get('release_sha256'), str),
                    'MAIN_RELEASED requires an actual immutable release object')
            authentic = context.verified_release(context.root / 'main_release.json', state['release']['release_sha256'])
            require(state['release'] == authentic, 'Release metadata was not derived from the actual reviewed main release')
    return True


def genesis_record(context):
    state = copy.deepcopy(context.source_state)
    if context.test_only:
        state['environment'] = 'TEST_ONLY'
    validate_domain(state, context)
    return dict(schema=SCHEMA, sequence=0, previous_sha256=None,
                event=dict(kind='GENESIS_FROM_FROZEN_C063', source_snapshot_sha256=context.source_snapshot_sha256,
                           execution_input_manifest_sha256=context.execution_input_manifest_sha256,
                           synthetic_test_history=context.test_only), state=state)


def apply_event(previous, event, context):
    """Reconstruct a successor; an arbitrary self-consistent successor is invalid."""
    validate_domain(previous, context)
    require(isinstance(event, dict), 'Typed event required')
    kind = event.get('kind')
    if kind == 'AUTHORIZE_ONCE':
        exact_keys(event, ('kind', 'release_path', 'release_sha256'))
        require(not context.test_only and previous['environment'] == 'PREPARE_ONLY' and
                previous['release'] is None and previous['batch']['native_reservations'] == 0 and
                previous['inflight'] is None, 'Authorization cannot reset/replace an existing batch')
        release = context.verified_release(event['release_path'], event['release_sha256'])
        state = context.base.authorize(previous, release)
    elif kind == 'RESERVE_NEXT_ONCE':
        exact_keys(event, ('kind',))
        state, _ = context.base.reserve(previous)
    elif kind == 'ATTACH_ORIGINAL_HANDLE_ONCE':
        exact_keys(event, ('kind', 'token', 'handle_path', 'handle_sha256'))
        handle = context.verified_handle(previous, event['token'], event['handle_path'], event['handle_sha256'])
        state = context.base.attach_handle(previous, event['token'], handle)
    elif kind == 'IMPORT_ORIGINAL_RESULT_ONCE':
        exact_keys(event, ('kind', 'token', 'verification_path', 'verification_sha256'))
        require(previous['inflight'] is not None and previous['inflight']['token'] == event['token'], 'No matching inflight result')
        node = previous['nodes'][previous['inflight']['node_id']]
        attempt = node['attempts'][-1]
        require(attempt['execution_handle'] is not None, 'Original tool/PID handle must be retained before import')
        observation = context.verified_observation(previous, event['token'], event['verification_path'], event['verification_sha256'])
        state = context.base.ingest(previous, event['token'], observation)
    elif kind == 'HALT_ORIGINAL_ATTEMPT_NO_RETRY':
        exact_keys(event, ('kind', 'token', 'failure_path', 'failure_sha256'))
        failure = context.verified_failure(previous, event['token'], event['failure_path'], event['failure_sha256'])
        state = context.base.halt_inflight(previous, event['token'], failure['reason'],
                                           event['failure_sha256'], native_launch_observed=failure['native_launch_observed'])
    else:
        raise ValueError('Unregistered event; reset/delete/requeue/reissue is not an operation')
    validate_domain(state, context)
    # These are redundant to deterministic replay, but make lifecycle obligations explicit.
    require(state['next_attempt_serial'] >= previous['next_attempt_serial'] and
            state['batch']['native_reservations'] >= previous['batch']['native_reservations'], 'Attempt budget regressed')
    require(set(previous['nodes']).issubset(state['nodes']), 'Old node disappeared')
    for node_id, old_node in previous['nodes'].items():
        new_node = state['nodes'][node_id]
        require(len(new_node['attempts']) >= len(old_node['attempts']), 'Attempt history disappeared')
        for old_attempt, new_attempt in zip(old_node['attempts'], new_node['attempts']):
            for key in ('token', 'input', 'input_sha256'):
                require(new_attempt[key] == old_attempt[key], 'Reserved input/token changed')
            if old_attempt['execution_handle'] is not None:
                require(new_attempt['execution_handle'] == old_attempt['execution_handle'], 'Original handle reset')
            if old_attempt['observation'] is not None:
                require(new_attempt['observation'] == old_attempt['observation'], 'Original imported observation changed')
        if old_node['status'] == 'passed':
            require(new_node == old_node, 'Passed leaf changed or was recomputed')
    if previous['halted'] is not None:
        raise ValueError('A halted batch has no further transitions; no implicit next batch')
    return state


def verify_records(records, context):
    require(bool(records), 'No genesis')
    expected_genesis = genesis_record(context)
    require(records[0] == expected_genesis, 'Genesis is not exactly derived from the sealed C063 snapshot')
    previous_hash = digest(records[0])
    previous = records[0]['state']
    for sequence, record in enumerate(records[1:], 1):
        exact_keys(record, ('schema', 'sequence', 'previous_sha256', 'event', 'state'))
        require(record['schema'] == SCHEMA and record['sequence'] == sequence and
                record['previous_sha256'] == previous_hash, 'Snapshot sequence/hash chain broken')
        expected = apply_event(previous, record['event'], context)
        require(record['state'] == expected, 'Successor is not the complete replayed legal transition')
        previous, previous_hash = expected, digest(record)
    return dict(sequence=len(records)-1, sha256=previous_hash, state=previous,
                recovery_action='OBSERVE_ORIGINAL_HANDLE_OR_HALT_UNKNOWN_NEVER_RELAUNCH' if previous['inflight'] else
                                ('HALTED_BATCH_NO_TRANSITION' if previous['halted'] else 'MAIN_MAY_INSPECT_OR_RESERVE_WITHIN_SAME_RELEASE'))


def recover(directory, context):
    directory = Path(directory)
    files = sorted(path for path in directory.iterdir() if re.fullmatch(r'snapshot_[0-9]{8}\.json', path.name))
    records = []
    for sequence, path in enumerate(files):
        require(path.name == f'snapshot_{sequence:08d}.json', 'Snapshot gap')
        raw = path.read_bytes()
        record = json.loads(raw)
        require(raw == canonical(record), 'Partial/noncanonical snapshot')
        records.append(record)
    result = verify_records(records, context)
    result['pending_files'] = sorted(path.name for path in directory.iterdir() if path.name.endswith('.pending'))
    return result


def publish(directory, record):
    """Called only while main's exclusive writer lock is held."""
    directory = Path(directory)
    path = directory / f'snapshot_{record["sequence"]:08d}.json'
    pending = directory / f'{path.name}.{os.getpid()}.pending'
    raw = canonical(record)
    with pending.open('xb') as stream:
        stream.write(raw); stream.flush(); os.fsync(stream.fileno())
    os.link(pending, path)
    pending.unlink()
    fd = os.open(directory, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)
    return hashlib.sha256(raw).hexdigest()


def initialize_preparation_genesis(directory, context):
    """Preparation-only one-time initialization; never exposed as a resume reset."""
    directory = Path(directory)
    directory.mkdir()  # exclusive; a missing/deleted history cannot be silently reset
    marker = directory.parent / 'state_initialized_once.json'
    with marker.open('xb') as stream:
        stream.write(canonical(dict(source_snapshot_sha256=context.source_snapshot_sha256,
                                    execution_input_manifest_sha256=context.execution_input_manifest_sha256,
                                    note='No reset or reinitialization is permitted')))
        stream.flush(); os.fsync(stream.fileno())
    with (directory / '.writer.lock').open('a+b') as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        return publish(directory, genesis_record(context))


def append_event(directory, context, expected_previous_sha256, event):
    directory = Path(directory)
    require((directory.parent / 'state_initialized_once.json').exists(), 'Missing once-only state origin marker')
    with (directory / '.writer.lock').open('a+b') as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        latest = recover(directory, context)
        require(latest['sha256'] == expected_previous_sha256, 'Stale writer; do not retry a dispatch')
        state = apply_event(latest['state'], event, context)
        record = dict(schema=SCHEMA, sequence=latest['sequence']+1, previous_sha256=latest['sha256'],
                      event=copy.deepcopy(event), state=state)
        return publish(directory, record)
