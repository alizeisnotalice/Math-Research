"""C063 exact tree/state/adapter core. No process launch or kernel is present.

All computational evidence is supplied externally and remains conditional on
main's hash/proof acceptance. Snapshots are immutable; only main may append.
"""
import copy
import fcntl
import hashlib
import json
import math
import os
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import re

MODEL = 'body32_unit_cauchy_v1'
MODEL_DEFINITION = 'V={(1,s2,...,s6):si=+-1}; f_t(z)=t/(pi*(t^2+z^2)); 32 independent unit-weight factors'
LADDER = ('1/32', '1/64', '1/128', '1/250', '1/1000')
TARGET = Q(1, 65000)
BOX_CHECKER_SHA256 = '891a16137bbb483099d577e1ef7199e7a632f141d2517f2540713a19ad27c279'
STATUSES = {'queued', 'inflight', 'passed', 'split', 'empty', 'blocked'}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def rational(value):
    require(not isinstance(value, bool) and isinstance(value, (str, int)), 'Exact rational input required')
    return Q(value)


def canonical(data):
    return (json.dumps(data, sort_keys=True, separators=(',', ':'), ensure_ascii=False) + '\n').encode()


def digest(data):
    return hashlib.sha256(canonical(data)).hexdigest()


def box(raw):
    require(isinstance(raw, list) and len(raw) == 5, 'Five intervals required')
    out = []
    for pair in raw:
        require(isinstance(pair, list) and len(pair) == 2, 'Two interval endpoints required')
        lo, hi = map(rational, pair)
        require(Q(0) <= lo <= hi <= Q(1), 'Closed raw face box inside [0,1]^5 required')
        out.append((lo, hi))
    return out


def strings(bounds):
    return [[str(lo), str(hi)] for lo, hi in bounds]


def tighten(raw):
    b = box(raw)
    lower = [max(b[j][0] for j in range(i, 5)) for i in range(5)]
    upper = [min(b[j][1] for j in range(i + 1)) for i in range(5)]
    conflicts = [i for i in range(5) if lower[i] > upper[i]]
    return dict(box=None if conflicts else strings(zip(lower, upper)),
                tightened_lower=list(map(str, lower)), tightened_upper=list(map(str, upper)),
                empty=bool(conflicts), conflicting_coordinates=conflicts)


def midpoint(raw):
    return ['1'] + [str((lo + hi) / 2) for lo, hi in box(raw)]


def encode_center(center, *, max_ray_abs, max_rate):
    require(len(center) == 6 and rational(center[0]) == 1, 'Face center required')
    q = list(map(rational, center))
    require(all(Q(0) <= value <= Q(1) for value in q), 'Center outside positive face')
    den = math.lcm(*(value.denominator for value in q))
    nums = [int(value * den) for value in q]
    numerator_bound = max_ray_abs * sum(abs(value) for value in nums)
    denominator_bound = max_rate * den
    require(max_ray_abs > 0 and max_rate > 0, 'Pinned fan bounds required')
    require(max(nums) < 2**63 and den < 2**63 and numerator_bound < 2**63 and denominator_bound < 2**63,
            'C016 int64 eligibility exhausted; stop without rounding')
    require(numerator_bound * denominator_bound < 2**127, 'C016 comparator range exhausted')
    return dict(denominator=den, numerators=nums, numerator_bound=numerator_bound,
                denominator_bound=denominator_bound)


def split_box(raw):
    b = box(raw)
    widths = [hi - lo for lo, hi in b]
    require(max(widths) > 0, 'Singleton cannot be split')
    axis = widths.index(max(widths))
    seam = sum(b[axis]) / 2
    left, right = copy.deepcopy(b), copy.deepcopy(b)
    left[axis] = (b[axis][0], seam)
    right[axis] = (seam, b[axis][1])
    return dict(axis=axis, seam=str(seam), lower_raw_box=strings(left), upper_raw_box=strings(right),
                lower=tighten(strings(left)), upper=tighten(strings(right)), closed_seam=True)


def make_node(node_id, parent, depth, raw_box, cutoff_index):
    t = tighten(raw_box)
    return dict(id=node_id, parent=parent, depth=depth, raw_box=strings(box(raw_box)),
                box=t['box'], tightening=t, cutoff_index=cutoff_index,
                status='empty' if t['empty'] else 'queued', attempts=[], children=[], proof=None)


def new_state(input_manifest_sha256, *, test_only=False):
    require(re.fullmatch('[0-9a-f]{64}', input_manifest_sha256) is not None, 'Input manifest SHA256 required')
    root = make_node('r', None, 0, [['0', '1']] * 5, 0)
    return dict(schema='C063-unique-cover-state-v1', environment='TEST_ONLY' if test_only else 'PREPARE_ONLY',
                input_manifest_sha256=input_manifest_sha256, model=MODEL, target='1/65000',
                nodes={'r': root}, queue=['r'], inflight=None, passed_leaves=[], empty_leaves=[],
                halted=None, release=None, next_attempt_serial=1,
                batch=dict(id='first-batch-proposed', max_new_nodes=8, max_native_reservations=8,
                           started_nodes=[], native_reservations=0),
                production_kernel_calls=0, global_g_certificate=False, new_C6_bound=False)


def authorize(state, release):
    s = copy.deepcopy(state)
    require(s['environment'] != 'TEST_ONLY', 'Test state cannot receive production release')
    require(release.get('main_authorized') is True and release.get('F061_gradient_gate_accepted') is True,
            'Explicit main release and accepted F061 gradient gate required')
    require(release.get('input_manifest_sha256') == s['input_manifest_sha256'], 'Release input pins mismatch')
    require(release.get('timeout_seconds') == 180, '180-second runtime policy required')
    require(isinstance(release.get('F061_gradient_gate_sha256'), str) and
            re.fullmatch('[0-9a-f]{64}', release['F061_gradient_gate_sha256']) is not None,
            'Pinned accepted F061 gradient evidence required')
    require(isinstance(release.get('release_sha256'), str) and re.fullmatch('[0-9a-f]{64}', release['release_sha256']) is not None,
            'External immutable release hash required')
    s['release'] = copy.deepcopy(release)
    s['environment'] = 'MAIN_RELEASED'
    return s


def queue_sort(state):
    state['queue'].sort(key=lambda node_id: (state['nodes'][node_id]['depth'], node_id))


def expected_input(state, node):
    center = midpoint(node['box'])
    return dict(model=MODEL, model_definition=MODEL_DEFINITION, center=center,
                window=dict(kind='B', a=LADDER[node['cutoff_index']], T='1/2'),
                face_box=node['box'], target='1/65000', quadratic='frame32',
                input_manifest_sha256=state['input_manifest_sha256'])


def reserve(state):
    s = copy.deepcopy(state)
    require(s['environment'] in ('TEST_ONLY', 'MAIN_RELEASED'), 'Preparation does not authorize dispatch')
    require(s['halted'] is None and s['inflight'] is None, 'Halted or existing inflight work: no dispatch')
    require(bool(s['queue']), 'No queued node')
    node_id = s['queue'][0]
    node = s['nodes'][node_id]
    require(node['status'] == 'queued', 'Only a queued leaf may be reserved')
    batch = s['batch']
    require(batch['native_reservations'] < batch['max_native_reservations'], 'Batch native-reservation budget exhausted')
    if node_id not in batch['started_nodes']:
        require(len(batch['started_nodes']) < batch['max_new_nodes'], 'Batch new-node budget exhausted')
        batch['started_nodes'].append(node_id)
    expected = expected_input(s, node)
    serial = s['next_attempt_serial']
    token = f'{node_id}.a{node["cutoff_index"]}.n{serial:06d}'
    attempt = dict(token=token, input=expected, input_sha256=digest(expected),
                   state='RESERVED_BEFORE_ANY_LAUNCH', execution_handle=None, observation=None,
                   native_launch_observed=False, launch_outcome_unknown=False)
    node['attempts'].append(attempt)
    node['status'] = 'inflight'
    s['queue'].pop(0)
    s['inflight'] = dict(node_id=node_id, token=token)
    s['next_attempt_serial'] += 1
    batch['native_reservations'] += 1
    validate_state(s)
    return s, copy.deepcopy(attempt)


def attach_handle(state, token, handle):
    s = copy.deepcopy(state)
    require(s['inflight'] is not None and s['inflight']['token'] == token, 'Wrong inflight token')
    attempt = s['nodes'][s['inflight']['node_id']]['attempts'][-1]
    require(attempt['execution_handle'] is None, 'Original handle is immutable')
    require(isinstance(handle, dict) and bool(handle), 'Observed original handle required')
    attempt['execution_handle'] = copy.deepcopy(handle)
    attempt['state'] = 'ORIGINAL_HANDLE_RECORDED'
    if s['environment'] != 'TEST_ONLY' and isinstance(handle.get('native_pid'), int) and handle['native_pid'] > 0:
        attempt['native_launch_observed'] = True
        s['production_kernel_calls'] += 1
    return s


def interval(raw):
    require(isinstance(raw, list) and len(raw) == 2, 'Two exact interval endpoints required')
    lo, hi = map(rational, raw)
    require(lo <= hi, 'Reversed interval')
    return lo, hi


def validate_jet(jet, expected, *, test_only=False):
    for key in ('model', 'model_definition', 'window', 'center'):
        require(jet.get(key) == expected[key], f'Jet {key} mismatch')
    prov = jet.get('provenance', {})
    require(prov.get('kind') == 'C016_NATIVE_SAME_WINDOW_POSITIVE_BAND', 'Artificial/full-G/different source rejected')
    require(prov.get('precision_bits') == 512 and prov.get('full_group_terms') == 967680 and
            prov.get('computed_sign_pairs') == 483840, 'Native precision/group metadata mismatch')
    require(prov.get('synthetic_fixture') is test_only, 'Fixture cannot be production evidence')
    for name in ('native_sha256', 'receipt_sha256', 'gradient_gate_sha256'):
        require(isinstance(prov.get(name), str) and re.fullmatch('[0-9a-f]{64}', prov[name]) is not None,
                f'External authenticated {name} required')
    value = interval(jet['F_interval'])
    require(value[0] > 0, 'Positive band value lower endpoint required')
    require(len(jet.get('face_gradient_intervals', [])) == 5, 'Five same-window physical face gradients required')
    gradient = [interval(raw) for raw in jet['face_gradient_intervals']]
    require(isinstance(jet.get('source_proof'), str) and jet['source_proof'], 'External proof source required')
    return value, gradient


def certificate_input(jet, expected, *, test_only=False):
    validate_jet(jet, expected, test_only=test_only)
    return {**{k: expected[k] for k in ('model', 'model_definition', 'window', 'center', 'face_box', 'target', 'quadratic')},
            'jet': copy.deepcopy(jet)}


def check_saved_box(jet, expected, saved):
    """Exact arithmetic and corner binding only; does not authenticate exp itself."""
    require(saved.get('model') == MODEL and saved.get('window') == expected['window'] and
            saved.get('quadratic') == 'frame32' and saved.get('target') == '1/65000', 'Box model/window/target mismatch')
    require(saved.get('conditional_on_jet_source') == jet['source_proof'], 'Box proof source mismatch')
    corners = saved.get('corners', [])
    require(len(corners) == 32, 'Exactly 32 labelled corners required, including duplicates')
    require({tuple(row['corner_bits']) for row in corners} == set(product((0, 1), repeat=5)), 'Corner labels incomplete')
    center = list(map(rational, expected['center']))
    b = box(expected['face_box'])
    lo, hi = interval(jet['F_interval'])
    grad = [interval(raw) for raw in jet['face_gradient_intervals']]
    a = rational(expected['window']['a'])
    for row in corners:
        d = [b[i][bit] - center[i + 1] for i, bit in enumerate(row['corner_bits'])]
        lower_dot = sum(min(gl * di, gu * di) for (gl, gu), di in zip(grad, d))
        f = min(hi, max(lo, lower_dot))
        q = sum(di * di for di in d) / 32
        require(list(map(rational, row['displacement'])) == d, 'Wrong corner geometry')
        require(rational(row['gradient_dot_lower']) == lower_dot and rational(row['minimizing_center_value']) == f,
                'Wrong shared-denominator minimization')
        require(rational(row['quadratic']) == q and rational(row['exponent']) == lower_dot / f - q / a**2,
                'Wrong fixed curvature expression')
        interval([row['certificate_lower'], row['certificate_upper']])
    lower = min(rational(row['certificate_lower']) for row in corners)
    require(lower == rational(saved['tangent_bound_lower']), 'Wrong whole-box minimum')
    passed = lower >= TARGET
    require(saved.get('status') == ('PASS_CONDITIONAL_ON_MATCHING_CENTER_JET' if passed else 'UNRESOLVED'),
            'Wrong box verdict')
    return passed


def ingest(state, token, observation):
    s = copy.deepcopy(state)
    require(s['inflight'] is not None and s['inflight']['token'] == token, 'Stale/duplicate/wrong result token')
    node = s['nodes'][s['inflight']['node_id']]
    attempt = node['attempts'][-1]
    require(observation.get('input_sha256') == attempt['input_sha256'], 'Observation input hash mismatch')
    require(observation.get('receipt_exit_code') == 0 and observation.get('stderr_empty') is True,
            'Native failure must use halt_inflight, never a mathematical transition')
    expected = attempt['input']
    jet = observation['jet']
    (lo, hi), _ = validate_jet(jet, expected, test_only=s['environment'] == 'TEST_ONLY')
    if s['environment'] != 'TEST_ONLY':
        require(jet['provenance']['gradient_gate_sha256'] == s['release']['F061_gradient_gate_sha256'],
                'Gradient proof must be the main-accepted release evidence')
    attempt['observation'] = copy.deepcopy(observation)
    attempt['state'] = 'IMPORTED_ONCE'
    s['inflight'] = None
    if lo > TARGET:
        require('box_certificate' in observation, 'Center pass requires the full same-box certificate')
        require(observation.get('box_checker_sha256') == BOX_CHECKER_SHA256 and
                observation.get('box_certificate_sha256') == digest(observation['box_certificate']),
                'Pinned box checker and canonical certificate content hash required')
        passed = check_saved_box(jet, expected, observation['box_certificate'])
        if passed:
            node['status'] = 'passed'
            node['proof'] = dict(token=token, input_sha256=attempt['input_sha256'],
                                 certificate_sha256=observation['box_certificate_sha256'],
                                 native_sha256=jet['provenance']['native_sha256'])
            s['passed_leaves'].append(node['id'])
        else:
            if all(lo_b == hi_b for lo_b, hi_b in box(node['box'])):
                node['status'] = 'blocked'
                s['halted'] = dict(reason='UNRESOLVED_SINGLETON_BOX', node_id=node['id'], token=token)
            else:
                split = split_box(node['box'])
                node['status'] = 'split'
                node['split'] = split
                for suffix, raw in (('0', split['lower_raw_box']), ('1', split['upper_raw_box'])):
                    child = make_node(node['id'] + suffix, node['id'], node['depth'] + 1, raw, node['cutoff_index'])
                    require(child['id'] not in s['nodes'], 'Child overwrite forbidden')
                    s['nodes'][child['id']] = child
                    node['children'].append(child['id'])
                    (s['empty_leaves'] if child['status'] == 'empty' else s['queue']).append(child['id'])
    elif hi < TARGET:
        if node['cutoff_index'] + 1 < len(LADDER):
            node['cutoff_index'] += 1
            node['status'] = 'queued'
            s['queue'].append(node['id'])
        else:
            node['status'] = 'blocked'
            s['halted'] = dict(reason='UNRESOLVED_CUTOFF_LADDER_EXHAUSTED_NOT_FULL_G_COUNTEREXAMPLE', node_id=node['id'], token=token)
    else:
        node['status'] = 'blocked'
        s['halted'] = dict(reason='UNRESOLVED_CENTER_OVERLAP_OR_EQUALITY', node_id=node['id'], token=token)
    if s['environment'] != 'TEST_ONLY' and not attempt['native_launch_observed']:
        attempt['native_launch_observed'] = True
        s['production_kernel_calls'] += 1
    queue_sort(s)
    validate_state(s)
    return s


def halt_inflight(state, token, reason, retained_receipt_sha256, *, native_launch_observed=None):
    s = copy.deepcopy(state)
    require(s['inflight'] is not None and s['inflight']['token'] == token, 'Wrong inflight halt token')
    node = s['nodes'][s['inflight']['node_id']]
    node['status'] = 'blocked'
    node['attempts'][-1]['state'] = 'FAILED_NO_RETRY'
    node['attempts'][-1]['retained_receipt_sha256'] = retained_receipt_sha256
    attempt = node['attempts'][-1]
    if s['environment'] != 'TEST_ONLY':
        if native_launch_observed is True and not attempt['native_launch_observed']:
            attempt['native_launch_observed'] = True
            s['production_kernel_calls'] += 1
        elif native_launch_observed is None and not attempt['native_launch_observed']:
            attempt['launch_outcome_unknown'] = True
    s['halted'] = dict(reason=reason, node_id=node['id'], token=token)
    s['inflight'] = None
    validate_state(s)
    return s


def validate_state(state):
    require(state['schema'] == 'C063-unique-cover-state-v1', 'Wrong state schema')
    nodes = state['nodes']
    require('r' in nodes and nodes['r']['parent'] is None, 'Root missing')
    reached, todo = set(), ['r']
    while todo:
        node_id = todo.pop()
        require(node_id not in reached, 'Cycle or duplicate child')
        reached.add(node_id)
        node = nodes[node_id]
        require(node['id'] == node_id and node['status'] in STATUSES, 'Bad node identity/status')
        require(isinstance(node['cutoff_index'], int) and 0 <= node['cutoff_index'] < len(LADDER), 'Bad cutoff index')
        require(node['tightening'] == tighten(node['raw_box']) and node['box'] == node['tightening']['box'], 'Geometry changed')
        if node['status'] == 'split':
            expected = split_box(node['box'])
            require(node['split'] == expected and node['children'] == [node_id + '0', node_id + '1'], 'Split or closed seam changed')
            for child_id, raw in zip(node['children'], (expected['lower_raw_box'], expected['upper_raw_box'])):
                child = nodes[child_id]
                require(child['parent'] == node_id and child['depth'] == node['depth'] + 1 and child['raw_box'] == raw,
                        'Child does not cover its registered parent half')
                todo.append(child_id)
        else:
            require(not node['children'], 'Leaf has children')
        if node['status'] == 'empty':
            require(node['tightening']['empty'], 'False emptiness')
        if node['status'] == 'passed':
            require(node['proof'] is not None, 'Passed leaf proof missing')
        for attempt in node['attempts']:
            require(attempt['input_sha256'] == digest(attempt['input']), 'Attempt exact input changed')
            old_index = LADDER.index(attempt['input']['window']['a'])
            old_node = {**node, 'cutoff_index': old_index}
            require(attempt['input'] == expected_input(state, old_node), 'Attempt geometry/model/window binding changed')
    require(reached == set(nodes), 'Unreachable node in snapshot')
    queued = sorted((n['id'] for n in nodes.values() if n['status'] == 'queued'), key=lambda k: (nodes[k]['depth'], k))
    require(state['queue'] == queued and len(set(state['queue'])) == len(state['queue']), 'Queue omitted/duplicated/reordered leaf')
    require(set(state['passed_leaves']) == {n['id'] for n in nodes.values() if n['status'] == 'passed'} and
            len(set(state['passed_leaves'])) == len(state['passed_leaves']), 'Passed ledger mismatch')
    require(set(state['empty_leaves']) == {n['id'] for n in nodes.values() if n['status'] == 'empty'}, 'Empty ledger mismatch')
    active = [n for n in nodes.values() if n['status'] == 'inflight']
    require(len(active) <= 1, 'Only main may have one prepared inflight slot in this version')
    require((not active and state['inflight'] is None) or
            (len(active) == 1 and state['inflight'] == dict(node_id=active[0]['id'], token=active[0]['attempts'][-1]['token'])),
            'Inflight identity mismatch')
    require(state['global_g_certificate'] is False and state['new_C6_bound'] is False, 'Main alone may promote a final theorem')
    batch = state['batch']
    require(1 <= batch['max_new_nodes'] <= 8 and 1 <= batch['max_native_reservations'] <= 8, 'Unregistered first-batch budget')
    require(batch['native_reservations'] == sum(len(node['attempts']) for node in nodes.values()) and
            batch['native_reservations'] <= batch['max_native_reservations'], 'Reservation ledger mismatch')
    require(len(set(batch['started_nodes'])) == len(batch['started_nodes']) <= batch['max_new_nodes'] and
            set(batch['started_nodes']) == {node['id'] for node in nodes.values() if node['attempts']}, 'Started-node ledger mismatch')
    return dict(nonempty_frontier=[n['id'] for n in nodes.values() if n['status'] not in ('split', 'empty')],
                all_nonempty_leaves_passed=all(n['status'] in ('split', 'empty', 'passed') for n in nodes.values()),
                root_partition_preserved=True, global_g_certificate=False)


def recover(directory):
    directory = Path(directory)
    files = sorted(path for path in directory.iterdir() if re.fullmatch(r'snapshot_[0-9]{8}\.json', path.name))
    require(bool(files), 'No authoritative snapshot')
    previous = None
    latest = None
    for sequence, path in enumerate(files):
        require(path.name == f'snapshot_{sequence:08d}.json', 'Snapshot sequence gap')
        raw = path.read_bytes()
        record = json.loads(raw)
        require(raw == canonical(record), 'Noncanonical or partial snapshot')
        require(record['sequence'] == sequence and record['previous_sha256'] == previous, 'Broken snapshot hash chain')
        validate_state(record['state'])
        if latest is not None:
            require(record['state']['input_manifest_sha256'] == latest['state']['input_manifest_sha256'], 'Immutable input pin changed')
            for node_id in latest['state']['passed_leaves']:
                require(record['state']['nodes'][node_id] == latest['state']['nodes'][node_id], 'Previously passed leaf changed')
        previous = hashlib.sha256(raw).hexdigest()
        latest = record
    return dict(sequence=len(files)-1, sha256=previous, state=latest['state'],
                recovery_action='OBSERVE_ORIGINAL_HANDLE_OR_HALT_UNKNOWN_NEVER_RELAUNCH' if latest['state']['inflight'] else 'CONTINUE_FROM_UNIQUE_QUEUE_AFTER_MAIN_RELEASE',
                uncommitted_pending_files=[path.name for path in directory.iterdir() if path.name.endswith('.pending')])


def append_snapshot(directory, state, expected_previous_sha256):
    """One main writer, exclusive append, atomic hard-link publication, no replay."""
    directory = Path(directory)
    directory.mkdir(exist_ok=True)
    validate_state(state)
    with (directory / '.writer.lock').open('a+b') as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        existing = any(re.fullmatch(r'snapshot_[0-9]{8}\.json', path.name) for path in directory.iterdir())
        current = recover(directory) if existing else None
        require((current['sha256'] if current else None) == expected_previous_sha256, 'Stale writer predecessor hash')
        if current:
            require(current['state']['input_manifest_sha256'] == state['input_manifest_sha256'], 'Input pin changed')
            for node_id in current['state']['passed_leaves']:
                require(state['nodes'][node_id] == current['state']['nodes'][node_id], 'Passed leaf recomputation/change forbidden')
        sequence = current['sequence'] + 1 if current else 0
        record = dict(sequence=sequence, previous_sha256=expected_previous_sha256, state=state)
        raw = canonical(record)
        final = directory / f'snapshot_{sequence:08d}.json'
        pending = directory / f'snapshot_{sequence:08d}.{os.getpid()}.pending'
        with pending.open('xb') as stream:
            stream.write(raw)
            stream.flush()
            os.fsync(stream.fileno())
        os.link(pending, final)  # atomic, exclusive; an existing snapshot is never overwritten
        pending.unlink()
        fd = os.open(directory, os.O_RDONLY)
        try:
            os.fsync(fd)
        finally:
            os.close(fd)
        return hashlib.sha256(raw).hexdigest()
