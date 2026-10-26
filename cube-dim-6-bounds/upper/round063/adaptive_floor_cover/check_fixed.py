"""Exactly twelve predeclared pure-Fraction/tree/adapter configurations.

One invocation only; hard 180 s alarm; checkpoints retained on any failure.
All fabricated jets are explicitly TEST_ONLY, never observations of body32.
"""
import sys
sys.dont_write_bytecode = True
import copy
import hashlib
import json
import os
from datetime import datetime, timezone
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import signal
import time
import cover_core as C

P = Path(__file__).resolve().parent
PIN = '0' * 64


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rejects(function, *args, **kwargs):
    try:
        function(*args, **kwargs)
    except (ValueError, KeyError):
        return True
    raise AssertionError('Expected explicit rejection')


def state():
    return C.new_state(PIN, test_only=True)


def fixture_jet(attempt, lower, upper=None):
    upper = lower if upper is None else upper
    expected = attempt['input']
    return {**{key: expected[key] for key in ('model', 'model_definition', 'window', 'center')},
            'F_interval': [str(lower), str(upper)], 'face_gradient_intervals': [['0', '0']] * 5,
            'source_proof': 'SYNTHETIC_TEST_ONLY:constant-positive-function-not-body32',
            'provenance': dict(kind='C016_NATIVE_SAME_WINDOW_POSITIVE_BAND', precision_bits=512,
                               full_group_terms=967680, computed_sign_pairs=483840, synthetic_fixture=True,
                               native_sha256=PIN, receipt_sha256=PIN, gradient_gate_sha256=PIN)}


def fixture_box(attempt, jet):
    """Rational conservative enclosures for synthetic zero-gradient tangent data.

    For x>=0, exp(-x)>=4^-ceil(x), from e<4; upper bound is 1.
    This does not evaluate a production kernel or invoke the A057 checker.
    """
    expected = attempt['input']
    center = list(map(Q, expected['center']))
    bounds = C.box(expected['face_box'])
    a, f = Q(expected['window']['a']), Q(jet['F_interval'][0])
    corners = []
    for bits in product((0, 1), repeat=5):
        d = [bounds[i][bit] - center[i + 1] for i, bit in enumerate(bits)]
        q = sum(di * di for di in d) / 32
        loss = q / a**2
        ceiling = -((-loss.numerator) // loss.denominator)
        lower = f / (Q(4) ** ceiling)
        corners.append(dict(corner_bits=list(bits), displacement=list(map(str, d)),
                            gradient_dot_lower='0', minimizing_center_value=str(f),
                            quadratic=str(q), exponent=str(-loss),
                            certificate_lower=str(lower), certificate_upper=str(f),
                            passes_target=lower >= C.TARGET))
    lower = min(Q(row['certificate_lower']) for row in corners)
    return dict(status='PASS_CONDITIONAL_ON_MATCHING_CENTER_JET' if lower >= C.TARGET else 'UNRESOLVED',
                model=C.MODEL, window=expected['window'], quadratic='frame32', target='1/65000',
                conditional_on_jet_source=jet['source_proof'], tangent_bound_lower=str(lower),
                tangent_bound_upper=str(f), corners=corners, kernel_calls=0, global_g_certificate=False)


def observation(attempt, lower, upper=None):
    jet = fixture_jet(attempt, lower, upper)
    data = dict(input_sha256=attempt['input_sha256'], receipt_exit_code=0, stderr_empty=True, jet=jet)
    if lower > C.TARGET:
        certificate = fixture_box(attempt, jet)
        data.update(box_certificate=certificate, box_certificate_sha256=C.digest(certificate),
                    box_checker_sha256=C.BOX_CHECKER_SHA256)
    return data


def main():
    result_path = P / 'preparation_results.json'
    receipt_path = P / 'preparation_receipt.json'
    if result_path.exists() or receipt_path.exists():
        raise SystemExit('Existing preparation receipt/checkpoint: no automatic replay')
    tic = time.monotonic()
    receipt = dict(status='RUNNING_PREPARATION_ONLY', pid=os.getpid(), started_utc=datetime.now(timezone.utc).isoformat(),
                   command=[sys.executable, '-B', str(Path(__file__).resolve())], timeout_seconds=180,
                   new_native_calls=0, production_calls=0, compiles=0, retries=0)
    result = dict(status='RUNNING', configurations=[], production_calls=0, native_calls=0,
                  scalar_kernel_calls=0, compiles=0, symbolic_packages_invoked=0,
                  arithmetic='Python fractions.Fraction only', synthetic_jets_are_kernel_evidence=False,
                  global_g_certificate=False, new_C6_bound=False)
    def save():
        result['elapsed_seconds'] = time.monotonic() - tic
        result_path.write_text(json.dumps(result, indent=2) + '\n')
        receipt['wall_seconds'] = time.monotonic() - tic
        receipt_path.write_text(json.dumps(receipt, indent=2) + '\n')
    def record(number, topic, **details):
        assert number == len(result['configurations']) + 1 and number <= 12
        result['configurations'].append(dict(configuration=number, topic=topic, passed=True, **details))
        save()
    def timeout(signum, frame):
        raise TimeoutError('180-second preparation program limit; stop without replay')
    signal.signal(signal.SIGALRM, timeout)
    signal.setitimer(signal.ITIMER_REAL, 180)
    save()
    try:
        pins = json.loads((P / 'input_manifest.json').read_text())
        for path, expected_hash in pins['readonly_sha256'].items():
            assert sha(Path(path)) == expected_hash, path
        for name, expected_hash in pins['owned_source_sha256'].items():
            assert sha(P / name) == expected_hash, name
        receipt['input_manifest_sha256'] = sha(P / 'input_manifest.json')
        grid = [Q(0), Q(1, 2), Q(1)]
        cases = [
            [['0', '1']] * 5,
            [['0', '1'], ['0', '1/2'], ['1/2', '1'], ['0', '1'], ['0', '1/2']],
            [['1/2', '1'], ['0', '1'], ['0', '1/2'], ['1/2', '1/2'], ['0', '1/2']],
            [['0', '0'], ['0', '1'], ['0', '1'], ['0', '1'], ['0', '1']],
            [['0', '1/2'], ['1', '1'], ['0', '1'], ['0', '1'], ['0', '1']]
        ]
        points_checked = 0
        for raw in cases:
            t = C.tighten(raw)
            if not t['empty']:
                assert C.tighten(t['box'])['box'] == t['box']
            for point in product(grid, repeat=5):
                ordered = all(point[i] >= point[i + 1] for i in range(4))
                original = ordered and all(lo <= v <= hi for (lo, hi), v in zip(C.box(raw), point))
                tightened = not t['empty'] and ordered and all(lo <= v <= hi for (lo, hi), v in zip(C.box(t['box']), point))
                assert original == tightened
                points_checked += 1
        record(1, 'order-tightening exact finite-grid identity and idempotence', boxes=len(cases), grid_points_checked=points_checked)

        empty = C.tighten([['0', '1/2'], ['3/4', '1']] + [['0', '1']] * 3)
        assert empty['empty'] and empty['conflicting_coordinates']
        tie = C.tighten([['1/2', '1/2']] * 5)
        zero = C.tighten([['0', '0']] * 5)
        assert not tie['empty'] and not zero['empty']
        assert rejects(C.tighten, [['1', '0']] * 5)
        record(2, 'strict emptiness; tie and coordinate-zero boundary retention', conflicting_coordinates=empty['conflicting_coordinates'])

        center = C.midpoint([['0', '1']] * 5)
        fan_bounds = pins['fan_metadata']
        encoded = C.encode_center(center, max_ray_abs=fan_bounds['max_ray_abs'], max_rate=fan_bounds['max_rate'])
        assert center == ['1'] + ['1/2'] * 5 and encoded['numerators'] == [2, 1, 1, 1, 1, 1] and encoded['denominator'] == 2
        assert rejects(C.rational, 0.5) and rejects(C.rational, True)
        assert rejects(C.encode_center, ['1'] + ['1/' + str(2**64)] * 5,
                       max_ray_abs=fan_bounds['max_ray_abs'], max_rate=fan_bounds['max_rate'])
        record(3, 'midpoint and exact native-input encoding; no S2 rejection', root_center=center, denominator=encoded['denominator'])

        split = C.split_box([['0', '1']] * 5)
        assert split['axis'] == 0 and split['seam'] == '1/2' and split['closed_seam']
        grid = [Q(0), Q(1, 4), Q(1, 2), Q(3, 4), Q(1)]
        for point in product(grid, repeat=5):
            ordered = all(point[i] >= point[i + 1] for i in range(4))
            contained = []
            for child in (split['lower'], split['upper']):
                contained.append(not child['empty'] and ordered and all(lo <= v <= hi for (lo, hi), v in zip(C.box(child['box']), point)))
            assert ordered == any(contained)
            if ordered and point[0] == Q(1, 2):
                assert all(contained)
        record(4, 'longest-axis deterministic split and closed seam union', grid_points_checked=len(grid)**5)

        partial = [['1', '1'], ['1/2', '1/2'], ['0', '1/2'], ['0', '0'], ['0', '0']]
        assert C.split_box(partial)['axis'] == 2
        assert rejects(C.split_box, [['1/2', '1/2']] * 5)
        record(5, 'zero-width coordinates and singleton non-split')

        s, attempt = C.reserve(state())
        assert s['inflight']['token'] == attempt['token'] and s['queue'] == []
        assert rejects(C.reserve, s)
        s = C.attach_handle(s, attempt['token'], {'test_only_tool_handle': 'not-a-real-session'})
        assert rejects(C.attach_handle, s, attempt['token'], {'test_only_tool_handle': 'other'})
        assert rejects(C.reserve, C.new_state(PIN))
        record(6, 'unique queue, one inflight reservation and immutable original handle')

        s = state()
        indices = []
        for index in range(5):
            s, a = C.reserve(s)
            indices.append(s['nodes']['r']['cutoff_index'])
            s = C.ingest(s, a['token'], observation(a, C.TARGET / 2))
        assert indices == list(range(5)) and s['nodes']['r']['status'] == 'blocked'
        assert s['halted']['reason'].startswith('UNRESOLVED_CUTOFF_LADDER_EXHAUSTED')
        for bounds in [(C.TARGET / 2, C.TARGET * 2), (C.TARGET, C.TARGET)]:
            s, a = C.reserve(state())
            s = C.ingest(s, a['token'], observation(a, *bounds))
            assert s['nodes']['r']['cutoff_index'] == 0 and s['halted']['reason'] == 'UNRESOLVED_CENTER_OVERLAP_OR_EQUALITY'
        record(7, 'strict center comparisons and cutoff ladder; overlap/equality/last-cutoff stop', exact_ladder=list(C.LADDER))

        s = state()
        for _ in range(2):
            s, a = C.reserve(s)
            s = C.ingest(s, a['token'], observation(a, C.TARGET / 2))
        s, a = C.reserve(s)
        s = C.ingest(s, a['token'], observation(a, C.TARGET * 2))
        assert s['nodes']['r']['status'] == 'split' and s['queue'] == ['r0', 'r1']
        assert s['nodes']['r0']['cutoff_index'] == s['nodes']['r1']['cutoff_index'] == 2
        split_state = copy.deepcopy(s)
        record(8, 'center pass/whole-box failure splits, children inherit current a', child_ids=s['queue'], cutoff_index=2)

        s, a = C.reserve(state())
        obs = observation(a, C.TARGET * 2)
        jet = obs['jet']
        data = C.certificate_input(jet, a['input'], test_only=True)
        assert data['face_box'] == a['input']['face_box'] and data['jet'] == jet
        assert rejects(C.validate_jet, jet, a['input'])
        for key, value in [('model', 'different-model'), ('window', {'kind':'B','a':'1/250','T':'1/2'}),
                           ('center', ['1'] * 6), ('face_gradient_intervals', [['0','0']] * 4)]:
            wrong = copy.deepcopy(jet); wrong[key] = value
            assert rejects(C.validate_jet, wrong, a['input'], test_only=True)
        artificial = copy.deepcopy(jet); artificial['provenance']['kind'] = '052_ARTIFICIAL_FULL_G_JET'
        assert rejects(C.validate_jet, artificial, a['input'], test_only=True)
        wrong_box = copy.deepcopy(obs['box_certificate']); wrong_box['corners'][0]['displacement'][0] = '0'
        assert rejects(C.check_saved_box, jet, a['input'], wrong_box)
        record(9, 'same-model/window/center/box adapter; artificial and synthetic production evidence rejected')

        passed, a = C.reserve(state())
        passing_observation = observation(a, (Q(2)**100) * C.TARGET)
        passed = C.ingest(passed, a['token'], passing_observation)
        assert passed['passed_leaves'] == ['r'] and passed['queue'] == []
        assert C.validate_state(passed)['all_nonempty_leaves_passed'] and not passed['global_g_certificate']
        assert rejects(C.reserve, passed)
        assert rejects(C.ingest, passed, a['token'], passing_observation)
        record(10, 'passed-leaf no-recompute, exact frontier accounting and no global claim from fixtures')

        recovery_dir = P / 'test_snapshots' / 'recovery'
        recovery_dir.parent.mkdir(exist_ok=True)
        first_hash = C.append_snapshot(recovery_dir, state(), None)
        inflight, a = C.reserve(state())
        inflight = C.attach_handle(inflight, a['token'], {'fixture_only': 'original-test-handle-not-launched'})
        second_hash = C.append_snapshot(recovery_dir, inflight, first_hash)
        (recovery_dir / 'orphan.pending').write_text('uncommitted fixture; never dispatch\n')
        restored = C.recover(recovery_dir)
        assert restored['sha256'] == second_hash and restored['state'] == inflight
        assert restored['recovery_action'] == 'OBSERVE_ORIGINAL_HANDLE_OR_HALT_UNKNOWN_NEVER_RELAUNCH'
        assert restored['uncommitted_pending_files'] == ['orphan.pending'] and rejects(C.reserve, restored['state'])
        record(11, 'immutable snapshot hash chain; inflight recovery observes original handle only', snapshots=2)

        assert rejects(C.append_snapshot, recovery_dir, inflight, first_hash)
        tamper_dir = P / 'test_snapshots' / 'tampered'
        tamper_dir.mkdir()
        original = json.loads((recovery_dir / 'snapshot_00000000.json').read_text())
        original['previous_sha256'] = 'f' * 64
        (tamper_dir / 'snapshot_00000000.json').write_bytes(C.canonical(original))
        assert rejects(C.recover, tamper_dir)
        corrupt = copy.deepcopy(inflight)
        corrupt['nodes']['r']['attempts'][0]['input']['center'][1] = '0'
        assert rejects(C.validate_state, corrupt)
        passed_dir = P / 'test_snapshots' / 'passed'
        passed_hash = C.append_snapshot(passed_dir, passed, None)
        changed = copy.deepcopy(passed); changed['nodes']['r']['proof']['native_sha256'] = '1' * 64
        assert rejects(C.append_snapshot, passed_dir, changed, passed_hash)
        budget = state(); budget['batch']['max_native_reservations'] = 1
        budget, a = C.reserve(budget)
        budget = C.ingest(budget, a['token'], observation(a, C.TARGET / 2))
        assert rejects(C.reserve, budget) and budget['queue'] == ['r']
        split_state['batch']['max_new_nodes'] = 1
        assert rejects(C.reserve, split_state)
        record(12, 'tamper/stale/duplicate rejection and finite reservation budgets', original_state_native_budget=8, new_node_budget=8)

        initial = C.new_state(sha(P / 'input_manifest.json'))
        initial_hash = C.append_snapshot(P / 'state', initial, None)
        result.update(status='PASS_TWELVE_FRACTION_TREE_ADAPTER_CONFIGURATIONS_ZERO_PRODUCTION',
                      initial_snapshot_sha256=initial_hash, initial_queue=['r'], initial_inflight=None,
                      old_061_seed_recalculations=0, Python=sys.version)
        receipt.update(status='PASS_PREPARATION_ONLY_NO_DISPATCH', exit_code=0,
                       ended_utc=datetime.now(timezone.utc).isoformat())
        save()
        print(result['status'])
    except BaseException as error:
        signal.setitimer(signal.ITIMER_REAL, 0)
        result.update(status='FAILED_STOP_NO_AUTOMATIC_RETRY', error=repr(error))
        receipt.update(status='FAILED_STOP_NO_AUTOMATIC_RETRY', exit_code=1, error=repr(error),
                       ended_utc=datetime.now(timezone.utc).isoformat())
        save()
        raise
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)


if __name__ == '__main__':
    main()
