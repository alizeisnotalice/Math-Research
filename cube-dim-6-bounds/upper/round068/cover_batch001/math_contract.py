"""Native/jet/whole-box reconstruction for a single reserved C068 input.

This module never launches a kernel. It reads immutable original artifacts,
checks all six physical fields and exact Euler, and recomputes the 32-corner
Arb certificate with the pinned A057 checker. Status strings alone are unused.
"""
import hashlib
import json
import math
from fractions import Fraction as Q
from pathlib import Path
import re


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def integer(text):
    require(isinstance(text, str) and re.fullmatch(r'0|-?[1-9][0-9]*', text) is not None, 'Canonical integer dyadic field required')
    return int(text)


def dyadic(raw):
    require(isinstance(raw, dict) and set(raw) == {'lower_mantissa', 'upper_mantissa', 'common_binary_exponent'},
            'Exact three-field native dyadic endpoint object required')
    lo, hi, exponent = (integer(raw[key]) for key in ('lower_mantissa', 'upper_mantissa', 'common_binary_exponent'))
    require(lo <= hi and abs(exponent) <= 100000, 'Invalid/unsupported dyadic interval')
    factor = Q(2) ** exponent
    return lo * factor, hi * factor


def native_command(context, expected):
    center = list(map(Q, expected['center']))
    encoded = context.base.encode_center(expected['center'], **context.fan_bounds)
    a, T = Q(expected['window']['a']), Q(expected['window']['T'])
    require(expected['window']['kind'] == 'B' and T == Q(1, 2) and str(a) in context.base.LADDER, 'Unregistered window')
    return [context.paths['binary'], 'full', context.paths['fan'],
            str(a.numerator), str(a.denominator), str(T.numerator), str(T.denominator),
            str(encoded['denominator']), *map(str, encoded['numerators'])]


def validate_native(context, expected, token, native, receipt, native_bytes_sha256, stderr_bytes_sha256):
    require(receipt.get('schema') == 'C068-native-receipt-v1', 'Wrong native receipt schema')
    require(receipt.get('token') == token and receipt.get('input_sha256') == context.base.digest(expected), 'Native attempt identity mismatch')
    require(receipt.get('command') == native_command(context, expected), 'Native command/parameters mismatch')
    for key in ('model', 'model_definition', 'window', 'center'):
        require(receipt.get(key) == expected[key], f'Native receipt {key} mismatch')
    require(receipt.get('synthetic_fixture') is context.test_only, 'Fixture cannot enter a production result')
    require(receipt.get('native_exit_code') == 0 and receipt.get('stderr_empty') is True, 'Native failure is not mathematical evidence')
    require(receipt.get('native_sha256') == native_bytes_sha256 and receipt.get('stderr_sha256') == stderr_bytes_sha256,
            'Original native output bytes changed')
    require(stderr_bytes_sha256 == hashlib.sha256(b'').hexdigest(), 'Native stderr must be empty')
    require(receipt.get('source_sha256') == context.external_sha(context.paths['source']) and
            receipt.get('binary_sha256') == context.external_sha(context.paths['binary']) and
            receipt.get('fan_sha256') == context.external_sha(context.paths['fan']), 'Source/binary/fan receipt pins mismatch')
    require(receipt.get('runtime_sha256') == context.runtime_pins, 'Native runtime receipt pins mismatch')
    require(receipt.get('timeout_seconds') == 180, 'Native wall-clock policy changed')
    require(native.get('status') == 'all_real_nodes_supported' and native.get('precision_bits') == 512,
            'Wrong native engine/domain/precision')
    require(native.get('full_group_terms') == 967680 and native.get('computed_sign_pairs') == 483840,
            'Incomplete body32 group/sign-pair enumeration')
    for field in ('confluent_terms', 'zero_node_terms'):
        value = native.get(field)
        require(type(value) is int and 0 <= value <= 967680 and value % 2 == 0, f'Invalid bound/zero/repeat counter: {field}')
    require(type(native.get('unique_node_values')) is int and 1 <= native['unique_node_values'] <= 6 * 483840,
            'Invalid distinct-node count')
    require(isinstance(native.get('kernel_seconds'), (int, float)) and not isinstance(native['kernel_seconds'], bool)
            and math.isfinite(native['kernel_seconds']) and 0 <= native['kernel_seconds'] <= 180,
            'Invalid native timing receipt')
    require(isinstance(native.get('gradient'), list) and len(native['gradient']) == 6, 'All six physical gradients required')
    F = dyadic(native['band'])
    gradients = [dyadic(raw) for raw in native['gradient']]
    pa, pT = dyadic(native['p_a']), dyadic(native['p_T'])
    native_euler = dyadic(native['Euler_xgrad_plus5B_plusTpT_minusapa'])
    require(F[0] > 0 and pa[0] > 0 and pT[0] > 0, 'Positive band/endpoints unresolved')
    require(native_euler[0] <= 0 <= native_euler[1], 'Native Euler interval misses zero')
    def scaled(bounds, factor):
        return tuple(sorted(factor * value for value in bounds))
    terms = [scaled(F, Q(5)), scaled(pT, Q(expected['window']['T'])), scaled(pa, -Q(expected['window']['a']))]
    terms.extend(scaled(bounds, Q(x)) for bounds, x in zip(gradients, expected['center']))
    euler = tuple(sum(interval[index] for interval in terms) for index in (0, 1))
    require(euler[0] <= 0 <= euler[1], 'Independent exact-rational Euler reconstruction misses zero')
    return dict(F=F, physical_gradients=gradients, p_a=pa, p_T=pT, independent_Euler=euler,
                native_Euler=native_euler,
                counts={name: native[name] for name in ('full_group_terms', 'computed_sign_pairs', 'confluent_terms',
                                                       'zero_node_terms', 'unique_node_values')})


def reconstruct(context, state, token, folder):
    """Re-read original bytes and derive the only admissible imported observation."""
    require(state['inflight'] is not None and state['inflight']['token'] == token, 'No matching reserved original attempt')
    node = state['nodes'][state['inflight']['node_id']]
    attempt = node['attempts'][-1]
    expected = attempt['input']
    folder = Path(folder)
    require(folder.resolve() == context.attempt_directory(token).resolve(), 'Attempt artifact directory mismatch')
    receipt_path, native_path, stderr_path = folder / 'native_receipt.json', folder / 'native.json', folder / 'stderr.txt'
    receipt, native = read(receipt_path), read(native_path)
    intent = read(folder / 'intent.json')
    require(intent['token'] == token and intent['input'] == expected and intent['input_sha256'] == attempt['input_sha256'],
            'Pre-launch immutable intent does not match reservation')
    require(receipt['intent_sha256'] == sha(folder / 'intent.json'), 'Original launch intent changed')
    if not context.test_only:
        require(receipt['release_sha256'] == state['release']['release_sha256'] and
                intent['release_sha256'] == state['release']['release_sha256'], 'Native release mismatch')
        handle = attempt['execution_handle']
        require(handle is not None and handle['native_pid'] == receipt['owned_native_pid'] and
                handle['wrapper_pid'] == receipt['owned_wrapper_pid'], 'Original observed process handle mismatch')
        require(receipt['original_process_sha256'] == sha(folder/'original_process.json'), 'Native process receipt bytes changed')
        original = read(folder/'original_process.json')
        require(original['token'] == token and original['input_sha256'] == attempt['input_sha256'] and
                original['owned_native_pid'] == receipt['owned_native_pid'] and
                original['owned_wrapper_pid'] == receipt['owned_wrapper_pid'], 'Native process/intent binding mismatch')
    checked = validate_native(context, expected, token, native, receipt, sha(native_path), sha(stderr_path))
    source_proof = f'C068_ORIGINAL_NATIVE:{native_path}#sha256={sha(native_path)};receipt={sha(receipt_path)}'
    jet = {**{key: expected[key] for key in ('model', 'model_definition', 'window', 'center')},
           'F_interval': list(map(str, checked['F'])),
           'face_gradient_intervals': [list(map(str, pair)) for pair in checked['physical_gradients'][1:]],
           'source_proof': source_proof,
           'provenance': dict(kind='C016_NATIVE_SAME_WINDOW_POSITIVE_BAND', precision_bits=512,
                              full_group_terms=967680, computed_sign_pairs=483840, synthetic_fixture=context.test_only,
                              native_sha256=sha(native_path), receipt_sha256=sha(receipt_path),
                              gradient_gate_sha256=context.f061_manifest_sha256,
                              independent_full_gradient_crosschecked_for_this_new_attempt=False)}
    context.base.validate_jet(jet, expected, test_only=context.test_only)
    observation = dict(input_sha256=attempt['input_sha256'], receipt_exit_code=0, stderr_empty=True, jet=jet)
    if checked['F'][0] > context.base.TARGET:
        data = context.base.certificate_input(jet, expected, test_only=context.test_only)
        certificate = context.checker.certify(data, precision_bits=384)
        context.base.check_saved_box(jet, expected, certificate)
        observation.update(box_certificate=certificate, box_certificate_sha256=context.base.digest(certificate),
                           box_checker_sha256=context.base.BOX_CHECKER_SHA256)
    verdict = ('CENTER_STRICTLY_ABOVE' if checked['F'][0] > context.base.TARGET else
               'CENTER_STRICTLY_BELOW' if checked['F'][1] < context.base.TARGET else 'CENTER_OVERLAP_OR_EQUALITY_STOP')
    return dict(schema='C068-reconstructed-native-and-box-v1', token=token, input_sha256=attempt['input_sha256'],
                native_artifact_sha256=sha(native_path), receipt_artifact_sha256=sha(receipt_path),
                stderr_artifact_sha256=sha(stderr_path), intent_artifact_sha256=sha(folder / 'intent.json'),
                source_sha256=receipt['source_sha256'], binary_sha256=receipt['binary_sha256'], fan_sha256=receipt['fan_sha256'],
                complete_physical_gradient_intervals=[list(map(str, pair)) for pair in checked['physical_gradients']],
                endpoint_intervals={'p_a':list(map(str, checked['p_a'])), 'p_T':list(map(str, checked['p_T']))},
                independent_Euler_interval=list(map(str, checked['independent_Euler'])),
                native_Euler_interval=list(map(str, checked['native_Euler'])), native_counts=checked['counts'],
                independent_zero_repeat_count_recomputation=False,
                center_classification=verdict, observation=observation,
                exp_reconstructed_with_pinned_A057_checker=checked['F'][0] > context.base.TARGET,
                synthetic_fixture=context.test_only, universal_gradient_crosscheck_claim=False,
                global_g_certificate=False, new_C6_bound=False)
