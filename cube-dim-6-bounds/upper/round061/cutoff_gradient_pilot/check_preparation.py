"""Twelve registered finite saved/symbolic checks. No production evaluator."""
import sys
sys.dont_write_bytecode = True
import copy
import hashlib
import json
import platform
import time
from fractions import Fraction as Q
from pathlib import Path
import flint
from flint import arb, ctx
import sympy as sp
from jet_adapter import MODEL, MODEL_DEFINITION, box_inputs, build_jet, dyadic, proposed_command, registered_case

P = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def contained(native, reference):
    lo, hi = dyadic(native)
    rl, ru = dyadic(reference)
    return lo <= rl <= ru <= hi


def five_distances(center):
    answer = []
    for k in range(1, 6):
        answer.append(sum((v - sum(center[:k]) / k) ** 2 for v in center[:k])
                      + sum((v - sum(center[k:]) / (6 - k)) ** 2 for v in center[k:]))
    return answer


def main():
    result_path = P / 'preparation_results.json'
    if result_path.exists():
        raise SystemExit('Refuse another formal batch; retain the existing receipt')
    tic = time.monotonic()
    registration = json.loads((P / 'registration.json').read_text())
    pins = json.loads((P / 'input_manifest.json').read_text())
    result = dict(status='RUNNING_REGISTERED_PREPARATION_ONLY', configurations=[],
                  new_production_kernel_calls=0, new_native_invocations=0, new_compiles=0,
                  automatic_compile_fixes=0, old_seed_recalculations=0, new_C6_bound=False,
                  global_floor_1_65000_certified=False, A057_independently_accepted=False)

    def save():
        result['elapsed_seconds'] = time.monotonic() - tic
        result_path.write_text(json.dumps(result, indent=2) + '\n')

    def record(number, **data):
        assert number == len(result['configurations']) + 1
        assert number <= registration['formal_configuration_cap'] == 12
        result['configurations'].append(dict(configuration=number, pass_check=True, **data))
        save()

    save()
    try:
        assert sha(P / 'registration.json') == pins['registration_sha256']
        for path, digest in pins['source_sha256'].items():
            assert sha(Path(path)) == digest, ('source pin', path)
        engine = Path(registration['source_path']).parent
        old_registration = json.loads((engine / 'registration.json').read_text())
        native_anchors = json.loads((engine / 'anchors_native.json').read_text())
        references = json.loads((engine / 'anchors_validation.json').read_text())
        assert len(native_anchors['anchors']) == len(references['cases']) == 6
        assert references['native_sha256'] == sha(engine / 'anchors_native.json')
        for index, (native, reference) in enumerate(zip(native_anchors['anchors'], references['cases'])):
            assert contained(native['value'], reference['band_DD5_reference'])
            assert all(contained(x, y) for x, y in zip(native['partials'], reference['partials_DD6_reference']))
            for label, exact in reference['density_exact_rationals'].items():
                lower, upper = dyadic(native[label])
                assert lower <= Q(exact) <= upper
            if index == 0:
                assert all(dyadic(raw) == (Q(0), Q(0)) for raw in native['partials'])
            record(index + 1, type='saved C016 scalar/Hermite evidence, not a new scalar run',
                   old_nodes=old_registration['anchors'][index], DD5_and_all_six_DD6_contained=True,
                   old_two_density_exact_rationals_contained=True)

        _, point, _, center_raw = registered_case('S3-primary')
        center = list(map(Q, center_raw))
        columns = [[1] * 6 for _ in range(3)]
        columns[1][4] = -1
        columns[2][5] = -1
        weights = [Q(3, 10), Q(1, 5), Q(1, 2)]
        assert [sum(weights[j] * columns[j][i] for j in range(3)) for i in range(6)] == center
        assert sp.Matrix.hstack(*map(sp.Matrix, columns)).rank() == 3
        distances = five_distances(center)
        assert min(distances) == Q(16, 125) > Q(1, 125000) ** 2
        record(7, type='prescribed genuine S3, outside paid S2 tube', center=center_raw,
               three_body_columns=columns, nonzero_coefficients=list(map(str, weights)), rank=3,
               five_S2_squared_distances=list(map(str, distances)), min_distance_squared=str(min(distances)))

        _, low, _, low_center_raw = registered_case('low-primary')
        low_center = list(map(Q, low_center_raw))
        old_low_path = next(Path(path) for path in pins['source_sha256'] if path.endswith('round038/green_floor_diagnostic/results.json'))
        old_low = json.loads(old_low_path.read_text())['points'][2]
        assert old_low['numerators'] == low['point_numerators'] and old_low['denominator'] == low['denominator']
        assert low_center[0] == 1 and all(Q(0) < low_center[i] < low_center[i - 1] for i in range(1, 6))
        distances = five_distances(low_center)
        assert min(distances) > Q(1, 125000) ** 2
        tokens = iter(map(int, Path(registration['fan_path']).read_text().split()))
        assert next(tokens) == 21
        rays = []
        for _ in range(21):
            assert next(tokens) > 0
            for _ in range(6):
                q = next(tokens)
                ray = [next(tokens) for _ in range(6)]
                rays.append((q, ray))
        assert next(tokens, None) is None
        max_ray = max(abs(v) for _, ray in rays for v in ray)
        max_q = max(q for q, _ in rays)
        largest_numerator = max_ray * max(sum(abs(v) for v in point['point_numerators']), sum(abs(v) for v in low['point_numerators']))
        largest_denominator = max_q * max(point['denominator'], low['denominator'])
        assert largest_numerator < 2 ** 63 and largest_denominator < 2 ** 63
        assert largest_numerator * largest_denominator < 2 ** 127
        record(8, type='ordinary low-point identity/domain and all-node integer-range bound', center=low_center_raw,
               source=str(old_low_path), five_S2_squared_distances=list(map(str, distances)),
               min_distance_squared=str(min(distances)), numerator_bound=largest_numerator,
               denominator_bound=largest_denominator, comparator_product_bound=largest_numerator * largest_denominator)

        x, a, z = sp.symbols('x a z', real=True)
        b = sp.symbols('b0:6', real=True)
        denominator = a * a + x * x
        g = [a / denominator]
        for k in range(1, 6):
            g.append(sp.factor(-(2 * x * g[k - 1] + (g[k - 2] if k >= 2 else 0)) / denominator))
        assert all(sp.simplify(sp.diff(sp.atan(x / a), x, k) / sp.factorial(k) - g[k - 1] / k) == 0 for k in range(1, 7))
        rational_kernel = 1 / sp.prod(z - bi for bi in b)
        assert sp.simplify(sp.diff(rational_kernel, b[0]) - rational_kernel / (z - b[0])) == 0
        record(9, type='symbolic normalized atan jets1..6 and labelled copied-node contour identity',
               derivative_orders=list(range(1, 7)), contour_inserts_one_copy_without_extra_factor=True)

        ctx.prec = 512
        for window in [registration['window'], registration['predeclared_fallback_window']]:
            a_value, T_value = Q(window['a']), Q(window['T'])
            assert 0 < a_value < T_value
            assert all(sp.simplify((g[k - 1] / k).subs(x, 0).subs(a, sp.Rational(a_value.numerator, a_value.denominator))
                                   - (0 if k % 2 == 0 else sp.Rational((-1) ** ((k - 1) // 2), k) / sp.Rational(a_value.numerator, a_value.denominator) ** k)) == 0 for k in range(1, 7))
        saved_time_path = next(Path(path) for path in pins['source_sha256'] if path.endswith('round007/g_time_moments/results.json'))
        rows = json.loads(saved_time_path.read_text())['time_nodes']
        old_a = arb(next(row['H_arb384'] for row in rows if row['tau'] == '1/250'))
        old_T = arb(next(row['H_arb384'] for row in rows if row['tau'] == '16/125'))
        old_control_band = old_a - old_T
        control_lower = Q(str(old_control_band.lower().fmpq()))
        assert control_lower > Q(1, 65000)
        record(10, type='positive windows/allzero normalized jets; saved old-control prerequisite only',
               windows=[registration['window'], registration['predeclared_fallback_window']],
               old_control_band_window={'kind': 'B', 'a': '1/250', 'T': '16/125'},
               old_control_band_lower_exact=str(control_lower), primary_old_control_above_target_by_positive_window_inclusion=True,
               pilot_low_point_prerequisite_uncomputed=True)

        def raw(m, exponent=0):
            return dict(lower_mantissa=str(m), upper_mantissa=str(m), common_binary_exponent=str(exponent))

        _, _, window, center = registered_case('low-primary')
        native = dict(status='all_real_nodes_supported', precision_bits=512, full_group_terms=967680,
                      computed_sign_pairs=483840, band=raw(2, -16), gradient=[raw(i, -20) for i in range(6)],
                      Euler_xgrad_plus5B_plusTpT_minusapa=raw(0))
        receipt = dict(model=MODEL, model_definition=MODEL_DEFINITION, window=window, center=center,
                       command=proposed_command('low-primary'), native_exit_code=0, stderr_empty=True, synthetic_fixture=True)
        jet = build_jet(native, receipt, 'low-primary', 'SYNTHETIC_TEST_ONLY:metadata-adapter', allow_synthetic=True)
        assert jet['face_gradient_intervals'] == [[str(Q(i, 2 ** 20))] * 2 for i in range(1, 6)]
        boxes = box_inputs(jet)
        assert len(boxes) == 8 and [Q(row['halfwidth_exact']) for row in boxes] == [Q(1, 250 * 2 ** k) for k in range(1, 9)]
        assert all(row['jet'] is jet and row['window'] == window for row in boxes)
        assert all(Q(pair[0]) <= Q(ci) <= Q(pair[1]) for row in boxes for pair, ci in zip(row['face_box'], center[1:]))
        record(11, type='synthetic adapter valid metadata/physical face extraction/eight dyadic boxes',
               synthetic_only=True, authentic_kernel_evidence=False, kernel_recalculations=0,
               halfwidths=[row['halfwidth_exact'] for row in boxes])

        wrong = copy.deepcopy(receipt)
        wrong['window']['a'] = '1/1000'
        rejected = False
        try:
            build_jet(native, wrong, 'low-primary', 'SYNTHETIC_TEST_ONLY:wrong-window', allow_synthetic=True)
        except ValueError:
            rejected = True
        assert rejected
        record(12, type='synthetic wrong-window rejection', synthetic_only=True, rejected=True)

        source = Path(registration['source_path']).read_text()
        build = json.loads((engine / 'build_record.json').read_text())
        assert sha(Path(registration['source_path'])) == build['source_sha256']
        assert sha(Path(registration['binary_path'])) == build['binary_sha256']
        assert 'sizeof(ulong)==8&&sizeof(fmpz)==8&&sizeof(arb_struct)==48' in source
        assert 'std::string(flint_version)!="3.6.0"' in source
        assert flint.__FLINT_VERSION__ == '3.6.0'
        result['inherited_ABI'] = dict(source_binary_and_runtime_pins_match=True,
                                       compile_time_sizes_from_unchanged_C016=dict(ulong=8, fmpz=8, arb_struct=48),
                                       native_runtime_version_gate='3.6.0', current_python_flint=flint.__version__,
                                       current_python_FLINT=flint.__FLINT_VERSION__, new_ABI_native_invocations=0,
                                       binary_sha256=build['binary_sha256'], source_sha256=build['source_sha256'])
        result.update(status='PASS_TWELVE_PREPARATION_CONFIGURATIONS_ZERO_PRODUCTION_KERNELS',
                      python=sys.version, sympy=sp.__version__, platform=platform.platform(),
                      registration_sha256=sha(P / 'registration.json'), input_manifest_sha256=sha(P / 'input_manifest.json'),
                      source_sha256={name: sha(P / name) for name in ['check_preparation.py', 'jet_adapter.py']})
        save()
    except BaseException as error:
        result.update(status='FAILED_STOP_NO_AUTOMATIC_RETRY', error=repr(error))
        save()
        raise
    print(result['status'])


if __name__ == '__main__':
    main()
