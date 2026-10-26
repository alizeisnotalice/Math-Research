"""C061 metadata/interval adapter only. No kernel, process launch or compilation."""
import json
from fractions import Fraction as Q
from pathlib import Path

P = Path(__file__).resolve().parent
MODEL = 'body32_unit_cauchy_v1'
MODEL_DEFINITION = 'V={(1,s2,...,s6):si=+-1}; f_t(z)=t/(pi*(t^2+z^2)); 32 independent unit-weight factors'
CASE_IDS = ('low-primary', 'low-fallback', 'S3-primary', 'S3-fallback')


def registered_case(case_id):
    if case_id not in CASE_IDS:
        raise ValueError('Unregistered point/window')
    registration = json.loads((P / 'registration.json').read_text())
    point_id = 'A038-low3' if case_id.startswith('low-') else 'S3'
    point = next(row for row in registration['future_points'] if row['id'] == point_id)
    window = registration['window'] if case_id.endswith('-primary') else registration['predeclared_fallback_window']
    center = [str(Q(v, point['denominator'])) for v in point['point_numerators']]
    return registration, point, window, center


def proposed_command(case_id):
    """Return reviewable argv; this function never executes it."""
    registration, point, window, _ = registered_case(case_id)
    a, T = Q(window['a']), Q(window['T'])
    return [registration['binary_path'], 'full', registration['fan_path'],
            str(a.numerator), str(a.denominator), str(T.numerator), str(T.denominator),
            str(point['denominator']), *map(str, point['point_numerators'])]


def dyadic(raw):
    if set(raw) != {'lower_mantissa', 'upper_mantissa', 'common_binary_exponent'}:
        raise ValueError('Exact native dyadic endpoints required')
    factor = Q(2) ** int(raw['common_binary_exponent'])
    bounds = (int(raw['lower_mantissa']) * factor, int(raw['upper_mantissa']) * factor)
    if bounds[0] > bounds[1]:
        raise ValueError('Reversed interval')
    return bounds


def build_jet(native, receipt, case_id, source_proof, *, allow_synthetic=False):
    """Convert a same-window saved native result; provenance is a separate gate."""
    _, _, window, center = registered_case(case_id)
    expected = dict(model=MODEL, model_definition=MODEL_DEFINITION, window=window, center=center)
    if any(receipt.get(key) != value for key, value in expected.items()):
        raise ValueError('Saved result metadata does not match this point/window/model')
    synthetic = receipt.get('synthetic_fixture', False)
    if synthetic and not allow_synthetic:
        raise ValueError('Synthetic jet is not a production result')
    if receipt.get('command') != proposed_command(case_id):
        raise ValueError('Native command does not match the registered parameters')
    if receipt.get('native_exit_code') != 0 or receipt.get('stderr_empty') is not True:
        raise ValueError('Failed native cannot supply a jet')
    if native.get('status') != 'all_real_nodes_supported' or native.get('precision_bits') != 512:
        raise ValueError('Wrong result type or precision')
    if native.get('full_group_terms') != 967680 or native.get('computed_sign_pairs') != 483840:
        raise ValueError('Incomplete body32 group')
    if len(native.get('gradient', [])) != 6:
        raise ValueError('Six physical gradients required before selecting the five face fields')
    value = dyadic(native['band'])
    gradient = [dyadic(raw) for raw in native['gradient']]
    euler = dyadic(native['Euler_xgrad_plus5B_plusTpT_minusapa'])
    if not euler[0] <= 0 <= euler[1]:
        raise ValueError('Saved Euler interval misses zero')
    if not isinstance(source_proof, str) or not source_proof:
        raise ValueError('An external verified result source is required')
    if synthetic and not source_proof.startswith('SYNTHETIC_TEST_ONLY:'):
        raise ValueError('Synthetic provenance must remain explicit')
    return dict(**expected, F_interval=list(map(str, value)),
                face_gradient_intervals=[list(map(str, interval)) for interval in gradient[1:]],
                source_proof=source_proof, synthetic_fixture=synthetic)


def classification(jet):
    lower, upper = map(Q, jet['F_interval'])
    target = Q(1, 65000)
    if lower > target:
        return 'CENTER_STRICTLY_ABOVE_TARGET'
    if upper < target:
        return 'CENTER_STRICTLY_BELOW_TARGET'
    return 'UNRESOLVED_STOP_NO_PRECISION_INCREASE'


def box_inputs(jet):
    """Eight predetermined nested face boxes, reusing this one matching jet."""
    if classification(jet) != 'CENTER_STRICTLY_ABOVE_TARGET':
        raise ValueError('Center must be strictly above the target before any box test')
    center = list(map(Q, jet['center']))
    a = Q(jet['window']['a'])
    if center[0] != 1 or any(v < 0 or v > 1 for v in center[1:]):
        raise ValueError('Registered positive face point required')
    return [dict(model=MODEL, model_definition=MODEL_DEFINITION, window=jet['window'],
                 center=jet['center'], target='1/65000', quadratic='frame32', jet=jet,
                 halfwidth_exact=str(a / (2 ** k)), halfwidth_index=k,
                 face_box=[[str(max(Q(0), ci - a / (2 ** k))),
                            str(min(Q(1), ci + a / (2 ** k)))] for ci in center[1:]])
            for k in range(1, 9)]
