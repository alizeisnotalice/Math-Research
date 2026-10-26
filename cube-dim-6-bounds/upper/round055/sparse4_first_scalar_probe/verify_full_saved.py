"""Read-only verification of C055's saved sole native result; never launches."""
import sys
sys.dont_write_bytecode = True
import hashlib
import json
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
from fractions import Fraction as Q
from pathlib import Path

P = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(name):
    return json.loads((P / name).read_text())


def bounds(value):
    factor = Q(2) ** int(value['common_binary_exponent'])
    result = (int(value['lower_mantissa']) * factor,
              int(value['upper_mantissa']) * factor)
    assert result[0] <= result[1]
    return result


def decimal_bound(value, rounding):
    with localcontext() as context:
        context.prec = 200
        context.rounding = rounding
        value = Decimal(value.numerator) / Decimal(value.denominator)
        return format(value.quantize(Decimal('1e-50'), rounding=rounding), 'f')


def main():
    assert sha(P / 'precheck_manifest.json') == 'e5a61613fa0af7d235c5bd8017e8a9b1171b3577d06509c31ccb75f74088892c'
    frozen = read('precheck_manifest.json')
    for name, digest in frozen['whitelist_sha256'].items():
        assert sha(P / name) == digest, ('frozen local', name)
    for name, digest in frozen['external_readonly_sha256'].items():
        assert sha(Path(name)) == digest, ('frozen external', name)
    assert sha(P / 'full_release.json') == '70ddaf10a385b49e9cb9c96db88ad02fd5303c2a8d041c77a158b54ffdd37491'
    release = read('full_release.json')
    assert release['status'] == 'MAIN_RELEASED_ONE_NATIVE_AFTER_B055'
    assert release['B055_review_accepted_by_main'] is True
    assert sha(Path(release['B055_manifest_path'])) == release['B055_manifest_sha256'] == 'c1171119bbefd07f9120633791d13d5c76dbc43a43034a5481462329f392e1ab'
    assert release['parameters'] == dict(beta='1/80', point_numerators=[50, 50, 50, 45, 27, 0], point_denominator=50, T='32/825', Arb_bits=512)
    assert release['max_native_attempts'] == 1 and release['native_timeout_seconds'] == 1800
    execution = read('full_execution.json')
    assert execution['attempt'] == 1 and execution['phase'] == 'EXITED' and execution['exit_code'] == 0
    assert execution['timeout_seconds'] == 1800 and execution['process_wall_seconds'] < 1800
    assert execution['main_release_sha256'] == sha(P / 'full_release.json')
    assert execution['precheck_manifest_sha256'] == sha(P / 'precheck_manifest.json')
    assert execution['source_sha256'] == sha(P / 'sparse4_probe.cpp')
    assert execution['input_sha256'] == sha(P / 'native_input.txt')
    assert execution['binary_sha256'] == sha(P / 'scratch/sparse4_probe') == frozen['binary_sha256']
    assert execution['stdout_sha256'] == sha(P / 'native_full.stdout.json')
    assert execution['stderr_sha256'] == sha(P / 'native_full.stderr.txt')
    assert (P / 'native_full.stderr.txt').read_bytes() == b''
    assert (P / 'full_token.txt').read_text() == 'C055_MAIN_RELEASED_FINAL_DEPENDENCIES_CHECKED\n'
    native = read('native_full.stdout.json')
    assert native['status'] == 'COMPLETE_ONE_JOINT_SCALAR_PROBE' and native['precision_bits'] == 512
    assert native['paired_image_cone_terms'] == native['G_DD5_evaluations'] == native['H_DD5_evaluations'] == 25399680
    assert native['G_logical_terms'] == native['H_logical_terms'] == 609592320
    assert native['total_DD5_evaluations'] == 50799360
    checkpoints = [json.loads(line) for line in (P / 'full_checkpoint.jsonl').read_text().splitlines()]
    assert len(checkpoints) == 61
    for i, row in enumerate(checkpoints[:-1], 1):
        assert row['phase'] == 'running' and row['complete'] is False and row['images_complete'] == 32 * i
        assert row['paired_image_cone_terms'] == row['G_DD5_evaluations'] == row['H_DD5_evaluations'] == 32 * i * 13229
        assert row['logical_terms_each_field'] == row['paired_image_cone_terms'] * 24
        assert row['total_DD5_evaluations'] == row['paired_image_cone_terms'] * 2
    last = checkpoints[-1]
    assert last['complete'] is True and last['phase'] == 'computed'
    assert last['paired_image_cone_terms'] == last['G_DD5_evaluations'] == last['H_DD5_evaluations'] == 25399680
    assert last['total_DD5_evaluations'] == 50799360
    assert last['seconds'] == native['kernel_seconds']
    G, H, K = [bounds(native[key]) for key in ['G', 'H_T', 'K_T']]
    assert K[0] <= G[0] - H[1] and K[1] >= G[1] - H[0]
    threshold = Q(32, 825) / 4096
    assert threshold == Q(1, 105600) and K[1] < threshold
    comparison = read('point_comparison.json')
    assert comparison['status'] == 'COMPLETE_ONE_FIXED_POINT_COMPARISON'
    assert comparison['classification'] == 'EXCLUDES_ONLY_THIS_FIXED_FRAME_TIME_NECESSARY_CONDITION'
    assert comparison['threshold_exact'] == str(threshold)
    assert comparison['zero_node_logical_terms_each_field'] == native['zero_node_paired_terms'] * 24 == 6238488
    assert comparison['confluent_logical_terms_each_field'] == native['confluent_paired_terms'] * 24 == 2039688
    for key in ['global_mT_certified', 'global_g_or_I_certified', 'new_Cn_bound']:
        assert comparison[key] is False
    decimal = {key: dict(lower=decimal_bound(a, ROUND_FLOOR), upper=decimal_bound(b, ROUND_CEILING))
               for key, (a, b) in dict(G=G, H_T=H, K_T=K, threshold_minus_K_T=(threshold - K[1], threshold - K[0])).items()}
    return dict(status='PASS_SAVED_FULL_RESULT_NO_NATIVE_REPLAY', frozen_local_files_checked=len(frozen['whitelist_sha256']),
                frozen_external_files_checked=len(frozen['external_readonly_sha256']), precheck_manifest_sha256=sha(P / 'precheck_manifest.json'),
                native_pid=execution['native_pid'], native_attempt=1, native_exit_code=0, process_wall_seconds=execution['process_wall_seconds'],
                checkpoint_lines=len(checkpoints), completed_images=1920, paired_terms=25399680, logical_terms_each_field=609592320,
                total_DD5_field_evaluations=50799360, interval_K_encloses_subtraction=True, threshold_exact=str(threshold),
                K_upper_strictly_below_threshold=True, classification=comparison['classification'],
                outward_decimal_enclosures_50_places=decimal, new_Cn_bound=False, global_mT_certified=False, global_g_or_I_certified=False)


if __name__ == '__main__':
    print(json.dumps(main(), indent=2))
