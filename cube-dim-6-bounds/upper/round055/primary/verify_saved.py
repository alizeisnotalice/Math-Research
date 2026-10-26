"""Read-only dyadic probe/count verification. No private paths or kernels."""
import json
from fractions import Fraction as Q
from pathlib import Path


def bounds(row):
    factor = Q(2)**int(row['common_binary_exponent'])
    return int(row['lower_mantissa'])*factor, int(row['upper_mantissa'])*factor


def main():
    p = Path(__file__).resolve().parents[1]/'sparse4_first_scalar_probe'
    native = json.loads((p/'native_full.stdout.json').read_text())
    receipt = json.loads((p/'full_execution.json').read_text())
    checkpoints = [json.loads(line) for line in (p/'full_checkpoint.jsonl').read_text().splitlines()]
    assert native['status'] == 'COMPLETE_ONE_JOINT_SCALAR_PROBE'
    assert receipt['attempt'] == 1 and receipt['phase'] == 'EXITED' and receipt['exit_code'] == 0
    assert (p/'native_full.stderr.txt').read_bytes() == b''
    assert native['precision_bits'] == 512
    assert native['paired_image_cone_terms'] == 13229*1920 == 25399680
    assert native['G_logical_terms'] == native['H_logical_terms'] == 609592320
    assert native['G_DD5_evaluations'] == native['H_DD5_evaluations'] == 25399680
    assert native['total_DD5_evaluations'] == 50799360
    assert len(checkpoints) == 61 and checkpoints[-1]['complete']
    for i, row in enumerate(checkpoints[:-1], 1):
        assert row['images_complete'] == 32*i
        assert row['paired_image_cone_terms'] == 32*i*13229
    g, h, k = [bounds(native[name]) for name in ['G', 'H_T', 'K_T']]
    assert all(lower <= upper for lower, upper in [g, h, k])
    assert k[0] <= g[0]-h[1] <= g[1]-h[0] <= k[1]
    assert 0 < k[0] <= k[1] < Q(1, 105600)
    print('PASS saved complete probe/counts/interval K<T/4096; fixed frame/time only;0 kernel')


if __name__ == '__main__':
    main()
