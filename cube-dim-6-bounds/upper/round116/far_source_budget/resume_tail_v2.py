#!/usr/bin/env python3
"""Explicit output-only repair: keep cfg1 checkpoint, invoke cfg_tail ONCE."""
import sys
sys.set_int_max_str_digits(0)
import json
import os
import resource
import time
from pathlib import Path
import run_prepared as original

HERE = Path(__file__).resolve().parent


def main():
    start, cpu = time.monotonic(), time.process_time()
    release = json.loads((HERE/'parent_release_v2.json').read_text())
    reg = json.loads((HERE/'resume_registration_v2.json').read_text())
    original.req(release.get('authorized') is True and release.get('round') == '116v2',
                 'explicit output-only continuation release')
    original.req(release['script_sha256'] == original.sha(__file__), 'v2 entry pin')
    original.req(release['registration_sha256'] == original.sha(HERE/'resume_registration_v2.json'),
                 'v2 registration pin')
    for name, expected in reg['input_sha256'].items():
        original.req(original.sha(HERE/name) == expected, 'fixed input '+name)
    for row in json.loads((HERE/'source_manifest.json').read_text())['sources']:
        original.req(original.sha(row['path']) == row['sha256'], 'unchanged source '+row['role'])
    checkpoint = json.loads((HERE/'checkpoint.json').read_text())
    old_receipt = json.loads((HERE/'run_receipt.json').read_text())
    original.req(checkpoint['count'] == 1 and len(checkpoint['completed_configurations']) == 1,
                 'one complete saved cfg1 only')
    original.req(old_receipt['status'] == 'FAILED' and old_receipt['completed'] == 1,
                 'original failure preserved')
    original.req('4300 digits' in old_receipt['failure'], 'only documented output limitation')
    original.req(not (HERE/'results.json').exists(), 'no complete result overwritten')
    with (HERE/'started_v2.json').open('x') as f:
        json.dump({'pid': os.getpid(), 'utc': original.utc(),
                   'entry_sha256': original.sha(__file__),
                   'release_sha256': original.sha(HERE/'parent_release_v2.json'),
                   'cfg1_checkpoint_sha256': original.sha(HERE/'checkpoint.json'),
                   'cfg_window_calls': 0, 'maximum_cfg_tail_calls': 1}, f, indent=2)
        f.write('\n')
    status, failure, calls = 'FAILED', None, 0
    try:
        inp = json.loads((HERE/'upstream.json').read_text())
        pi = (original.endpoint(inp['pi_saved036'], 'lower'),
              original.endpoint(inp['pi_saved036'], 'upper'))
        calls += 1
        tail = original.cfg_tail(inp, pi)
        both = checkpoint['completed_configurations'] + [tail]
        status = 'COMPLETED_2_FIXED_CONFIGURATIONS_AFTER_OUTPUT_ONLY_V2'
        data = {'status': status, 'configurations': both,
                'cfg1_original_checkpoint_sha256': original.sha(HERE/'checkpoint.json'),
                'original_failure_receipt_sha256': original.sha(HERE/'run_receipt.json'),
                'v2_release_sha256': original.sha(HERE/'parent_release_v2.json'),
                'cfg_window_calls_v2': 0, 'cfg_tail_calls_v2': calls,
                'new_C6_bound': False, 'full_far_paid': False}
        original.atomic('checkpoint_v2.json', data)
        original.atomic('results.json', data)
        print(json.dumps({'status': status, 'pid': os.getpid(),
                          'cfg1_reused_not_rerun': True,
                          'window_I_lower_display': both[0]['I_envelope_lower_display'],
                          'tail_norm_display': tail['tail_norm_display']}), flush=True)
    except BaseException as e:
        failure = type(e).__name__ + ': ' + str(e)
        raise
    finally:
        ru = resource.getrusage(resource.RUSAGE_SELF)
        original.atomic('run_receipt_v2.json', {
            'status': status, 'failure': failure, 'pid': os.getpid(), 'utc': original.utc(),
            'cfg_window_calls': 0, 'cfg_tail_calls': calls,
            'CPU_seconds': time.process_time()-cpu,
            'wall_seconds_sample_before_final_receipt': time.monotonic()-start,
            'timing_scope': 'child start before release checks, sample before final receipt/exit; external observer separately covers complete child exit',
            'max_rss_bytes_macos': ru.ru_maxrss, 'children': 0,
            'runtime_limit': None, 'CPU_limit': None,
            'output_only_repair': 'sys.set_int_max_str_digits(0); mathematical cfg_tail unchanged',
            'first_run_reclassified_as_success': False,
        })


if __name__ == '__main__':
    main()
