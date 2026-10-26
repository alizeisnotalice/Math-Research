"""Run the pre-existing 12-check preparation batch once, with a 180 s watchdog.

This wrapper cannot launch a native kernel or compile. It never overwrites the
registration, and refuses any previous receipt or preparation-result file.
"""
import sys
sys.dont_write_bytecode = True
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import time
from datetime import datetime, timezone

P = Path(__file__).resolve().parent
LIMIT = 180


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def utc():
    return datetime.now(timezone.utc).isoformat()


def main():
    receipt_path = P / 'preparation_execution_receipt.json'
    result_path = P / 'preparation_results.json'
    output_path = P / 'preparation_stdout.txt'
    error_path = P / 'preparation_stderr.txt'
    for path in (receipt_path, result_path, output_path, error_path):
        if path.exists():
            raise SystemExit(f'Refuse replay or overwrite: {path}')
    policy = json.loads((P / 'policy_amendment_180s.json').read_text())
    assert policy['per_program_wall_timeout_seconds'] == LIMIT
    assert policy['production_kernel_calls_authorized_now'] == 0
    assert sha(P / 'registration.json') == policy['original_registration_sha256']
    assert sha(Path(policy['A057_sealed_manifest_path'])) == policy['A057_sealed_manifest_sha256']
    command = [sys.executable, '-B', str(P / 'check_preparation.py')]
    tic = time.monotonic()
    receipt = dict(schema='C061-preparation-execution-v1', status='REGISTERED_NOT_YET_LAUNCHED',
                   started_utc=utc(), command=command, cwd=str(P), timeout_seconds=LIMIT,
                   policy_amendment_sha256=sha(P / 'policy_amendment_180s.json'),
                   registration_sha256=sha(P / 'registration.json'),
                   input_manifest_sha256=sha(P / 'input_manifest.json'),
                   code_sha256={name: sha(P / name) for name in
                                ('prepare_once.py', 'check_preparation.py', 'jet_adapter.py')},
                   new_native_invocations=0, new_production_kernel_calls=0, new_compiles=0,
                   old_seed_recalculations=0, C055_reruns=0, automatic_retries=0,
                   global_g_certificate=False, new_C6_bound=False)
    with receipt_path.open('x') as stream:
        json.dump(receipt, stream, indent=2)
        stream.write('\n')

    def save():
        receipt['wall_seconds'] = time.monotonic() - tic
        temp = receipt_path.with_suffix('.json.tmp')
        temp.write_text(json.dumps(receipt, indent=2) + '\n')
        temp.replace(receipt_path)

    child = None
    try:
        with output_path.open('xb') as out, error_path.open('xb') as err:
            child = subprocess.Popen(command, cwd=P, stdout=out, stderr=err,
                                     start_new_session=True)
            receipt.update(status='RUNNING_PREPARATION_ONLY', owned_child_pid=child.pid,
                           owned_process_group_id=child.pid, launched_utc=utc())
            save()
            try:
                returncode = child.wait(timeout=max(0.001, LIMIT - (time.monotonic() - tic)))
            except subprocess.TimeoutExpired:
                receipt.update(status='TIMEOUT_STOP_NO_RETRY', timeout_utc=utc())
                save()
                # start_new_session=True makes this group exclusive to this wrapper.
                os.killpg(child.pid, signal.SIGKILL)
                child.wait()
                receipt.update(native_exit_code=None, preparation_exit_code=child.returncode,
                               killed_only_owned_process_group=True)
                save()
                raise SystemExit('Preparation reached 180 s; owned group reaped; no replay')
        receipt.update(preparation_exit_code=returncode, ended_utc=utc(),
                       stdout_sha256=sha(output_path), stderr_sha256=sha(error_path),
                       stderr_empty=error_path.stat().st_size == 0)
        if result_path.exists():
            result = json.loads(result_path.read_text())
            receipt.update(result_sha256=sha(result_path), result_status=result['status'],
                           completed_configurations=len(result.get('configurations', [])))
        else:
            result = {}
        passed = (returncode == 0 and receipt['stderr_empty'] and
                  result.get('status') == 'PASS_TWELVE_PREPARATION_CONFIGURATIONS_ZERO_PRODUCTION_KERNELS' and
                  len(result.get('configurations', [])) == 12)
        receipt['status'] = 'PASS_PREPARATION_ONLY_NO_PRODUCTION_RELEASE' if passed else 'FAILED_STOP_NO_RETRY'
        save()
        if not passed:
            raise SystemExit('Preparation failed; receipt/checkpoint retained; no replay')
        print(receipt['status'])
    except BaseException as error:
        if child is not None and child.poll() is None:
            os.killpg(child.pid, signal.SIGKILL)
            child.wait()
            receipt['killed_only_owned_process_group'] = True
        receipt['ended_utc'] = utc()
        receipt['error'] = repr(error)
        if receipt['status'] not in ('TIMEOUT_STOP_NO_RETRY', 'FAILED_STOP_NO_RETRY'):
            receipt['status'] = 'FAILED_STOP_NO_RETRY'
        if result_path.exists():
            receipt['retained_checkpoint_sha256'] = sha(result_path)
        save()
        raise


if __name__ == '__main__':
    main()
