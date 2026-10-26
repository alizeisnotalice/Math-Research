"""One bounded main-scheduler operation. No native launch and no batch loop.

Commands: status | authorize | reserve | attach TOKEN TOOL_SESSION_ID |
          verify TOKEN | import TOKEN | halt TOKEN
Each program has a hard 180-second limit. A failed operation is not retried.
"""
import sys
sys.dont_write_bytecode = True
import json
import signal
import time
import os
from datetime import datetime, timezone
from pathlib import Path
from runtime_support import C068, require_worker, PRODUCTION, bounded_reap
sys.path.insert(0, str(C068))
import core_v2 as V
from batch_context import Context, P, sha, read, write_new, require
import math_contract as M


def main(guard):
    require_worker(guard)
    tic = time.monotonic()
    # The enclosing worker/watchdog own both deadlines; never cancel them here.
    receipt_dir = P / 'operation_receipts'
    receipt_dir.mkdir(exist_ok=True)
    receipt_path = receipt_dir / f'{time.time_ns()}-{os.getpid()}.json'
    receipt = dict(schema='C068-bounded-control-program-v1', pid=os.getpid(), command=sys.argv,
                   status='RUNNING', start_utc=datetime.now(timezone.utc).isoformat(), timeout_seconds=180,
                   native_calls=0, automatic_retries=0)
    write_new(receipt_path, receipt)
    def save_receipt():
        receipt['wall_seconds'] = time.monotonic()-tic
        temporary = receipt_path.with_suffix('.pending')
        temporary.write_text(json.dumps(receipt, indent=2)+'\n')
        temporary.replace(receipt_path)
    try:
        require(len(sys.argv) >= 2, 'One explicit scheduler operation required')
        context = Context()
        latest = V.recover(P / 'state', context)
        state, command = latest['state'], sys.argv[1]
        receipt['before_snapshot_sha256'] = latest['sha256']
        event = None
        if command == 'status':
            require(len(sys.argv) == 2, 'status takes no arguments')
        elif command == 'authorize':
            require(len(sys.argv) == 2, 'authorize uses the unique main_release.json')
            event = dict(kind='AUTHORIZE_ONCE', release_path=str(P/'main_release.json'), release_sha256=sha(P/'main_release.json'))
        elif command == 'reserve':
            require(len(sys.argv) == 2, 'reserve takes no arguments and never launches')
            budget = state['batch']
            if budget['native_reservations'] >= 8 or (state['queue'] and state['queue'][0] not in budget['started_nodes'] and len(budget['started_nodes']) >= 8):
                receipt['bounded_stop'] = 'BATCH_BUDGET_EXHAUSTED_QUEUE_RETAINED'
            elif not state['queue'] and state['inflight'] is None:
                receipt['bounded_stop'] = 'NO_QUEUED_WORK_MAIN_MUST_AUDIT_FRONTIER'
            else:
                event = dict(kind='RESERVE_NEXT_ONCE')
        elif command in ('attach', 'verify', 'import', 'halt'):
            require(len(sys.argv) == (4 if command == 'attach' else 3), 'Wrong bounded-operation arguments')
            token = sys.argv[2]
            require(state['inflight'] is not None and state['inflight']['token'] == token, 'No matching original inflight token')
            node = state['nodes'][state['inflight']['node_id']]
            attempt = node['attempts'][-1]
            folder = context.attempt_directory(token)
            if command == 'attach':
                process = read(folder/'original_process.json')
                intent = read(folder/'intent.json')
                handle = dict(tool_session_id=int(sys.argv[3]), native_pid=process['owned_native_pid'], wrapper_pid=intent['owned_wrapper_pid'])
                data = dict(token=token, input_sha256=attempt['input_sha256'], handle=handle,
                            intent_sha256=sha(folder/'intent.json'), observed_original_process_record=process,
                            original_process_sha256=sha(folder/'original_process.json'), synthetic_fixture=False)
                path = folder/'original_handle.json'
                if path.exists():
                    require(read(path) == data, 'Existing original handle differs; never replace it')
                else:
                    write_new(path, data)
                event = dict(kind='ATTACH_ORIGINAL_HANDLE_ONCE', token=token, handle_path=str(path), handle_sha256=sha(path))
            elif command == 'verify':
                require(attempt['execution_handle'] is not None, 'Record original tool/PID handle before verification')
                rebuilt = M.reconstruct(context, state, token, folder)
                path = folder/'verification.json'
                if path.exists():
                    require(read(path) == rebuilt, 'Original verification differs from reconstructed native/Arb result')
                else:
                    write_new(path, rebuilt)
                print(json.dumps(dict(token=token, verification_path=str(path), sha256=sha(path),
                                      center=rebuilt['center_classification'], native_calls_this_program=0)))
            elif command == 'import':
                path = folder/'verification.json'
                event = dict(kind='IMPORT_ORIGINAL_RESULT_ONCE', token=token, verification_path=str(path), verification_sha256=sha(path))
            else:
                path = folder/'failure.json'
                event = dict(kind='HALT_ORIGINAL_ATTEMPT_NO_RETRY', token=token, failure_path=str(path), failure_sha256=sha(path))
        else:
            raise ValueError('Unregistered operation; no reset/retry/automatic-next-batch exists')
        if event is not None:
            V.append_event(P/'state', context, latest['sha256'], event)
            latest = V.recover(P/'state', context)
            state = latest['state']
        print(json.dumps(dict(operation=command, snapshot_sha256=latest['sha256'], sequence=latest['sequence'],
                              queue=state['queue'], inflight=state['inflight'], passed_leaves=state['passed_leaves'],
                              halted=state['halted'], reservations=state['batch']['native_reservations'],
                              started_nodes=state['batch']['started_nodes'], observed_native_calls=state['production_kernel_calls'],
                              native_calls_this_program=0, global_g_certificate=False,
                              program_wall_seconds=time.monotonic()-tic)))
        receipt.update(status='COMPLETED', exit_code=0, after_snapshot_sha256=latest['sha256'],
                       end_utc=datetime.now(timezone.utc).isoformat())
        save_receipt()
    except BaseException as error:
        receipt.update(status='STOP_NO_RETRY', exit_code=1, error=repr(error),
                       end_utc=datetime.now(timezone.utc).isoformat())
        save_receipt()
        raise

if __name__ == '__main__':
    raise SystemExit('Use the reviewed runtime envelope; standalone execution is forbidden')
