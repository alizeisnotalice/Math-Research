"""Exactly one author batch: four new scenarios, first failure stops, no body."""
import signal
import time
START=time.monotonic()
signal.signal(signal.SIGALRM,signal.SIG_DFL)
signal.setitimer(signal.ITIMER_REAL,179.8)
import json
import os
import subprocess
import sys
sys.dont_write_bytecode=True
from runtime_support import P,NAMES,FIXTURE,require,read,write_new,sha,verify_pins,bounded_reap


def observe(pid):
    proc=subprocess.run(['/bin/ps','-p',str(pid),'-o','pid=,ppid=,pgid=,stat='],
                        capture_output=True,text=True,timeout=.2)
    raw=proc.stdout.strip()
    return dict(pid=pid,ps_exit_code=proc.returncode,stdout=raw,
                running=bool(raw) and 'Z' not in raw.split()[-1])


def main():
    verify_pins()
    require(not (P/'run_started.json').exists() and not (P/'fixtures').exists(),'No repeat or continuation of original run')
    write_new(P/'run_started.json',dict(pid=os.getpid(),start_monotonic=START,registration_sha256=sha(P/'registration.json'),
              input_manifest_sha256=sha(P/'input_manifest.json'),native_calls=0,automatic_retries=0))
    results=[]
    for name in NAMES:
        sentinel=None
        supervisor=None
        tic=time.monotonic()
        record=dict(name=name,passed=False,body_kernel_calls=0,native_calls=0,automatic_retries=0)
        try:
            if name=='early_return_descendants':
                sentinel=subprocess.Popen([sys.executable,'-B',str(P/'fixture_stubs.py'),'sentinel',name],start_new_session=True)
                ready_end=time.monotonic()+.2
                while not (P/'sentinel_started.json').exists() and time.monotonic()<ready_end:
                    time.sleep(.002)
                require((P/'sentinel_started.json').exists(),'Unrelated sentinel readiness missing')
                record['sentinel']=read(P/'sentinel_started.json')
            supervisor=subprocess.Popen([sys.executable,'-B',str(P/'envelope.py'),'fixture',name],
                                         stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            record['supervisor_pid']=supervisor.pid
            stdout,stderr=supervisor.communicate(timeout=1.95)
            record.update(supervisor_exit=supervisor.returncode,stdout=stdout.decode(),stderr=stderr.decode(),
                          supervisor_wall_seconds=time.monotonic()-tic)
            folder=P/'fixtures'/name
            handles=read(folder/'handles.json')
            record['handles']=handles
            before=read(folder/'before_receipt.json')
            record['before_receipt']=before
            require(record['supervisor_wall_seconds']<1.95,'Short supervisor envelope exceeded')
            require(handles['worker_pid']==handles['worker_pgid'],'Own identity anchor not a group leader')
            if name=='receipt_retained_reader_hang':
                require(supervisor.returncode==-signal.SIGKILL,'Receipt blocked supervisor was not independently killed')
                require((folder/'receipt_full_pipe.json').exists() and not (folder/'receipt.json').exists(),
                        'True receipt blocking/checkpoint evidence missing')
                require(not (folder/'receipt_unexpected_return.json').exists(),'Receipt blocking write returned')
                record['receipt_blocked']=read(folder/'receipt_full_pipe.json')
                record['watchdog_reaped_by_supervisor_claim']=False
            else:
                record['completion']=read(folder/'completion.json')
                require(supervisor.returncode==(0 if name=='early_return_descendants' else 1),'Unexpected supervisor status')
                if name=='receipt_error_gate_unacked_eperm':
                    record['receipt_failure']=read(folder/'receipt_failure.json')
                    require(before['cleanup_error']['errno']==1 and not before['group_kill_acknowledged'] and
                            not before['group_sigkill_sent_by_supervisor'],'Unacknowledged EPERM was treated as success')
                    require(before['original_worker_observation']['identity_anchor_unreaped'],'Lost original identity observation')
                    require(record['receipt_failure']['status']=='STOP_RECEIPT_WRITE_ERROR_NO_RETRY' and
                            not record['completion']['final_receipt_written'],'Receipt error falsely marked completed')
                    require(before['task_result']['detail']['actual_absent_release_rejected'],'Missing release not dynamically rejected')
                else:
                    record['receipt']=read(folder/'receipt.json')
                    require(record['completion']['reaps']['worker']['reaped'] and
                            record['completion']['reaps']['watchdog']['reaped'],'Original direct processes not reaped')
            if name=='retained_reader_cleanup':
                record['cleanup_blocked']=read(folder/'cleanup_full_pipe.json')
                require(record['cleanup_blocked']['retained_reader'] and record['cleanup_blocked']['filled_bytes']>0,
                        'No retained-reader full-pipe evidence')
                require(not (folder/'cleanup_unexpected_return.json').exists(),'Cleanup full-pipe returned')
                require(any(e['kind']=='GROUP_SIGKILL_SENT' for e in before['events']) and
                        not before['group_sigkill_sent_by_supervisor'],'Redundant post-ACK group kill not eliminated')
            if name=='early_return_descendants':
                require(sentinel.poll() is None,'Unrelated sentinel was terminated by tested group cleanup')
                require(record['sentinel']['pgid']!=handles['worker_pgid'],'Sentinel is not an unrelated own PG')
                record['sentinel_alive_after_workload_cleanup']=True
                for role in ('child','grandchild'):
                    child=read(folder/f'{role}.json')
                    require(child['pgid']==handles['worker_pgid'],'Descendant escaped registered PG')
                    record[role]=child
                    record[role+'_observation']=observe(child['pid'])
                    require(not record[role+'_observation']['running'],'Known descendant still running')
            record['passed']=True
        except BaseException as error:
            record.update(error=repr(error),failed=True)
            if supervisor is not None and supervisor.poll() is None:
                supervisor.kill()  # this harness's own original child only
                record['emergency_supervisor_reap']=bounded_reap(supervisor,time.monotonic()+.2)
            raise
        finally:
            if sentinel is not None:
                record['sentinel_original_cleanup_sent']=sentinel.poll() is None
                if sentinel.poll() is None:
                    sentinel.terminate()
                record['sentinel_reap']=bounded_reap(sentinel,time.monotonic()+.2)
            record['configuration_wall_seconds']=time.monotonic()-tic
            write_new(P/(name+'_observation.json'),record)
        results.append(record)
        print(json.dumps(dict(name=name,status='PASS',supervisor_pid=supervisor.pid,
                              wall_seconds=record['configuration_wall_seconds'])),flush=True)
    write_new(P/'results.json',dict(status='PASS_FOUR_NEW_RUNTIME_CONFIGURATIONS',harness_pid=os.getpid(),
              wall_seconds=time.monotonic()-START,configurations=results,native_calls=0,body_kernel_calls=0,
              old_math_fixture_reruns=0,old_normal_soft_reruns=0,automatic_retries=0,
              production_gate='CLOSED_PENDING_INDEPENDENT_REVIEW_AND_MAIN_RELEASE'))


if __name__=='__main__':
    main()
