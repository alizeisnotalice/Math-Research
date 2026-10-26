"""One registered run, six configurations, no mathematical fixture or native call."""
import signal
import time
START = time.monotonic()
signal.signal(signal.SIGALRM, signal.SIG_DFL)
signal.setitimer(signal.ITIMER_REAL, 179.8)
import json
import os
from pathlib import Path
import subprocess
import sys
sys.dont_write_bytecode=True
from runtime_support import P, NAMES, FIXTURE, PRODUCTION, require, read, write_new, sha, verify_pins


def state_of(pid):
    # Only the known PID from our own Popen receipt, never a name/process scan.
    result = subprocess.run(['/bin/ps', '-p', str(pid), '-o', 'pid=,ppid=,pgid=,stat='],
                            capture_output=True, text=True, timeout=0.3)
    text = result.stdout.strip()
    return dict(pid=pid, ps_exit_code=result.returncode, observed=text,
                running=bool(text) and 'Z' not in text.split()[-1])


def main():
    verify_pins()
    registration=read(P/'registration.json')
    require(registration['fixtures']==list(NAMES) and registration['max_configurations']==6 and
            registration['native_calls']==0, 'Registration mismatch')
    require(not (P/'fixtures').exists() and not (P/'fixture_run_started.json').exists(),
            'Original fixture run already started; no rerun or overwrite')
    write_new(P/'fixture_run_started.json', dict(pid=os.getpid(), started_monotonic=START,
              registration_sha256=sha(P/'registration.json'), input_manifest_sha256=sha(P/'input_manifest.json'),
              max_configurations=6, native_calls=0, automatic_retries=0))
    results=[]
    for name in NAMES:
        tic=time.monotonic()
        child=subprocess.Popen([sys.executable, '-B', str(P/'envelope.py'), 'fixture', name],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        try:
            stdout,stderr=child.communicate(timeout=1.95)
        except subprocess.TimeoutExpired:
            child.kill()
            try:
                stdout,stderr=child.communicate(timeout=0.2)
            except subprocess.TimeoutExpired:
                stdout,stderr=b'',b'HARNESS_REAP_UNRESOLVED'
            write_new(P/'fixture_failure.json', dict(name=name, supervisor_pid=child.pid,
                      reason='Supervisor exceeded registered short envelope; stop without rerun',
                      returncode=child.returncode, wall_seconds=time.monotonic()-tic))
            raise
        folder=P/'fixtures'/name
        observed=dict(name=name, supervisor_pid=child.pid, supervisor_returncode=child.returncode,
                      wall_seconds=time.monotonic()-tic, stdout=stdout.decode(), stderr=stderr.decode(),
                      native_calls=0, automatic_retries=0)
        handles=read(folder/'handles.json')
        observed['handles']=handles
        observed['pid_observations']=[]
        if name=='receipt_write_hang':
            require(child.returncode == -signal.SIGKILL and (folder/'before_receipt.json').exists() and
                    not (folder/'receipt.json').exists(), 'Receipt-block watchdog kill/partial record missing')
            receipt=read(folder/'before_receipt.json')
            observed['receipt_write_blocked_and_supervisor_killed']=True
            observed['watchdog_reaped_by_supervisor']=False
        else:
            receipt=read(folder/'receipt.json')
            completion=read(folder/'completion.json')
            observed['completion']=completion
            require(completion['reaps']['worker']['reaped'] and completion['reaps']['watchdog']['reaped'],
                    'Direct owned children not reaped in ordinary completion')
            expected=0 if name in ('normal_child','gate_and_deadline_isolation') else 1
            require(child.returncode==expected, 'Unexpected supervisor return code')
        observed['receipt']=receipt
        require(observed['wall_seconds']<1.95 and receipt['deadline_policy']==FIXTURE, 'Short deadline violated')
        require(handles['worker_pid']==handles['worker_pgid'] and
                handles['worker_pid']!=read(folder/'start.json')['supervisor_pgid'], 'Worker not isolated')
        if name=='soft_timeout':
            require(read(folder/'stub_cleanup.json')['reaped'], 'Soft timeout original child not reaped')
        if name in ('cleanup_hang','descendant_group'):
            require(any(e['kind']=='GROUP_SIGKILL_SENT' for e in receipt['events']), 'Independent group fuse not observed')
        for filename in ('child.json','tree_child.json','grandchild.json'):
            path=folder/filename
            if path.exists():
                item=read(path)
                require(item['pgid']==handles['worker_pgid'], 'Stub descendant escaped own group')
                observation=state_of(item['pid'])
                require(not observation['running'], 'Known stub still running after group cleanup')
                observed['pid_observations'].append(observation)
        if name=='descendant_group':
            require((folder/'grandchild.json').exists(), 'Actual grandchild fixture not reached')
        write_new(folder/'harness_observation.json', observed)
        results.append(observed)
        print(json.dumps(dict(name=name,status='PASS',wall_seconds=observed['wall_seconds'],
                              supervisor_pid=child.pid,supervisor_returncode=child.returncode)),flush=True)
    write_new(P/'results.json', dict(schema='C069-six-runtime-fixtures-v1', status='PASS_SIX_RUNTIME_CONFIGURATIONS',
              harness_pid=os.getpid(), wall_seconds=time.monotonic()-START, configurations=results,
              native_calls=0, body_kernel_calls=0, mathematical_fixtures_repeated=0, automatic_retries=0,
              production_deadline_policy=PRODUCTION, production_gate='CLOSED_PENDING_REVIEW_AND_MAIN_RELEASE'))


if __name__=='__main__':
    main()
