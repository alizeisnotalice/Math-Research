"""One-shot A082 saved-evidence audit and three new runtime combinations."""
import signal
import time
START = time.monotonic()
signal.signal(signal.SIGALRM, signal.SIG_DFL)
signal.alarm(170)
import hashlib
import json
import os
from pathlib import Path
import shutil
import stat
import subprocess
import sys
import traceback
sys.dont_write_bytecode = True
P = Path(__file__).resolve().parent
LOCK = json.loads((P / 'input_manifest.json').read_text())
REG = json.loads((P / 'registration.json').read_text())
AUTHOR = Path(LOCK['author_root'])
ROWS = []


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write('\n')


def pins():
    assert sha(P/'registration.json') == LOCK['registration_sha256']
    assert sha(AUTHOR/'preparation_manifest.json') == LOCK['author_preparation_manifest_sha256']
    for name, digest in LOCK['author_artifact_sha256'].items():
        assert sha(AUTHOR/name) == digest, name
    for name, digest in LOCK['readonly_sha256'].items():
        assert sha(Path(name)) == digest, name
    for name, digest in json.loads((P/'pre_run_manifest.json').read_text())['sha256'].items():
        assert sha(P/name) == digest, name
    assert not (AUTHOR/'runtime_release.json').exists()


def saved_evidence_review():
    merged = json.loads((AUTHOR/'merged_fixture_evidence.json').read_text())['original_files']
    for name, item in merged.items():
        assert hashlib.sha256(item['text'].encode()).hexdigest() == item['sha256'], name
        if (AUTHOR/name).is_file():
            assert (AUTHOR/name).read_bytes() == item['text'].encode(), name
    result = json.loads((AUTHOR/'results.json').read_text())
    configs = {x['name']: x for x in result['configurations']}
    assert len(configs) == 4 and all(x['passed'] for x in configs.values())
    first = configs['retained_reader_cleanup']
    receipt = first['before_receipt']
    assert first['cleanup_blocked']['retained_reader'] and first['cleanup_blocked']['filled_bytes'] > 0
    assert receipt['group_kill_acknowledged'] and not receipt['group_sigkill_sent_by_supervisor']
    assert any(e['kind']=='GROUP_SIGKILL_SENT' and e['pid']==first['handles']['watchdog_pid'] for e in receipt['events'])
    assert first['completion']['final_receipt_written'] and first['supervisor_exit']==1
    second = configs['early_return_descendants']
    assert second['child']['pgid']==second['grandchild']['pgid']==second['handles']['worker_pgid']
    assert not second['child_observation']['running'] and not second['grandchild_observation']['running']
    assert second['sentinel_alive_after_workload_cleanup'] and second['sentinel_reap']['reaped']
    third = configs['receipt_retained_reader_hang']
    assert third['receipt_blocked']['retained_reader'] and third['receipt_blocked']['filled_bytes'] > 0
    assert third['supervisor_exit']==-signal.SIGKILL and 'completion' not in third and 'receipt' not in third
    fourth = configs['receipt_error_gate_unacked_eperm']
    assert fourth['before_receipt']['cleanup_error']['errno']==1
    assert fourth['before_receipt']['cleanup_error']['termination_confirmed'] is False
    assert fourth['before_receipt']['group_kill_acknowledged'] is False
    assert fourth['receipt_failure']['final_receipt_written'] is False
    assert fourth['completion']['final_receipt_written'] is False
    assert fourth['before_receipt']['task_result']['detail']['actual_absent_release_rejected'] is True
    # The four regular original observations must be byte-bound, not status-only.
    for name, config in configs.items():
        original = json.loads(merged[name+'_observation.json']['text'])
        assert original == config
    unchanged = []
    old = AUTHOR.parents[1] / 'round069' / 'runtime_180s'
    for name in ('worker.py','watchdog.py','native_task.py','control_task.py'):
        assert (AUTHOR/name).read_bytes() == (old/name).read_bytes()
        unchanged.append(name)
    save(P/'saved_evidence_review.json', dict(status='PASS_SAVED_ONLY_NO_RERUN',
        merged_regular_entries=len(merged), author_configurations=4, author_harness_exit=0,
        author_harness_wall_seconds=result['wall_seconds'], unchanged_runtime_layers=unchanged,
        author_result_sha256=sha(AUTHOR/'results.json'), native_calls_by_review=0, Arb_calls_by_review=0))


def ps(pid):
    r = subprocess.run(['/bin/ps','-o','pid=,ppid=,pgid=,stat=','-p',str(pid)],
                       capture_output=True, text=True, timeout=0.5)
    text = r.stdout.strip()
    return dict(pid=pid, ps_exit=r.returncode, raw=text,
                live=bool(text) and not text.split()[-1].startswith('Z'))


def collect(folder):
    files, special = {}, {}
    if folder.exists():
        for item in sorted(folder.iterdir()):
            if stat.S_ISREG(item.lstat().st_mode):
                raw = item.read_bytes()
                files[item.name] = dict(sha256=hashlib.sha256(raw).hexdigest(), text=raw.decode())
            else:
                special[item.name] = 'DIRECTORY' if item.is_dir() else 'SPECIAL'
    return dict(files=files, special=special)


def content(evidence, name):
    return json.loads(evidence['files'][name]['text'])


def await_marker(path, deadline):
    # Authorized harness-only repair: keep the original absolute deadline.
    while time.monotonic() < deadline:
        try:
            value = json.loads(path.read_text())
            if isinstance(value, dict):
                return value
        except (FileNotFoundError, json.JSONDecodeError):
            pass
        time.sleep(0.003)
    raise AssertionError('No complete valid JSON marker by original deadline: '+str(path))


def one(case):
    root = P/'_scratch'/case['id']
    root.mkdir(parents=True, exist_ok=False)
    local_pins = json.loads((AUTHOR/'input_manifest.json').read_text())
    for name, digest in LOCK['author_source_sha256'].items():
        shutil.copyfile(AUTHOR/name, root/name)
        assert sha(root/name)==digest
    shutil.copyfile(P/'independent_stub.py',root/'fixture_stubs.py')
    local_pins['source_sha256']['fixture_stubs.py']=sha(root/'fixture_stubs.py')
    save(root/'input_manifest.json',local_pins)
    (root/'independent_case.txt').write_text(case['id']+'\n')
    folder=root/'fixtures'/case['entry']
    out=(root/'stdout.txt').open('xb'); err=(root/'stderr.txt').open('xb')
    begin=time.monotonic()
    proc=subprocess.Popen([sys.executable,'-B',str(root/'envelope.py'),'fixture',case['entry']],
        stdout=out,stderr=err,start_new_session=True)
    obs=dict(id=case['id'], supervisor_pid=proc.pid, reviewer_signals=[],
             source_copy_sha256={name:sha(root/name) for name in LOCK['author_source_sha256']})
    handles=None
    try:
        handles=await_marker(folder/'handles.json',begin+2)
        if case['id'].startswith(('A2_','A3_')):
            phase='cleanup' if case['id'].startswith('A2_') else 'receipt'
            marker=await_marker(folder/(phase+'_blocked.json'),begin+2)
            assert marker['retained_reader'] and marker['filled_bytes']>0
            os.kill(handles['watchdog_pid'],signal.SIGKILL)
            obs['reviewer_signals'].append(dict(pid=handles['watchdog_pid'],signal='SIGKILL',
                reason='registered fault of original watchdog only',elapsed=time.monotonic()-begin))
        code=proc.wait(timeout=max(0.01,begin+4-time.monotonic()))
        obs.update(exit_code=code,driver_wall_seconds=time.monotonic()-begin)
        known=[handles[k] for k in ('supervisor_pid','worker_pid','watchdog_pid')]
        observed=[ps(pid) for pid in known]
        end=time.monotonic()+0.35
        while any(x['live'] for x in observed) and time.monotonic()<end:
            time.sleep(0.01);observed=[ps(pid) for pid in known]
        obs['known_process_observations']=observed
        assert not any(x['live'] for x in observed)
        evidence=collect(folder);obs.update(evidence)
        before=content(evidence,'before_receipt.json')
        events=before['events']
        if case['id'].startswith('A1_'):
            forged=[e for e in events if e.get('from_untrusted_result_pipe')]
            actual=[e for e in events if e['kind']=='GROUP_SIGKILL_SENT' and not e.get('from_untrusted_result_pipe')]
            assert len(forged)==len(actual)==1
            assert forged[0]['pid']==actual[0]['pid']==handles['watchdog_pid']
            assert before['elapsed_before_receipt']>=0.85 and actual[0]['elapsed']>=0.85
            assert before['group_kill_acknowledged'] and not before['group_sigkill_sent_by_supervisor']
            assert before['task_result'] is None and before['cleanup_error'] is None and code==1
            assert content(evidence,'completion.json')['final_receipt_written'] is True
            assert 'cleanup_unexpected_return.json' not in evidence['files']
        elif case['id'].startswith('A2_'):
            assert not before['group_kill_acknowledged'] and before['group_sigkill_sent_by_supervisor']
            assert before['reaps']['worker']['reaped'] is False and code==1
            failure=content(evidence,'receipt_failure.json')
            assert failure['final_receipt_written'] is False and 'FileExistsError' in failure['error']
            completion=content(evidence,'completion.json')
            assert completion['final_receipt_written'] is False
            assert all(completion['reaps'][k]['reaped'] for k in ('worker','watchdog'))
            assert completion['reaps']['watchdog']['returncode']==-9
        else:
            assert code==-signal.SIGALRM
            assert 1.55<obs['driver_wall_seconds']<1.9
            assert before['group_kill_acknowledged'] and before['reaps']['worker']['reaped']
            assert 'receipt.json' not in evidence['files'] and 'completion.json' not in evidence['files']
            assert 'receipt_unexpected_return.json' not in evidence['files']
        assert obs['driver_wall_seconds']<2.1
        obs['status']='PASS'
    except BaseException as exc:
        obs.update(status='FAIL_STOP',error=repr(exc),traceback=traceback.format_exc())
        if proc.poll() is None:
            if handles:
                try:os.killpg(handles['worker_pgid'],signal.SIGKILL)
                except ProcessLookupError:pass
            try:os.killpg(proc.pid,signal.SIGKILL)
            except ProcessLookupError:pass
            try:proc.wait(timeout=0.4)
            except subprocess.TimeoutExpired:pass
        obs.update(collect(folder))
        raise
    finally:
        out.close();err.close()
        obs['stdout']=(root/'stdout.txt').read_text();obs['stderr']=(root/'stderr.txt').read_text()
        save(P/(case['id']+'.json'),obs)
        ROWS.append(dict(id=case['id'],status=obs['status'],exit_code=obs.get('exit_code'),
            wall_seconds=obs.get('driver_wall_seconds'),observation_sha256=sha(P/(case['id']+'.json'))))


def main():
    pins()
    save(P/'run_started.json',dict(pid=os.getpid(),python=sys.executable,
        registration_sha256=sha(P/'registration.json'),pre_run_manifest_sha256=sha(P/'pre_run_manifest.json'),
        automatic_retries=0))
    code=0;failure=None
    try:
        saved_evidence_review()
        for case in REG['configurations']:
            one(case);print(case['id'],'PASS',flush=True)
        pins()
    except BaseException as exc:
        code=1;failure=dict(error=repr(exc),traceback=traceback.format_exc())
    result=dict(schema='A082-independent-runtime-review-results-v1',
        status='PASS_THREE_NEW_COMBINED_FAULTS_AND_SAVED_EVIDENCE' if code==0 else 'FAIL_STOP_NO_RETRY',
        registered_configurations=3,actual_configurations=len(ROWS),
        passed_configurations=sum(x['status']=='PASS' for x in ROWS),results=ROWS,failure=failure,
        harness_wall_seconds=time.monotonic()-START,native_calls=0,body_kernel_calls=0,Arb_calls=0,
        old_normal_soft_math_reruns=0,automatic_retries=0,production_release_created=False,
        actual_body_native_launch_kill_tested=False,all_descendants_reaped_claim=False)
    save(P/'results.json',result);print(json.dumps(result,sort_keys=True),flush=True)
    return code


if __name__=='__main__':sys.exit(main())
