#!/usr/bin/env python3
"""Second A090 saved-only group: original runtime receipts and exact comparison.

No production imports, oracle calls, native launch, signal, or source writes.
Only ps reads the explicitly recorded task PIDs. It never scans by name.
"""
import sys
sys.dont_write_bytecode = True
import signal
signal.signal(signal.SIGALRM, signal.SIG_DFL)
signal.setitimer(signal.ITIMER_REAL, 179.5)
import resource
resource.setrlimit(resource.RLIMIT_CPU,(20,20))
import os
import hashlib
import json
from pathlib import Path
import subprocess
import time
from fractions import Fraction as Q
from datetime import datetime, timezone

START=time.monotonic()
CPU=time.process_time()
OUT=Path(__file__).resolve().parent
ROOT=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds')
C=ROOT/'upper/round068/cover_batch001'
R=ROOT/'upper/round082/runtime_revision'
A=ROOT/'upper/round090/cover_batch001_actual'
HASHES={}
POLICY=dict(soft=170.0,group=178.0,parent=179.0,watchdog=179.5,supervisor=179.8,worker_fallback=178.0,reap=1.0)


def need(ok,text):
    if not ok:
        raise ValueError(text)


def sha(p):
    p=Path(p)
    digest=hashlib.sha256(p.read_bytes()).hexdigest()
    old=HASHES.setdefault(str(p),digest)
    need(old==digest,'Changing input '+str(p))
    return digest


def read(p):
    sha(p)
    return json.loads(Path(p).read_text())


def write_new(name,data):
    with (OUT/name).open('x') as f:
        json.dump(data,f,indent=2,sort_keys=True)
        f.write('\n')
        f.flush()
        os.fsync(f.fileno())


def main():
    core=read(OUT/'core_results.json')
    seal=read(core['final_seal_path'])
    need(sha(core['final_seal_path'])==core['final_seal_sha256'],'Same final scientific seal')
    observed=read(core['observed_original_sessions_path'])
    by_token={x['token']:x for x in core['attempts']}
    observed_by_token={x['token']:x for x in observed['attempts']}
    release_sha=sha(R/'runtime_release.json')
    operations=[]
    owned_pids=set()
    native_ops={}
    controls=[]
    bundle=read(A/'runtime_receipts.json')
    bundled_by_folder={x['original_directory']:x for x in bundle['operations']}
    need(len(bundled_by_folder)==42,'42 actual runtime operations in seal bundle')
    for folder in sorted((R/'operations').iterdir()):
        need(folder.is_dir(),'Unexpected non-directory in production operations')
        need(str(folder) in bundled_by_folder,'Actual operation absent in frozen bundle')
        bundled=bundled_by_folder[str(folder)]
        need(set(bundled['files'])=={'start.json','handles.json','before_receipt.json','receipt.json','completion.json'},
             'Five original runtime files in bundle')
        for name,item in bundled['files'].items():
            raw=item['original_utf8'].encode('utf-8')
            need(hashlib.sha256(raw).hexdigest()==item['sha256']==sha(folder/name) and
                 raw==(folder/name).read_bytes(),'Frozen runtime UTF8 exact bytes binding')
        start,handles,before,receipt,completion=(read(folder/name) for name in
            ('start.json','handles.json','before_receipt.json','receipt.json','completion.json'))
        need(start['mode']=='production' and start['runtime_release_sha256']==release_sha,
             'Operation release/mode')
        need(start['deadline_policy']==POLICY and receipt['deadline_policy']==POLICY,'Fixed deadline policy')
        need(start['args']==receipt['args'] and receipt['mode']=='production','Operation args/mode')
        need(before==receipt,'Before/final receipt consistency')
        need(receipt['status']=='TASK_RESULT_OBSERVED' and receipt['task_result']['logical_exit_code']==0,
             'Actual production operation error/timeout')
        need(receipt['task_result']['detail']['operation']==start['args'] and
             receipt['task_result']['detail']['runtime_release_sha256']==release_sha,'Task result binding')
        need(receipt['cleanup_error'] is None and completion['cleanup_error'] is None and
             completion['final_receipt_written'] is True and
             receipt['group_kill_acknowledged'] is True,'Cleanup/receipt complete')
        need(receipt['automatic_retries']==completion['automatic_retries']==start['automatic_retries']==0,
             'No automatic retry')
        need(receipt['all_descendants_reaped_claim'] is False and completion['all_descendants_reaped_claim'] is False,
             'Preserve reap scope')
        need(0<=receipt['elapsed_before_receipt']<=completion['elapsed']<=180,'Cleanup-inclusive elapsed<=180')
        need(completion['reaps']['worker']['reaped'] is True and completion['reaps']['watchdog']['reaped'] is True,
             'Original direct worker/watchdog reaped')
        need(completion['reaps']['worker']['returncode']==-9 and completion['reaps']['watchdog']['returncode']==0,
             'Normal workload-anchor and watchdog return codes')
        need(handles['supervisor_pid']==start['supervisor_pid']==receipt['supervisor_pid'] and
             handles['worker_pid']==handles['worker_pgid']==handles['inherited_native_pgid']==receipt['worker_pid'] and
             handles['watchdog_pid']==receipt['watchdog_pid'],'Runtime original handles')
        events=receipt['events']
        need(any(e['kind']=='WATCHDOG_ARMED' and e['pid']==handles['watchdog_pid'] and
                 e['worker_pgid']==handles['worker_pid'] for e in events),'Original watchdog armed')
        need(any(e['kind']=='WORKER_STARTED' and e['pid']==handles['worker_pid'] and
                 e['pgid']==handles['worker_pid'] for e in events),'Original worker started')
        need(not any(e['kind'] in ('SUPERVISOR_ERROR','SOFT_WORK_CUTOFF_SENT','GROUP_KILL_ERROR_UNRESOLVED') for e in events),
             'No timeout/cleanup error events')
        need(not (folder/'receipt_failure.json').exists(),'Receipt failure artifact')
        for name in ('worker.stderr','watchdog.stderr'):
            need((folder/name).read_bytes()==b'','Runtime stderr not empty')
        for path in folder.iterdir():
            need(path.is_file(),'Unexpected operation artifact type')
            sha(path)
        args=start['args']
        if args[0]=='native':
            token=args[1]
            need(token in by_token and token not in native_ops,'Unique original native operation')
            handle=by_token[token]['handle']
            need(handle['wrapper_pid']==handles['worker_pid'],'Original native wrapper operation binding')
            original=read(C/'attempts'/token/'original_process.json')
            need(original['owned_native_pid']==handle['native_pid'],'Native PID original process')
            owned_pids.add(handle['native_pid'])
            native_ops[token]=str(folder)
            observation=observed_by_token[token]
            if 'runtime_operation' in observation:
                need(Path(observation['runtime_operation'])==folder,'Tool observation original operation')
        else:
            need(args[0]=='control','Only original native/control operation')
            controls.append(args[1:])
        for key in ('supervisor_pid','worker_pid','watchdog_pid'):
            owned_pids.add(handles[key])
        operations.append(dict(folder=str(folder),args=args,handles=handles,
             status=receipt['status'],elapsed_including_cleanup=completion['elapsed'],
             final_receipt_written=completion['final_receipt_written'],reaps=completion['reaps'],
             all_descendants_reaped_claim=False))
    need(len(native_ops)==8 and set(native_ops)==set(by_token),'Exactly eight original native operation folders')
    need(len(operations)==42 and len(controls)==34,'42 operations:8 native/34 controls')
    # Main's final bounded reserve must observe exhausted budget and retain queue.
    control_receipts=[]
    for path in sorted((C/'operation_receipts').iterdir()):
        need(path.is_file() and path.suffix=='.json','No partial control receipts')
        data=read(path)
        need(data['status']=='COMPLETED' and data['exit_code']==0 and data['native_calls']==0 and
             data['automatic_retries']==0 and 0<=data['wall_seconds']<=180,'Control successful finite wall')
        owned_pids.add(data['pid'])
        control_receipts.append(data)
    budget=[x for x in control_receipts if x.get('bounded_stop')=='BATCH_BUDGET_EXHAUSTED_QUEUE_RETAINED']
    need(len(budget)==1 and budget[0]['before_snapshot_sha256']==budget[0]['after_snapshot_sha256']==core['snapshot_sha256'],
         'Real exhausted reserve retained immutable final queue')
    need(len(control_receipts)==len(controls),'Control receipts/actual operation cardinality')
    environment=read(A/'preflight_environment_failure.json')
    schema_failure=read(A/'preflight_schema_failure.json')
    need(environment['actual_native_calls']==0 and environment['actual_state_transitions']==0 and
         environment['error']=='ModuleNotFoundError: flint','First environment setup error preserved')
    need(schema_failure['native_calls']==0 and schema_failure['production_transitions']==0 and
         schema_failure['scientific_retry'] is False and 'artifact_sha256' in schema_failure['error'],
         'Second manifest-key setup error preserved separately')
    # No process discovery by name, kill or unbounded wait. Empty output means the
    # recorded PIDs are absent now, not that supervisor reaped all descendants.
    process=subprocess.run(['/bin/ps','-p',','.join(map(str,sorted(owned_pids))),
        '-o','pid=,ppid=,pgid=,stat=,lstart=,command='],capture_output=True,text=True,timeout=5)
    need(process.returncode in (0,1),'ps error')
    need(process.stdout.strip()=='','A recorded original task PID is still present; inspect without signaling')
    ledger=read(ROOT/'results/bounds_ledger.json')
    dual=read(ROOT/'coordination/phase_20261025/dual_breakthrough_checkpoint.json')
    frozen=read(ROOT/'coordination/phase_20261025/frozen_thresholds.json')
    current=ledger['lower']
    lo,hi=map(Q,current['exact_interval'])
    threshold=Q('2073015624830/1000000000000')
    need(current['id']=='R6_rectangular_exact_event_round075' and lo<=hi and lo>threshold and hi<Q(64),
         'Current complete075 interval and exact upper comparison')
    need(ledger['upper']['exact']=='64' and ledger['public_best_verified']['upper']['exact']=='64',
         'Current/public retained upper64')
    need(dual['dual_checkpoint'] is None and dual['post_dual_completed_rounds']==0 and
         dual['upper_since_freeze_certified'] is False,'Dual null/count0')
    lower_cert=ROOT/current['certificate']
    need(sha(lower_cert)==current['certificate_sha256'],'Current lower certificate identity only; no rerun')
    frozen_upper=Q(frozen['lower']['exact_interval'][1])
    need(lo>frozen_upper,'Complete current075 interval exceeds frozen exact old upper')
    comparison=dict(current_lower_id=current['id'],current_lower_relation=current['relation'],
        current_lower_complete_exact_interval=current['exact_interval'],
        current_lower_certificate_sha256=sha(lower_cert),current_lower_explicit_strict_threshold=str(threshold),
        current075_complete_lower_above_frozen_complete_upper=True,
        current075_complete_upper_below_64=True,retained_upper='64',
        public_upper='64',new090_C6_bound=False,new090_global_g_certificate=False,
        dual_checkpoint=None,post_dual_completed_rounds=0)
    for p,digest in list(HASHES.items()):
        need(hashlib.sha256(Path(p).read_bytes()).hexdigest()==digest,'Changed runtime evidence: '+p)
    result=dict(status='PASS_SAVED_RUNTIME_HANDLES_BUDGET_AND_COMPARISON',
        core_results_sha256=sha(OUT/'core_results.json'),operations=operations,
        operation_count=len(operations),native_operation_count=len(native_ops),control_operation_count=len(controls),
        control_receipt_count=len(control_receipts),budget_exhausted_receipt=budget[0],
        all_programs_cleanup_inclusive_within_180_seconds=True,
        maximum_operation_elapsed=max(x['elapsed_including_cleanup'] for x in operations),
        known_task_pid_observation=dict(pids=sorted(owned_pids),ps_exit_code=process.returncode,
             stdout=process.stdout,stderr=process.stderr,observation_utc=datetime.now(timezone.utc).isoformat(),
             assertion='Recorded task PIDs absent in this observation; no all-descendant reaping claim',signals_sent=0),
        comparison=comparison,source_hashes=HASHES,native_calls=0,Arb_calls=0,source_writes=0,
        all_descendants_reaped_claim=False)
    write_new('runtime_results.json',result)
    write_new('runtime_run_receipt.json',dict(status='PASS',pid=os.getpid(),exit_code=0,
        utc=datetime.now(timezone.utc).isoformat(),python=sys.executable,
        wall_seconds=time.monotonic()-START,cpu_seconds=time.process_time()-CPU,
        cpu_limit_seconds=20,wall_limit_seconds=180,cleanup_and_receipt_included=True,
        native_calls=0,Arb_calls=0,source_writes=0,automatic_retries=0,
        result_sha256=sha(OUT/'runtime_results.json')))
    print(json.dumps(dict(status=result['status'],operations=len(operations),native_operations=len(native_ops),
        max_elapsed=result['maximum_operation_elapsed'],known_pids_absent=len(owned_pids))))


if __name__=='__main__':
    try:
        main()
    except BaseException as error:
        write_new('runtime_failure_receipt.json',dict(status='STOP_NO_RETRY',pid=os.getpid(),error=repr(error),
             wall_seconds=time.monotonic()-START,cpu_seconds=time.process_time()-CPU,native_calls=0,
             Arb_calls=0,source_writes=0))
        raise
