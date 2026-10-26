#!/usr/bin/env python3
"""Authorized group2: finish source bindings; no recover or mathematical replay.

First-run exact arithmetic and event replay are evidenced only by the original
trace reaching line225 after lines191-222. This program does not rerun them.
"""
import sys
sys.dont_write_bytecode=True
import signal
signal.signal(signal.SIGALRM,signal.SIG_DFL)
signal.setitimer(signal.ITIMER_REAL,179.5)
import resource
resource.setrlimit(resource.RLIMIT_CPU,(20,20))
import hashlib
import json
import os
from pathlib import Path
import time
from datetime import datetime,timezone

START=time.monotonic()
CPU=time.process_time()
OUT=Path(__file__).resolve().parent
ROOT=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds')
C=ROOT/'upper/round068/cover_batch001'
R=ROOT/'upper/round082/runtime_revision'
A=ROOT/'upper/round090/cover_batch001_actual'
HASHES={}


def need(ok,text):
    if not ok:
        raise ValueError(text)


def sha(p):
    p=Path(p)
    value=hashlib.sha256(p.read_bytes()).hexdigest()
    old=HASHES.setdefault(str(p),value)
    need(old==value,'Changed input '+str(p))
    return value


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
    reg=read(OUT/'continuation_registration.json')
    for name,key in [('check_saved.py','first_script_sha256'),('core_failure_receipt.json','first_failure_receipt_sha256'),
                     ('first_run_tool_receipt.json','first_tool_observation_sha256')]:
        need(sha(OUT/name)==reg[key],'Frozen original reviewer artifact')
    failure=read(OUT/'core_failure_receipt.json')
    need(failure['native_calls']==failure['Arb_corner_calls']==failure['source_writes']==0,'Retained harness failure scope')
    seal_path=A/'science_manifest.json'
    need(sha(seal_path)=='5b9f004437d34fc0fd58be76ef762f1a7f975fa7cc4ce0e61dac2a89cc02f82d','Final scientific seal')
    seal=read(seal_path)
    for item in seal['files']:
        path=ROOT/item['path']
        need(path.is_relative_to(ROOT) and '..' not in Path(item['path']).parts,'Seal path scope')
        need(path.stat().st_size==item['bytes'] and sha(path)==item['sha256'],'Sealed original changed')
    # Source/release/reference checks that were pending when the first harness
    # JSON-parsed markdown. All prose files are authenticated only as bytes.
    for name in ('core_acceptance.json','runtime_acceptance.json'):
        acceptance=read(A/name)
        for item in acceptance['independent_review_artifacts']:
            need(sha(item['path'])==item['sha256'],'Independent review reference binding')
    for name,digest in read(R/'preparation_manifest.json')['artifact_sha256'].items():
        need(sha(R/name)==digest,'Runtime sealed source binding')
    pins=read(C/'input_manifest.json')
    for path,digest in pins['readonly_sha256'].items():
        need(sha(path)==digest,'Core readonly pin')
    for name,digest in pins['owned_source_sha256'].items():
        need(sha(C/name)==digest,'Core source pin')
    for name in ('input_manifest.json','preparation_manifest.json'):
        sha(C/name)
    runtime_pins=read(R/'input_manifest.json')
    for name,digest in runtime_pins['source_sha256'].items():
        need(sha(R/name)==digest,'Runtime source pin')
    for path,digest in runtime_pins['readonly_sha256'].items():
        need(sha(path)==digest,'Runtime readonly pin')
    sha(C/'state_initialized_once.json')
    # Read, rather than recover, already sealed final record and original
    # verification summaries. This does not add a second event replay/corner run.
    last_path=C/'state/snapshot_00000025.json'
    last=read(last_path)
    state=last['state']
    observed_path=A/'observed_original_sessions.json'
    observed=read(observed_path)
    summaries=[]
    for item in observed['attempts']:
        folder=C/'attempts'/item['token']
        saved=read(folder/'verification.json')
        expected=read(folder/'reserved_input.json')['input']
        receipt=read(folder/'native_receipt.json')
        native=read(folder/'native.json')
        box=saved['observation']['box_certificate']
        summaries.append(dict(token=item['token'],center=expected['center'],window=expected['window'],
            face_box=expected['face_box'],input_sha256=saved['input_sha256'],handle=read(folder/'original_handle.json')['handle'],
            center_classification=saved['center_classification'],F_interval=saved['observation']['jet']['F_interval'],
            whole_box_pass=box['status']=='PASS_CONDITIONAL_ON_MATCHING_CENTER_JET',
            whole_box_lower=box['tangent_bound_lower'],native_counts=saved['native_counts'],
            kernel_seconds=native['kernel_seconds'],native_wrapper_wall_seconds=receipt['wall_seconds'],
            verification_sha256=sha(folder/'verification.json'),receipt_sha256=sha(folder/'native_receipt.json')))
    for p,digest in list(HASHES.items()):
        need(hashlib.sha256(Path(p).read_bytes()).hexdigest()==digest,'Changed source binding '+p)
    result=dict(status='PASS_SOURCE_BINDING_CONTINUATION_AFTER_RETAINED_HARNESS_ERROR',
        final_seal_path=str(seal_path),final_seal_sha256=sha(seal_path),
        observed_original_sessions_path=str(observed_path),observed_original_sessions_sha256=sha(observed_path),
        sequence=last['sequence'],snapshot_sha256=sha(last_path),queue=state['queue'],
        inflight=state['inflight'],halted=state['halted'],passed_leaves=state['passed_leaves'],empty_leaves=state['empty_leaves'],
        reservations=state['batch']['native_reservations'],distinct_started_nodes=state['batch']['started_nodes'],
        original_native_calls_observed=state['production_kernel_calls'],attempts=summaries,
        summary_fields='Read from frozen originals; exact native/corner validation occurred before retained first-run exception',
        pure_core_recover_calls_this_program=0,pure_core_recover_calls_total=1,
        first_run_completed_scope=dict(script_sha256=reg['first_script_sha256'],tool_chunk='06194a',
            evidence='Observed traceback reached main line225; original lines191-222 returned and assertions passed',
            core_recover_returned=True,original_eight_saved_imports_exactly_checked=True,exact_corner_rows_checked=256,
            final_result_output_was_not_written=True,overall_first_harness_exit_code=1),
        native_calls=0,Arb_corner_calls=0,source_writes=0,source_hashes=HASHES,
        independent_exponential_recalculation=False,independent_full_gradient_recalculation=False,
        independent_zero_repeat_enumeration=False,global_g_certificate=False,new_C6_bound=False)
    write_new('core_results.json',result)
    write_new('continuation_run_receipt.json',dict(status='PASS',pid=os.getpid(),exit_code=0,
        utc=datetime.now(timezone.utc).isoformat(),python=sys.executable,
        wall_seconds=time.monotonic()-START,cpu_seconds=time.process_time()-CPU,
        cpu_limit_seconds=20,wall_limit_seconds=180,cleanup_and_receipt_included=True,
        pure_core_recover_calls_this_program=0,native_calls=0,Arb_corner_calls=0,
        source_writes=0,automatic_retries=0,result_sha256=sha(OUT/'core_results.json')))
    print(json.dumps(dict(status=result['status'],sequence=result['sequence'],queue=result['queue'],
        passed_leaves=result['passed_leaves'],recover_calls_this_program=0)))


if __name__=='__main__':
    try:
        main()
    except BaseException as error:
        write_new('continuation_failure_receipt.json',dict(status='STOP_NO_RETRY',pid=os.getpid(),
            error=repr(error),wall_seconds=time.monotonic()-START,cpu_seconds=time.process_time()-CPU,
            pure_core_recover_calls=0,native_calls=0,Arb_calls=0,source_writes=0))
        raise
