import sys, json, hashlib, signal, time, os
from pathlib import Path
from datetime import datetime, timezone
sys.dont_write_bytecode=True
signal.signal(signal.SIGALRM, signal.SIG_DFL)
signal.setitimer(signal.ITIMER_REAL, 175)
start=time.monotonic()
P=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-6/cube-dim-6-bounds')
C=P/'upper/round068/cover_batch001'; R=P/'upper/round082/runtime_revision'; O=P/'upper/round090/cover_batch001_actual'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def write(p,v):
    with Path(p).open('x') as f: json.dump(v,f,indent=2,sort_keys=True);f.write('\n');f.flush();os.fsync(f.fileno())
def ref(p):return dict(path=str(p),sha256=sha(p))
def sealed(p,key):
    m=read(p)
    for name,d in m[key].items(): assert sha(p.parent/name)==d,(p,name)
    return ref(p)
assert (O/'preflight_environment_failure.json').is_file()
assert not (O/'registration.json').exists() and not (C/'main_release.json').exists() and not (R/'runtime_release.json').exists()
core=sealed(C/'preparation_manifest.json','artifact_sha256')
runtime=sealed(R/'preparation_manifest.json','artifact_sha256')
ac=P/'audit/round068/cover_revision_review'; ar=P/'audit/round082/runtime_revision_review'
sealed(ac/'manifest.json','whitelist_sha256');sealed(ar/'manifest.json','artifact_sha256')
sys.path.insert(0,str(C))
from batch_context import Context
import core_v2 as V
ctx=Context(); before=V.recover(C/'state',ctx); state=before['state']
assert state['environment']=='PREPARE_ONLY' and state['queue']==['r'] and state['inflight'] is None and state['batch']['native_reservations']==0
write(O/'registration.json',dict(schema='main090-finite-eight-node-pilot-v1',time=datetime.now(timezone.utc).isoformat(),owner='main',math_python=sys.executable,
    input_sha256=sha(C/'input_manifest.json'),preparation_manifest=core,runtime_preparation=runtime,starting_state=before,
    max_new_nodes=8,max_total_native_reservations=8,target_g='1/65000',T='1/2',cutoff_ladder=['1/32','1/64','1/128','1/250','1/1000'],
    total_program_wall_seconds=180,resource='one guarded native at a time; main serial; existing mathematical venv',
    stop_conditions=['8 nodes or 8 total reservations','native or verification failure','timeout','unknown orphan','unresolved cutoff ladder'],
    no_automatic_next_batch=True,old_seed_replays=0,source_pin_verification_only=True,global_g_certificate=False,new_C6_bound=False,
    current_lower='2073015624830/1000000000000',current_upper='64',environment_failure_preserved=True))
write(O/'core_acceptance.json',dict(status='PASS_C068_REVISION_ACCEPTED_BY_MAIN',owner='main',
    C068_preparation_manifest_sha256=core['sha256'],independent_review_artifacts=[ref(ac/'manifest.json'),ref(ac/'report.md'),ref(ac/'proof_review.md')],
    future_native_jets_accepted_under_engine_proof=True,F061_is_not_a_universal_new_point_crosscheck=True,finite_caps='8 nodes/8 total native, no renewal',
    global_g_certificate=False,new_C6_bound=False))
v=read(C/'main_release.template.json');v.update(production_authorized=True,template_only=False,
    input_manifest_sha256=sha(C/'input_manifest.json'),preparation_manifest_sha256=core['sha256'],
    future_native_jets_accepted_under_engine_proof=True,revision_acceptance=ref(O/'core_acceptance.json'))
write(C/'main_release.json',v)
write(O/'runtime_acceptance.json',dict(status='PASS_RUNTIME_180S_ACCEPTED_BY_MAIN',owner='main',
    runtime_preparation_manifest_sha256=runtime['sha256'],independent_review_artifacts=[ref(ar/'manifest.json'),ref(ar/'report.md')],
    prior071_rejected_preserved=True,actual_body_deadline_kill_not_pretested=True,total_wall_cap_including_cleanup_seconds=180,
    no_automatic_retries=True,global_g_certificate=False,new_C6_bound=False))
v=read(R/'runtime_release.template.json');v.update(production_authorized=True,
    C068_main_release_sha256=sha(C/'main_release.json'),preparation_manifest_sha256=runtime['sha256'],runtime_acceptance=ref(O/'runtime_acceptance.json'))
write(R/'runtime_release.json',v)
ctx.verified_release(C/'main_release.json',sha(C/'main_release.json'))
sys.path.insert(0,str(R));import runtime_support
assert runtime_support.production_gate(['control','status'])==sha(R/'runtime_release.json')
receipt=dict(status='PASS_PREFLIGHT_RELEASED_NOT_YET_NATIVE',pid=os.getpid(),python=sys.executable,time=datetime.now(timezone.utc).isoformat(),wall_seconds=time.monotonic()-start,
    core_release=ref(C/'main_release.json'),runtime_release=ref(R/'runtime_release.json'),initial_state=before,
    native_calls=0,production_transitions=0,exit_code=0)
write(O/'preflight_receipt.json',receipt)
print(json.dumps(receipt,sort_keys=True))
