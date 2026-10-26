"""Eight predeclared C068 attack/repair fixtures, one run, zero native launches."""
import sys
sys.dont_write_bytecode = True
import copy
import hashlib
import json
import os
import signal
import time
from datetime import datetime, timezone
from fractions import Fraction as Q
from pathlib import Path
import core_v2 as V
import math_contract as M
from batch_context import Context, P, read, sha, write_new

ARB_CALLS = 0
CONTEXT_SERIAL = 0


def rejects(call, *args, **kwargs):
    try:
        call(*args, **kwargs)
    except (ValueError, KeyError, FileNotFoundError, FileExistsError):
        return True
    raise AssertionError('Expected rejection of registered attack')


def context(label):
    global CONTEXT_SERIAL
    CONTEXT_SERIAL += 1
    root = P/'fixtures'/f'{CONTEXT_SERIAL:02d}-{label}'
    root.mkdir(parents=True)
    c = Context(test_only=True, fixture_root=root)
    original = c.checker.certify
    def counted(*args, **kwargs):
        global ARB_CALLS
        ARB_CALLS += 1
        return original(*args, **kwargs)
    c.checker.certify = counted
    return c


def records(c):
    return [V.genesis_record(c)]


def event(history, c, e):
    previous = history[-1]
    state = V.apply_event(previous['state'], e, c)
    history.append(dict(schema=V.SCHEMA, sequence=len(history), previous_sha256=V.digest(previous), event=e, state=state))
    return state


def reserve(history, c):
    state = event(history, c, dict(kind='RESERVE_NEXT_ONCE'))
    node = state['nodes'][state['inflight']['node_id']]
    return state, node['attempts'][-1]


def raw(q):
    q = Q(q)
    require_power = q.denominator & (q.denominator-1) == 0
    assert require_power
    return dict(lower_mantissa=str(q.numerator), upper_mantissa=str(q.numerator),
                common_binary_exponent=str(-(q.denominator.bit_length()-1)))


def make_native_fixture(c, attempt, value=Q(1,32768)):
    token, expected = attempt['token'], attempt['input']
    folder = c.attempt_directory(token)
    folder.mkdir(parents=True)
    a, T, pt = Q(expected['window']['a']), Q(1,2), Q(1,2**20)
    pa = (5*value+T*pt)/a
    native = dict(status='all_real_nodes_supported', precision_bits=512, full_group_terms=967680,
                  computed_sign_pairs=483840, confluent_terms=0, zero_node_terms=0, unique_node_values=1,
                  kernel_seconds=0, band=raw(value), gradient=[raw(0) for _ in range(6)],
                  p_a=raw(pa), p_T=raw(pt), Euler_xgrad_plus5B_plusTpT_minusapa=raw(0))
    intent = dict(token=token, input=expected, input_sha256=attempt['input_sha256'],
                  release_sha256='0'*64, owned_wrapper_pid=900001, synthetic_fixture=True)
    write_new(folder/'intent.json', intent)
    write_new(folder/'native.json', native)
    (folder/'stderr.txt').write_bytes(b'')
    receipt = dict(schema='C068-native-receipt-v1', token=token, input_sha256=attempt['input_sha256'],
                   command=M.native_command(c, expected), **{k:expected[k] for k in ('model','model_definition','window','center')},
                   source_sha256=c.external_sha(c.paths['source']), binary_sha256=c.external_sha(c.paths['binary']),
                   fan_sha256=c.external_sha(c.paths['fan']), runtime_sha256=c.runtime_pins,
                   intent_sha256=sha(folder/'intent.json'), release_sha256='0'*64, native_exit_code=0,
                   stderr_empty=True, native_sha256=sha(folder/'native.json'), stderr_sha256=sha(folder/'stderr.txt'),
                   owned_wrapper_pid=900001, owned_native_pid=900002, timeout_seconds=180, synthetic_fixture=True)
    write_new(folder/'native_receipt.json', receipt)
    handle = dict(token=token, input_sha256=attempt['input_sha256'], synthetic_fixture=True,
                  handle={'fixture_original_handle':'NEVER_LAUNCHED'}, intent_sha256=sha(folder/'intent.json'))
    write_new(folder/'original_handle.json', handle)
    return folder, native, receipt


def attach(history, c, attempt):
    path = c.attempt_directory(attempt['token'])/'original_handle.json'
    return event(history, c, dict(kind='ATTACH_ORIGINAL_HANDLE_ONCE', token=attempt['token'], handle_path=str(path), handle_sha256=sha(path)))


def complete_fixture(history, c, value=Q(1,32768)):
    state, attempt = reserve(history,c)
    folder, _, _ = make_native_fixture(c,attempt,value)
    state = attach(history,c,attempt)
    rebuilt = M.reconstruct(c,state,attempt['token'],folder)
    path = folder/'verification.json'
    write_new(path,rebuilt)
    state = event(history,c,dict(kind='IMPORT_ORIGINAL_RESULT_ONCE',token=attempt['token'],verification_path=str(path),verification_sha256=sha(path)))
    return state, attempt, rebuilt


def forged_successor(history, e, state):
    return history+[dict(schema=V.SCHEMA,sequence=len(history),previous_sha256=V.digest(history[-1]),event=e,state=state)]


def main():
    output, receipt_path = P/'integration_results.json', P/'integration_receipt.json'
    if output.exists() or receipt_path.exists():
        raise SystemExit('Existing integration receipt/checkpoint: no automatic replay')
    tic=time.monotonic()
    def timeout(signum,frame):
        raise TimeoutError('180-second integration fixture limit')
    signal.signal(signal.SIGALRM,timeout); signal.setitimer(signal.ITIMER_REAL,180)
    receipt=dict(status='RUNNING',pid=os.getpid(),started_utc=datetime.now(timezone.utc).isoformat(),
                 command=[sys.executable,'-B',str(Path(__file__).resolve())],timeout_seconds=180,
                 native_calls=0,compiles=0,old_C061_recalculations=0,automatic_retries=0)
    result=dict(status='RUNNING',fixtures=[],native_calls=0,production_calls=0,
                original_C063_formal_replays=0,old_C061_recalculations=0,synthetic_jets_are_kernel_evidence=False,
                global_g_certificate=False,new_C6_bound=False)
    def save():
        result.update(elapsed_seconds=time.monotonic()-tic,actual_Arb32_corner_checker_calls=ARB_CALLS)
        output.write_text(json.dumps(result,indent=2)+'\n')
        receipt['wall_seconds']=time.monotonic()-tic
        receipt_path.write_text(json.dumps(receipt,indent=2)+'\n')
    def record(n,topic,**data):
        assert n==len(result['fixtures'])+1 and n<=8
        result['fixtures'].append(dict(fixture=n,topic=topic,passed=True,**data));save()
    save()
    try:
        c=context('R1-fixed-root'); h=records(c)
        shrink=copy.deepcopy(h[0]); shrink['state']['nodes']['r']=c.base.make_node('r',None,0,[['0','1/2']]*5,0)
        assert rejects(V.verify_records,[shrink],c)
        for field,value in [('target','1/64000'),('model','different_model')]:
            changed=copy.deepcopy(h[0]['state']);changed[field]=value
            assert rejects(V.validate_domain,changed,c)
        changed=copy.deepcopy(h[0]['state']);changed['nodes']['r']['depth']=1
        assert rejects(V.validate_domain,changed,c)
        record(1,'R1 root shrink, target/model and root depth attacks rejected',omitted_witness_in_old_bug=['1']*5)

        c=context('R2-no-refund');h=records(c);initial=copy.deepcopy(h[0]['state']);s,a=reserve(h,c)
        attack=forged_successor(h,dict(kind='RESERVE_NEXT_ONCE'),initial)
        assert rejects(V.verify_records,attack,c)
        assert rejects(V.apply_event,s,dict(kind='RESET_TO_INITIAL'),c)
        assert rejects(V.apply_event,s,dict(kind='RESERVE_NEXT_ONCE',state=initial),c)
        assert s['batch']['native_reservations']==1 and s['inflight']['token']==a['token']
        record(2,'R2 old inflight/attempt/reservation cannot disappear or be reissued',retained_token=a['token'])

        c=context('R3-no-promotion');h=records(c)
        promoted=copy.deepcopy(h[0]['state']);promoted['environment']='MAIN_RELEASED';promoted['release']=None
        assert rejects(V.validate_domain,promoted,c)
        assert rejects(V.verify_records,forged_successor(h,dict(kind='RESERVE_NEXT_ONCE'),promoted),c)
        production=Context();bad=copy.deepcopy(production.source_state);bad['environment']='MAIN_RELEASED'
        assert rejects(V.validate_domain,bad,production)
        assert rejects(V.apply_event,h[0]['state'],dict(kind='AUTHORIZE_ONCE',release_path='unused',release_sha256='0'*64),c)
        record(3,'R3 test-to-production and MAIN_RELEASED without actual reviewed release rejected')

        c=context('lifecycle-topology');h=records(c);s,a,rebuilt=complete_fixture(h,c)
        assert s['nodes']['r']['status']=='split' and s['queue']==['r0','r1']
        damaged=copy.deepcopy(s);damaged['queue'].pop();del damaged['nodes']['r1'];damaged['nodes']['r']['children']=['r0']
        assert rejects(V.verify_records,forged_successor(h,dict(kind='RESERVE_NEXT_ONCE'),damaged),c)
        orphan=context('orphan-stop');oh=records(orphan);os_,oa=reserve(oh,orphan)
        folder=orphan.attempt_directory(oa['token']);folder.mkdir(parents=True)
        failure=dict(token=oa['token'],input_sha256=oa['input_sha256'],synthetic_fixture=True,
                     reason='UNKNOWN_ORPHAN_STOP',native_launch_observed=None,retained_artifacts=[])
        path=folder/'failure.json';write_new(path,failure)
        os_=event(oh,orphan,dict(kind='HALT_ORIGINAL_ATTEMPT_NO_RETRY',token=oa['token'],failure_path=str(path),failure_sha256=sha(path)))
        assert os_['nodes']['r']['status']=='blocked' and os_['batch']['native_reservations']==1
        assert rejects(V.apply_event,os_,dict(kind='RESERVE_NEXT_ONCE'),orphan)
        record(4,'child/queue topology and unknown-orphan retained stop; no reset/requeue')

        c=context('native-contract');h=records(c);s,a=reserve(h,c);folder,native,nr=make_native_fixture(c,a)
        args=(c,a['input'],a['token'],native,nr,sha(folder/'native.json'),sha(folder/'stderr.txt'))
        checked=M.validate_native(*args);assert len(checked['physical_gradients'])==6
        failures=0
        for key,value in [('gradient',native['gradient'][:5]),('precision_bits',384),('full_group_terms',967678),
                          ('confluent_terms',967682),('zero_node_terms',1)]:
            wrong=copy.deepcopy(native);wrong[key]=value
            assert rejects(M.validate_native,c,a['input'],a['token'],wrong,nr,sha(folder/'native.json'),sha(folder/'stderr.txt'));failures+=1
        wrong=copy.deepcopy(native);wrong['gradient'][0]=raw(Q(1,2**20))
        assert rejects(M.validate_native,c,a['input'],a['token'],wrong,nr,sha(folder/'native.json'),sha(folder/'stderr.txt'));failures+=1
        for key,value in [('source_sha256','f'*64),('command',nr['command'][:-1]+['123']),('window',{'kind':'B','a':'1/64','T':'1/2'})]:
            wrong=copy.deepcopy(nr);wrong[key]=value
            assert rejects(M.validate_native,c,a['input'],a['token'],native,wrong,sha(folder/'native.json'),sha(folder/'stderr.txt'));failures+=1
        assert rejects(M.dyadic,dict(lower_mantissa='nan',upper_mantissa='1',common_binary_exponent='0'));failures+=1
        record(5,'complete native all6grad/param/source/dyadic/count/zero-repeat/exact-Euler binding',negative_mutations_rejected=failures)

        c=context('forged-math-pass');h=records(c);s,a=reserve(h,c);folder,_,_=make_native_fixture(c,a)
        s=attach(h,c,a);rebuilt=M.reconstruct(c,s,a['token'],folder)
        assert rebuilt['observation']['box_certificate']['status']=='UNRESOLVED'
        forged=copy.deepcopy(rebuilt);cert=forged['observation']['box_certificate']
        cert['status']='PASS_CONDITIONAL_ON_MATCHING_CENTER_JET'
        cert['tangent_bound_lower']=str(2*c.base.TARGET)
        for row in cert['corners']:
            row['certificate_lower']=str(2*c.base.TARGET);row['certificate_upper']=str(2*c.base.TARGET);row['passes_target']=True
        forged['observation']['box_certificate_sha256']=c.base.digest(cert)
        path=folder/'verification.json';write_new(path,forged)
        assert rejects(c.verified_observation,s,a['token'],path,sha(path))
        assert rejects(V.apply_event,s,dict(kind='IMPORT_ORIGINAL_RESULT_ONCE',token=a['token'],verification_path=str(path),verification_sha256=sha(path)),c)
        record(6,'forged PASS and self-consistent fake endpoint hashes rejected by real Arb reconstruction',real_checker_precision_bits=384,corners_per_call=32)

        c=context('eight-budget');h=records(c)
        for _ in range(8):
            s,a,rb=complete_fixture(h,c)
        assert len(s['batch']['started_nodes'])==8 and s['batch']['native_reservations']==8
        assert rejects(V.apply_event,s,dict(kind='RESERVE_NEXT_ONCE'),c)
        assert V.verify_records(h,c)['state']==s
        changed=copy.deepcopy(s);changed['batch']['native_reservations']=0
        assert rejects(V.verify_records,forged_successor(h,dict(kind='RESERVE_NEXT_ONCE'),changed),c)
        passed=list(s['passed_leaves'])
        assert passed
        for node_id in passed:
            assert s['nodes'][node_id]['proof'] is not None
        record(7,'eight distinct nodes/eight cumulative reservations and passed-leaf proof history',started_nodes=s['batch']['started_nodes'],passed_fixture_leaves=passed)

        c=context('durable-original');directory=c.root/'state'
        g=V.initialize_preparation_genesis(directory,c)
        r=V.append_event(directory,c,g,dict(kind='RESERVE_NEXT_ONCE'))
        restored=V.recover(directory,c);s=restored['state'];a=s['nodes']['r']['attempts'][-1]
        assert restored['recovery_action']=='OBSERVE_ORIGINAL_HANDLE_OR_HALT_UNKNOWN_NEVER_RELAUNCH'
        assert rejects(V.append_event,directory,c,r,dict(kind='RESERVE_NEXT_ONCE'))
        folder,_,_=make_native_fixture(c,a)
        path=folder/'original_handle.json'
        r=V.append_event(directory,c,r,dict(kind='ATTACH_ORIGINAL_HANDLE_ONCE',token=a['token'],handle_path=str(path),handle_sha256=sha(path)))
        s=V.recover(directory,c)['state'];rebuilt=M.reconstruct(c,s,a['token'],folder)
        vp=folder/'verification.json';write_new(vp,rebuilt)
        # The native fixture already exists; recovery imports it, with no second producer call.
        r=V.append_event(directory,c,r,dict(kind='IMPORT_ORIGINAL_RESULT_ONCE',token=a['token'],verification_path=str(vp),verification_sha256=sha(vp)))
        final=V.recover(directory,c);assert final['state']['inflight'] is None
        assert rejects(V.append_event,directory,c,r,dict(kind='IMPORT_ORIGINAL_RESULT_ONCE',token=a['token'],verification_path=str(vp),verification_sha256=sha(vp)))
        (directory/'unfinished.pending').write_text('not a snapshot\n')
        assert V.recover(directory,c)['pending_files']==['unfinished.pending']
        assert rejects(V.initialize_preparation_genesis,directory,c)
        record(8,'durable typed-event recovery imports original once; duplicate, reset and orphan redispatch rejected',published_sequence=final['sequence'])

        production=Context()
        genesis_sha=V.initialize_preparation_genesis(P/'state',production)
        result.update(status='PASS_EIGHT_TARGETED_INTEGRATION_FIXTURES_ZERO_NATIVE',
                      production_genesis_sha256=genesis_sha,production_genesis_source_sha256=production.source_snapshot_sha256,
                      input_manifest_sha256=production.execution_input_manifest_sha256,
                      initial_queue=['r'],initial_inflight=None,initial_environment='PREPARE_ONLY')
        receipt.update(status='PASS_PREPARATION_ONLY',exit_code=0,ended_utc=datetime.now(timezone.utc).isoformat())
        save();print(result['status'])
    except BaseException as error:
        signal.setitimer(signal.ITIMER_REAL,0)
        result.update(status='FAILED_STOP_NO_AUTOMATIC_RETRY',error=repr(error))
        receipt.update(status='FAILED_STOP_NO_AUTOMATIC_RETRY',exit_code=1,error=repr(error),ended_utc=datetime.now(timezone.utc).isoformat())
        save();raise
    finally:
        signal.setitimer(signal.ITIMER_REAL,0)


if __name__=='__main__':
    main()
