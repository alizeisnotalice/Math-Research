#!/usr/bin/env python3
"""Nine independently constructed revision checks; no native or CLI runs."""
import signal
signal.alarm(180)
import sys
sys.dont_write_bytecode=True
from pathlib import Path
from fractions import Fraction as Q
from tempfile import TemporaryDirectory
from datetime import datetime,timezone
import copy,hashlib,json,os,time

ROOT=Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def reject(call):
    try:
        call()
    except (ValueError,KeyError,FileNotFoundError,FileExistsError) as error:
        return type(error).__name__+': '+str(error)
    raise AssertionError('Forbidden operation was accepted')


def dyadic(lo,hi=None):
    lo,hi=Q(lo),Q(lo if hi is None else hi)
    denominator=max(lo.denominator,hi.denominator)
    assert denominator & (denominator-1)==0
    return {'lower_mantissa':str(int(lo*denominator)), 'upper_mantissa':str(int(hi*denominator)),
            'common_binary_exponent':str(-(denominator.bit_length()-1))}


def main():
    tic=time.monotonic(); start=datetime.now(timezone.utc).isoformat()
    assert not (ROOT/'results.json').exists(), 'Do not overwrite an accepted audit run'
    pins=json.loads((ROOT/'input_manifest.json').read_text())
    for path,digest in pins['source_sha256'].items(): assert sha(path)==digest,path
    author=Path(pins['author_root'])
    sys.path.insert(0,str(author))
    import core_v2 as V
    import math_contract as M
    from batch_context import Context,write_new,read
    checks=[]; arb_calls=0

    with TemporaryDirectory(prefix='independent_fixtures_',dir=ROOT) as temporary:
        fixture_root=Path(temporary)
        def context(label):
            nonlocal arb_calls
            root=fixture_root/label; root.mkdir()
            c=Context(test_only=True,fixture_root=root)
            original=c.checker.certify
            def counted(*args,**kwargs):
                nonlocal arb_calls
                arb_calls+=1
                return original(*args,**kwargs)
            c.checker.certify=counted
            return c

        def advance(history,c,event):
            state=V.apply_event(history[-1]['state'],event,c)
            history.append({'schema':V.SCHEMA,'sequence':len(history),'previous_sha256':V.digest(history[-1]),
                            'event':copy.deepcopy(event),'state':state})
            return state

        def reserve(history,c):
            state=advance(history,c,{'kind':'RESERVE_NEXT_ONCE'})
            node=state['nodes'][state['inflight']['node_id']]
            return state,node['attempts'][-1]

        def make_originals(c,attempt,value=Q(1,2**14)):
            expected=attempt['input']; token=attempt['token']; folder=c.attempt_directory(token)
            folder.mkdir(parents=True)
            # Mixed-sign, interval-valued physical gradient, including nonzero x1 field.
            gradients=[(Q(n,2**24)-Q(1,2**28),Q(n,2**24)+Q(1,2**28)) for n in (-3,2,-1,4,-2,1)]
            x=list(map(Q,expected['center'])); a=Q(expected['window']['a']); T=Q(1,2); pT=Q(1,2**22)
            pa=tuple((5*value+T*pT+sum(xi*g[j] for xi,g in zip(x,gradients)))/a for j in (0,1))
            assert pa[0]>0
            native={'status':'all_real_nodes_supported','precision_bits':512,'full_group_terms':967680,
                    'computed_sign_pairs':483840,'confluent_terms':2,'zero_node_terms':4,'unique_node_values':7,
                    'kernel_seconds':0,'band':dyadic(value),'gradient':[dyadic(*g) for g in gradients],
                    'p_a':dyadic(*pa),'p_T':dyadic(pT),'Euler_xgrad_plus5B_plusTpT_minusapa':dyadic(-Q(1,2**16),Q(1,2**16))}
            intent={'token':token,'input':expected,'input_sha256':attempt['input_sha256'],
                    'release_sha256':'0'*64,'owned_wrapper_pid':700001,'synthetic_fixture':True}
            write_new(folder/'intent.json',intent);write_new(folder/'native.json',native)
            (folder/'stderr.txt').write_bytes(b'')
            receipt={'schema':'C068-native-receipt-v1','token':token,'input_sha256':attempt['input_sha256'],
                     'command':M.native_command(c,expected),**{k:expected[k] for k in ('model','model_definition','window','center')},
                     'source_sha256':c.external_sha(c.paths['source']),'binary_sha256':c.external_sha(c.paths['binary']),
                     'fan_sha256':c.external_sha(c.paths['fan']),'runtime_sha256':c.runtime_pins,
                     'intent_sha256':sha(folder/'intent.json'),'release_sha256':'0'*64,'native_exit_code':0,
                     'stderr_empty':True,'native_sha256':sha(folder/'native.json'),'stderr_sha256':sha(folder/'stderr.txt'),
                     'owned_wrapper_pid':700001,'owned_native_pid':700002,'timeout_seconds':180,'synthetic_fixture':True}
            write_new(folder/'native_receipt.json',receipt)
            handle={'token':token,'input_sha256':attempt['input_sha256'],'synthetic_fixture':True,
                    'handle':{'fixture':'NO_PROCESS_WAS_LAUNCHED'},'intent_sha256':sha(folder/'intent.json')}
            write_new(folder/'original_handle.json',handle)
            return folder,native,receipt,gradients

        def attach(history,c,attempt):
            path=c.attempt_directory(attempt['token'])/'original_handle.json'
            return advance(history,c,{'kind':'ATTACH_ORIGINAL_HANDLE_ONCE','token':attempt['token'],
                                      'handle_path':str(path),'handle_sha256':sha(path)})

        def complete(history,c,value=Q(1,2**14)):
            state,attempt=reserve(history,c)
            folder,native,receipt,gradients=make_originals(c,attempt,value)
            state=attach(history,c,attempt)
            rebuilt=M.reconstruct(c,state,attempt['token'],folder)
            path=folder/'verification.json';write_new(path,rebuilt)
            state=advance(history,c,{'kind':'IMPORT_ORIGINAL_RESULT_ONCE','token':attempt['token'],
                                     'verification_path':str(path),'verification_sha256':sha(path)})
            return state,attempt,rebuilt

        def rewrite(path,data):
            # Only the newly-created attack fixture is changed, never an author artifact.
            Path(path).write_bytes(V.canonical(data))

        # 1. Actual production Context, read-only unexecuted genesis; no release is created.
        real=Context()
        actual=V.recover(author/'state',real)
        state=actual['state']
        assert actual['sequence']==0 and state['environment']=='PREPARE_ONLY' and state['release'] is None
        assert state['batch']['native_reservations']==state['production_kernel_calls']==0
        assert state['queue']==['r'] and state['inflight'] is None and state['nodes']['r']['raw_box']==V.ROOT_BOX
        no_release=reject(lambda:V.apply_event(state,{'kind':'RESERVE_NEXT_ONCE'},real))
        checks.append({'configuration':1,'status':'PASS_REAL_GENESIS_READONLY','snapshot_sha256':actual['sha256'],
                       'pre_release_reservation_rejection':no_release,'real_release_created':False})

        # 2. Typed-event successor injection for R1/R3, beyond merely testing a standalone validator.
        c=context('state-injection');h=[V.genesis_record(c)];s,a=reserve(h,c)
        failures={}
        for label,mutation in [('shrunk_root',lambda st:st['nodes'].__setitem__('r',c.base.make_node('r',None,0,[['0','1/3']]*5,0))),
                               ('promoted_test_environment',lambda st:st.update(environment='MAIN_RELEASED',release=None))]:
            altered=copy.deepcopy(h);mutation(altered[-1]['state'])
            failures[label]=reject(lambda:V.verify_records(altered,c))
        failures['arbitrary_successor_field']=reject(lambda:V.apply_event(h[0]['state'],{'kind':'RESERVE_NEXT_ONCE','state':s},c))
        failures['test_authorize']=reject(lambda:V.apply_event(h[0]['state'],{'kind':'AUTHORIZE_ONCE','release_path':'unused','release_sha256':'0'*64},c))
        checks.append({'configuration':2,'status':'PASS_R1_R3_REPLAY_AND_EVENT_API','rejections':failures})

        # 3. Actual new persistence API, stale writer and retained orphan stop.
        c=context('durable-orphan');directory=c.root/'state'
        h0=V.initialize_preparation_genesis(directory,c)
        h1=V.append_event(directory,c,h0,{'kind':'RESERVE_NEXT_ONCE'})
        latest=V.recover(directory,c);s=latest['state'];token=s['inflight']['token'];attempt=s['nodes']['r']['attempts'][-1]
        failures={'stale_writer':reject(lambda:V.append_event(directory,c,h0,{'kind':'RESERVE_NEXT_ONCE'})),
                  'reset_event':reject(lambda:V.append_event(directory,c,h1,{'kind':'RESET','state':c.source_state})),
                  'inflight_reissue':reject(lambda:V.append_event(directory,c,h1,{'kind':'RESERVE_NEXT_ONCE'})),
                  'reinitialize':reject(lambda:V.initialize_preparation_genesis(directory,c))}
        folder=c.attempt_directory(token);folder.mkdir(parents=True)
        failure=folder/'failure.json'
        write_new(failure,{'token':token,'input_sha256':attempt['input_sha256'],'synthetic_fixture':True,
                           'reason':'UNKNOWN_ORPHAN_STOP','native_launch_observed':None})
        h2=V.append_event(directory,c,h1,{'kind':'HALT_ORIGINAL_ATTEMPT_NO_RETRY','token':token,
                                         'failure_path':str(failure),'failure_sha256':sha(failure)})
        stopped=V.recover(directory,c)
        assert stopped['state']['batch']['native_reservations']==1 and stopped['state']['nodes']['r']['attempts'][0]['token']==token
        failures['post_orphan_dispatch']=reject(lambda:V.append_event(directory,c,h2,{'kind':'RESERVE_NEXT_ONCE'}))
        checks.append({'configuration':3,'status':'PASS_R2_DURABLE_ORPHAN_MONOTONICITY','rejections':failures,
                       'retained_original_token':token,'retained_reservations':1,'recovery_action':stopped['recovery_action']})

        # 4. Same node, two genuine reconstruction/import operations at distinct synthetic cutoffs.
        c=context('cutoff-binding');h=[V.genesis_record(c)]
        s,a0,r0=complete(h,c,Q(1,2**18));s,a1,r1=complete(h,c,Q(1,2**18))
        assert a0['input']['window']['a']=='1/32' and a1['input']['window']['a']=='1/64'
        assert a0['token']!=a1['token'] and a0['input_sha256']!=a1['input_sha256']
        assert len(s['nodes']['r']['attempts'])==2 and s['nodes']['r']['cutoff_index']==2
        folder=c.attempt_directory(a1['token']);native=read(folder/'native.json');receipt=read(folder/'native_receipt.json')
        bad=copy.deepcopy(receipt);bad['command']=M.native_command(c,a0['input'])
        e_command=reject(lambda:M.validate_native(c,a1['input'],a1['token'],native,bad,sha(folder/'native.json'),sha(folder/'stderr.txt')))
        bad=copy.deepcopy(receipt);bad['window']=a0['input']['window']
        e_window=reject(lambda:M.validate_native(c,a1['input'],a1['token'],native,bad,sha(folder/'native.json'),sha(folder/'stderr.txt')))
        checks.append({'configuration':4,'status':'PASS_NODE_WINDOW_UNIQUENESS','tokens':[a0['token'],a1['token']],
                       'old_command_rejection':e_command,'old_window_rejection':e_window,'normal_gradient_retained_in_both':True})

        # 5. Eight cumulative attempts, including two cutoff steps at the same node.
        c=context('mixed-budget');h=[V.genesis_record(c)]
        ids=[]
        for i in range(8):
            s,a,r=complete(h,c,Q(1,2**18) if i<2 else Q(1,2**14));ids.append(a['token'].split('.')[0])
        assert s['batch']['native_reservations']==8 and len(set(ids))==6 and len(s['batch']['started_nodes'])==6
        error=reject(lambda:V.apply_event(s,{'kind':'RESERVE_NEXT_ONCE'},c))
        checked=V.verify_records(h,c)
        assert checked['state']==s
        checks.append({'configuration':5,'status':'PASS_MIXED_CUTOFF_TOTAL_ATTEMPT_CAP','attempt_node_ids':ids,
                       'reservations':8,'distinct_nodes':6,'ninth_rejection':error,'full_history_replayed':True})

        # 6. Topology corruption and a real-Arb passed synthetic child.
        c=context('topology-passed');h=[V.genesis_record(c)]
        s,a,r=complete(h,c)
        assert s['nodes']['r']['status']=='split'
        altered=copy.deepcopy(h);del altered[-1]['state']['nodes']['r1'];altered[-1]['state']['queue'].remove('r1')
        missing=reject(lambda:V.verify_records(altered,c))
        altered=copy.deepcopy(h);node=altered[-1]['state']['nodes']['r0'];node['status']='passed';node['proof']={'fake':True}
        altered[-1]['state']['queue'].remove('r0');altered[-1]['state']['passed_leaves'].append('r0')
        counterfeit=reject(lambda:V.verify_records(altered,c))
        s,passed_attempt,passed_result=complete(h,c,Q(2)**80)
        passed_id=passed_attempt['token'].split('.')[0]
        assert s['nodes'][passed_id]['status']=='passed' and not s['global_g_certificate']
        altered=copy.deepcopy(h);altered[-1]['state']['nodes'][passed_id]['proof']['token']='changed'
        changed=reject(lambda:V.verify_records(altered,c))
        duplicate=reject(lambda:V.apply_event(s,h[-1]['event'],c))
        checks.append({'configuration':6,'status':'PASS_TOPOLOGY_AND_PASSED_HISTORY','missing_sibling_rejection':missing,
                       'counterfeit_pass_rejection':counterfeit,'passed_node_change_rejection':changed,
                       'duplicate_import_rejection':duplicate,'passed_synthetic_child':passed_id,'global_g_certificate':False})

        # 7. Independent interval gradient fixture and exact rational corner checks.
        c=context('six-gradients');h=[V.genesis_record(c)];s,a=reserve(h,c)
        folder,native,receipt,grad=make_originals(c,a);s=attach(h,c,a)
        rebuilt=M.reconstruct(c,s,a['token'],folder)
        assert rebuilt['complete_physical_gradient_intervals']==[list(map(str,g)) for g in grad]
        assert rebuilt['observation']['jet']['face_gradient_intervals']==[list(map(str,g)) for g in grad[1:]]
        certificate=rebuilt['observation']['box_certificate'];assert certificate['precision_bits']==384 and len(certificate['corners'])==32
        for row in certificate['corners']:
            delta=list(map(Q,row['displacement']));dot=sum(min(g[0]*d,g[1]*d) for g,d in zip(grad[1:],delta))
            f=Q(1,2**14);quadratic=sum(d*d for d in delta)/32
            assert Q(row['gradient_dot_lower'])==dot and Q(row['minimizing_center_value'])==f
            assert Q(row['quadratic'])==quadratic and Q(row['exponent'])==dot/f-quadratic/Q(a['input']['window']['a'])**2
        broken=copy.deepcopy(native);broken['gradient'][0]=dyadic(1)
        rewrite(folder/'native.json',broken);changed_receipt=copy.deepcopy(receipt);changed_receipt['native_sha256']=sha(folder/'native.json')
        rewrite(folder/'native_receipt.json',changed_receipt)
        error=reject(lambda:M.reconstruct(c,s,a['token'],folder))
        assert 'Euler reconstruction' in error
        checks.append({'configuration':7,'status':'PASS_ALL_SIX_INTERVAL_FIELDS_AND_REAL_ARB_CORNERS',
                       'complete_physical_gradient_intervals':[list(map(str,g)) for g in grad],
                       'actual_corner_count':32,'Arb_bits':384,'normal_field_corruption_rejection':error,
                       'native_Euler_still_contains_zero':True})

        # 8. Original artifact replacement after verification and a rehashed false box.
        c=context('originals-and-exp');h=[V.genesis_record(c)];s,a=reserve(h,c)
        folder,native,receipt,grad=make_originals(c,a);s=attach(h,c,a)
        rebuilt=M.reconstruct(c,s,a['token'],folder);path=folder/'verification.json';write_new(path,rebuilt)
        original_bytes={name:(folder/name).read_bytes() for name in ['native.json','native_receipt.json','verification.json']}
        changed=copy.deepcopy(native);changed['confluent_terms']=6;rewrite(folder/'native.json',changed)
        changed_receipt=copy.deepcopy(receipt);changed_receipt['native_sha256']=sha(folder/'native.json');rewrite(folder/'native_receipt.json',changed_receipt)
        changed_original=reject(lambda:c.verified_observation(s,a['token'],str(path),sha(path)))
        for name,raw in original_bytes.items(): (folder/name).write_bytes(raw)
        fake=copy.deepcopy(rebuilt);cert=fake['observation']['box_certificate']
        for row in cert['corners']:
            row['certificate_lower']='1';row['certificate_upper']='2';row['passes_target']=True
        cert['tangent_bound_lower']='1';cert['tangent_bound_upper']='2';cert['status']='PASS_CONDITIONAL_ON_MATCHING_CENTER_JET'
        fake['observation']['box_certificate_sha256']=c.base.digest(cert);rewrite(path,fake)
        false_box=reject(lambda:c.verified_observation(s,a['token'],str(path),sha(path)))
        checks.append({'configuration':8,'status':'PASS_LATE_ORIGINAL_AND_REHASHED_EXP_FORGERY_REJECTION',
                       'changed_original_rejection':changed_original,'fake_box_rejection':false_box,
                       'mutated_files_are_new_review_fixtures_only':True})

        # 9. Unit-exercise production-only mathematical checks without authorizing any real Context.
        c=context('production-contract-unit');h=[V.genesis_record(c)];s,a=reserve(h,c)
        folder,native,receipt,grad=make_originals(c,a)
        view=copy.copy(c);view.test_only=False  # unit branch only, never passed to V.authorize or a runner
        production_state=copy.deepcopy(s);release='4'*64
        production_state['release']={'release_sha256':release}
        production_state['nodes']['r']['attempts'][-1]['execution_handle']={'native_pid':700002,'wrapper_pid':700001}
        intent=read(folder/'intent.json');intent['release_sha256']=release;intent['synthetic_fixture']=False;rewrite(folder/'intent.json',intent)
        process={'token':a['token'],'input_sha256':a['input_sha256'],'owned_native_pid':700002,'owned_wrapper_pid':700001}
        write_new(folder/'original_process.json',process)
        receipt.update(synthetic_fixture=False,release_sha256=release,intent_sha256=sha(folder/'intent.json'),
                       original_process_sha256=sha(folder/'original_process.json'));rewrite(folder/'native_receipt.json',receipt)
        baseline=M.reconstruct(view,production_state,a['token'],folder)
        assert baseline['universal_gradient_crosscheck_claim'] is False
        errors={}
        for label,key,value in [('release','release_sha256','5'*64),('source','source_sha256','0'*64),
                                ('window','window',{'kind':'B','a':'1/64','T':'1/2'})]:
            bad=copy.deepcopy(receipt);bad[key]=value;rewrite(folder/'native_receipt.json',bad)
            errors[label]=reject(lambda:M.reconstruct(view,production_state,a['token'],folder))
        rewrite(folder/'native_receipt.json',receipt)
        bad_process=copy.deepcopy(process);bad_process['owned_native_pid']=700003;rewrite(folder/'original_process.json',bad_process)
        bad=copy.deepcopy(receipt);bad['original_process_sha256']=sha(folder/'original_process.json');rewrite(folder/'native_receipt.json',bad)
        errors['original_pid']=reject(lambda:M.reconstruct(view,production_state,a['token'],folder))
        checks.append({'configuration':9,'status':'PASS_PRODUCTION_MATH_CONTRACT_UNIT_BRANCH','rejections':errors,
                       'fictional_PIDs_never_observed':True,'actual_main_release_created':False,
                       'actual_production_Context_authorized':False,'actual_native_calls':0,
                       'unit_fixture_not_a_real_native_receipt':True})

    assert len(checks)==9
    for path,digest in pins['source_sha256'].items(): assert sha(path)==digest,path
    output={'status':'PASS_C068_REVISION_CORE_AND_MATH_CONTRACT_WITH_MAIN_RELEASE_PENDING',
            'actual_configurations':9,'checks':checks,'actual_Arb384_checker_calls':arb_calls,
            'actual_Arb_corner_evaluations':32*arb_calls,'native_calls':0,'production_calls':0,
            'old_C063_writer_runs':0,'author_formal_program_runs':0,'old_formal_replays':0,
            'matlab_calls':0,'subagents':0,'author_files_unchanged':True,'actual_release_created':False,
            'real_launch_path_tested':False,'real_timeout_kill_path_tested':False,
            'global_g_certificate':False,'actual_body32_cover':False,'new_C6_bound':False,
            'pid':os.getpid(),'python':sys.version,'start_utc':start,'end_utc':datetime.now(timezone.utc).isoformat(),
            'elapsed_seconds':time.monotonic()-tic,'hard_process_limit_seconds':180,
            'source_sha256':{name:sha(ROOT/name) for name in ['registration.md','input_manifest.json','verify_revision.py']}}
    (ROOT/'results.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps({key:output[key] for key in ['status','actual_configurations','actual_Arb384_checker_calls','elapsed_seconds','native_calls']},indent=2))


if __name__=='__main__':
    main()
