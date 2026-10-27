#!/usr/bin/env python3
"""One dual-gated15-worker raw-H certificate and exact true-cost budget."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations_with_replacement,islice
from math import comb
from concurrent.futures import ProcessPoolExecutor,as_completed
import gzip,hashlib,importlib.util,json,multiprocessing,os,sys,time,traceback
sys.set_int_max_str_digits(0);sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def require(ok,msg):
    if not ok:raise ArithmeticError(msg)
def canon(x):return json.dumps(x,sort_keys=True,separators=(',',':')).encode()
def write_new(n,x):
    with(P/n).open('x')as f:json.dump(x,f,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
def load_joint():
    manifest=read(P/'source_manifest.json');path=Path(manifest['joint_core_path'])
    expected=next(r['sha256']for r in manifest['sources']if r['path']==str(path));require(sha(path)==expected,'PINNED128_JOINT_SOURCE')
    spec=importlib.util.spec_from_file_location('c146_reused_joint',path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m.Joint
def release_requirements():
    return {'round':'146','authorized':True,'scope':'ONE_FIXED144_RAW_CANDIDATE_H_AND_TRUE_COST','workers':15,'arithmetic_bits':192,'prefix_states':735471,'tail_rows':16,'logical_root_requests':735487,'raw_a_quantized':False,'nice_minimum':10,'time_or_CPU_limit':None,'optimizer_kernel_erfc_native_calls_permitted':0,'run_sha256':sha(P/'driver.py'),'proof_sha256':sha(P/'PROOF.md'),'registration_sha256':sha(P/'registration.json'),'inputs_sha256':sha(P/'inputs.json'),'jobs_sha256':sha(P/'jobs.json'),'source_manifest_sha256':sha(P/'source_manifest.json')}
def preflight():
    reg=read(P/'registration.json');inp=read(P/'inputs.json');manifest=read(P/'source_manifest.json');jobs=read(P/'jobs.json')
    for field,name in [('run_sha256','driver.py'),('proof_sha256','PROOF.md'),('inputs_sha256','inputs.json'),('jobs_sha256','jobs.json'),('source_manifest_sha256','source_manifest.json')]:require(reg[field]==sha(P/name),'FIXED_PREPARED_'+name)
    for row in manifest['sources']:require(sha(row['path'])==row['sha256'],'FROZEN_SOURCE:'+row['path'])
    x=[F(v,int(inp['candidate']['denominator']))for v in inp['candidate']['numerators']];m=inp['raw_model_input'];raw=[F(m['beta0'])]+[F(v)for row in m['beta_by_h']for v in row];a=[v/b for v,b in zip(x,raw)]
    require(len(x)==49 and sum(x)==1 and x[0]>0 and min(x)>=0,'EXACT_NEW_SIMPLEX')
    require(list(map(str,a))==inp['raw_coefficients_exact'],'UNQUANTIZED_ACTUAL_RAW_A')
    caps=list(map(F,inp['valid_cap_vector_exact']));C=sum(v*b for v,b in zip(a,caps))
    require(C==F(inp['current_cost_exact'])==1-x[32]+x[32]*F(inp['current_valid_beta2_upper'])/F(inp['bL_definition_only']),'ACTUAL49_TERM_COST')
    require(len(jobs)==15 and all(j['job']==i for i,j in enumerate(jobs))and jobs[0]['start']==0 and jobs[-1]['end_exclusive']==735471 and all(j['states']==j['end_exclusive']-j['start']for j in jobs)and all(jobs[i]['end_exclusive']==jobs[i+1]['start']for i in range(14)),'DISJOINT_COMPLETE_JOBS')
    return inp,jobs
def worker(job):
    require(os.getpriority(os.PRIO_PROCESS,0)>=10,'WORKER_NICE')
    require(sha(P/'inputs.json')==read(P/'registration.json')['inputs_sha256'],'WORKER_FROZEN_INPUTS')
    import flint
    require(flint.__version__=='0.9.0'and flint.__FLINT_VERSION__=='3.6.0','WORKER_PINNED_FLINT')
    inp=read(P/'inputs.json');Joint=load_joint();joint=Joint(inp['raw_model_input'],inp['candidate'])
    records=[];h=hashlib.sha256();counts=[0]*17;mass=[0,0];low=upper=0;first=last=None
    for state in islice(combinations_with_replacement(range(17),8),job['start'],job['end_exclusive']):
        if first is None:first=state
        last=state;h.update(bytes(state));wl,wu=joint.w(state);vl,vu=joint.root(wl,wu)
        pl,pu=joint.probnum(state,0),joint.probnum(state,1);low+=pu*vl;upper+=pl*vu;mass[0]+=pl;mass[1]+=pu;counts[state[-1]]+=1;records.append([str(vl),str(vu)])
    require(len(records)==job['states'],'EXACT_JOB_STATE_COUNT')
    name='chunks/roots_'+str(job['job']).zfill(2)+'.json.gz';path=P/name
    with path.open('xb')as raw:
        with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0)as f:f.write(canon(records))
        raw.flush();os.fsync(raw.fileno())
    require(path.stat().st_size<100_000_000,'ROOT_FILE_SIZE_CAP')
    result={**job,'first':list(first),'last':list(last),'state_bytes_sha256':h.hexdigest(),'root_file':name,'root_file_sha256':sha(path),'prefix_count_by_K':counts,'lower_at_pU_numerator':str(low),'upper_at_pL_numerator':str(upper),'probability_numerators_by_law':list(map(str,mass)),'integer_tables_sha256':hashlib.sha256(canon(joint.tables())).hexdigest(),'root_cache_info':str(joint.root.cache_info())}
    with(P/'chunks'/('header_'+str(job['job']).zfill(2)+'.json')).open('x')as f:json.dump(result,f,indent=2);f.write('\n')
    return result
def run():
    release=read(P/'parent_release.json')
    for k,v in release_requirements().items():require(release.get(k)==v,'RELEASE:'+k)
    require(release.get('main_approved')is True and release.get('independent_review_approved')is True,'EXACT_HASH_DUAL_GATE')
    require(not(P/'started.json').exists(),'NO_AUTOMATIC_RESTART');os.nice(max(0,10-os.getpriority(os.PRIO_PROCESS,0)))
    require(all(os.environ.get(k)=='1'for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']),'ONE_THREAD_ENV')
    inp,jobs=preflight();write_new('started.json',{'round':'146','pid':os.getpid(),'workers':15,'nice':os.getpriority(os.PRIO_PROCESS,0),'release_sha256':sha(P/'parent_release.json'),'time_limit':None});(P/'chunks').mkdir(exist_ok=False);start=time.monotonic();completed=[]
    try:
        with ProcessPoolExecutor(max_workers=15,mp_context=multiprocessing.get_context('spawn'))as pool:
            futures=[pool.submit(worker,j)for j in jobs]
            for future in as_completed(futures):completed.append(future.result())
        completed.sort(key=lambda r:r['job']);require(len(completed)==15,'ALL_WORKERS_EXITED')
        Joint=load_joint();joint=Joint(inp['raw_model_input'],inp['candidate']);tables=joint.tables();th=hashlib.sha256(canon(tables)).hexdigest();require(all(r['integer_tables_sha256']==th for r in completed),'IDENTICAL_UNQUANTIZED_TABLES')
        write_new('integer_tables.json',tables)
        digest=hashlib.sha256();countK=[0]*17
        for state in combinations_with_replacement(range(17),8):digest.update(bytes(state));countK[state[-1]]+=1
        require(digest.hexdigest()==inp['expected_state_hash']and countK==[comb(k+7,7)for k in range(17)],'COMPLETE_STATE_ORDER')
        require(all(sum(r['prefix_count_by_K'][k]for r in completed)==countK[k]for k in range(17)),'COMPLETE_PREFIX_COUNT')
        totalL=sum(int(r['lower_at_pU_numerator'])for r in completed);totalU=sum(int(r['upper_at_pL_numerator'])for r in completed);mass=[sum(int(r['probability_numerators_by_law'][i])for r in completed)for i in range(2)]
        tail=[];tailU=0;tailMass=[0,0]
        for K in range(17,33):
            wl=joint.common[0]+joint.table[-1][K][0];wu=joint.common[1]+joint.table[-1][K][1];vl,vu=joint.root(wl,wu);p=joint.p[0][-1];pp=joint.p[1][-1]
            pl=comb(32,K)*p**K*(joint.S-p)**(32-K);pu=comb(32,K)*pp**K*(joint.S-pp)**(32-K);tailU+=pl*vu;tailMass[0]+=pl;tailMass[1]+=pu
            tail.append({'K':K,'W_units':[str(wl),str(wu)],'v_units':[str(vl),str(vu)],'probability_numerators_by_law':[str(pl),str(pu)]})
        require(all(mass[i]+tailMass[i]==joint.S**32 for i in range(2)),'BOTH_WHOLE_LAWS_UNIT_MASS')
        HL,HU=F(totalL,joint.D),F(totalU+tailU,joint.D);require(0<HL<=HU,'COMPLETE_NEW_H_INTERVAL')
        L5,U5=HL**5,HU**5;x=F(inp['x_changed_exact']);bL=F(inp['bL_definition_only']);C=F(inp['current_cost_exact']);theta=F(inp['threshold_g0_power6_exact']);old6=F(read(inp['saved128_certificate_path'])['M_power6_upper'])
        def threshold(target,Q):return bL*(target/Q-(1-x))/x
        result={'status':'COMPLETE_ONE_NEW_RAW_H_AND_TRUE_COST_CERTIFICATE','candidate':inp['candidate'],'raw_a_quantized':False,'H_interval':[str(HL),str(HU)],'H_power5_interval':[str(L5),str(U5)],'current_true_cost_upper_exact':str(C),'current_valid_beta2_upper':inp['current_valid_beta2_upper'],'current_valid_M6_power6_upper':str(C*U5),'current_cap_beats_g0_threshold':C*U5<theta,'current_cap_improves_saved128':C*U5<old6,
          'true_beta2_upper_threshold_for_strict_g0':{'sufficient_exact':str(threshold(theta,U5)),'necessary_for_this_raw_candidate_exact':str(threshold(theta,L5))},'true_beta2_upper_threshold_for_strict_saved128_improvement':{'sufficient_exact':str(threshold(old6,U5)),'necessary_for_this_raw_candidate_exact':str(threshold(old6,L5))},
          'prefix_states':735471,'prefix_count_by_K':countK,'state_bytes_sha256':digest.hexdigest(),'chunks':completed,'tail_records':tail,'prefix_probability_numerators_by_law':list(map(str,mass)),'tail_probability_numerators_by_law':list(map(str,tailMass)),'probability_denominator_power2':6144,'weighted_sum_denominator_power2':6336,'both_complete_endpoint_laws_unit_mass':True,'logical_root_requests':735487,'workers':15,'arithmetic_bits':192,'tail_discarded':False,'probability_renormalized':False,'optimizer_kernel_erfc_native_calls':0,'old_baseline_replays':0,'new_C6_bound':False,'global_g0_proved':False,'release_sha256':sha(P/'parent_release.json'),'run_sha256':sha(Path(__file__)),'wall_before_receipt_seconds':time.monotonic()-start}
        result['full_I_upper_power5_from123']=str(F(3,16)**6*C*U5/384**5)
        write_new('candidate_certificate.json',result);write_new('exit_receipt.json',{'status':result['status'],'certificate_sha256':sha(P/'candidate_certificate.json'),'workers_completed':15,'prefix_states':735471,'tail_rows':16,'logical_root_requests':735487,'optimizer_kernel_erfc_native_calls':0,'new_C6_bound':False})
        print(json.dumps({'status':result['status'],'current_cap_beats_g0_threshold':result['current_cap_beats_g0_threshold'],'current_cap_improves_saved128':result['current_cap_improves_saved128'],'logical_root_requests':735487}))
    except BaseException:
        write_new('exit_receipt.json',{'status':'FAILED_STOP_NO_RETRY','failure':traceback.format_exc(),'completed_worker_results':len(completed),'optimizer_kernel_erfc_native_calls':0,'new_C6_bound':False});raise
if __name__=='__main__':
    if sys.argv[1:]==['--saved-only-preflight']:
        inp,jobs=preflight();print(json.dumps({'status':'SAVED_ONLY_PREFLIGHT_PASS','workers':15,'raw_a_quantized':False,'actual_cost_matches49_term_sum':True,'new_H_or_root_calls':0}))
    elif not sys.argv[1:]:run()
    else:raise SystemExit('usage: driver.py [--saved-only-preflight]')
