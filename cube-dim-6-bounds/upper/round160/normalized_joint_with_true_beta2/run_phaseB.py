#!/usr/bin/env python3
"""Separately gated fifteen-worker new normalized H certificate; zero optimizer."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations_with_replacement,islice
from math import comb
from concurrent.futures import ProcessPoolExecutor,as_completed
import gzip,hashlib,importlib.util,json,multiprocessing,os,sys,time,traceback
sys.dont_write_bytecode=True
from common import P,read,sha,req,new,canon,coeffs,gated,source_roles
D=P/'phaseB'
def load_joint():
    s=source_roles();path=s['128_joint'];pins=read(P/'STATIC_SOURCE_PINS.json')['sources'];expected=next(r['sha256']for r in pins if r['role']=='128_joint')
    req(sha(path)==expected,'PINNED128_JOINT_SOURCE')
    spec=importlib.util.spec_from_file_location('c160_reused_joint',path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m.Joint
def candidate_preflight(inp):
    reg=read(D/'registration.json');candidate=read(D/'candidate.json');co=read(D/'coefficients.json');sel=read(D/'selection.json')
    for k,n in [('candidate_sha256','candidate.json'),('coefficients_sha256','coefficients.json'),('selection_sha256','selection.json'),('jobs_sha256','jobs.json')]:req(reg[k]==sha(D/n),'FIXED_PHASE_B_'+n)
    result=read(P/'phaseA/result.json');req(sha(P/'phaseA/result.json')==reg['phaseA_result_sha256']==sel['phaseA_result_sha256'],'FROZEN_ACTUAL_PHASE_A')
    req(result['candidate']==candidate and result['candidate_canonical_sha256']==sel['candidate_canonical_sha256']==hashlib.sha256(canon(candidate)).hexdigest(),'SELECTED_SAME_NEW_CANDIDATE')
    req(sel['main_selected']is True and sel['phaseB_prepare_authorized']is True and result['H_power5_below_saved128_float_diagnostic_only']is True,'EVIDENCED_PRIMARY_SELECTION')
    x,caps,a=coeffs(inp['model_input'],candidate);req(co['a_exact']==list(map(str,a))and co['cap_vector_exact']==list(map(str,caps))and co['sum_a_cap_exact']=='1'and co['raw_a_quantized']is False,'EXACT49_ACTUAL_A_AND_COST')
    req(co['canonical_a_sha256']==hashlib.sha256(canon(list(map(str,a)))).hexdigest(),'NEW_A_HASH')
    jobs=read(D/'jobs.json');req(len(jobs)==15 and all(j['job']==i for i,j in enumerate(jobs))and jobs[0]['start']==0 and jobs[-1]['end_exclusive']==735471 and all(j['states']==j['end_exclusive']-j['start']>0 for j in jobs)and all(jobs[i]['end_exclusive']==jobs[i+1]['start']for i in range(14)),'COMPLETE_DISJOINT_PREFIX_JOBS')
    return candidate,co,jobs
def worker(job):
    req(os.getpriority(os.PRIO_PROCESS,0)>=10,'WORKER_NICE')
    reg=read(D/'registration.json');req(sha(P/'inputs.json')==reg['inputs_sha256']and sha(D/'candidate.json')==reg['candidate_sha256'],'WORKER_FROZEN_INPUTS')
    import flint
    req(flint.__version__=='0.9.0'and flint.__FLINT_VERSION__=='3.6.0','WORKER_PINNED_FLINT')
    inp=read(P/'inputs.json');Joint=load_joint();joint=Joint(inp['model_input'],read(D/'candidate.json'))
    records=[];h=hashlib.sha256();counts=[0]*17;mass=[0,0];low=upper=0;first=last=None
    for state in islice(combinations_with_replacement(range(17),8),job['start'],job['end_exclusive']):
        if first is None:first=state
        last=state;h.update(bytes(state));wl,wu=joint.w(state);vl,vu=joint.root(wl,wu)
        pl,pu=joint.probnum(state,0),joint.probnum(state,1);low+=pu*vl;upper+=pl*vu;mass[0]+=pl;mass[1]+=pu;counts[state[-1]]+=1;records.append([str(vl),str(vu)])
    req(len(records)==job['states'],'EXACT_JOB_STATE_COUNT');name='chunks/roots_'+str(job['job']).zfill(2)+'.json.gz';path=D/name
    with path.open('xb')as raw:
        with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0)as f:f.write(canon(records))
        raw.flush();os.fsync(raw.fileno())
    req(path.stat().st_size<100_000_000,'ROOT_FILE_SIZE_CAP')
    result={**job,'first':list(first),'last':list(last),'state_bytes_sha256':h.hexdigest(),'root_file':name,'root_file_sha256':sha(path),'prefix_count_by_K':counts,'lower_at_pU_numerator':str(low),'upper_at_pL_numerator':str(upper),'probability_numerators_by_law':list(map(str,mass)),'integer_tables_sha256':hashlib.sha256(canon(joint.tables())).hexdigest(),'root_cache_info':str(joint.root.cache_info())}
    new(D/'chunks','header_'+str(job['job']).zfill(2)+'.json',result);return result
def run():
    inp,s,reg,release=gated('phaseB');candidate,co,jobs=candidate_preflight(inp)
    new(D,'started.json',{'round':'160','phase':'phaseB','pid':os.getpid(),'workers':15,'nice':os.getpriority(os.PRIO_PROCESS,0),'release_sha256':sha(D/'parent_release.json'),'time_limit':None});(D/'chunks').mkdir(exist_ok=False);start=time.monotonic();cpu=time.process_time();completed=[]
    try:
        with ProcessPoolExecutor(max_workers=15,mp_context=multiprocessing.get_context('spawn'))as pool:
            futures=[pool.submit(worker,j)for j in jobs]
            for future in as_completed(futures):completed.append(future.result())
        completed.sort(key=lambda r:r['job']);req(len(completed)==15,'ALL_WORKERS_EXITED')
        Joint=load_joint();joint=Joint(inp['model_input'],candidate);tables=joint.tables();th=hashlib.sha256(canon(tables)).hexdigest();req(all(r['integer_tables_sha256']==th for r in completed),'SAME_NEW_A_TABLES_IN_ALL_WORKERS');new(D,'integer_tables.json',tables)
        actualA=[tables['a0']]+[v for row in tables['a_by_h']for v in row];req(actualA==co['a_exact'],'CERTIFIER_ACTUAL_A_EQUALS_FROZEN_A')
        digest=hashlib.sha256();countK=[0]*17
        for state in combinations_with_replacement(range(17),8):digest.update(bytes(state));countK[state[-1]]+=1
        req(digest.hexdigest()==inp['expected_state_hash']and countK==[comb(k+7,7)for k in range(17)],'COMPLETE_STATE_ORDER')
        req(all(sum(r['prefix_count_by_K'][k]for r in completed)==countK[k]for k in range(17)),'COMPLETE_PREFIX_COUNTS')
        totalL=sum(int(r['lower_at_pU_numerator'])for r in completed);totalU=sum(int(r['upper_at_pL_numerator'])for r in completed);mass=[sum(int(r['probability_numerators_by_law'][i])for r in completed)for i in range(2)]
        tail=[];tailU=0;tailMass=[0,0]
        for K in range(17,33):
            wl=joint.common[0]+joint.table[-1][K][0];wu=joint.common[1]+joint.table[-1][K][1];vl,vu=joint.root(wl,wu);p=joint.p[0][-1];pp=joint.p[1][-1]
            pl=comb(32,K)*p**K*(joint.S-p)**(32-K);pu=comb(32,K)*pp**K*(joint.S-pp)**(32-K);tailU+=pl*vu;tailMass[0]+=pl;tailMass[1]+=pu
            tail.append({'K':K,'W_units':[str(wl),str(wu)],'v_units':[str(vl),str(vu)],'probability_numerators_by_law':[str(pl),str(pu)]})
        req(all(mass[i]+tailMass[i]==joint.S**32 for i in range(2)),'BOTH_COMPLETE_ENDPOINT_LAWS_UNIT_MASS')
        HL,HU=F(totalL,joint.D),F(totalU+tailU,joint.D);req(0<HL<=HU,'COMPLETE_NEW_H_INTERVAL')
        L5,U5=HL**5,HU**5;theta0=F(inp['model_input']['threshold_M_power6_exact']);thetas=F(inp['threshold_gstar_power6_exact']);old6=F(read(s['128_certificate'])['M_power6_upper'])
        result={'status':'COMPLETE_ONE_NEW_NORMALIZED_FULL_H_CERTIFICATE','candidate':candidate,'raw_a_quantized':False,'canonical_a_sha256':co['canonical_a_sha256'],'normalizing_cap_cost_exact':'1','true_cost_upper_at_most_one':True,'new_beta2_cap_exact':inp['new_beta2_cap_exact'],'source156_result_sha256':release['source156_result_sha256'],'source156_saved_audit':release['source156_saved_audit'],'H_interval':[str(HL),str(HU)],'H_power5_interval':[str(L5),str(U5)],'M_power6_upper':str(U5),'old128_M_power6_upper':str(old6),'threshold_g0_power6_exact':str(theta0),'threshold_gstar_power6_exact':str(thetas),'improves_saved128':U5<old6,'beats_g0_M_slot':U5<theta0,'beats_gstar_M_slot':U5<thetas,'exact_margins_target_minus_upper':{'old128':str(old6-U5),'g0':str(theta0-U5),'gstar':str(thetas-U5)},
          'prefix_states':735471,'prefix_count_by_K':countK,'state_bytes_sha256':digest.hexdigest(),'chunks':completed,'tail_records':tail,'prefix_probability_numerators_by_law':list(map(str,mass)),'tail_probability_numerators_by_law':list(map(str,tailMass)),'probability_denominator_power2':6144,'weighted_sum_denominator_power2':6336,'both_complete_endpoint_laws_unit_mass':True,'logical_root_requests':735487,'workers':15,'arithmetic_bits':192,'tail_discarded':False,'probability_renormalized':False,'optimizer_kernel_erfc_native_density_calls':0,'old_baseline_replays':0,'candidate_M_upper_requires_independent_saved_audit':True,'actual_M_lower':False,'new_C6_bound':False,'global_g_proved':False,'release_sha256':sha(D/'parent_release.json'),'run_sha256':sha(Path(__file__)),'wall_before_receipt_seconds':time.monotonic()-start,'parent_CPU_before_receipt_seconds':time.process_time()-cpu}
        result['full_I_upper_power5_from123']=str(F(3,16)**6*U5/384**5)
        new(D,'candidate_certificate.json',result);new(D,'exit_receipt.json',{'status':result['status'],'certificate_sha256':sha(D/'candidate_certificate.json'),'workers_completed':15,'prefix_states':735471,'tail_rows':16,'logical_root_requests':735487,'optimizer_kernel_erfc_native_density_calls':0,'new_C6_bound':False})
        print(json.dumps({k:result[k]for k in ['status','improves_saved128','beats_g0_M_slot','beats_gstar_M_slot','logical_root_requests','new_C6_bound']}))
    except BaseException:
        new(D,'exit_receipt.json',{'status':'FAILED_STOP_NO_RETRY','failure':traceback.format_exc(),'completed_worker_results':len(completed),'new_C6_bound':False});raise
if __name__=='__main__':
    req(not sys.argv[1:],'usage: run_phaseB.py (after its separate dual release only)');run()
