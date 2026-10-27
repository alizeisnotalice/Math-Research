"""One fixed15-worker audit, no wallcap, no automatic restart or new points."""
from pathlib import Path
import argparse,json,multiprocessing as mp,os,queue,resource,time,traceback
from common import read,sha,atomic,locks,fault_cleanup
from worker import execute
from structure import finish
H=Path(__file__).resolve().parent
def stage(phase,specs,plan,O,resume):
    ctx=mp.get_context('spawn');ready=ctx.Queue();results=ctx.Queue();start=ctx.Event();release=ctx.Event();stop=ctx.Event();pilot=ctx.Value('i',-1);ps=[];infos={};out={};fault=None;cleanup=None
    try:
        for s in specs:
            p=ctx.Process(target=execute,args=(phase,s,plan,str(O),ready,results,start,release,pilot,stop,resume),name=f'174-{phase}{s["audit_shard"]:02d}');p.start();ps.append(p)
        while len(infos)<15:
            if stop.is_set():raise RuntimeError('fault before pilot')
            try:r=ready.get(timeout=.1);assert r['sid'] not in infos;infos[r['sid']]=r
            except queue.Empty:pass
            if any(p.exitcode is not None for p in ps):raise RuntimeError('child exited before registration')
        reg={'master_pid':os.getpid(),'phase':phase,'workers':[infos[i] for i in range(15)],'no_math_before_pilot_gate':True};atomic(O/f'{phase}_pool_registration.json',reg);print(json.dumps(reg),flush=True)
        pending=[i for i in range(15) if infos[i]['next_index']<infos[i]['total']];pilot.value=pending[0] if pending else -1
        if not pending:release.set()
        start.set()
        while len(out)<15:
            if stop.is_set():raise RuntimeError('fault stop enters cleanup immediately')
            try:
                r=results.get(timeout=.1)
                if r['type']=='pilot_pass':
                    assert r['sid']==pilot.value and r['next_index']==infos[pilot.value]['next_index']+1 and not release.is_set();atomic(O/f'{phase}_pilot.json',r);print(json.dumps(r),flush=True);release.set()
                elif r['type']=='summary':assert r['sid'] not in out and r['status']=='complete';out[r['sid']]=r
                else:raise RuntimeError(str(r))
            except queue.Empty:
                if all(p.exitcode is not None for p in ps):raise RuntimeError('all exited with missing result')
            if stop.is_set() or any(p.exitcode not in (None,0) for p in ps):raise RuntimeError('worker fault enters cleanup immediately')
        while any(p.is_alive() for p in ps):
            if stop.is_set() or any(p.exitcode not in (None,0) for p in ps):raise RuntimeError('post-summary fault')
            for p in ps:p.join(timeout=.02)
        for p in ps:p.join();assert p.exitcode==0
    except BaseException as e:
        fault={'exception':repr(e),'traceback':traceback.format_exc(),'no_auto_restart':True};cleanup=fault_cleanup(ps,stop,start,release,O/f'{phase}_fault_cleanup.json')
    finally:
        ready.cancel_join_thread();results.cancel_join_thread();ready.close();results.close()
    r={'phase':phase,'complete':fault is None and len(out)==15,'workers':[{'pid':p.pid,'exitcode':p.exitcode,'alive':p.is_alive()} for p in ps],'fault':fault,'cleanup':cleanup};atomic(O/f'{phase}_exit.json',r)
    if any(p.is_alive() for p in ps):os._exit(3)
    assert r['complete'],r
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan-sha',required=True);ap.add_argument('--approval',required=True);ap.add_argument('--authorization',required=True);ap.add_argument('--authorization-sha',required=True);ap.add_argument('--resume-reviewed',action='store_true');a=ap.parse_args();assert a.approval=='audit_bandcurvature174-numeric020-approved-after-review';assert sha(H/'plan.json')==a.plan_sha
    p=read(H/'plan.json');p['self_sha256']=a.plan_sha
    for n,h in p['code_locks'].items():assert sha(H/n)==h
    locks(p['input_locks']);locks(p['protected_locks']);assert sha(a.authorization)==a.authorization_sha;auth=read(a.authorization)
    for k,v in [('plan_sha256',a.plan_sha),('sources',455),('band_products',22295),('pair_events',35872),('actual_corner_cases',34170),('workers',15)]:assert auth[k]==v
    O=Path(p['output_path']);assert O.exists()==a.resume_reviewed
    if a.resume_reviewed:assert auth.get('resume_reviewed') is True and auth.get('successful_prefix_locks');locks(auth['successful_prefix_locks'])
    else:O.mkdir(parents=True)
    assert not (O/'result.json').exists();began=time.monotonic();fault=None;result=None
    atomic(O/'master_registration.json',{'pid':os.getpid(),'plan_sha256':a.plan_sha,'authorization_sha256':a.authorization_sha,'normal_wall_cap':None,'resume':a.resume_reviewed,'finite_counts':[455,22295,35872,34170,15]})
    try:
        stage('mass',p['mass_shards'],p,O,a.resume_reviewed);locks(p['input_locks']);locks(p['protected_locks'])
        for n,h in p['code_locks'].items():assert sha(H/n)==h
        q={}
        for s in p['mass_shards']:
            for cid,r in read(O/'mass'/f'shard{s["audit_shard"]:02d}'/'q_answers.json').items():assert cid not in q;q[cid]=r
        assert len(q)==455
        qp=O/'independent_q_registry.json'
        if qp.exists():assert read(qp)==q
        else:atomic(qp,q)
        stage('corners',p['corner_shards'],p,O,a.resume_reviewed);locks(p['input_locks'])
        for n,h in p['code_locks'].items():assert sha(H/n)==h
        result=finish(p,O);locks(p['input_locks']);locks(p['protected_locks'])
        for n,h in p['code_locks'].items():assert sha(H/n)==h
        atomic(O/'result.json',result)
    except BaseException as e:fault={'exception':repr(e),'traceback':traceback.format_exc(),'prefix_preserved':True,'automatic_restart':False};atomic(O/f'fault{os.getpid()}.json',fault)
    u=resource.getrusage(resource.RUSAGE_SELF);c=resource.getrusage(resource.RUSAGE_CHILDREN);r={'status':'PASS' if result is not None and fault is None else 'FAULT_prefix_preserved','master_pid':os.getpid(),'seconds':time.monotonic()-began,'cpu_user':u.ru_utime,'cpu_system':u.ru_stime,'children_cpu_user':c.ru_utime,'children_cpu_system':c.ru_stime,'peak_rss_bytes':u.ru_maxrss,'fault':fault,'new_H_calls':0,'no_old_success_math_replay':True};atomic(O/f'actual_exit{os.getpid()}.json',r);print(json.dumps(r),flush=True);return 0 if r['status']=='PASS' else 2
if __name__=='__main__':raise SystemExit(main())
