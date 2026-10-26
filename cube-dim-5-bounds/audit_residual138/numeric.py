#!/usr/bin/env python3
"""Finite direct-PF EPS128 audit; explicit review SHA and main authorization."""
import sys
sys.dont_write_bytecode=True
import argparse,gzip,json,multiprocessing as mp,os,queue,time,traceback,resource
from pathlib import Path
from fractions import Fraction as F
import common as C
import pf128 as P
from bind import raw_gate

def worker(wid,tasks,results,ack,stop):
    try:
        os.nice(5);assert C.sha(C.FAN)==C.FAN_SHA;fan=json.loads(gzip.decompress(C.FAN.read_bytes()))
        results.put(dict(kind='ready',worker=wid,pid=os.getpid()))
        while not stop.is_set():
            job=tasks.get()
            if job is None:return
            bid,rows=job
            for item in rows:
                if stop.is_set():return
                entry,task=item['entry'],item['task'];raw_gate(entry,task);began=time.monotonic()
                balls,geo=P.values(entry['point'],fan);saved=[]
                for ball,original in zip(balls,entry['record']['values']):
                    lo,hi=map(F,original);assert lo<ball[0]<=ball[1]<hi,'direct PF enclosure outside production ball'
                    outward=P.dyadic(ball);assert lo<F(outward[0])<=F(outward[1])<hi;saved.append(outward)
                assert len(saved)==5 and F(saved[0][0])>0
                proof=dict(index=item['index'],global_ordinal=entry['global_ordinal'],point=entry['point'],pair_cache_key=entry['record']['pair_cache_key'],
                    raw_sha256=entry['point_checkpoint_sha256'],independent_values=saved,all_five_strictly_contained=True,geometry=geo,worker_pid=os.getpid(),seconds=time.monotonic()-began,
                    route='direct EPS128 integer rational poles;768bit log/pi,512bit outward output',peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
                results.put(dict(kind='point',worker=wid,block=bid,proof=proof))
                if ack.get() is not True:return
            results.put(dict(kind='block_done',worker=wid,block=bid))
    except BaseException:
        results.put(dict(kind='fault',worker=wid,pid=os.getpid(),error=traceback.format_exc()));stop.set()

def release(pool,tasks,acks,stop,fault):
    """The same idle/ACK/fault cleanup is exercised by fault_probe.py."""
    if fault:stop.set()
    for ack in acks:ack.put(False)
    for _ in pool:tasks.put(None)
    started=[p for p in pool if p.pid is not None]
    if fault:
        for p in started:p.join(timeout=2)
        for p in started:
            if p.is_alive():p.terminate()
    for p in started:p.join()
    return [dict(pid=p.pid,exitcode=p.exitcode,alive=p.is_alive()) for p in started]

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--registration',required=True);ap.add_argument('--reviewed-registration-sha',required=True);ap.add_argument('--main-authorization',required=True);a=ap.parse_args()
    assert a.main_authorization=='actualaudit138-approved-after-source-review','main must explicitly authorize actual values'
    rp=Path(a.registration);assert C.sha(rp)==a.reviewed_registration_sha;r=C.read(rp);lock=C.read(C.HERE/'source_lock.json')
    assert C.sha(C.HERE/'source_lock.json')==r['source_lock_sha256'];C.check_locks(lock['code_locks']);C.check_locks(lock['input_locks']);C.check_locks(r['input_locks'])
    assert r['workers']==15;C.assert_pids_absent(r['production_pids']);began=time.monotonic();tasks=C.read(r['tasks']);pm=C.read(r['point_manifest']);entries=pm['records'];assert len(entries)==len(tasks)==r['point_count']
    assert len({tuple(e['point']) for e in entries})==len(entries)
    out=C.SCRATCH/('numeric_batch%04d'%r['batch']);out.mkdir(exist_ok=True)
    locks={**r['input_locks'],str(rp):a.reviewed_registration_sha,**lock['code_locks']}
    plans=[dict(block=0,start=0,end=1)]+[dict(block=i,start=j,end=min(j+128,len(entries))) for i,j in enumerate(range(1,len(entries),128),1)]
    states={};pending=[]
    for plan in plans:
        bid=plan['block'];p=out/('block%04d.json'%bid)
        state=C.read(p) if p.exists() else dict(**plan,points=[],input_locks=locks,status='verified_prefix_only')
        assert state['input_locks']==locks and all(state[k]==plan[k] for k in ('block','start','end'))
        assert len(state['points'])<=plan['end']-plan['start']
        for j,proof in enumerate(state['points']):
            e=entries[plan['start']+j];assert proof['index']==plan['start']+j and proof['point']==e['point'] and proof['raw_sha256']==e['point_checkpoint_sha256']
            assert len(proof['independent_values'])==5
            for ball,original in zip(proof['independent_values'],e['record']['values']):
                lo,hi=map(F,original);aa,bb=map(F,ball);assert lo<aa<=bb<hi
        states[bid]=state
        rows=[dict(index=i,entry=entries[i],task=tasks[i]) for i in range(plan['start']+len(state['points']),plan['end'])]
        if rows:pending.append((bid,rows))
    if not pending:print(json.dumps(dict(status='already_all_verified',new_evaluations=0)));return
    index=dict(task='audit_residual138_numeric',batch=r['batch'],status='running',master_pid=os.getpid(),registration_sha256=a.reviewed_registration_sha,
        main_authorization=a.main_authorization,input_locks=locks,plans=plans,workers=[],old_values_recomputed=False,successful_saved_points_skipped=True,wall_cap=None)
    def save():
        prefix=0
        for plan in plans:
            prefix+=len(states[plan['block']]['points'])
            if len(states[plan['block']]['points'])<plan['end']-plan['start']:break
        index.update(verified_points=sum(len(s['points']) for s in states.values()),global_prefix=prefix,
            blocks=[dict(block=bid,verified=len(s['points']),sha256=C.sha(out/('block%04d.json'%bid)) if (out/('block%04d.json'%bid)).exists() else None) for bid,s in sorted(states.items())])
        C.write(out/'index.json',index)
    ctx=mp.get_context('spawn');tq=ctx.Queue();rq=ctx.Queue();stop=ctx.Event();acks=[ctx.Queue() for _ in range(r['workers'])]
    pool=[ctx.Process(target=worker,args=(i,tq,rq,acks[i],stop),name='PF138-%02d'%i) for i in range(r['workers'])];fault=None
    save();remaining={bid for bid,rows in pending};sent=set();pilot=next((x for x in pending if x[0]==0),None)
    def submit_rest():
        for job in pending:
            if job[0] not in sent:tq.put(job);sent.add(job[0])
    try:
        for p in pool:p.start()
        while len(index['workers'])<r['workers']:
            try:msg=rq.get(timeout=5)
            except queue.Empty:
                assert all(p.is_alive() for p in pool),'worker died during initialization';continue
            assert msg['kind']=='ready',msg;index['workers'].append(msg)
        save();print(json.dumps(dict(master_pid=os.getpid(),worker_pids=[p.pid for p in pool],new_points_remaining=sum(len(rows) for _,rows in pending))),flush=True)
        if pilot:tq.put(pilot);sent.add(0)
        else:submit_rest()
        while remaining:
            try:msg=rq.get(timeout=5)
            except queue.Empty:
                assert all(p.is_alive() for p in pool),'worker exited before completion';continue
            if msg['kind']=='point':
                bid=msg['block'];s=states[bid];proof=msg['proof'];assert proof['index']==s['start']+len(s['points'])
                s['points'].append(proof);s['status']='complete' if len(s['points'])==s['end']-s['start'] else 'verified_prefix_only'
                C.write(out/('block%04d.json'%bid),s);save();acks[msg['worker']].put(True)
            elif msg['kind']=='block_done':
                remaining.remove(msg['block']);print(json.dumps(dict(batch=r['batch'],block=msg['block'],verified=index['verified_points'])),flush=True)
                if msg['block']==0:index['pilot_five_components_passed']=True;save();submit_rest()
            else:raise RuntimeError(msg)
    except BaseException:
        fault=traceback.format_exc();index.update(status='fault_saved_prefix',fault=fault);save()
    finally:
        # No fault can leave idle workers blocked forever on tasks.get or ACK.
        index['pool_exits']=release(pool,tq,acks,stop,fault)
        index.update(status='all_five_components_strictly_verified' if fault is None else 'fault_saved_prefix',seconds=time.monotonic()-began)
        save()
    assert fault is None,fault
    C.check_locks(lock['code_locks']);C.check_locks(lock['input_locks']);C.check_locks(locks)
    assert index['verified_points']==index['global_prefix']==len(entries) and all(p.exitcode==0 for p in pool)
    print(json.dumps(dict(status=index['status'],verified=index['verified_points'],seconds=index['seconds'])),flush=True)

if __name__=='__main__':main()
