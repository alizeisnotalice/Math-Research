#!/usr/bin/env python3
"""Finite EPS128 independent values via EPS64 at2x; durable coordinator only."""
import sys
sys.dont_write_bytecode=True
import argparse,json,gzip,os,time,multiprocessing as mp,traceback
from pathlib import Path
from fractions import Fraction as F
from interfaces import REG,HERE,common as C,arith as A


def atomic(path,value):
    p=Path(path);q=Path(str(p)+'.tmp')
    with q.open('wb') as f:f.write(json.dumps(value,sort_keys=True,indent=2).encode()+b'\n');f.flush();os.fsync(f.fileno())
    q.replace(p)


def worker(wid,tasks,results,ack,stop):
    try:
        assert C.file_sha(REG['fan'])==REG['fan_sha256'];fan=json.loads(gzip.decompress(Path(REG['fan']).read_bytes()))
        results.put({'kind':'ready','worker':wid,'pid':os.getpid()})
        while not stop.is_set():
            job=tasks.get()
            if job is None:break
            bid,rows=job
            for row in rows:
                if stop.is_set():return
                began=time.monotonic();assert C.file_sha(row['path'])==row['sha256']
                raw=json.loads(Path(row['path']).read_text());record=row['record'];pair=['1/128',','.join(row['point'])]
                assert raw['stage'] in ('complete','H128_target_obstruction') and raw['partial'] is None
                assert raw['epsilon']==record['epsilon']=='1/128' and raw['pair_cache_key']==pair
                assert raw['record']==record and raw['point']==record['point']==row['point']
                assert raw['kernel_sha256']==record['kernel_sha256']==row['kernel_sha256']
                assert raw['target']==record['qualification_target']==REG['target'] and record['pair_cache_key']==pair and record['model']==row['model']
                doubled=list(map(lambda x:str(2*F(x)),row['point']));assert F(doubled[0])==2
                intervals,geometry=A.direct_values(doubled,fan);saved=[]
                for j,(interval,production) in enumerate(zip(intervals,record['values'])):
                    factor=16 if j==0 else 32;scaled=(interval[0]*factor,interval[1]*factor);L,U=map(F,production)
                    assert L<scaled[0]<=scaled[1]<U,'scaled independent enclosure outside H128 production ball'
                    out=A.dyadic_outward(scaled);assert L<F(out[0])<=F(out[1])<U;saved.append(list(out))
                assert len(saved)==5 and F(saved[0][0])>0
                proof={'index':row['index'],'point':row['point'],'pair_cache_key':pair,'independent_EPS64_input':doubled,
                    'H_multiplier':16,'gradient_multiplier':32,'raw_checkpoint_sha256':row['sha256'],'all5_strict_containment':True,
                    'independent_H128_values_512_dyadic':saved,'independent_geometry':geometry,'worker_pid':os.getpid(),'seconds':time.monotonic()-began}
                results.put({'kind':'point','worker':wid,'block':bid,'proof':proof})
                if ack.get() is not True:return
            results.put({'kind':'block_done','worker':wid,'block':bid})
    except BaseException:
        results.put({'kind':'fault','worker':wid,'pid':os.getpid(),'error':traceback.format_exc()});stop.set()


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--registration',required=True);args=parser.parse_args();began=time.monotonic()
    rp=HERE/args.registration;r=json.loads(rp.read_text());outdir=HERE/r['output_directory'];outdir.mkdir(exist_ok=True)
    locks={**r['input_locks'],str(rp):C.file_sha(rp)}
    for p,h in locks.items():assert C.file_sha(p)==h,('input changed',p)
    for p,h in REG['protected_old_locks'].items():assert C.file_sha(p)==h
    scale=json.loads(Path(r['scaling_receipt']).read_text());assert scale['status']=='EPS128_to_EPS64_doubled_input_exact_binding_qualified_no_values'
    manifest=json.loads(Path(r['point_manifest']).read_text())
    assert manifest['target']==REG['target']=='11/100000' and manifest['epsilon']=='1/128' and manifest['kernel_sha256']==r['kernel_sha256']
    rows=[{**row,'path':row['point_checkpoint'],'sha256':row['point_checkpoint_sha256']} for row in manifest['records']]
    assert len(rows)==r['new_points']==4096
    assert len({','.join(row['point']) for row in rows})==4096 and all(F(row['point'][0])==1 for row in rows)
    assert not {','.join(row['point']) for row in rows}.intersection(r['old_H128_keys'])
    assert {tuple(row['point']) for row in rows}=={tuple(row['point']) for row in json.loads(Path(r['production_tasks']).read_text())}
    plans=[{'block':0,'start':0,'stop_exclusive':1}]
    for bid,offset in enumerate(range(1,len(rows),r['block_size']),1):plans.append({'block':bid,'start':offset,'stop_exclusive':min(offset+r['block_size'],len(rows))})
    states={};pending=[]
    for plan in plans:
        bid=plan['block'];path=outdir/('block%04d.json'%bid)
        state=json.loads(path.read_text()) if path.exists() else {**plan,'points':[],'status':'verified_prefix_only','input_locks':locks}
        assert state['input_locks']==locks
        for offset,proof in enumerate(state['points']):
            i=plan['start']+offset;row=rows[i];assert proof['index']==i and proof['point']==row['point'] and proof['raw_checkpoint_sha256']==row['sha256']
            for ball,production in zip(proof['independent_H128_values_512_dyadic'],row['record']['values']):
                lo,hi=map(F,ball);L,U=map(F,production);assert L<lo<=hi<U
        assert len(state['points'])<=plan['stop_exclusive']-plan['start'];states[bid]=state
        if len(state['points'])<plan['stop_exclusive']-plan['start']:
            pending.append((bid,[{'index':i,'kernel_sha256':r['kernel_sha256'],'model':manifest['model'],**rows[i]} for i in range(plan['start']+len(state['points']),plan['stop_exclusive'])]))
    if not pending:print(json.dumps({'already_verified_all':4096,'new_evaluations':0}),flush=True);return
    index={'task':'audit_g110116_numeric128','status':'pending_values','master_pid':os.getpid(),'input_registration_sha256':C.file_sha(rp),
           'point_manifest_sha256':r['point_manifest_sha256'],'input_locks':locks,'plans':plans,'workers':[],'old_values_recomputed':False,'successful_prefixes_never_repeated':True}
    def save():
        prefix=0
        for plan in plans:
            prefix+=len(states[plan['block']]['points'])
            if len(states[plan['block']]['points'])<plan['stop_exclusive']-plan['start']:break
        index.update(verified_points=sum(len(s['points']) for s in states.values()),global_verified_prefix=prefix,
            block_progress=[{'block':b,'verified':len(s['points']),'sha256':C.file_sha(outdir/('block%04d.json'%b)) if (outdir/('block%04d.json'%b)).exists() else None} for b,s in sorted(states.items())])
        atomic(outdir/'index.json',index)
    ctx=mp.get_context('spawn');tasks=ctx.Queue();results=ctx.Queue();stop=ctx.Event();acks=[ctx.Queue() for _ in range(4)]
    pool=[ctx.Process(target=worker,args=(i,tasks,results,acks[i],stop)) for i in range(4)]
    for process in pool:process.start()
    while len(index['workers'])<4:
        message=results.get();assert message['kind']=='ready',message;index['workers'].append(message)
    save();print(json.dumps({'master_pid':os.getpid(),'worker_pids':[p.pid for p in pool],'finite_points_remaining':sum(len(rows) for _,rows in pending)}),flush=True)
    pilot=next((job for job in pending if job[0]==0),None);sent=set();remaining={b for b,_ in pending};fault=None
    if pilot:tasks.put(pilot);sent.add(0)
    else:
        for job in pending:tasks.put(job);sent.add(job[0])
    while remaining:
        message=results.get();kind=message['kind']
        if kind=='point':
            bid=message['block'];state=states[bid];proof=message['proof'];assert proof['index']==state['start']+len(state['points'])
            state['points'].append(proof);state['status']='complete' if len(state['points'])==state['stop_exclusive']-state['start'] else 'verified_prefix_only'
            atomic(outdir/('block%04d.json'%bid),state);save();acks[message['worker']].put(True)
        elif kind=='block_done':
            remaining.remove(message['block']);print(json.dumps({'block':message['block'],'verified_points':index['verified_points']}),flush=True)
            if message['block']==0:
                index['first_point_quality_passed']=True;save()
                for job in pending:
                    if job[0] not in sent:tasks.put(job);sent.add(job[0])
        else:fault=message;stop.set();break
    if fault:
        index.update(status='fault_preserved_prefix',fault=fault);save()
        for ack in acks:ack.put(False)
    else:
        for _ in pool:tasks.put(None)
    for process in pool:process.join()
    index['pool_exits']=[{'pid':p.pid,'exitcode':p.exitcode,'is_alive':p.is_alive()} for p in pool]
    index.update(status='all4096_H128_plus4grad_strictly_verified' if fault is None else 'fault_preserved_prefix',wall_seconds=time.monotonic()-began)
    for p,h in {**locks,**REG['protected_old_locks']}.items():assert C.file_sha(p)==h,('input changed after audit',p)
    save();assert fault is None,fault;assert index['verified_points']==4096 and all(p.exitcode==0 for p in pool)
    print(json.dumps({k:index[k] for k in ['master_pid','status','verified_points','global_verified_prefix','wall_seconds']}),flush=True)


if __name__=='__main__':main()
