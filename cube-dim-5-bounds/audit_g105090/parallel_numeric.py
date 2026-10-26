#!/usr/bin/env python3
"""12 finite workers, fixed blocks, one durable coordinator, no old values."""
import sys
sys.dont_write_bytecode=True
import argparse,json,gzip,os,time,resource,multiprocessing as mp,traceback
from pathlib import Path
from fractions import Fraction as F
from queue import Empty
from interfaces import REG,HERE,common,arith


def atomic(path,obj):
    path=Path(path);temp=Path(str(path)+'.tmp')
    with temp.open('wb') as f:
        f.write(json.dumps(obj,sort_keys=True,indent=2).encode()+b'\n');f.flush();os.fsync(f.fileno())
    temp.replace(path)


def worker(worker_id,tasks,results,ack,stop):
    try:
        assert common.file_sha(REG['fan'])==REG['fan_sha256']
        fan=json.loads(gzip.decompress(Path(REG['fan']).read_bytes()))
        results.put({'kind':'ready','worker':worker_id,'pid':os.getpid()})
        while not stop.is_set():
            task=tasks.get()
            if task is None:break
            bid,jobs=task
            for job in jobs:
                if stop.is_set():return
                began=time.monotonic()
                assert common.file_sha(job['path'])==job['sha256'],'raw point SHA mismatch'
                raw=json.loads(Path(job['path']).read_text())
                assert raw['stage']=='complete' and raw['partial'] is None and raw['point']==job['point'] and raw['record']==job['record']
                values,geometry=arith.direct_values(job['point'],fan);saved=[]
                for interval,production in zip(values,job['record']['values']):
                    L,U=map(F,production);assert L<interval[0]<=interval[1]<U,'independent enclosure outside production ball'
                    out=arith.dyadic_outward(interval);assert L<F(out[0])<=F(out[1])<U,'outward512 saved ball not strictly contained';saved.append(list(out))
                item={'index':job['index'],'point_manifest_index':job['index'],'point':job['point'],'raw_checkpoint_sha256':job['sha256'],
                      'independent_values_512_dyadic':saved,'all5_strict_containment':True,'geometry':geometry,'worker_pid':os.getpid(),'seconds':time.monotonic()-began}
                results.put({'kind':'point','worker':worker_id,'block':bid,'item':item})
                if ack.get() is not True:return
            results.put({'kind':'block_done','worker':worker_id,'block':bid})
    except BaseException:
        results.put({'kind':'fault','worker':worker_id,'pid':os.getpid(),'error':traceback.format_exc()});stop.set()


def old_keys():
    assert common.file_sha(REG['source086'])==REG['source086_sha256']
    reader=common.Stream(REG['source086'])
    while True:
        key=reader.key();assert key is not None
        if key=='cache':
            keys={k for k,_ in reader.items(True)};reader.f.close();return keys
        reader.value()


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--batch',type=int,required=True);args=parser.parse_args()
    began=time.monotonic();outdir=HERE/('numeric_batch%04d'%args.batch);br=json.loads((outdir/'input_registration.json').read_text())
    assert common.file_sha(br['point_manifest'])==br['point_manifest_sha256'];manifest=json.loads(Path(br['point_manifest']).read_text())
    assert manifest['source086_sha256']==REG['source086_sha256'] and manifest['target']==REG['target'] and manifest['epsilon']==REG['epsilon']
    points=manifest['points'];assert len(points)==manifest['new_points']==br['new_points']
    source_keys=old_keys();assert len(source_keys)==180040
    assert len({tuple(p['point']) for p in points})==len(points) and all(','.join(p['point']) not in source_keys for p in points),'would recompute accepted old anchor'
    code_locks={str((HERE/n).resolve()):common.file_sha(HERE/n) for n in ['registration.json','interfaces.py','parallel_numeric.py']}
    code_locks.update({REG['accepted081_arith']:REG['accepted081_arith_sha256'],REG['accepted081_common']:REG['accepted081_common_sha256'],str((outdir/'input_registration.json').resolve()):common.file_sha(outdir/'input_registration.json')})
    plans=[];states={};pending=[];blocksize=br['block_size']
    for bid,start in enumerate(range(0,len(points),blocksize)):
        end=min(start+blocksize,len(points));plans.append({'block':bid,'start':start,'stop_exclusive':end})
        path=outdir/('block%04d.json'%bid)
        if path.exists():
            state=json.loads(path.read_text());assert state['code_locks']==code_locks and state['point_manifest_sha256']==br['point_manifest_sha256']
        else:state={'block':bid,'start':start,'stop_exclusive':end,'code_locks':code_locks,'point_manifest_sha256':br['point_manifest_sha256'],'points':[],'status':'verified_prefix_only'}
        assert len(state['points'])<=end-start
        for offset,proof in enumerate(state['points']):
            idx=start+offset;assert proof['index']==idx and proof['point']==points[idx]['point'] and proof['raw_checkpoint_sha256']==points[idx]['sha256']
            for interval,production in zip(proof['independent_values_512_dyadic'],points[idx]['record']['values']):
                a,b=map(F,interval);L,U=map(F,production);assert L<a<=b<U
        states[bid]=state
        if len(state['points'])<end-start:pending.append((bid,[{'index':i,**points[i]} for i in range(start+len(state['points']),end)]))
    registration={'task':'audit_g105090_numeric','batch':args.batch,'new_points':len(points),'old180040_recomputed':False,'workers':12,'block_size':blocksize,'blocks':plans,'code_locks':code_locks,'point_manifest_sha256':br['point_manifest_sha256'],'wall_deadline':None,'alarm':False,'writer':'single coordinator fsync atomic; per-point ACK; block and global success prefixes'}
    rp=outdir/'registration.json'
    if rp.exists():assert json.loads(rp.read_text())==registration
    else:atomic(rp,registration)
    if not pending:
        print(json.dumps({'batch':args.batch,'already_verified_all':len(points),'new_evaluations':0}),flush=True);return
    index={'task':'audit_g105090_numeric','batch':args.batch,'master_pid':os.getpid(),'status':'pending_values','new_points':len(points),'blocks':plans,'workers':[],'code_locks':code_locks,'point_manifest_sha256':br['point_manifest_sha256'],'successful_prefixes_never_repeated':True,'old180040_recomputed':False}
    def save_index():
        verified=sum(len(s['points']) for s in states.values());prefix=0
        for plan in plans:
            s=states[plan['block']];prefix+=len(s['points'])
            if len(s['points'])<plan['stop_exclusive']-plan['start']:break
        index.update({'verified_points':verified,'global_contiguous_verified_prefix':prefix,'block_progress':[{'block':bid,'verified':len(s['points']),'sha256':common.file_sha(outdir/('block%04d.json'%bid)) if (outdir/('block%04d.json'%bid)).exists() else None} for bid,s in sorted(states.items())]});atomic(outdir/'index.json',index)
    context=mp.get_context('spawn');tasks=context.Queue();results=context.Queue();stop=context.Event();acks=[context.Queue() for _ in range(12)]
    workers=[context.Process(target=worker,args=(i,tasks,results,acks[i],stop)) for i in range(12)]
    for process in workers:process.start()
    ready=set()
    while len(ready)<12:
        message=results.get();assert message['kind']=='ready',message;ready.add(message['worker']);index['workers'].append(message)
    save_index();print(json.dumps({'master_pid':os.getpid(),'workers':[p.pid for p in workers],'finite_blocks':len(pending),'new_points_remaining':sum(len(j) for _,j in pending)}),flush=True)
    for task in pending:tasks.put(task)
    done=0;fault=None
    while done<len(pending):
        message=results.get();kind=message['kind']
        if kind=='point':
            bid=message['block'];state=states[bid];proof=message['item'];assert proof['index']==state['start']+len(state['points'])
            state['points'].append(proof);state['status']='complete' if len(state['points'])==state['stop_exclusive']-state['start'] else 'verified_prefix_only'
            atomic(outdir/('block%04d.json'%bid),state);save_index();acks[message['worker']].put(True)
        elif kind=='block_done':
            done+=1;print(json.dumps({'batch':args.batch,'complete_block':message['block'],'verified_points':index['verified_points'],'global_prefix':index['global_contiguous_verified_prefix']}),flush=True)
        else:fault=message;stop.set();break
    if fault:
        index['status']='fault_preserved_prefix';index['fault']=fault;save_index()
        for a in acks:a.put(False)
    else:
        for _ in workers:tasks.put(None)
    for process in workers:process.join()
    index['pool_exits']=[{'pid':p.pid,'exitcode':p.exitcode,'is_alive':p.is_alive()} for p in workers]
    index['status']='numeric_complete_all_new_H_plus4grad' if fault is None else 'fault_preserved_prefix'
    index['wall_seconds']=time.monotonic()-began;index['peak_master_rss_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    for path,h in REG['protected_old_locks'].items():assert common.file_sha(path)==h,'protected old file changed'
    save_index();assert fault is None,fault;assert all(p.exitcode==0 for p in workers)
    print(json.dumps({k:index[k] for k in ['batch','master_pid','status','verified_points','global_contiguous_verified_prefix','wall_seconds']}),flush=True)


if __name__=='__main__':main()
