#!/usr/bin/env python3
"""Prepared15-explicit-worker transport export,not yet authorized or run.

One fixed shard per process. No task generation or mathematical checker call.
Every completed shard is decoded and compared to its write-time projection.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse,gzip,hashlib,json,os,shutil,subprocess,sys,time
import codec as C
H=Path(__file__).resolve().parent
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def read(p):
    p=Path(p)
    with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:return json.load(f)
def put(p,x):
    with Path(p).open('x') as f:json.dump(x,f,sort_keys=True,separators=(',',':'));f.write('\n');f.flush();os.fsync(f.fileno())
def lines(p):
    with gzip.open(p,'rt') as f:
        for s in f:yield json.loads(s)
def checked(plan_sha):
    assert sha(H/'plan_v2.json')==plan_sha;p=read(H/'plan_v2.json')
    for n,h in p['code_locks'].items():assert sha(H/n)==h
    assert sha(H/'node_index.jsonl.gz')==p['index_sha256']
    return p
def worker(out,sid,plan_sha,resume):
    plan=checked(plan_sha);start=time.monotonic();dest=out/f'worker{sid:02d}';dest.mkdir(exist_ok=resume)
    if not resume:put(dest/'registration.json',{'pid':os.getpid(),'sid':sid,'plan_sha256':plan_sha,'count':plan['shard_counts'][sid],'stage':'registered_before_source_reads'})
    reg=read(dest/'registration.json');assert reg['sid']==sid and reg['plan_sha256']==plan_sha
    assert not (dest/'result.json').exists(),'completed worker must not be restarted'
    rows=[r for r in lines(H/'node_index.jsonl.gz') if r['shard']==sid];assert len(rows)==plan['shard_counts'][sid]
    chunks=dest/'chunks';chunks.mkdir(exist_ok=resume);parts=[];projection=[]
    try:
        for j in range(0,len(rows),256):
            p=chunks/f'chunk{j//256:04d}.gz';rp=p.with_suffix('.json');batch=rows[j:j+256]
            if rp.exists():
                receipt=read(rp);assert receipt['first_id']==batch[0]['id'] and receipt['count']==len(batch) and sha(p)==receipt['sha256']
            else:
                assert not p.exists();tmp=p.with_suffix('.tmp');proof=[]
                with tmp.open('xb') as f:
                    with gzip.GzipFile(fileobj=f,mode='wb',mtime=0,filename='') as z:
                        if j==0:z.write(C.MAGIC+sid.to_bytes(4,'big')+len(rows).to_bytes(4,'big'))
                        for row in batch:
                            src=Path(row['packet_path']);assert sha(src)==row['packet_sha256'];data,ps,origin=C.encode(row,read(src),plan['times']);z.write(data)
                            proof.append({'id':row['id'],'source_gzip_sha256':row['packet_sha256'],'original_raw_sha256':origin,'normalized_projection_sha256':ps})
                    f.flush();os.fsync(f.fileno())
                tmp.rename(p);receipt={'first_id':batch[0]['id'],'count':len(batch),'sha256':sha(p),'projections':proof};put(rp,receipt)
            parts.append(p);projection.extend(receipt['projections'])
        assert len(projection)==len(rows)
        target=out/f'endpoints{sid:02d}.bin.gz'
        if not target.exists():
            with target.with_suffix('.tmp').open('xb') as f:
                for p in parts:
                    assert sha(p)==read(p.with_suffix('.json'))['sha256']
                    with p.open('rb') as src:shutil.copyfileobj(src,f)
                f.flush();os.fsync(f.fileno())
            target.with_suffix('.tmp').rename(target)
        assert target.stat().st_size<plan['shard_hard_limit_bytes'];before_sha=sha(target)
        # Actual transport roundtrip for this newly completed shard only.
        with gzip.open(target,'rb') as f:
            assert C.take(f,8)==C.MAGIC and int.from_bytes(C.take(f,4),'big')==sid and int.from_bytes(C.take(f,4),'big')==len(rows)
            for row,written in zip(rows,projection):
                _,ps,origin=C.decode(f,row,plan['times'])
                assert written=={'id':row['id'],'source_gzip_sha256':row['packet_sha256'],'original_raw_sha256':origin,'normalized_projection_sha256':ps}
            assert f.read(1)==b''
        assert sha(target)==before_sha
        for row in rows:assert sha(row['packet_path'])==row['packet_sha256']
        checked(plan_sha)
        put(dest/'projection_map.json',projection)
        transport={'status':'PASS_pure_transport_decode_projection_roundtrip','pid':os.getpid(),'sid':sid,'records':len(rows),'endpoint_intervals':73*len(rows),
          'compressed_bytes':target.stat().st_size,'shard_sha256':before_sha,'source_before_after_sha':True,'write_projection_equals_read_projection':True,
          'projection_map_sha256':sha(dest/'projection_map.json'),'new_H_R_factor_mean_range_chord_oldpower_calls':0,'mathematical_checker_invoked':False}
        put(dest/'transport.json',transport)
        put(dest/'result.json',{**transport,'plan_sha256':plan_sha,'chunks':len(parts),'seconds':time.monotonic()-start})
    except BaseException as ex:
        put(dest/f'fault_{os.getpid()}.json',{'pid':os.getpid(),'sid':sid,'exception':repr(ex),'complete_chunks_preserved':True,'automatic_retry':False});raise
def stop(procs):
    for p in procs:
        if p.poll() is None:p.terminate()
    for p in procs:
        try:p.wait(timeout=5)
        except subprocess.TimeoutExpired:p.kill();p.wait(timeout=5)
def run(out,plan_sha,approval):
    assert approval=='integral167-parallel15-transport-export-approved';plan=checked(plan_sha)
    for p,h in plan['source_locks'].items():assert sha(p)==h
    out.mkdir(exist_ok=False);start=time.monotonic();procs=[];logs=[]
    put(out/'registration.json',{'task':'167_exact15_transport_export','pid':os.getpid(),'utc':datetime.now(timezone.utc).isoformat(),'plan_sha256':plan_sha,
      'workers':15,'packets':97971,'shard_counts':plan['shard_counts'],'wall_cap':None,'new_H_R_mean_oldpower_calls':0})
    try:
        for sid in range(15):
            log=(out/f'worker{sid:02d}.log').open('xb');logs.append(log)
            p=subprocess.Popen([sys.executable,'-B',__file__,'worker','--output',str(out),'--plan-sha',plan_sha,'--sid',str(sid)],stdout=log,stderr=subprocess.STDOUT);procs.append(p)
        for f in logs:f.close()
        put(out/'pool_registration.json',{'master_pid':os.getpid(),'worker_pids':[p.pid for p in procs],'sids':list(range(15)),'no_automatic_restart':True})
        while any(p.poll() is None for p in procs):
            if any(p.poll() not in (None,0) for p in procs):raise RuntimeError('fixed shard worker failed')
            time.sleep(.1)
        assert all(p.returncode==0 for p in procs)
        results=[read(out/f'worker{sid:02d}/result.json') for sid in range(15)]
        assert [r['sid'] for r in results]==list(range(15)) and sum(r['records'] for r in results)==97971
        for sid,r in enumerate(results):
            assert r['plan_sha256']==plan_sha and r['mathematical_checker_invoked'] is False
            assert sha(out/f'endpoints{sid:02d}.bin.gz')==r['shard_sha256']
        # Master copies only necessary metadata,then assembles the projection index.
        for relative,spec in plan['portable_files'].items():
            assert sha(spec['path'])==spec['sha256'];dest=out/relative;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(spec['path'],dest);assert sha(dest)==spec['sha256']
        shutil.copyfile(H/'plan_v2.json',out/'export_plan_v2.json');assert sha(out/'export_plan_v2.json')==plan_sha
        with (out/'projection_map.jsonl.gz').open('xb') as f:
            with gzip.GzipFile(fileobj=f,mode='wb',mtime=0,filename='') as z:
                for sid in range(15):
                    for r in read(out/f'worker{sid:02d}/projection_map.json'):z.write((json.dumps(r,sort_keys=True,separators=(',',':'))+'\n').encode())
        for p,h in plan['source_locks'].items():assert sha(p)==h
        checked(plan_sha)
        whitelist=['registration.json','pool_registration.json','export_plan_v2.json','projection_map.jsonl.gz',*plan['portable_files'],*[f'endpoints{s:02d}.bin.gz' for s in range(15)],*[f'worker{s:02d}/transport.json' for s in range(15)],*[f'worker{s:02d}/result.json' for s in range(15)]]
        files={n:sha(out/n) for n in whitelist}
        put(out/'completion.json',{'status':'exact15_payload_and_transport_complete_no_math_replay','master_pid':os.getpid(),'worker_pids':[p.pid for p in procs],
          'worker_exitcodes':[p.returncode for p in procs],'all_workers_joined':True,'plan_sha256':plan_sha,'packets':97971,'shards':15,
          'payload_bytes':sum(r['compressed_bytes'] for r in results),'all_shard_transport_checks_pass':True,'files':files,'new_H_R_factor_mean_oldpower_calls':0,'seconds':time.monotonic()-start})
    except BaseException as ex:
        stop(procs)
        for f in logs:
            if not f.closed:f.close()
        put(out/f'fault_{os.getpid()}.json',{'exception':repr(ex),'all_workers_stopped':all(p.poll() is not None for p in procs),'worker_exitcodes':[p.returncode for p in procs],
          'complete_chunks_and_successful_shards_preserved':True,'action':'stop;no automatic restart;primary reviews any source-identical continuation'});raise
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['run','worker']);ap.add_argument('--plan-sha',required=True);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--approval');ap.add_argument('--sid',type=int);ap.add_argument('--resume',action='store_true');a=ap.parse_args()
    if a.mode=='run':run(a.output,a.plan_sha,a.approval)
    else:assert 0<=a.sid<15;worker(a.output,a.sid,a.plan_sha,a.resume)
