#!/usr/bin/env python3
"""Prepared fixed15-shard exact serialization. Not authorized/run yet.

Each <=256-record chunk is transactional; completed chunks are verified/copied
without re-reading source packets on explicit resume. No arithmetic evaluator.
"""
from pathlib import Path
import argparse,gzip,hashlib,json,os,shutil
import codec as C
H=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):
    p=Path(p)
    with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:return json.load(f)
def put(p,x):
    with Path(p).open('x') as f:json.dump(x,f,sort_keys=True,separators=(',',':'));f.write('\n');f.flush();os.fsync(f.fileno())
def lines(p):
    with gzip.open(p,'rt') as f:
        for line in f:yield json.loads(line)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plan-sha',required=True);ap.add_argument('--approval',required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--resume',action='store_true');a=ap.parse_args()
    assert a.approval=='integral167-exact-export-approved' and sha(H/'export_plan.json')==a.plan_sha
    plan=read(H/'export_plan.json')
    for n,h in plan['code_locks'].items():assert sha(H/n)==h
    for p,h in plan['source_locks'].items():assert sha(p)==h
    assert sha(H/'node_index.jsonl.gz')==plan['index_sha256']
    if a.resume:assert a.output.is_dir() and not (a.output/'completion.json').exists()
    else:a.output.mkdir(exist_ok=False);put(a.output/'registration.json',{'pid':os.getpid(),'plan_sha256':a.plan_sha,'packets':97971,'shards':15,'new_H_R_oldpower_calls':0})
    assert read(a.output/'registration.json')['plan_sha256']==a.plan_sha
    chunks=a.output/'scratch_chunks';chunks.mkdir(exist_ok=True);groups=[[] for _ in range(15)]
    for r in lines(H/'node_index.jsonl.gz'):groups[r['shard']].append(r)
    try:
        for sid,rows in enumerate(groups):
            assert len(rows)==plan['shard_counts'][sid];parts=[]
            for j in range(0,len(rows),256):
                p=chunks/f'{sid:02d}_{j//256:04d}.gz';receipt=p.with_suffix('.json');batch=rows[j:j+256]
                if receipt.exists():
                    x=read(receipt);assert x['first_id']==batch[0]['id'] and x['count']==len(batch) and sha(p)==x['sha256']
                else:
                    assert not p.exists();tmp=p.with_suffix('.tmp');projection_rows=[]
                    with tmp.open('xb') as dest:
                        with gzip.GzipFile(fileobj=dest,mode='wb',mtime=0,filename='') as z:
                            if j==0:z.write(C.MAGIC+sid.to_bytes(4,'big')+len(rows).to_bytes(4,'big'))
                            for row in batch:
                                raw=Path(row['packet_path']);assert sha(raw)==row['packet_sha256'];data,proj,original=C.encode(row,read(raw),plan['times']);z.write(data)
                                projection_rows.append({'id':row['id'],'source_gzip_sha256':row['packet_sha256'],'original_raw_sha256':original,'normalized_projection_sha256':proj})
                        dest.flush();os.fsync(dest.fileno())
                    tmp.rename(p);put(receipt,{'first_id':batch[0]['id'],'count':len(batch),'sha256':sha(p),'projections':projection_rows})
                parts.append(p)
            target=a.output/f'endpoints{sid:02d}.bin.gz'
            if not target.exists():
                with target.with_suffix('.tmp').open('xb') as f:
                    for p in parts:
                        with p.open('rb') as s:shutil.copyfileobj(s,f)
                    f.flush();os.fsync(f.fileno())
                target.with_suffix('.tmp').rename(target)
            assert target.stat().st_size<plan['shard_hard_limit_bytes']
        for p,h in plan['source_locks'].items():assert sha(p)==h
        # Portable metadata and proof sources are byte copies,not recomputations.
        for relative,spec in plan['portable_files'].items():
            src=Path(spec['path']);assert sha(src)==spec['sha256'];dest=a.output/relative;dest.parent.mkdir(parents=True,exist_ok=True)
            if not dest.exists():shutil.copyfile(src,dest)
            assert sha(dest)==spec['sha256']
        if not (a.output/'export_plan.json').exists():shutil.copyfile(H/'export_plan.json',a.output/'export_plan.json')
        assert sha(a.output/'export_plan.json')==a.plan_sha
        if not (a.output/'projection_map.jsonl.gz').exists():
            with (a.output/'projection_map.jsonl.gz.tmp').open('xb') as f:
                with gzip.GzipFile(filename='',mode='wb',mtime=0,fileobj=f) as z:
                    for p in sorted(chunks.glob('*.json')):
                        for row in read(p)['projections']:z.write((json.dumps(row,sort_keys=True,separators=(',',':'))+'\n').encode())
            (a.output/'projection_map.jsonl.gz.tmp').rename(a.output/'projection_map.jsonl.gz')
        files={str(p.relative_to(a.output)):sha(p) for p in a.output.rglob('*') if p.is_file() and chunks not in p.parents}
        put(a.output/'completion.json',{'status':'exact_payload_exported_no_math_replay','pid':os.getpid(),'plan_sha256':a.plan_sha,'files':files,'new_H_R_factor_mean_oldpower_calls':0})
    except BaseException as ex:
        fault=a.output/f'fault_{os.getpid()}.json';put(fault,{'error':repr(ex),'completed_chunks_preserved':True,'action':'stop;explicit source-identical resume only'});raise
if __name__=='__main__':main()
