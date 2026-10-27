#!/usr/bin/env python3
"""179A exact byte transport. No mathematical producer/checker is imported.

Projection means copying saved JSON fields, not computing any numerical value.
"""
from pathlib import Path
import argparse,gzip,hashlib,json,os,tarfile,time

HERE=Path(__file__).resolve().parent
PLAN_SHA='c037e917cf87506e0a7e3e0b0969d9854704ced4e47e36c2de201368ee5e5504'
GROUP_CAP=33554432
FILE_CAP=50000000

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def read(p):
    p=Path(p)
    with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:return json.load(f)

def put(p,x):
    with Path(p).open('x') as f:json.dump(x,f,sort_keys=True,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())

def check_inputs(sel,refs,repo):
    for h,b in sel['blobs'].items():
        p=Path(b['source_path']);assert p.is_file() and not p.is_symlink() and p.stat().st_size==b['original_bytes'] and sha(p)==h
    for rel,x in refs['files'].items():
        p=repo/rel;assert p.resolve().is_relative_to(repo.resolve()) and p.is_file() and p.stat().st_size==x['bytes'] and sha(p)==x['sha256']

def tar_size(group):
    # Only container-byte accounting;this is not an integral/DFS sum.
    n=1024
    for b in group:n+=512+b['bytes']+(-b['bytes'])%512
    return n+(-n)%10240

def main(out,repo,regpath,regsha):
    assert sha(HERE/'plan.json')==PLAN_SHA and sha(regpath)==regsha
    reg=read(regpath);assert reg['plan_sha256']==PLAN_SHA
    assert sha(__file__)==reg['encode_sha256'] and sha(HERE/'decode179.py')==reg['decode_sha256']
    plan=read(HERE/'plan.json')
    for name,key in [('selection_final.json','selection_final_sha256'),('field_recipes.json','field_recipes_sha256'),('repository_refs.json','repository_refs_sha256'),('projection_plan.json','projection_plan_sha256')]:assert sha(HERE/name)==plan[key]
    sel=read(HERE/'selection_final.json');refs=read(HERE/'repository_refs.json');pp=read(HERE/'projection_plan.json')
    assert len(sel['paths'])==1525 and len(refs['files'])==95
    check_inputs(sel,refs,repo);out.mkdir(exist_ok=False);work=out/'staging';work.mkdir();start=time.monotonic()
    put(out/'encode_registration.json',{'pid':os.getpid(),'plan_sha256':PLAN_SHA,'preregistration_sha256':regsha,'encode_sha256':sha(__file__),
        'decode_sha256':reg['decode_sha256'],'profile':'A','literal_blobs':302,'projection_payloads':32,'virtual_math_targets':1173,'new_math_calls':0})
    physical={};projections={};progress=[]
    for h,b in sel['blobs'].items():
        if b['storage']=='literal':physical[h]={'sha256':h,'bytes':b['original_bytes'],'path':b['source_path'],'kind':'original_literal'}
    assert len(physical)==302
    def source(logical):
        h=sel['paths'][logical]['blob_sha256'];return Path(sel['blobs'][h]['source_path'])
    def project(name,rows):
        p=work/(name+'.jsonl');count=0
        with p.open('xb') as f:
            for row in rows:
                f.write(json.dumps(row,separators=(',',':'),ensure_ascii=True).encode('utf-8')+b'\n');count+=1
            f.flush();os.fsync(f.fileno())
        assert p.stat().st_size==pp['projected_raw_bytes_by_payload'][name]
        h=sha(p);assert h not in physical
        physical[h]={'sha256':h,'bytes':p.stat().st_size,'path':str(p),'kind':'saved_field_projection'}
        projections[name]={'blob_sha256':h,'bytes':p.stat().st_size,'records':count}
        progress.append({'projection':name,**projections[name]});put(out/f'projection_checkpoint{len(progress):02d}.json',progress[-1])
    def original_rows(case,worker):
        selected=[]
        for h,b in sel['blobs'].items():
            r=b.get('recipe',{})
            if r.get('kind')=='producer_math_rows' and r['case']==case and r['worker']==worker:selected.append((r['chunk'],b['source_path']))
        assert len(selected)==26
        for ch,p in sorted(selected):yield from read(p)
    cols164=pp['schema']['producer164_math_rows']['extra_columns']
    for i in range(15):
        project(f'extra162_gain_worker{i:02d}',([r['key'],r['total_gain_relative112_units']] for r in original_rows('162',i)))
        project(f'extra164_witness_worker{i:02d}',([r[c] for c in cols164] for r in original_rows('164',i)))
    eff=read(source('producer162/actual/effective_changes.json'))
    project('old156_effective_rows',([r['key'],r['old156_effU']] for r in eff))
    project('global_effective_order',([r['key']] for r in eff))
    assert len(projections)==32 and len(physical)==334
    control=out/'control';control.mkdir();controls={}
    for name in ['plan.json','selection_final.json','field_recipes.json','repository_refs.json','projection_plan.json']:
        source=HERE/name;dest=control/name
        with source.open('rb') as f,dest.open('xb') as g:
            for b in iter(lambda:f.read(1048576),b''):g.write(b)
        assert sha(source)==sha(dest);controls['control/'+name]={'bytes':dest.stat().st_size,'sha256':sha(dest)}
    for name,source in [('preregistration.json',Path(regpath))]+[(n,HERE/n) for n in ['encode179.py','decode179.py','run179.py']]:
        dest=control/name
        with source.open('rb') as f,dest.open('xb') as g:
            for b in iter(lambda:f.read(1048576),b''):g.write(b)
        assert sha(source)==sha(dest);controls['control/'+name]={'bytes':dest.stat().st_size,'sha256':sha(dest)}
    # Container groups are fixed by content and byte sizes,never by numerical values.
    groups=[]
    for b in sorted(physical.values(),key=lambda b:(-b['bytes'],b['sha256'])):
        assert tar_size([b])<=GROUP_CAP
        for g in groups:
            if tar_size(g+[b])<=GROUP_CAP:g.append(b);break
        else:groups.append([b])
    packs=[]
    for i,g in enumerate(groups):
        plain=work/f'payload{i:02d}.tar'
        with tarfile.open(plain,'w',format=tarfile.USTAR_FORMAT) as t:
            for b in sorted(g,key=lambda b:b['sha256']):
                info=tarfile.TarInfo('blobs/'+b['sha256']);info.size=b['bytes'];info.mode=0o444;info.uid=info.gid=0;info.mtime=0;info.uname=info.gname=''
                with Path(b['path']).open('rb') as f:t.addfile(info,f)
        assert plain.stat().st_size==tar_size(g)<=GROUP_CAP
        pack=out/f'payload{i:02d}.tar.gz'
        with pack.open('xb') as dest:
            with gzip.GzipFile(fileobj=dest,mode='wb',filename='',mtime=0,compresslevel=9) as z:
                with plain.open('rb') as f:
                    for b in iter(lambda:f.read(1048576),b''):z.write(b)
            dest.flush();os.fsync(dest.fileno())
        assert pack.stat().st_size<FILE_CAP
        x={'name':pack.name,'bytes':pack.stat().st_size,'sha256':sha(pack),'uncompressed_tar_bytes':plain.stat().st_size,
            'members':[b['sha256'] for b in sorted(g,key=lambda b:b['sha256'])]};packs.append(x);put(out/f'pack_checkpoint{i:02d}.json',x)
    # Strip original scratch source paths from the portable control manifest.
    targets={}
    for h,b in sel['blobs'].items():
        targets[h]={'original_bytes':b['original_bytes'],'storage':b['storage']}
        if b['storage']=='projection_recipe':targets[h]['recipe']=b['recipe']
    manifest={'format':'ARC179A1','profile':'A','plan_sha256':PLAN_SHA,'preregistration_sha256':regsha,
        'encode_sha256':reg['encode_sha256'],'decode_sha256':reg['decode_sha256'],'target_paths':sel['paths'],'target_blobs':targets,
        'physical_blobs':{h:{'bytes':b['bytes'],'kind':b['kind']} for h,b in physical.items()},'projections':projections,'extra164_columns':cols164,
        'repository_refs':refs['files'],'control_files':controls,'packs':packs,'group_cap_bytes':GROUP_CAP,'physical_file_hard_gate_bytes':FILE_CAP,
        'restorable_original_math_JSON_files':1173,'prepared_cache_gzip_not_claimed_byte_restored':True,'new_math_calls':0}
    put(out/'archive_manifest.json',manifest)
    check_inputs(sel,refs,repo)
    assert sha(__file__)==reg['encode_sha256'] and sha(HERE/'decode179.py')==reg['decode_sha256'] and sha(HERE/'plan.json')==PLAN_SHA
    receipt={'status':'encoded_once_pending_independent_byte_decode','pid':os.getpid(),'seconds':time.monotonic()-start,'profile':'A',
        'plan_sha256':PLAN_SHA,'preregistration_sha256':regsha,'archive_manifest_sha256':sha(out/'archive_manifest.json'),'literal_blobs':302,'projection_payloads':32,
        'physical_blobs':len(physical),'packs':packs,'packed_bytes':sum(p['bytes'] for p in packs),'all_sources_before_after_unchanged':True,
        'H_binary_decodes':0,'mathematical_algorithm_calls':0,'fallback_or_extra_raw_sources':False}
    put(out/'encode_receipt.json',receipt);print(json.dumps(receipt,sort_keys=True))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);ap.add_argument('--repository',type=Path,required=True)
    ap.add_argument('--preregistration',type=Path,required=True);ap.add_argument('--preregistration-sha',required=True);a=ap.parse_args()
    try:main(a.output,a.repository,a.preregistration,a.preregistration_sha)
    except BaseException as ex:
        if a.output.exists() and not (a.output/'encode_fault.json').exists():put(a.output/'encode_fault.json',{'pid':os.getpid(),'exception':repr(ex),'prefix_preserved':True,'action':'stop,no automatic retry,profileB or extra raw source'})
        raise
