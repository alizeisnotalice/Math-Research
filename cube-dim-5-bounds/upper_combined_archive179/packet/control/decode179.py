#!/usr/bin/env python3
"""Independent179A byte decoder;input:new packs and tracked repository only.

Never opens original producer scratch files or imports numerical code.
"""
from pathlib import Path,PurePosixPath
import argparse,gzip,hashlib,json,os,shutil,tarfile,time

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def read(p):
    p=Path(p)
    with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:return json.load(f)
def lines(p):
    with Path(p).open() as f:
        for s in f:yield json.loads(s)
def gzlines(p):
    with gzip.open(p,'rt') as f:
        for s in f:yield json.loads(s)
def put(p,x):
    with Path(p).open('x') as f:json.dump(x,f,sort_keys=True,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
def safe(rel):
    p=PurePosixPath(rel);assert not p.is_absolute() and '..' not in p.parts and str(p)==rel
    return p

def main(archive,msha,repo,out):
    mf=archive/'archive_manifest.json';assert sha(mf)==msha;m=read(mf);assert m['format']=='ARC179A1' and m['profile']=='A'
    assert m['plan_sha256']=='c037e917cf87506e0a7e3e0b0969d9854704ced4e47e36c2de201368ee5e5504'
    for rel,x in m['control_files'].items():
        f=archive/str(safe(rel));assert f.stat().st_size==x['bytes'] and sha(f)==x['sha256']
    control=archive/'control';assert sha(control/'plan.json')==m['plan_sha256'];plan=read(control/'plan.json')
    assert sha(control/'preregistration.json')==m['preregistration_sha256']
    preg=read(control/'preregistration.json');assert preg['encode_sha256']==m['encode_sha256'] and preg['decode_sha256']==m['decode_sha256']
    assert sha(control/'encode179.py')==m['encode_sha256'] and sha(control/'decode179.py')==m['decode_sha256'] and sha(control/'run179.py')==preg['controller_sha256']
    for name,key in [('selection_final.json','selection_final_sha256'),('field_recipes.json','field_recipes_sha256'),('repository_refs.json','repository_refs_sha256'),('projection_plan.json','projection_plan_sha256')]:assert sha(control/name)==plan[key]
    selected=read(control/'selection_final.json');assert m['target_paths']==selected['paths']
    expected={h:{'original_bytes':x['original_bytes'],'storage':x['storage'],**({'recipe':x['recipe']} if x['storage']=='projection_recipe' else {})} for h,x in selected['blobs'].items()}
    assert m['target_blobs']==expected and m['repository_refs']==read(control/'repository_refs.json')['files']
    assert m['extra164_columns']==read(control/'field_recipes.json')['extra164_projection_columns']
    expected_projections=read(control/'projection_plan.json')['projected_raw_bytes_by_payload']
    assert set(m['projections'])==set(expected_projections) and len(m['projections'])==32
    assert all(m['projections'][k]['bytes']==v for k,v in expected_projections.items())
    assert sha(__file__)==m['decode_sha256'];out.mkdir(exist_ok=False);cache=out/'blob_cache';cache.mkdir();work=out/'decode_work';work.mkdir();restored=out/'restored';restored.mkdir();start=time.monotonic()
    put(out/'decode_registration.json',{'pid':os.getpid(),'decode_sha256':sha(__file__),'archive_manifest_sha256':msha,'source_rule':'newpacks+trackedrepository only;no original producer scratch','new_math_calls':0})
    # Repository bytes are hashed;thermal endpoint binary is never opened as gzip.
    for rel,x in m['repository_refs'].items():
        p=repo/str(safe(rel));assert p.resolve().is_relative_to(repo.resolve()) and p.is_file() and p.stat().st_size==x['bytes'] and sha(p)==x['sha256']
    assert len(m['repository_refs'])==95
    seen=set();pack_checks=[]
    for i,p in enumerate(m['packs']):
        assert str(safe(p['name']))==Path(p['name']).name
        packed=archive/p['name'];assert packed.stat().st_size==p['bytes']<50000000 and sha(packed)==p['sha256']
        plain=work/f'payload{i:02d}.tar';n=0
        with gzip.open(packed,'rb') as f,plain.open('xb') as dest:
            for b in iter(lambda:f.read(1048576),b''):
                n+=len(b);assert n<=33554432;dest.write(b)
            dest.flush();os.fsync(dest.fileno())
        assert n==p['uncompressed_tar_bytes'];members=[]
        with tarfile.open(plain,'r:') as t:
            for member in t:
                assert member.isfile() and not member.issym() and not member.islnk()
                parts=safe(member.name).parts;assert len(parts)==2 and parts[0]=='blobs';h=parts[1]
                assert len(h)==64 and all(c in '0123456789abcdef' for c in h) and h in m['physical_blobs'] and h not in seen
                assert member.size==m['physical_blobs'][h]['bytes'];target=cache/h
                f=t.extractfile(member);assert f is not None
                with target.open('xb') as dest:shutil.copyfileobj(f,dest,1048576)
                assert target.stat().st_size==member.size and sha(target)==h;seen.add(h);members.append(h)
        assert members==p['members'];pack_checks.append({'name':p['name'],'sha256':p['sha256'],'uncompressed_bytes':n,'verified_members':len(members)})
    assert seen==set(m['physical_blobs']) and len(seen)==334
    # All ordinary aliases use the same verified content blob;no mathematical file is read from scratch.
    for rel,x in m['target_paths'].items():
        safe(rel);h=x['blob_sha256']
        if m['target_blobs'][h]['storage']=='literal':
            dest=restored/rel;dest.parent.mkdir(parents=True,exist_ok=True)
            with (cache/h).open('rb') as f,dest.open('xb') as g:shutil.copyfileobj(f,g,1048576)
    def logical(rel):return cache/m['target_paths'][rel]['blob_sha256']
    def projection(name):return cache/m['projections'][name]['blob_sha256']
    math_checks={}
    def restore_math(h,rows):
        assert m['target_blobs'][h]['storage']=='projection_recipe';dest=cache/h
        with dest.open('x') as f:json.dump(rows,f,sort_keys=True,separators=(',',':'),ensure_ascii=True);f.write('\n');f.flush();os.fsync(f.fileno())
        assert dest.stat().st_size==m['target_blobs'][h]['original_bytes'] and sha(dest)==h
        math_checks[h]={'sha256':h,'bytes':dest.stat().st_size,'original_bytes_identical':True}
    base={r['key']:r for r in gzlines(repo/'cube-dim-5-bounds/upper_integral167/packet/node_index.jsonl.gz')};assert len(base)==97971
    combined_fields=['key','source112_mean_i','cached162_mean_i','cached164_slope_numerator','cached164_slope_denominator','source112_selfU',
        'candidate162_selfU','candidate164_selfU','combined_total_gain_relative112_units','new_selfU']
    for i in range(15):
        extra162=iter(lines(projection(f'extra162_gain_worker{i:02d}')));extra164=iter(lines(projection(f'extra164_witness_worker{i:02d}')))
        recipes={};buffers={};positions={};counts={}
        for case in ['162','164','166']:
            entries=[]
            for h,x in m['target_blobs'].items():
                r=x.get('recipe',{})
                if r.get('kind')=='producer_math_rows' and r['case']==case and r['worker']==i:entries.append((r['chunk'],h))
            recipes[case]=[h for ch,h in sorted(entries)];assert len(entries)==26
            # Frozen result metadata preserves exact original chunk boundaries.
            wr=read(logical(f'producer{case}/actual/worker{i:02d}/result.json'));counts[case]=[c['count'] for c in wr['chunks']]
            assert [c['sha256'] for c in wr['chunks']]==recipes[case]
            buffers[case]=[];positions[case]=0
        count=0
        for g in gzlines(logical(f'shared169/worker{i:02d}/gain_checks.jsonl.gz')):
            e162,e164=next(extra162),next(extra164);k=g['key'];assert e162[0]==e164[0]==k;b=base[k]
            r162={'key':k,'source112_mean_i':g['source112_mean_i'],'retained156_mean_i':b['new_mean_i'],'new_mean_i':g['cached162_mean_i'],
                'source112_selfU':g['source112_selfU'],'retained156_selfU':b['new_selfU'],'new_selfU':g['candidate162_selfU'],'total_gain_relative112_units':e162[1],'packet_sha256':b['packet_sha256']}
            r164={'key':k,'source112_mean_i':g['source112_mean_i'],'retained156_mean_i':b['new_mean_i'],'source112_selfU':g['source112_selfU'],
                'retained156_selfU':b['new_selfU'],'new_selfU':g['candidate164_selfU'],'slope_lower_numerator':g['cached164_slope_numerator'],'slope_lower_denominator':g['cached164_slope_denominator']}
            r164.update(dict(zip(m['extra164_columns'][1:],e164[1:])))
            rows={'162':r162,'164':r164,'166':{f:g[f] for f in combined_fields}}
            for case,row in rows.items():
                buffers[case].append(row);pos=positions[case]
                if len(buffers[case])==counts[case][pos]:restore_math(recipes[case][pos],buffers[case]);buffers[case]=[];positions[case]+=1
            count+=1
        assert next(extra162,None) is None and next(extra164,None) is None
        assert all(positions[c]==26 and not buffers[c] for c in positions)
        put(out/f'worker_byte_restore{i:02d}.json',{'worker':i,'rows':count,'original_JSON_chunks':78,'all_original_SHA_length_checks':True,'new_math_calls':0})
    del base
    old156={k:v for k,v in lines(projection('old156_effective_rows'))};order=[x[0] for x in lines(projection('global_effective_order'))]
    tree=read(logical('shared169/tree_checks.json'))['effective_nodes'];assert len(old156)==len(order)==len(tree)==97972 and len(set(order))==97972
    for case in ['162','164','166']:
        rows=[]
        for k in order:
            t=tree[k]
            if case=='166':rows.append(dict(t))
            else:rows.append({'key':k,'old156_effU':old156[k],'new_effU':t['candidate'+case+'_effU']})
        targets=[h for h,x in m['target_blobs'].items() if x.get('recipe',{}).get('kind')=='producer_effective_rows' and x['recipe']['case']==case]
        assert len(targets)==1;restore_math(targets[0],rows)
    assert len(math_checks)==1173
    # Materialize every logical alias and independently hash every restored path.
    checks={}
    for rel,x in m['target_paths'].items():
        h=x['blob_sha256'];dest=restored/str(safe(rel))
        if not dest.exists():
            dest.parent.mkdir(parents=True,exist_ok=True)
            with (cache/h).open('rb') as f,dest.open('xb') as g:shutil.copyfileobj(f,g,1048576)
        assert dest.stat().st_size==x['original_bytes'] and sha(dest)==h
        checks[rel]={'sha256':h,'bytes':dest.stat().st_size,'original_bytes_identical':True}
    assert len(checks)==1525
    put(out/'restoration_checks.json',{'archive_manifest_sha256':msha,'paths':checks,'original_math_blobs':math_checks,'new_math_calls':0})
    assert sha(__file__)==m['decode_sha256'] and sha(mf)==msha
    for rel,x in m['repository_refs'].items():
        f=repo/str(safe(rel));assert f.resolve().is_relative_to(repo.resolve()) and f.stat().st_size==x['bytes'] and sha(f)==x['sha256']
    receipt={'status':'independent_byte_decode_all_original_targets_pass','pid':os.getpid(),'seconds':time.monotonic()-start,'archive_manifest_sha256':msha,
        'decode_sha256':sha(__file__),'physical_blobs':len(seen),'projections':len(m['projections']),'original_mathematical_JSON':len(math_checks),'logical_paths':len(checks),
        'packs':pack_checks,'repository_references':95,'all_original_lengths_SHA_match':True,'restoration_checks_sha256':sha(out/'restoration_checks.json'),
        'original_scratch_files_opened':0,'H_binary_decodes':0,'mathematical_algorithm_calls':0,'mathematical_qualification_claim':False}
    put(out/'decode_receipt.json',receipt);print(json.dumps(receipt,sort_keys=True))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--archive',type=Path,required=True);ap.add_argument('--manifest-sha',required=True)
    ap.add_argument('--repository',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    try:main(a.archive,a.manifest_sha,a.repository,a.output)
    except BaseException as ex:
        if a.output.exists() and not (a.output/'decode_fault.json').exists():put(a.output/'decode_fault.json',{'pid':os.getpid(),'exception':repr(ex),'prefix_preserved':True,'action':'stop,no automatic retry,fallback or original scratch lookup'})
        raise
