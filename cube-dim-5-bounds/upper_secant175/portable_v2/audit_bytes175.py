#!/usr/bin/env python3
"""Independent pair decoder + full field comparison; never imports writer/view recipe."""
from pathlib import Path
import argparse,gzip,hashlib,itertools,json,os,time
from io175 import sha,read,put,checked,lines,canonical
HERE=Path(__file__).resolve().parent

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--archive',type=Path,required=True);ap.add_argument('--index-sha',required=True)
    ap.add_argument('--plan-sha',required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    assert sha(HERE/'archive_plan.json')==a.plan_sha;plan=read(HERE/'archive_plan.json')
    assert sha(a.archive/'archive_index.json')==a.index_sha;index=read(a.archive/'archive_index.json')
    for name,h in plan['code_hashes'].items():assert sha(HERE/name)==h
    for r in plan['source_files'].values():checked(r)
    a.output.mkdir(exist_ok=False);start=time.monotonic();total=0;fieldtotal=0;checks=[]
    put(a.output/'registration.json',{'pid':os.getpid(),'code_sha256':sha(__file__),'plan_sha256':a.plan_sha,'archive_index_sha256':a.index_sha,'new_math_calls':0})
    try:
        root=Path(plan['original175_output_root']);pre=read(plan['actual175_preflight']['path'])
        for i in range(15):
            pr=index['pairs'][i];checked(pr)
            def originals():
                for name in plan['175_chunk_order'][i]:yield from read(root/name)
            pairs=lines(pr['path']);n=0
            for pair,r in itertools.zip_longest(pairs,originals()):
                assert pair is not None and r is not None and type(pair) is list and len(pair)==2
                assert type(pair[0]) is int and type(pair[1]) is int and pair==[r['Y'],r['Z']];n+=1
            assert n==pr['count']==pre['shards'][i]['count'];total+=n
            new=index['new_parameter_views']['shards'][i];checked(new);old=pre['shards'][i];checked(old)
            count=0;oldhash=hashlib.sha256();newhash=hashlib.sha256();fields=set()
            for s,t in itertools.zip_longest(lines(new['path']),lines(old['path'])):
                assert s is not None and t is not None;sb=canonical(s);tb=canonical(t);assert sb==tb
                assert s['key']==t['key'];oldhash.update(tb+b'\n');newhash.update(sb+b'\n');fields.update(s);fieldtotal+=len(s);count+=1
            assert count==n and oldhash.digest()==newhash.digest()
            checks.append({'worker':i,'pairs':n,'all_originalY_Z_fields_equal':True,'all_parameter_fields_equal':True,'field_names':sorted(fields),
                'canonical_parameters_sha256':oldhash.hexdigest(),'old_prepared_gzip_sha256':old['sha256'],'new_view_gzip_sha256':new['sha256']})
        assert total==97971
        tr=index['new_parameter_views']['topology'];checked(tr);old=pre['topology'];checked(old)
        nt,ot=read(tr['path']),read(old['path']);assert canonical(nt)==canonical(ot)
        assert len(nt['nodes'])==97972 and len(nt['roots'])==52360
        for r in plan['source_files'].values():checked(r)
        put(a.output/'checks.json',{'shards':checks,'root_pairs':total,'parameter_field_comparisons':fieldtotal,
            'all_topology_fields_order_and_roots_equal':True,'canonical_topology_sha256':hashlib.sha256(canonical(nt)).hexdigest(),
            'old_topology_gzip_sha256':old['sha256'],'new_topology_gzip_sha256':tr['sha256'],
            'new_math_calls':0,'full175output_JSON_byte_restoration_claim':False,'old_prepared_gzip_SHA_claim':False})
        put(a.output/'receipt.json',{'pid':os.getpid(),'status':'all97971pairs_and_all15parameter_shards_and_fulltopology_field_equivalence_pass',
            'seconds':time.monotonic()-start,'archive_index_sha256':a.index_sha,'checks_sha256':sha(a.output/'checks.json'),
            'all_sources_unchanged':True,'new_math_calls':0,'H_decodes':0,'portable_mathematical_verifier_executed':False})
    except BaseException as e:put(a.output/'fault.json',{'pid':os.getpid(),'exception':repr(e),'prefix_preserved':True,'automatic_retry':False});raise
