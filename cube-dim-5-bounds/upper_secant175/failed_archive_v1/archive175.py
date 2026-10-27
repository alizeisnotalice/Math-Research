#!/usr/bin/env python3
"""One authorized root-pair byte projection and new parameter-view generation."""
from pathlib import Path
import argparse,datetime,json,os,time
from io175 import sha,read,put,ref,checked,gzwrite
from view_recipe import build_views
HERE=Path(__file__).resolve().parent

def locks(plan):
    out={}
    for k,r in plan['source_files'].items():checked(r);out[k]=r
    for n,h in plan['code_hashes'].items():assert sha(HERE/n)==h
    return out

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--plan-sha',required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    assert sha(HERE/'archive_plan.json')==a.plan_sha;plan=read(HERE/'archive_plan.json');a.output.mkdir(exist_ok=False);start=time.monotonic()
    put(a.output/'registration.json',{'pid':os.getpid(),'task':'175_rootpair_byte_projection_only','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'plan_sha256':a.plan_sha,'code_hashes':plan['code_hashes'],'root_pairs':97971,'new_math_calls':0,'automatic_retry':False})
    try:
        before=locks(plan);put(a.output/'sources_before.json',before);dest=a.output/'pairs';dest.mkdir()
        root=Path(plan['original175_output_root']);pairs=[];total=0;rawtotal=0
        for i,names in enumerate(plan['175_chunk_order']):
            rows=[]
            for name in names:
                for r in read(root/name):
                    assert type(r['Y']) is int and type(r['Z']) is int;rows.append([r['Y'],r['Z']])
            raw=b''.join(json.dumps(r,separators=(',',':')).encode()+b'\n' for r in rows)
            p=dest/f'roots{i:02d}.jsonl.gz';gzwrite(p,raw);pairs.append({'worker':i,'count':len(rows),'uncompressed_bytes':len(raw),**ref(p)});total+=len(rows);rawtotal+=len(raw)
        assert total==97971 and rawtotal==4212061
        views=build_views(plan['restored179_root'],plan['repository_root'],a.output/'parameter_views',plan['view_bindings'])
        after=locks(plan);assert before==after;put(a.output/'sources_after.json',after)
        put(a.output/'archive_index.json',{'task':'175_root_pair_transport','pairs':pairs,'total_pairs':total,'raw_pair_bytes':rawtotal,
            'compressed_pair_bytes':sum(r['bytes'] for r in pairs),'new_parameter_views':views,'new_math_calls':0,
            'full175output_JSON_byte_restoration_claim':False,'original_prepared_gzip_SHA_claim':False})
        put(a.output/'encode_receipt.json',{'pid':os.getpid(),'status':'byte_projection_and_new_views_written_pending_independent_byte_audit',
            'seconds':time.monotonic()-start,'archive_index_sha256':sha(a.output/'archive_index.json'),'all_sources_before_after_unchanged':True,
            'sources_before_sha256':sha(a.output/'sources_before.json'),'sources_after_sha256':sha(a.output/'sources_after.json'),'new_math_calls':0})
    except BaseException as e:put(a.output/'fault.json',{'pid':os.getpid(),'exception':repr(e),'stop_preserve_no_retry':True});raise
