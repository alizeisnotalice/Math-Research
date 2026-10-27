#!/usr/bin/env python3
"""Authorized byte-scope continuation: reuse15sealedpairs;newviews only."""
from pathlib import Path
import argparse,datetime,os,time
from io175 import sha,read,put,ref,checked
from view_recipe import build_views
HERE=Path(__file__).resolve().parent

def locks(plan):
    for r in plan['source_files'].values():checked(r)
    for name,h in plan['code_hashes'].items():assert sha(HERE/name)==h
    return plan['source_files']

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--plan-sha',required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    assert sha(HERE/'archive_plan.json')==a.plan_sha;plan=read(HERE/'archive_plan.json');a.output.mkdir(exist_ok=False);start=time.monotonic()
    put(a.output/'registration.json',{'pid':os.getpid(),'task':'175_authorized_byte_only_remaining_views_and_index','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'plan_sha256':a.plan_sha,'code_hashes':plan['code_hashes'],'reused_sealed_pair_files':15,'new_pair_compressions':0,'new_math_calls':0})
    try:
        before=locks(plan);put(a.output/'sources_before.json',before);saved=read(checked(plan['preserved_pairs_index']));pairs=saved['pairs']
        assert len(pairs)==15 and sum(r['count'] for r in pairs)==97971
        for r in pairs:checked(r)
        views=build_views(plan['restored179_root'],plan['repository_root'],a.output/'views002',plan['view_bindings'])
        after=locks(plan);assert before==after;put(a.output/'sources_after.json',after)
        put(a.output/'archive_index.json',{'task':'175_root_pair_transport','pairs':pairs,'total_pairs':97971,'raw_pair_bytes':4212061,
            'compressed_pair_bytes':sum(r['bytes'] for r in pairs),'new_parameter_views':views,'new_math_calls':0,
            'full175output_JSON_byte_restoration_claim':False,'original_prepared_gzip_SHA_claim':False,'new_pair_compressions':0})
        put(a.output/'encode_receipt.json',{'pid':os.getpid(),'status':'authorized_remaining_views_written_reusing_all15_original_pairgzip',
            'seconds':time.monotonic()-start,'archive_index_sha256':sha(a.output/'archive_index.json'),'all_sources_before_after_unchanged':True,
            'sources_before_sha256':sha(a.output/'sources_before.json'),'sources_after_sha256':sha(a.output/'sources_after.json'),
            'new_pair_compressions':0,'new_math_calls':0})
    except BaseException as e:put(a.output/'fault.json',{'pid':os.getpid(),'exception':repr(e),'stop_preserve_no_retry':True});raise
