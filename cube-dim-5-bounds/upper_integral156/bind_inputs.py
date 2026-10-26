#!/usr/bin/env python3
"""Semantic binding of serialized shards/topology to frozen original fields.

Only equality/eligibility/source metadata; no new factors/products/means/gains.
"""
from pathlib import Path
import argparse,json,os,time
from common import sha,read,put,lines,lockcheck

HERE=Path(__file__).resolve().parent

def main(prepath,output):
    src=read(HERE/'sources.json')['source_locks'];lockcheck(src);pre=read(prepath)
    put(output.with_suffix('.registration.json'),{'pid':os.getpid(),'code_sha256':sha(__file__),'common_sha256':sha(HERE/'common.py'),
        'preflight_results':{'path':str(prepath),'sha256':sha(prepath)},'source_original_nodes':97972,'source_shard_rows':97907,
        'new_factor_product_mean_gain_UH_evaluations':0})
    start=time.monotonic()
    for x in pre['external_files'].values():assert sha(x['path'])==x['sha256']
    reg=read(pre['external_files']['registration.json']['path']);ini=read(pre['external_files']['initial_registration.json']['path'])
    assert ini['prepare_sha256']==sha(HERE/'prepare.py') and ini['common_sha256']==sha(HERE/'common.py')
    assert ini['sources_sha256']==sha(HERE/'sources.json')
    rows={r['key']:r for r in read(src['112_nodes']['path'])};certs={c['key']:c for c in read(src['112_certificates']['path'])}
    original=read(src['112_tree']['path']);nodes=original['nodes'];tree=read(reg['topology']['path'])
    assert len(rows)==len(certs)==len(nodes)==len(tree['nodes'])==97972
    assert tree['roots']==original['roots'] and tree['old_complete_UH']==original['budget']['UH']
    for n in tree['nodes']:
        old=nodes[n['key']]
        assert n=={**{k:old[k] for k in ['key','level','parent','selfL','selfU','effU']},'children':old['children'] or []}
        assert old['weight']==rows[n['key']]['weight']
    reuse=read(reg['reuse150']['path']);assert reuse==read(src['150_nodes']['path']);reused={r['key'] for r in reuse}
    assert len(reused)==64
    seen=set();counts=[]
    for idx in reg['shards']:
        assert sha(idx['path'])==idx['sha256'];count=0
        for item in lines(idx['path']):
            k=item['key'];assert k not in seen and k not in reused;seen.add(k);count+=1
            r,c,n=rows[k],certs[k],nodes[k]
            assert n['selfU']==c['new_selfU']<c['old_selfU'] and r['old_U_i']>=1<<144 and r['old_L_i']<r['old_U_i']
            expected={k:r[k] for k in ['key','geometry','level','weight','point','old_L_i','old_U_i','packet_path','packet_sha256']}
            group='range092' if '/upper_range092/' in r['packet_path'] else 'adaptive095' if '/upper_adaptive095/' in r['packet_path'] else 'adaptive103'
            expected.update({'old_mean_i':c['new_M_i'],'old_selfL':n['selfL'],'old_selfU':n['selfU'],'source_group':group})
            assert item==expected
        assert count==idx['count'];counts.append(count)
    expected={k for k,r in rows.items() if nodes[k]['selfU']==certs[k]['new_selfU']<certs[k]['old_selfU'] and r['old_U_i']>=1<<144 and r['old_L_i']<r['old_U_i']}-reused
    assert seen==expected and len(seen)==97907
    assert set(nodes)-seen-reused=={'5:31,31,31,31'}
    assert read(reg['unchanged']['path'])==[{'key':'5:31,31,31,31','reasons':['112strict_chord_not_active']}]
    lockcheck(src)
    result={'task':'upper_integral156_semantic_source_binding','pid':os.getpid(),'seconds':time.monotonic()-start,
        'preflight_results_sha256':sha(prepath),'sources_sha256':sha(HERE/'sources.json'),'source_original_nodes':97972,'source_shard_rows':97907,
        'shard_counts':counts,'copied150_roots':64,'unchanged_self_nodes':1,'all_shards_and_topology_semantically_bound':True,
        'all_frozen_sources_unchanged':True,'new_factor_product_mean_gain_UH_evaluations':0,'global_Cn_claim':False}
    put(output,result);print(json.dumps(result,sort_keys=True))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--preflight',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();main(a.preflight,a.output)
