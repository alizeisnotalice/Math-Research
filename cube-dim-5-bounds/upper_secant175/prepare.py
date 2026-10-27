#!/usr/bin/env python3
"""175 metadata copy/binding only: no root, slope, gain or tree arithmetic."""
from pathlib import Path
import argparse,gzip,json,os,time
from common import sha,read,put,lines,inputcheck
HERE=Path(__file__).resolve().parent

def prepare(out):
    src=read(HERE/'sources.json');inputcheck(src);out.mkdir(exist_ok=False);start=time.monotonic()
    put(out/'registration.json',{'task':'175_metadata_only','pid':os.getpid(),'sources_sha256':sha(HERE/'sources.json'),
        'prepare_sha256':sha(__file__),'common_sha256':sha(HERE/'common.py'),'new_root_slope_gain_DFS_cases':0})
    locks=src['source_locks'];pre=read(locks['metadata166']['path']);oldtree=read(pre['topology']['path'])
    root=Path(locks['completion166']['path']).parent;done=read(root/'completion.json');accepted=read(locks['accepted166']['path'])
    assert accepted['status']=='primary_accepts_complete_I166_and_C5_upper' and accepted['new_UH']==done['new_UH']
    assert accepted['sources']['producer166']['sha256']==locks['completion166']['sha256']
    assert accepted['sources']['audit169']['sha256']==locks['audit169']['sha256']
    eff=read(root/'effective_nodes.json');emap={r['key']:r for r in eff};assert len(emap)==97972
    nodes={n['key']:n for n in oldtree['nodes']};assert len(nodes)==97972 and len(oldtree['roots'])==52360
    original112={r['key']:r for r in read(locks['112_nodes']['path'])}
    cert112={r['key']:r for r in read(locks['112_certificates']['path'])}
    assert set(original112)==set(cert112)==set(nodes)
    indices=[];seen=set()
    for i,idx in enumerate(pre['shards']):
        w=read(root/f'worker{i:02d}/result.json');assert w['status']=='complete' and w['worker']==i
        def original_rows():
            for c in w['chunks']:
                assert sha(c['path'])==c['sha256'];yield from read(c['path'])
        old=iter(original_rows());p=out/f'index{i:02d}.jsonl.gz';count=0
        with gzip.open(p,'xt') as f:
            for s in lines(idx['path']):
                r=next(old);k=s['key'];assert k==r['key'] and k not in seen;seen.add(k);count+=1
                assert r['cached162_mean_i']==s['mean162_i'] and r['source112_mean_i']==s['mean112_i']
                assert r['cached164_slope_numerator']==s['slope164_numerator'] and r['cached164_slope_denominator']==s['slope164_denominator']
                assert r['source112_selfU']==s['self112']==nodes[k]['selfU']
                assert r['candidate162_selfU']==s['self162'] and r['candidate164_selfU']==s['self164']
                assert r['new_selfU']==emap[k]['new_selfU'] and s['selfL']<=r['new_selfU']
                o,c=original112[k],cert112[k]
                assert c['packet_sha256']==o['packet_sha256'] and c['old_selfU']==o['old_selfU']
                assert c['new_selfU']==s['self112'] and c['new_selfU']<c['old_selfU']
                assert o['old_L_i']==s['range_L_i'] and o['old_U_i']==s['range_U_i']
                assert c['new_M_i']==s['mean112_i'] and o['old_selfL']==s['selfL']
                assert o['weight']==s['weight'] and o['level']==s['level']
                assert 0<=s['range_L_i']<s['range_U_i'] and s['range_U_i']>=1<<144
                assert s['range_L_i']<=s['mean162_i']<=s['mean112_i']<=s['range_U_i']
                v=dict(s);v.update({'self166':r['new_selfU'],'strict112_active':True,
                    'certificate112_old_selfU':c['old_selfU'],'source112_packet_sha256':c['packet_sha256']})
                f.write(json.dumps(v,sort_keys=True,separators=(',',':'))+'\n')
            assert next(old,None) is None
        assert count==idx['count'];indices.append({'worker':i,'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size,'count':count})
    assert len(seen)==97971 and set(nodes)-seen=={'5:31,31,31,31'}
    selected=[]
    for n in oldtree['nodes']:
        k=n['key'];r=emap[k];assert r['candidate162_effU']==n['candidate162_effU'] and r['candidate164_effU']==n['candidate164_effU']
        selected.append({'key':k,'level':n['level'],'parent':n['parent'],'children':n['children'],
            'self112':n['selfU'],'selfL':n['selfL'],'self166':r['new_selfU'],'eff166':r['new_effU']})
    tree={'nodes':selected,'roots':oldtree['roots'],'baseline166_UH':done['new_UH']}
    p=out/'topology.json.gz'
    with gzip.open(p,'xt') as f:json.dump(tree,f,sort_keys=True,separators=(',',':'))
    inputcheck(src)
    result={'task':'175_metadata_only','status':'all_exact_cached_fields_bound','pid':os.getpid(),'seconds':time.monotonic()-start,
        'sources_sha256':sha(HERE/'sources.json'),'shards':indices,'topology':{'path':str(p),'sha256':sha(p),'bytes':p.stat().st_size},
        'active_nodes':97971,'nodes':97972,'roots':52360,'nonactive_key':'5:31,31,31,31','new_root_slope_gain_DFS_cases':0,
        'mean_source':'qualified162 copied from bound166 metadata','slope_source':'qualified164 copied from bound166 metadata',
        'baseline166_completion_sha256':locks['completion166']['sha256'],'all_inputs_before_after_unchanged':True}
    result['strict_original112_active_chord_bound_by_source_certificate']=True
    result['range_domain_note']='This fixed strictactive subset hasL<U;the secant lemma also allowsL=U by continuity,not a global rejection.'
    put(out/'results.json',result);print({k:v for k,v in result.items() if k not in ['shards','topology']})

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--output',type=Path,required=True);prepare(a.parse_args().output)
