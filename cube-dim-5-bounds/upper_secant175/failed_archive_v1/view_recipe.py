"""Field-copy recipe from restored179+tracked112/167. No root/gain/DFS."""
from pathlib import Path
from io175 import read,lines,checked,canonical,gzwrite,ref

def build_views(restored179,repository,output,bindings):
    restored179=Path(restored179);repository=Path(repository);output=Path(output);output.mkdir(exist_ok=False)
    for name,r in bindings['repository'].items():
        checked({'path':str(repository/name),'bytes':r['bytes'],'sha256':r['sha256']})
    for name,r in bindings['restored179'].items():
        checked({'path':str(restored179/name),'bytes':r['bytes'],'sha256':r['sha256']})
    project=repository/'cube-dim-5-bounds'
    index={r['key']:r for r in lines(project/'upper_integral167/packet/node_index.jsonl.gz')}
    orig={r['key']:r for r in read(project/'upper_path112/nodes.json.gz')}
    cert={r['key']:r for r in read(project/'upper_path112/node_certificates.json.gz')}
    tree=read(project/'upper_path112/numeric_hierarchy.json.gz');nodes=tree['nodes']
    old166={r['key']:r for r in read(restored179/'producer166/actual/effective_nodes.json')}
    assert set(index)==set(orig)==set(cert)==set(nodes)==set(old166)
    shards=[];seen=set()
    for i in range(15):
        rows=[]
        for name in bindings['166_chunk_order'][i]:
            for r in read(restored179/name):
                k=r['key'];assert k not in seen;seen.add(k);q=index[k];o=orig[k];c=cert[k];n=nodes[k]
                assert q['L_i']==o['old_L_i'] and q['U_i']==o['old_U_i']
                assert c['packet_sha256']==o['packet_sha256']==q['packet_sha256']
                assert c['old_selfU']==o['old_selfU'] and c['new_selfU']==r['source112_selfU']==q['self112U']==n['selfU']
                assert c['new_selfU']<c['old_selfU'] and c['new_M_i']==r['source112_mean_i']==q['m112_i']
                assert q['level']==o['level']==n['level'] and q['weight']==o['weight']
                s={'key':k,'level':q['level'],'weight':q['weight'],'range_L_i':q['L_i'],'range_U_i':q['U_i'],
                    'mean112_i':r['source112_mean_i'],'mean156_i':q['new_mean_i'],'mean162_i':r['cached162_mean_i'],
                    'self112':r['source112_selfU'],'selfL':n['selfL'],'self162':r['candidate162_selfU'],'self164':r['candidate164_selfU'],
                    'slope164_numerator':r['cached164_slope_numerator'],'slope164_denominator':r['cached164_slope_denominator'],
                    'self166':r['new_selfU'],'strict112_active':True,'certificate112_old_selfU':c['old_selfU'],
                    'source112_packet_sha256':c['packet_sha256']}
                assert s['self166']==old166[k]['new_selfU'];rows.append(s)
        p=output/f'view{i:02d}.jsonl.gz';gzwrite(p,b''.join(canonical(r)+b'\n' for r in rows))
        shards.append({'worker':i,'count':len(rows),**ref(p)})
    assert len(seen)==97971 and set(nodes)-seen=={'5:31,31,31,31'}
    newnodes=[]
    for k,n in nodes.items():
        r=old166[k];newnodes.append({'key':k,'level':n['level'],'parent':n['parent'],'children':n['children'] or [],
            'self112':n['selfU'],'selfL':n['selfL'],'self166':r['new_selfU'],'eff166':r['new_effU']})
    t={'nodes':newnodes,'roots':tree['roots'],'baseline166_UH':bindings['baseline166_UH']}
    p=output/'topology.json.gz';gzwrite(p,canonical(t))
    return {'shards':shards,'topology':ref(p),'H_decodes':0,'root_slope_gain_DFS_calls':0,'original_prepared_gzip_SHA_claim':False}
