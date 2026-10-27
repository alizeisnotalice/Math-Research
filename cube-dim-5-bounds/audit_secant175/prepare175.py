"""Metadata only: independently bind175 inputs to certified169 source and output."""
from pathlib import Path
from datetime import datetime,timezone
import argparse,gzip,json,os,time
from io175 import sha,read,put,lines,lock
HERE=Path(__file__).resolve().parent
ROOT=Path('/Users/zhengzhihao/Desktop/T')
AUTHOR=ROOT/'Math-Research-cube5-curvature128-sol/cube-dim-5-bounds/upper_secant175'
OLD=HERE.parent/'audit_combine169'
DEST=ROOT/'cube-finite-dim-scratch/audit_secant175/prepared001'

def main(psha):
    start=time.monotonic();DEST.mkdir(parents=True,exist_ok=False)
    assert sha(AUTHOR/'plan.json')==psha;plan=read(AUTHOR/'plan.json');src=read(AUTHOR/'sources.json')
    lock(src['source_locks']);b=read(OLD/'source_binding.json');c=read(OLD/'candidate_binding.json')
    lock(b['source_locks']);lock(c['candidate_locks'])
    assert sha(OLD/'run001/results.json')=='75ac91b7c743dad1e705aa4b60dfce39b2ea8d033eb6becd011afb74e0dabe59'
    own=read(OLD/'run001/results.json');assert own['new_gain_cases']==97971 and own['nodes']==97972 and own['worker_exitcodes']==[0]*15
    accepted=read(src['source_locks']['accepted166']['path'])
    assert accepted['status']=='primary_accepts_complete_I166_and_C5_upper' and accepted['new_UH']==own['new_UH']
    assert accepted['sources']['audit169']['sha256']==sha(OLD/'run001/results.json')
    locks={**{'prior169_'+k:v for k,v in b['source_locks'].items()},**{'prior169candidate_'+k:v for k,v in c['candidate_locks'].items()},
           **{'producer175_'+k:v for k,v in src['source_locks'].items()}}
    def add(name,path):locks[name]={'path':str(path),'sha256':sha(path)}
    for name in ['plan.json',*plan['code_and_source_hashes']]:
        add('producer175_'+name,AUTHOR/name)
    for name in ['source_binding.json','candidate_binding.json','execution_manifest.json','run001/results.json','run001/tree_checks.json']:
        add('qualified169_'+name,OLD/name)
    m=ROOT/'Math-Research-cube-dim-5/cube-dim-5-bounds'
    for name in ['round020/secant_precision.md','round020/secant_precision_registration.json','round020/secant_precision_result.json']:
        add('precision_'+name,m/name)
    pre=read(plan['preflight']['path']);assert sha(plan['preflight']['path'])==plan['preflight']['sha256']
    assert pre['new_root_slope_gain_DFS_cases']==0 and pre['status']=='all_exact_cached_fields_bound'
    assert pre['strict_original112_active_chord_bound_by_source_certificate'] is True
    add('producer175_preflight',plan['preflight']['path']);add('producer175_preflight_registration',Path(plan['preflight']['path']).parent/'registration.json')
    treeold=read(b['tree']['path']);assert sha(b['tree']['path'])==b['tree']['sha256']
    actual=read(OLD/'run001/tree_checks.json')['effective_nodes'];assert len(actual)==97972
    rows112={r['key']:r for r in read(src['source_locks']['112_nodes']['path'])}
    cert112={r['key']:r for r in read(src['source_locks']['112_certificates']['path'])}
    shards=[];seen=set()
    for i,idx in enumerate(pre['shards']):
        assert idx['worker']==i and sha(idx['path'])==idx['sha256']
        old=b['shards'][i];assert sha(old['path'])==old['sha256']
        wr=read(OLD/f'run001/worker{i:02d}/result.json');assert wr['status']=='complete'
        assert sha(OLD/f'run001/worker{i:02d}/updates.json')==wr['updates_sha256']
        updates=iter(read(OLD/f'run001/worker{i:02d}/updates.json'));given=iter(lines(idx['path']))
        path=DEST/f'input{i:02d}.jsonl.gz';count=0
        with gzip.open(path,'xt') as out:
            for s in lines(old['path']):
                key,value=next(updates);assert key==s['key'] and key not in seen and value==actual[key]['new_selfU']
                o,c112=rows112[key],cert112[key]
                assert c112['new_selfU']==s['self112']<c112['old_selfU']==o['old_selfU']
                assert c112['packet_sha256']==o['packet_sha256'] and c112['new_M_i']==s['mean112_i']
                assert o['old_L_i']==s['range_L_i'] and o['old_U_i']==s['range_U_i'] and o['old_selfL']==s['selfL']
                assert o['level']==s['level'] and o['weight']==s['weight']
                expected={**s,'self166':value,'strict112_active':True,
                    'certificate112_old_selfU':c112['old_selfU'],'source112_packet_sha256':c112['packet_sha256']}
                assert next(given)==expected
                assert 0<=s['range_L_i']<s['range_U_i'] and s['range_U_i']>=1<<144
                assert s['range_L_i']<=s['mean162_i']<=s['mean112_i']<=s['range_U_i']
                assert s['slope164_numerator']>0 and s['slope164_denominator']>0 and 16*s['slope164_numerator']>=s['slope164_denominator']
                out.write(json.dumps(expected,sort_keys=True,separators=(',',':'))+'\n');seen.add(key);count+=1
        assert next(updates,None) is None and next(given,None) is None and count==idx['count']
        shards.append({'worker':i,'count':count,'path':str(path),'sha256':sha(path),'producer_shard':idx})
    assert len(seen)==97971 and set(treeold['nodes'])-seen=={'5:31,31,31,31'}
    st=read(pre['topology']['path']);assert sha(pre['topology']['path'])==pre['topology']['sha256']
    supplied={n['key']:n for n in st['nodes']};assert set(supplied)==set(treeold['nodes'])
    slim={}
    for key,n in treeold['nodes'].items():
        r=actual[key];expected={k:n[k] for k in ['key','level','parent','children']}
        expected.update({'self112':n['selfU'],'selfL':n['selfL'],'self166':r['new_selfU'],'eff166':r['new_effU']})
        assert supplied[key]==expected;slim[key]=expected
    assert st['roots']==treeold['roots'] and st['baseline166_UH']==own['new_UH']
    path=DEST/'tree.json.gz'
    with gzip.open(path,'xt') as out:json.dump({'nodes':slim,'roots':st['roots'],'baseline166_UH':own['new_UH']},out,separators=(',',':'))
    lock(locks)
    put(HERE/'source_binding.json',{'task':'175_independent_metadata_only_source_binding','producer_plan_sha256':psha,'producer_directory':str(AUTHOR),
        'producer_preflight':plan['preflight'],'producer_sources':src['source_locks'],'source_locks':locks,'shards':shards,
        'tree':{'path':str(path),'sha256':sha(path)},'new_roots_slopes_gains_DFS':0,'old_H_mean_slope_witness_calls':0,
        'actual_producer175_bound':False})
    put(HERE/'preparation_receipt.json',{'pid':os.getpid(),'utc':datetime.now(timezone.utc).isoformat(),'seconds':time.monotonic()-start,
        'code_sha256':sha(__file__),'source_binding_sha256':sha(HERE/'source_binding.json'),'new_roots_slopes_gains_DFS':0})
    print('PASS metadata only; source fields copied; no real root/slope/gain/DFS')

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--producer-plan-sha',required=True);main(a.parse_args().producer_plan_sha)
