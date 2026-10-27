"""Bind frozen producer175 actual bytes/identities only; no mathematical cases."""
from pathlib import Path,PurePosixPath
from datetime import datetime,timezone
import os,time
from io175 import sha,read,put,lock
HERE=Path(__file__).resolve().parent
ROOT=Path('/Users/zhengzhihao/Desktop/T')
INDEX=ROOT/'cube-finite-dim-scratch/upper_secant175/dispatch001/output_manifest.json'

def main():
    start=time.monotonic();b=read(HERE/'source_binding.json');lock(b['source_locks'])
    assert sha(INDEX)=='bec0b6b7761c8e0d990279d01a22ff69bd1f45e0964fec768c83435e7c9a462c'
    index=read(INDEX);assert len(index['files'])==473 and index['bytes']==177937642
    out=ROOT/'cube-finite-dim-scratch/upper_secant175/run001';locks={}
    for name,row in index['files'].items():
        rel=PurePosixPath(name);assert not rel.is_absolute() and '..' not in rel.parts
        locks['output_'+name]={'path':str(out/name),'sha256':row['sha256']}
    lock(locks)
    def add(name,path,h=None):
        actual=sha(path)
        if h is not None:assert h==actual
        locks[name]={'path':str(path),'sha256':actual}
    add('actual_output_index',INDEX)
    author=Path(b['producer_directory']);p=read(author/'plan.json');assert sha(author/'plan.json')==b['producer_plan_sha256']
    add('producer_freeze',author/'manifest.json','b76868291983c0ef0c62e25cd8e600d426a8e2d27cd41aeacc536c00a9f5dba8')
    add('producer_publication',author/'candidate_binding.json','30dd76ffe0061d099f7d904b68433ec43cee77ee77e69f9518c28b20502b2b25')
    add('producer_execution_manifest',author/'execution_manifest.json','8cf809d8e2209a42b458840b67776336bde8f37c8a977a372f215f557d6b1af2')
    reg,r,done=read(out/'registration.json'),read(out/'results.json'),read(out/'completion.json')
    ident={'plan_sha256':b['producer_plan_sha256'],'preflight_sha256':p['preflight']['sha256'],
        'consumer_sha256':p['code_and_source_hashes']['consumer.py'],'sources_sha256':p['code_and_source_hashes']['sources.json']}
    assert reg['identity']==ident and all(done[k]==v for k,v in r.items())
    assert done['all_inputs_before_after_unchanged'] and done['registration_sha256']==sha(out/'registration.json')
    assert r['pid']==reg['pid']==60969 and r['new_fourth_roots']==195942 and r['new_gains']==97971
    assert r['nodes']==97972 and r['roots']==52360 and r['global_Cn_claim'] is False
    codes=read(out/'worker_exitcodes.json');assert codes['exitcodes']==[0]*15 and codes['pids']==r['worker_pids']
    for i in range(15):
        wr=read(out/f'worker{i:02d}/result.json');assert wr['status']=='complete' and wr['worker']==i and wr['identity']==ident
        assert wr['completed']==b['shards'][i]['count'] and wr['pid']==codes['pids'][i]
        assert read(out/f'worker{i:02d}/registration.json')['source_index']==b['shards'][i]['producer_shard']
    lock(locks);lock(b['source_locks'])
    put(HERE/'candidate_binding.json',{'task':'175_frozen_actual_producer_identity_binding','execution_ready':True,
        'source_binding_sha256':sha(HERE/'source_binding.json'),'producer_plan_sha256':b['producer_plan_sha256'],
        'candidate_dir':str(out),'candidate_locks':locks,'producer_identity':ident,
        'new_roots_slopes_gains_DFS':0,'source_ready_only_separate_primary_math_authorization_required':True})
    put(HERE/'binding_receipt.json',{'pid':os.getpid(),'utc':datetime.now(timezone.utc).isoformat(),'seconds':time.monotonic()-start,
        'code_sha256':sha(__file__),'candidate_binding_sha256':sha(HERE/'candidate_binding.json'),'new_roots_slopes_gains_DFS':0})
    print('PASS frozen producer bytes only; independent mathematical audit not run')

if __name__=='__main__':main()
