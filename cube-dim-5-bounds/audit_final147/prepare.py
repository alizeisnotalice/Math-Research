#!/usr/bin/env python3
"""Frozen-count/path contract and real input binding;zero newPF/log/q/corners."""
import sys
sys.dont_write_bytecode=True
import ast,json,gzip,os,difflib
from pathlib import Path
import common as C
from bind import raw_gate

def main():
    assert not (C.HERE/'source_lock.json').exists(),'inspect existing143;no duplicate setup'
    cfg=C.CONFIG;N=cfg['fixed_points'];assert type(N) is int and N>0 and cfg['batches']==[1] and cfg['workers']==15
    assert type(cfg['block_size']) is int and 1<=cfg['block_size']<=128
    pp=Path(cfg['prepared_manifest']);assert C.sha(pp)==cfg['prepared_manifest_sha256'];m=C.read(pp)
    assert len(m['finite_batches'])==1 and m['maximum_new_H128']==N and m['fulljet_hits']==0
    locks={str(pp):cfg['prepared_manifest_sha256'],str(C.HERE/'config.json'):C.sha(C.HERE/'config.json'),**m['source_locks'],**m['external_files'],**m['legacy_code_locks'],**{str(C.AUTHOR/n):h for n,h in m['repository_files'].items()}}
    for path,sha in [(cfg['source_result_manifest'],cfg['source_result_manifest_sha256']),(cfg['actual_result_manifest'],cfg['actual_result_manifest_sha256'])]:assert C.sha(path)==sha;locks[path]=sha
    old=Path(cfg['predecessor_dir']);assert C.sha(old/'manifest.json')==cfg['predecessor_manifest_sha256']
    oldmanifest=C.read(old/'manifest.json');oldlock=C.read(old/'source_lock.json')
    protected={str(old/'manifest.json'):C.sha(old/'manifest.json'),**{str(old/n):h for n,h in oldmanifest['files'].items()},**oldlock['protected_old_locks']}
    C.check_locks(locks);C.check_locks(protected)
    taskpath=Path(cfg['task_file']);assert C.sha(taskpath)==cfg['task_sha256'];tasks=C.read(taskpath);assert len(tasks)==N
    assert tasks==C.read(C.PROD/'tasks.json') and all(not t['cached'] for t in tasks)
    keys={t['pair_cache_key'][1] for t in tasks};assert len(keys)==N
    oldkeys=set()
    for item in cfg['completed_key_sources']:
        p=Path(item['path']);locks[str(p)]=C.sha(p);value=C.read(p)
        if item['schema']=='dictkeys':oldkeys.update(value)
        else:assert item['schema']=='taskkeys';oldkeys.update(t['pair_cache_key'][1] for t in value)
    assert len(oldkeys)==cfg['expected_old128'] and not keys.intersection(oldkeys)
    for path,h in list(locks.items()):
        if path.endswith('jobs.json.gz'):
            assert C.sha(path)==h;jobs=json.loads(gzip.decompress(Path(path).read_bytes()));assert not keys.intersection(','.join(t['point']) for t in jobs)
    recipe=cfg['source_geometry_recipe'];basepath=Path(recipe['base']);assert C.sha(basepath)==oldlock['input_locks'][str(basepath)]
    base=C.read(basepath);overlaypath=Path(recipe['overlay']);assert C.sha(overlaypath)==oldlock['input_locks'][str(overlaypath)];overlay=C.read(overlaypath)
    promotionspath=Path(recipe['promotions']);promotions=C.read(promotionspath)
    for p in (basepath,overlaypath,promotionspath):locks[str(p)]=C.sha(p)
    roots=base['roots'];assert roots==[[a,b,c,d] for a in range(32) for b in range(a+1) for c in range(b+1) for d in range(c+1)]
    passed=base['passed']+[r['node'] for r in overlay['new_passed']]+[r['child_node'] for r in promotions['qualified_children'].values()]
    splits=base['splits']+[r['node'] for r in overlay['new_splits']];queue=promotions['pending_children'];counts=recipe['expected_counts']
    assert len(passed)==counts['passed'] and len(splits)==counts['splits'] and len(queue)==counts['pending']==N
    assert set(map(C.node,queue))=={C.node(t['node']) for t in tasks};tree=C.tree(roots,passed+queue,splits)
    # Reuse the exact genuine EPS64 packet-member receipts qualified by138.
    ancestors={C.node(a['ancestor_packet']['node']):a['ancestor_packet'] for a in base['ancestry'].values()};ancestry={}
    fan=json.loads(gzip.decompress(C.FAN.read_bytes()))
    for i,t in enumerate(tasks):
        assert t['ordinal']==t['batch_ordinal']==i and t['batch']==1 and t['epsilon']=='1/128' and t['target']=='11/100000'
        point=list(map(C.F,t['point']));key=','.join(map(str,point));assert point[0]==1 and t['pair_cache_key']==['1/128',key] and t['same_point_H64_guard'] is None
        lo,hi=C.box(roots,t['node']);assert all(a<x<b for a,x,b in zip(lo,point[1:],hi)) and t['box_lower']==list(map(str,lo)) and t['box_upper']==list(map(str,hi))
        n=C.node(t['node']);proof=t['H64_g105_ancestor'];anc=C.node(proof['node']);assert anc in ancestors and proof['lower']==ancestors[anc]['qualified_lower']=='21/200000' and ancestors[anc]['epsilon']=='1/64'
        assert anc[0]==n[0] and n[1].startswith(anc[1]);al,ah=C.box(roots,list(anc));assert all(a<=x<=b for a,x,b in zip(al,point[1:],ah)) and C.geometry(point,fan)['eligible']
        transfer=t['transfer_unexecuted_backup2350_reservation'];assert transfer['point']==t['point'] and C.sha(transfer['backup_path'])==transfer['backup_sha256']
        ancestry[str(i)]=dict(task_node=t['node'],point=t['point'],ancestor_packet=ancestors[anc],proof=proof)
    C.write(C.SCRATCH/'source_geometry.json',dict(roots=roots,passed=passed,splits=splits,queue=queue,ancestry=ancestry));locks[str(C.SCRATCH/'source_geometry.json')]=C.sha(C.SCRATCH/'source_geometry.json')
    b=C.PROD/'future_batch0001';assert C.sha(b/'point_manifest.json')==cfg['actual_point_manifest_sha256'];pm=C.read(b/'point_manifest.json');manifest=C.read(b/'manifest.json');summary=C.read(b/'summary.json')
    assert len(pm['records'])==N and pm['epsilon']=='1/128' and pm['target']=='11/100000' and pm['kernel_sha256']==C.KERNEL_SHA
    assert summary['status']=='finite_original_prefix_complete' and summary['faults']==summary['obstructions']==summary['inflight']==0
    local={str(b/'manifest.json'):C.sha(b/'manifest.json')}
    for n,h in manifest['files'].items():assert C.sha(b/n)==h;local[str(b/n)]=h
    for e,t in zip(pm['records'],tasks):raw_gate(e,t);local[e['point_checkpoint']]=e['point_checkpoint_sha256']
    pids=[summary['master_pid']]+[r['pid'] for r in summary['pool_evidence']];assert len(pids)==16 and all(r['exitcode']==0 and not r['alive'] for r in summary['pool_evidence']);C.assert_pids_absent(pids)
    locks.update(local)
    for p in C.HERE.glob('*.py'):ast.parse(p.read_text())
    C.check_locks(locks);C.check_locks(protected);code={str(p):C.sha(p) for p in C.HERE.glob('*.py')}
    lock=dict(task=C.HERE.name,status='prepared_real_fixed_inputs_noPF',input_locks=locks,protected_old_locks=protected,code_locks=code,finite_points=N,finite_components=5*N,finite_batches=[N],workers=15,block_max=128,pilot_each_batch=1,wall_cap=None,new_PF=0,new_q=0,new_corners=0,source_tree=tree,source_recipe=recipe)
    C.write(C.HERE/'source_lock.json',lock);lsha=C.sha(C.HERE/'source_lock.json')
    reg=dict(audit_id=C.HERE.name,task=C.HERE.name+'_new_values',batch=1,point_count=N,component_count=5*N,point_manifest=str(b/'point_manifest.json'),tasks=str(taskpath),input_locks=local,production_pids=pids,workers=15,block_size=cfg['block_size'],pilot_count=1,wall_cap=None,actual_value_authorization_still_required=True,old_values_recomputed=False,source_lock_sha256=lsha)
    rp=C.SCRATCH/'batch0001_registration.json';C.write(rp,reg)
    C.write(C.HERE/'all_plan.json',dict(audit_id=C.HERE.name,source_lock_sha256=lsha,finite_points=N,finite_components=5*N,workers=15,blocks_are_not_total_budget=True,wall_cap=None,ordered_batches=[dict(batch=1,registration=str(rp),sha256=C.sha(rp))]))
    for name in ('numeric.py','structure.py','run_all.py'):
        text=''.join(difflib.unified_diff((old/name).read_text().splitlines(True),(C.HERE/name).read_text().splitlines(True),fromfile='qualified_predecessor/'+name,tofile=C.HERE.name+'/'+name));(C.HERE/(name+'.diff')).write_text(text)
    print(json.dumps(dict(source_lock_sha256=lsha,all_plan_sha256=C.sha(C.HERE/'all_plan.json'),registration_sha256=C.sha(rp),fixed_points=N,fixed_components=5*N,new_PF=0,new_q=0,new_corners=0,workers=15)),flush=True)

if __name__=='__main__':main()
