#!/usr/bin/env python3
"""Zero-PF lock of real4274 fresh jets, ancestry and fixed15-worker plan."""
import sys
sys.dont_write_bytecode=True
import ast,gzip,json,os,difflib
from pathlib import Path
import common as C
from bind import raw_gate

PREPARED='99e932b0b86c9d89efbebd6a187cd26f977679254b4b244fba82b7b63cd56260'
MAIN=C.DESK/'Math-Research-cube-dim-5/cube-dim-5-bounds'

def main():
    assert not (C.HERE/'source_lock.json').exists(),'inspect existing138 first'
    p=C.AUTHOR/'residual_prepared_manifest.json';assert C.sha(p)==PREPARED;m=C.read(p)
    locks={str(p):PREPARED,**m['source_locks'],**m['external_files'],**m['legacy_code_locks'],**{str(C.AUTHOR/n):h for n,h in m['repository_files'].items()}}
    protected={}
    for directory,expected in [('audit_g123_sol','a961e0fd00c580cb495e951835aab2fa7d15408e754ff39670dc1cfcb9704c17'),('upper_join132','fbc30d041ad0d48480501633274d8883e4ee79868f4e19973a33eb48920745da')]:
        d=C.HERE.parent/directory;assert C.sha(d/'manifest.json')==expected;protected[str(d/'manifest.json')]=expected
        for n,h in C.read(d/'manifest.json')['files'].items():protected[str(d/n)]=h
    C.check_locks(locks);C.check_locks(protected)
    tasks=C.read(C.PROD/'tasks.json');assert len(tasks)==4274 and not any(t['cached'] for t in tasks)
    keys={t['pair_cache_key'][1] for t in tasks};assert len(keys)==4274
    old128=C.read(C.PROD/'completed_fulljet_reuse_registry.json');assert len(old128)==10527 and not keys.intersection(old128)
    # Frozen point-time registries are source data,not old H or corner reruns.
    prior=C.read(C.DESK/'cube-finite-dim-scratch/upper_continue123/finite_plan.json')
    for ref in prior['I_scalar_planned_registry_refs'].values():
        jobs=json.loads(gzip.decompress(Path(ref['jobs']).read_bytes()));assert not keys.intersection(','.join(j['point']) for j in jobs)
    for path in [C.DESK/'cube-finite-dim-scratch/upper_g110107/planned_requests.json',C.DESK/'cube-finite-dim-scratch/upper_continue123/planned_tasks.json']:
        assert not keys.intersection(','.join(t['point']) for t in C.read(path))
    # Use genuine EPS64 leaves of the primary-accepted g105 packet;no old q/log/corners.
    packet=MAIN/'upper_portable096/run098_v3';pm=C.read(packet/'manifest.json');pc=C.read(packet/'checker_summary.json');accept=MAIN/'round013/primary_globalg105_acceptance.json';acc=C.read(accept)
    assert acc['status']=='primary_accepts_complete_same_model_global_g105' and pm['g']==pc['g']==acc['qualified_g']=='21/200000'
    assert pc['status']=='PASS_portable_compact_closed_root_packet' and pc['packet_manifest_sha256']==C.sha(packet/'manifest.json')
    for p in [packet/'manifest.json',packet/'checker_summary.json',accept]:locks[str(p)]=C.sha(p)
    wanted={C.node(t['H64_g105_ancestor']['node']) for t in tasks};leaves={}
    for chunk in pm['chunks']:
        path=packet/chunk['file'];assert C.sha(path)==chunk['sha256'];locks[str(path)]=chunk['sha256']
        for leaf in json.loads(gzip.decompress(path.read_bytes()))['leaves']:
            n=C.node(leaf['node'])
            if n in wanted:
                assert leaf['epsilon']==leaf['cache_epsilon']=='1/64' and leaf['curvature_mode']=='full_frame_16'
                leaves[n]=dict(node=list(n),packet=str(path),packet_sha256=chunk['sha256'],leaf_sha256=C.digest(leaf),epsilon='1/64',qualified_lower='21/200000',source_mechanism=leaf['source_mechanism'])
    assert set(leaves)==wanted
    geop=C.DESK/'cube-finite-dim-scratch/audit_g123_sol/source_geometry.json';oldlock=C.read(C.HERE.parent/'audit_g123_sol/source_lock.json')
    assert C.sha(geop)==oldlock['input_locks'][str(geop)];locks[str(geop)]=C.sha(geop)
    oldgeo=C.read(geop);last=C.read(C.DESK/'cube-finite-dim-scratch/upper_continue123/batch0002/state.json')
    roots=oldgeo['roots'];assert roots==[[a,b,c,d] for a in range(32) for b in range(a+1) for c in range(b+1) for d in range(c+1)]
    passed=oldgeo['passed']+[r['node'] for r in last['new_passed']];splits=oldgeo['splits']+[r['node'] for r in last['new_splits']]
    source128=C.DESK/'cube-finite-dim-scratch/upper_curvature128/evaluation001/candidate_results.jsonl'
    passed.extend(r['child_node'] for r in C.lines(source128) if r['passed'])
    old137=C.read(C.DESK/'cube-finite-dim-scratch/upper_parallel137/run001/master_attempt0001_state.json')
    passed.extend(c['child_node'] for c in old137['qualified_children'].values());queue=old137['pending_children']
    assert len(passed)==105635 and len(splits)==57549 and len(queue)==4274
    tree=C.tree(roots,passed+queue,splits);assert set(map(C.node,queue))=={C.node(t['node']) for t in tasks}
    fan=json.loads(gzip.decompress(C.FAN.read_bytes()));ancestry={}
    for i,t in enumerate(tasks):
        assert t['ordinal']==i and t['batch']==i//4096+1 and t['batch_ordinal']==i%4096
        lo,hi=C.box(roots,t['node']);point=list(map(C.F,t['point']));n=C.node(t['node']);anc=C.node(t['H64_g105_ancestor']['node'])
        assert anc in leaves and n[0]==anc[0] and n[1].startswith(anc[1]) and t['H64_g105_ancestor']['lower']=='21/200000'
        assert point[0]==1 and all(a<x<b for a,x,b in zip(lo,point[1:],hi)) and t['box_lower']==list(map(str,lo)) and t['box_upper']==list(map(str,hi))
        al,ah=C.box(roots,list(anc));assert all(a<=x<=b for a,x,b in zip(al,point[1:],ah))
        assert t['epsilon']=='1/128' and t['target']=='11/100000' and t['same_point_H64_guard'] is None
        assert t['pair_cache_key']==['1/128',','.join(map(str,point))] and C.geometry(point,fan)['eligible']
        ancestry[str(i)]=dict(task_node=t['node'],point=t['point'],ancestor_packet=leaves[anc],proof=t['H64_g105_ancestor'])
    C.write(C.SCRATCH/'source_geometry.json',dict(roots=roots,passed=passed,splits=splits,queue=queue,ancestry=ancestry))
    locks[str(C.SCRATCH/'source_geometry.json')]=C.sha(C.SCRATCH/'source_geometry.json')
    registers=[]
    for batch,count in [(1,4096),(2,178)]:
        b=C.PROD/('future_batch%04d'%batch);summary=C.read(b/'summary.json');manifest=C.read(b/'manifest.json');pts=C.read(b/'point_manifest.json');selected=C.read(C.PROD/('input_batch%04d.json'%batch))
        assert summary['status']=='finite_original_prefix_complete' and summary['faults']==summary['obstructions']==summary['inflight']==0
        assert len(pts['records'])==len(selected)==count and pts['epsilon']=='1/128' and pts['target']=='11/100000' and pts['kernel_sha256']==C.KERNEL_SHA
        local={str(b/'manifest.json'):C.sha(b/'manifest.json')}
        for name,h in manifest['files'].items():assert C.sha(b/name)==h;local[str(b/name)]=h
        for e,t in zip(pts['records'],selected):raw_gate(e,t);local[e['point_checkpoint']]=e['point_checkpoint_sha256']
        pids=[summary['master_pid']]+[p['pid'] for p in summary['pool_evidence']];assert len(pids)==16 and all(p['exitcode']==0 and not p['alive'] for p in summary['pool_evidence']);C.assert_pids_absent(pids)
        registers.append(dict(batch=batch,point_count=count,component_count=5*count,point_manifest=str(b/'point_manifest.json'),tasks=str(C.PROD/('input_batch%04d.json'%batch)),input_locks=local,production_pids=pids,workers=15,block_size=128,pilot_count=1,wall_cap=None,actual_value_authorization_still_required=True,old_values_recomputed=False))
        locks.update(local)
    for p in C.HERE.glob('*.py'):ast.parse(p.read_text())
    C.check_locks(locks);C.check_locks(protected)
    code={str(p):C.sha(p) for p in C.HERE.glob('*.py')}
    lock=dict(task='audit_residual138',status='prepared_all4274_real_new_jets_no_PF_yet',input_locks=locks,protected_old_locks=protected,code_locks=code,finite_batches=[4096,178],finite_points=4274,finite_components=21370,workers=15,block_max=128,pilot_each_batch=1,wall_cap=None,automatic_all_frozen_batches_after_review=True,no_newborn_children=True,new_PF=0,new_q=0,new_corners=0,source_tree=tree)
    C.write(C.HERE/'source_lock.json',lock);lsha=C.sha(C.HERE/'source_lock.json')
    for r in registers:
        r.update(source_lock_sha256=lsha,task='audit_residual138_new_values');C.write(C.SCRATCH/('batch%04d_registration.json'%r['batch']),r)
    C.write(C.HERE/'all_plan.json',dict(source_lock_sha256=lsha,finite_points=4274,finite_components=21370,workers=15,blocks_are_not_total_budget=True,wall_cap=None,
        ordered_batches=[dict(batch=r['batch'],registration=str(C.SCRATCH/('batch%04d_registration.json'%r['batch'])),sha256=C.sha(C.SCRATCH/('batch%04d_registration.json'%r['batch']))) for r in registers]))
    deltas={}
    for name in ('numeric.py','structure.py','bind.py'):
        a=(C.HERE.parent/'audit_g123_sol'/name).read_text().splitlines(True);b=(C.HERE/name).read_text().splitlines(True)
        text=''.join(difflib.unified_diff(a,b,fromfile='qualified123/'+name,tofile='audit_residual138/'+name));(C.HERE/(name+'.diff')).write_text(text);deltas[name]=C.sha(C.HERE/(name+'.diff'))
    print(json.dumps(dict(source_lock_sha256=lsha,registrations=[dict(batch=r['batch'],sha256=C.sha(C.SCRATCH/('batch%04d_registration.json'%r['batch']))) for r in registers],new_PF=0,new_q=0,new_corners=0,workers=15,old71812_corner_recomputed=False)),flush=True)

if __name__=='__main__':main()
