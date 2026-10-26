#!/usr/bin/env python3
"""Freeze123 source and batches, exact old-object association; zero new H/log."""
import sys
sys.dont_write_bytecode=True
import ast,gzip,os,time,resource
from pathlib import Path
from fractions import Fraction as F
import common as C

PREPARED='1edbc7912626737020fec146c80247f3f5876603732d97b0d6ecdd45aa46165f'
SEED='52fe7f9f4ecd8d4c16c1c4d48423f78bc09821770fb8e3c5d8899a9dfa4f691b'

def main():
    began=time.monotonic();assert not (C.HERE/'source_lock.json').exists(),'inspect existing preparation first'
    p=C.AUTHOR/'production_prepared_manifest.json';assert C.sha(p)==PREPARED;m=C.read(p)
    locks={str(p):PREPARED,**m['source_locks'],**m['external_files'],**{str(C.AUTHOR/name):h for name,h in m['repository_files'].items()}}
    for name,h in [('audit_g110116','111840a3858a2a8b616d6fe731c16a28f8c04a18515dc9ac8fa39e879a484f2d'),('audit_reuse122','a20ee932eedb8d83e2644bf7fdf1e74ec0aa8f93d310f25ab8ee7e48368d5e40')]:
        p=C.AUDITS/name/'manifest.json';assert C.sha(p)==h;locks[str(p)]=h
        # Consume frozen audit artifacts as oldsource evidence, never rerun values.
        old=C.read(p)
        for relative,value in old['files'].items():
            q=Path(relative) if Path(relative).is_absolute() else p.parent/relative
            locks[str(q)]=value
    for p in [C.AUDITS/'audit_kernel081/common.py',C.AUDITS/'audit_kernel081/arith.py']:
        locks[str(p)]=C.sha(p)
    C.check_locks(locks)
    plan=C.read(C.PROD/'finite_plan.json');tasks=C.read(C.PROD/'planned_tasks.json')
    assert plan['seed_sha256']==SEED and C.sha(C.PROD/'seed.json')==SEED and not (C.PROD/'seed.json').stat().st_mode&0o222
    assert len(tasks)==6237 and [r['count'] for r in plan['batches']]==[4096,2141] and not any(t['cached'] for t in tasks)
    assert C.digest([t['node'] for t in tasks[:4096]])==plan['batches'][0]['selected_nodes_sha256']
    assert C.digest([t['node'] for t in tasks[4096:]])==plan['batches'][1]['selected_nodes_sha256']
    keys={','.join(map(str,map(F,t['point']))) for t in tasks};assert len(keys)==6237
    seed=C.snapshot(C.PROD/'seed.json',keys)
    assert len(seed['passed'])==98007 and len(seed['splits'])==51884 and len(seed['meta']['queue'])==6237
    assert len(seed['cache']['1/64'])==184796 and len(seed['cache']['1/128'])==4290
    assert seed['meta']['target']=='11/100000' and seed['meta']['globalg110_qualified'] is False
    assert not seed['meta']['inflight'] and not seed['meta']['resume_jobs']
    roots=seed['meta']['roots'];assert roots==[[a,b,c,d] for a in range(32) for b in range(a+1) for c in range(b+1) for d in range(c+1)]
    source_tree=C.tree(roots,[n for n,h in seed['passed']]+seed['meta']['queue'],[n for n,h in seed['splits']])
    S=C.DESK/'cube-finite-dim-scratch';old098=C.snapshot(S/'upper_mixed098/final_state.json')
    old105=C.snapshot(S/'upper_g110105/retarget_state.json');old107=C.read(S/'upper_g110107/batch0001/state.json');old118=C.read(S/'upper_reuse118/state.json')
    assert seed['passed']==old105['passed']+[(C.node(r['node']),C.digest(r)) for r in old107['new_passed']+old118['new_passed']]
    assert seed['splits']==old098['splits']+[(C.node(r['node']),C.digest(r)) for r in old107['new_splits']]
    assert roots==old098['meta']['roots'] and seed['meta']['queue']==old118['queue']
    assert seed['cache']['1/64']==old098['cache']['1/64']
    expected128={**old098['cache']['1/128'],**{k:C.digest(v) for k,v in old107['new_cache_by_epsilon']['1/128'].items()}}
    assert len(expected128)==194+4096 and seed['cache']['1/128']==expected128 and not keys.intersection(expected128)
    accepted=C.read(C.PROD/'accepted_mapper.json');assert len(accepted)==98007
    for i,row in enumerate(accepted):
        n,h=seed['passed'][i];assert row['seed_passed_index']==i and C.node(row['node'])==n and row['leaf_sha256']==h
        e,p=row['pair_cache_key'];assert e==row['epsilon'] and row['record_sha256']==seed['cache'][e][p]
    fan=C.read(C.FAN) if C.FAN.suffix!='.gz' else __import__('json').loads(gzip.decompress(C.FAN.read_bytes()))
    assert C.sha(C.FAN)==C.FAN_SHA
    planned={}
    for label,ref in plan['I_scalar_planned_registry_refs'].items():
        assert C.sha(ref['jobs'])==ref['jobs_sha256'] and C.sha(ref['registration'])==ref['registration_sha256']
        jobs=__import__('json').loads(gzip.decompress(Path(ref['jobs']).read_bytes()));times=C.read(ref['registration'])['time_endpoints'];assert F(1,128) in set(map(F,times))
        planned[label]={','.join(map(str,map(F,j['point']))) for j in jobs}
    assert all(not keys.intersection(ps) for ps in planned.values())
    seen=set();guard_count=0
    for i,t in enumerate(tasks):
        lo,hi=C.box(roots,t['node']);p=list(map(F,t['point']));key=','.join(map(str,p));assert p[0]==1 and all(a<x<b for a,x,b in zip(lo,p[1:],hi))
        assert t['ordinal']==i and t['batch']==i//4096+1 and t['batch_ordinal']==i%4096
        assert t['target']=='11/100000' and t['epsilon']=='1/128' and t['pair_cache_key']==['1/128',key] and not t['cached']
        assert t['box_lower']==list(map(str,lo)) and t['box_upper']==list(map(str,hi)) and key not in seen;seen.add(key)
        generic=C.geometry(p,fan);assert generic['eligible'] and t['geometry_tests'][-1]=={'candidate':t['candidate'],**generic}
        g=t['same_point_H64_guard']
        if g is not None:
            assert t['candidate']==t['original_H64_candidate'] and t['point']==t['original_H64_point']==g['point'] and g['epsilon']=='1/64'
            assert g['H_interval']==seed['records']['1/64'][key]['values'][0];guard_count+=1
        elif t['candidate']!=t['original_H64_candidate']:assert g is None
    assert {C.node(t['node']) for t in tasks}==set(map(C.node,seed['meta']['queue']))
    assert [t['node'] for t in tasks[:877]]==seed['meta']['queue'][:877]
    # Static parse never imports numerical modules, and never executes their functions.
    for p in C.HERE.glob('*.py'):ast.parse(p.read_text())
    for name in ('point_manifest.json','state.json','journal.jsonl','initial_overlay.json','registration.json','pool_registration.json','summary.json','manifest.json'):
        p=C.PROD/'batch0001'/name;assert p.exists();locks[str(p)]=C.sha(p)
    assert locks[str(C.PROD/'batch0001/point_manifest.json')]=='40342fbe1d0bf596dc3a31aac741237ad926b17539c8f96a5988f1e3e331b656'
    assert locks[str(C.PROD/'batch0001/state.json')]=='75aaf1ef83f62b2e72c2b524a18af9749ab1ff50f5c276a0610be2ba6834e4a9'
    pm=C.read(C.PROD/'batch0001/point_manifest.json');assert len(pm['records'])==4096 and pm['epsilon']=='1/128' and pm['target']=='11/100000' and pm['model']==C.MODEL and pm['kernel_sha256']==C.KERNEL_SHA
    for entry,task in zip(pm['records'],tasks):
        assert entry['ordinal']==task['batch_ordinal'] and entry['global_ordinal']==task['ordinal'] and entry['node']==task['node'] and entry['point']==task['point'] and not entry['cached']
        path=Path(entry['point_checkpoint']);assert C.sha(path)==entry['point_checkpoint_sha256'];raw=C.read(path)
        assert raw['record']==entry['record'] and raw['partial'] is None and raw['stage']=='complete';locks[str(path)]=entry['point_checkpoint_sha256']
    C.write(C.SCRATCH/'source_geometry.json',dict(roots=roots,passed=[list(n) for n,h in seed['passed']],splits=[list(n) for n,h in seed['splits']],queue=seed['meta']['queue']))
    locks[str(C.SCRATCH/'source_geometry.json')]=C.sha(C.SCRATCH/'source_geometry.json')
    own={str(p):C.sha(p) for p in C.HERE.glob('*.py')};C.check_locks(locks)
    C.write(C.HERE/'source_lock.json',dict(task='audit_g123_sol',status='prepared_no_actual_new_values_authorized',input_locks=locks,code_locks=own,epsilon='1/128',target='11/100000',model=C.MODEL,fan=str(C.FAN),fan_sha256=C.FAN_SHA,finite_batches=[4096,2141],maximum_new_components=31185,workers_proposed=4,wall_cap=None,actual_numeric_authorization_required=True,scope='old source consumed by immutable references;new123 values must be independently enclosed'))
    result=dict(status='source_seed_exact_oldobject_association_and6237_pointtime_preflight_passed',pid=os.getpid(),seconds=time.monotonic()-began,source_tree=source_tree,old_typed_objects=189086,old_passed=98007,old_values_recomputed=0,new_H=0,new_gradient=0,new_log=0,planned_original_points=6237,finite_batches=[4096,2141],same_point_H64_guards=guard_count,all_new_pointtime_keys_unique=True,scalar_registries_disjoint=True,old128_disjoint=True,source_locks=len(locks),peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,global_g110_qualified=False,S_credit=0)
    C.write(C.HERE/'preflight.json',result);print(__import__('json').dumps(result),flush=True)

if __name__=='__main__':main()
