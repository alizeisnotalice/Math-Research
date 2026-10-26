#!/usr/bin/env python3
"""Zero-H bind of one complete original production batch before value audit."""
import sys
sys.dont_write_bytecode=True
import argparse,time,os
from pathlib import Path
import common as C

def raw_gate(entry,task):
    assert entry['ordinal']==task['batch_ordinal'] and entry['global_ordinal']==task['ordinal']
    for field in ('node','point','pair_cache_key','cached','epsilon','target','point_checkpoint'):assert entry[field]==task[field]
    assert not task['cached'] and C.sha(entry['point_checkpoint'])==entry['point_checkpoint_sha256']
    raw=C.read(entry['point_checkpoint']);rec=entry['record']
    assert raw['stage']==entry['stage']=='complete' and raw['partial'] is None and raw['record']==rec
    assert raw['epsilon']==rec['epsilon']=='1/128' and raw['target']==rec['qualification_target']=='11/100000'
    assert raw['point']==rec['point']==task['point'] and raw['pair_cache_key']==rec['pair_cache_key']==task['pair_cache_key']
    assert raw['kernel_sha256']==rec['kernel_sha256']==C.KERNEL_SHA and raw['model']==rec['model']==C.MODEL
    assert rec['candidate']==task['candidate'] and rec['geometry_tests']==task['geometry_tests'] and raw['pid']==entry['worker_pid']
    assert raw['attempts']==rec['attempts'] and rec['bits']==raw['attempts'][-1]['bits']
    assert len(rec['values'])==5 and len(raw['attempts'])>=1 and raw['attempts'][-1]['numeric_quality'] is True
    assert raw['attempts'][-1]['all_cones']==3504 and raw['precision_index'] in (0,1,2) and raw['calls']>=1
    for a,b in map(lambda pair:map(C.F,pair),rec['values']):assert a<=b and (b-a)*2**40<=max(abs(a),abs(b),C.F(1,10**6))
    assert C.F(rec['values'][0][0])>C.TARGET
    guard=task['same_point_H64_guard']
    if guard is not None:assert C.F(rec['values'][0][1])>=C.F(guard['H_interval'][0])
    if task['candidate']!=task['original_H64_candidate']:assert guard is None
    return raw

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--batch',type=int,choices=(1,2),required=True);ap.add_argument('--reviewed-source-lock-sha',required=True);ap.add_argument('--production-manifest-sha',required=True);a=ap.parse_args()
    assert C.sha(C.HERE/'source_lock.json')==a.reviewed_source_lock_sha
    lock=C.read(C.HERE/'source_lock.json');C.check_locks(lock['code_locks']);C.check_locks(lock['input_locks'])
    b=C.PROD/('batch%04d'%a.batch);m=b/'manifest.json';assert C.sha(m)==a.production_manifest_sha
    manifest=C.read(m);locks={str(C.HERE/'source_lock.json'):a.reviewed_source_lock_sha,str(m):a.production_manifest_sha}
    for name,h in manifest['files'].items():assert C.sha(b/name)==h;locks[str(b/name)]=h
    summary=C.read(b/'summary.json');assert summary['status']=='finite_original_prefix_complete' and summary['faults']==summary['obstructions_H128']==summary['inflight']==0
    assert summary['no_children_dispatch_in_this_plan'] and summary['S_credit']==0 and summary['globalg110_qualified'] is False
    tasks=C.read(C.PROD/('input_batch%04d.json'%a.batch));pm=C.read(b/'point_manifest.json')
    assert len(pm['records'])==len(tasks)==(4096 if a.batch==1 else 2141)
    assert pm['epsilon']=='1/128' and pm['target']=='11/100000' and pm['model']==C.MODEL and pm['kernel_sha256']==C.KERNEL_SHA
    for entry,task in zip(pm['records'],tasks):raw_gate(entry,task);locks[entry['point_checkpoint']]=entry['point_checkpoint_sha256']
    p=C.SCRATCH/('batch%04d_registration.json'%a.batch);assert not p.exists(),'inspect existing registration; no blind rebind'
    result=dict(task='audit_g123_sol_new_values',batch=a.batch,pid=os.getpid(),registered_unix=time.time(),point_count=len(tasks),component_count=len(tasks)*5,
        source_lock_sha256=a.reviewed_source_lock_sha,input_locks=locks,production_manifest_sha256=a.production_manifest_sha,point_manifest=str(b/'point_manifest.json'),tasks=str(C.PROD/('input_batch%04d.json'%a.batch)),
        workers=4,block_size=128,pilot_count=1,wall_cap=None,actual_value_authorization_still_required=True,old_values_recomputed=False)
    C.write(p,result);print(__import__('json').dumps(dict(registration=str(p),sha256=C.sha(p),new_points=len(tasks),new_H=0,new_gradient=0,new_log=0)))

if __name__=='__main__':main()
