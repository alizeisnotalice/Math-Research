#!/usr/bin/env python3
"""Only frozen fresh137 checkpoints;no H/log;both finite inputs preregistered."""
import sys
sys.dont_write_bytecode=True
import common as C

def raw_gate(entry,task):
    assert entry['ordinal']==task['batch_ordinal'] and entry['global_ordinal']==task['ordinal']
    for f in ('node','point','epsilon','target','point_checkpoint'):assert entry[f]==task[f]
    assert not task['cached'] and task['cache_reuse_source'] is None and task['same_point_H64_guard'] is None
    assert C.sha(entry['point_checkpoint'])==entry['point_checkpoint_sha256']
    raw=C.read(entry['point_checkpoint']);rec=entry['record']
    assert raw['stage']==entry['stage']=='complete' and raw['partial'] is entry['partial'] is None and raw['record']==rec
    assert raw['epsilon']==rec['epsilon']=='1/128' and raw['target']==rec['qualification_target']=='11/100000'
    assert raw['point']==rec['point']==task['point'] and raw['pair_cache_key']==rec['pair_cache_key']==task['pair_cache_key']
    assert raw['model']==rec['model']==C.MODEL and raw['kernel_sha256']==rec['kernel_sha256']==C.KERNEL_SHA
    assert rec['candidate']==task['candidate'] and rec['geometry_tests']==task['geometry_tests']
    assert raw['pid']==entry['worker_pid'] and raw['attempts']==rec['attempts'] and rec['bits']==raw['attempts'][-1]['bits']
    assert len(rec['values'])==5 and raw['attempts'][-1]['numeric_quality'] is True and raw['attempts'][-1]['all_cones']==3504
    assert raw['precision_index'] in (0,1,2) and raw['calls']>=1
    for pair in rec['values']:
        a,b=map(C.F,pair);assert a<=b and (b-a)*2**40<=max(abs(a),abs(b),C.F(1,10**6))
    assert C.F(rec['values'][0][0])>C.TARGET
    return raw
