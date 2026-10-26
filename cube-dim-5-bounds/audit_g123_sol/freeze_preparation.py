#!/usr/bin/env python3
"""Finalize pre-review code/zero-H probe locks; never reruns old values."""
import sys
sys.dont_write_bytecode=True
import ast,os
import common as C

def main():
    p=C.HERE/'source_lock.json';lock=C.read(p);assert lock['status']=='prepared_no_actual_new_values_authorized'
    assert not list(C.SCRATCH.glob('numeric_batch*')),'code cannot change after actual audit starts'
    C.check_locks(lock['input_locks']);pre=C.read(C.HERE/'preflight.json');fault=C.read(C.HERE/'fault_probe.json');symbolic=C.read(C.HERE/'symbolic_probe.json')
    assert pre['new_H']==pre['new_gradient']==pre['new_log']==0 and pre['planned_original_points']==6237
    assert fault['new_H']==fault['new_gradient']==fault['new_log']==0 and fault['shutdown_code_sha256']==C.sha(C.HERE/'numeric.py')
    assert all(not x['alive'] and x['exitcode']==0 for x in fault['pool_exits'])
    assert symbolic['residual']=='0' and symbolic['new_H']==symbolic['new_numerical_log']==0
    for q in C.HERE.glob('*.py'):ast.parse(q.read_text())
    lock['code_locks']={str(q):C.sha(q) for q in C.HERE.glob('*.py')}
    for q in [C.SCRATCH/'source_geometry.json',C.HERE/'preflight.json',C.HERE/'fault_probe.json',C.HERE/'symbolic_probe.json']:
        lock['input_locks'][str(q)]=C.sha(q)
    lock['final_preparation_pid']=os.getpid();lock['zero_H_probes_passed']=True
    C.write(p,lock)
    C.write(C.HERE/'preparation_manifest.json',dict(task='audit_g123_sol_preparation',status='awaiting_main_review_and_actualaudit_authorization',source_lock_sha256=C.sha(p),
        files={q.name:C.sha(q) for q in sorted(C.HERE.iterdir()) if q.is_file() and q.name!='preparation_manifest.json'},new_H=0,new_gradient=0,new_log=0,actual_value_batches_authorized=False))
    print(__import__('json').dumps(dict(source_lock_sha256=C.sha(p),manifest_sha256=C.sha(C.HERE/'preparation_manifest.json'),code_files=len(lock['code_locks']),locks=len(lock['input_locks']),new_H=0)))

if __name__=='__main__':main()
