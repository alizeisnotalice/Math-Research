#!/usr/bin/env python3
"""A038 bounded continuation from B032 saved complete endpoints; no seed call."""
import sys
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
import hashlib,importlib.util,json,os,subprocess,time,traceback
from fractions import Fraction as Q
from pathlib import Path
HERE=Path(__file__).resolve().parent
B030=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-upper-astra/cube-dim-6-bounds/upper/round030/green_descent_diagnostic')
spec=importlib.util.spec_from_file_location('frozen_B030_rules',B030/'diagnose.py')
rule=importlib.util.module_from_spec(spec);spec.loader.exec_module(rule)
# Only the diagnostic threshold changes in memory; all frozen source files stay read-only.
rule.G0=Q(1,65000)

def read(path):return json.loads(Path(path).read_text())
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def save(name,obj):
    p=HERE/(name+'.tmp');p.write_text(json.dumps(obj,indent=2,ensure_ascii=False)+'\n');os.replace(p,HERE/name)

def run():
    assert not(HERE/'run_receipt.json').exists(),'one batch already started'
    reg=read(HERE/'registration.json')
    assert time.time()<reg['deadline_epoch'] and reg['g0']=='1/65000'
    for p,h in reg['source_sha256'].items():assert sha(p)==h
    assert sha(HERE/'run_batch.py')==reg['driver_sha256']
    _,rays=rule.prep.fan_data();forms=rule.prep.all_zero_forms(rays)
    assert len(forms)==9308
    assert hashlib.sha256(json.dumps(forms,separators=(',',':')).encode()).hexdigest()==reg['zero_forms_sha256']
    seen=set(tuple(row)for row in reg['old_G_exact_point_keys'])
    nums,den=reg['starting_numerators'],reg['starting_denominator']
    source_native=read(reg['starting_native_path'])
    assert rule.point_key(nums,den)in seen
    source_check=rule.validate(source_native,nums,den)
    assert source_check['strictly_above_g0']
    start=time.monotonic()
    receipt=dict(status='RUNNING',driver_pid=os.getpid(),start_epoch=time.time(),native_calls=[],
                 registration_sha256=sha(HERE/'registration.json'),starting_point_recomputed=False,
                 max_native=3,max_wall_seconds=600,new_H_K_p=0,new_C6_bound=False)
    results=dict(status='INCOMPLETE',points=[],starting_point_recomputed=False,new_C6_bound=False)
    save('run_receipt.json',receipt)
    (HERE/'native_outputs').mkdir(exist_ok=False)
    code=1
    try:
        for index in (1,2,3):
            assert time.time()<reg['deadline_epoch']
            candidate=rule.proposal(source_native,nums,den,forms,rays)
            candidate['step_index']=index
            candidate['starting_numerators']=nums;candidate['starting_denominator']=den
            prior_path=Path(reg['starting_native_path']) if index==1 else HERE/f'native_outputs/{index-1:02d}.stdout.json'
            candidate['starting_native_path']=str(prior_path)
            candidate['starting_native_sha256']=sha(prior_path)
            newnums=candidate['final_numerators'];newden=10**14
            if index==1:assert newnums==reg['first_candidate_numerators']
            key=rule.point_key(newnums,newden)
            assert key not in seen,'exact duplicate candidate: stop, no replacement'
            candidate['not_in_prior_G_registry']=True
            save(f'candidate_{index:02d}.json',candidate)
            for p,h in reg['source_sha256'].items():assert sha(p)==h
            assert candidate['final_preflight']['feasible'] and candidate['final_preflight']['int64_safe']
            assert candidate['final_preflight']['zero_form_count']==0
            cmd=[reg['binary_path'],'full',reg['fan_path'],str(newden),*map(str,newnums)]
            rec=dict(index=index,numerators=newnums,denominator=newden,exact_key=list(key),
                     candidate_path=f'candidate_{index:02d}.json',candidate_sha256=sha(HERE/f'candidate_{index:02d}.json'),
                     command=cmd,status='REGISTERED_BEFORE_NATIVE',stdout=f'native_outputs/{index:02d}.stdout.json',
                     stderr=f'native_outputs/{index:02d}.stderr.txt')
            receipt['native_calls'].append(rec);save('run_receipt.json',receipt)
            tic=time.monotonic()
            with(HERE/rec['stdout']).open('w')as stdout,(HERE/rec['stderr']).open('w')as stderr:
                child=subprocess.Popen(cmd,stdout=stdout,stderr=stderr,cwd=HERE)
                rec['native_pid']=child.pid;rec['status']='RUNNING';save('run_receipt.json',receipt)
                try:rc=child.wait(timeout=min(120,reg['deadline_epoch']-time.time()))
                except subprocess.TimeoutExpired:
                    child.kill();rc=child.wait();rec['timed_out']=True
                except BaseException:
                    if child.poll()is None:child.kill()
                    rec['exit_code']=child.wait();raise
            rec.update(exit_code=rc,elapsed_seconds=time.monotonic()-tic,
                       stdout_sha256=sha(HERE/rec['stdout']),stderr_sha256=sha(HERE/rec['stderr']))
            save('run_receipt.json',receipt)
            assert rc==0 and not rec.get('timed_out'),'native failure; no retry'
            out=read(HERE/rec['stdout']);check=rule.validate(out,newnums,newden)
            previous_g=rule.dyadic(source_native['G']);new_g=rule.dyadic(out['G'])
            record=dict(index=index,numerators=newnums,denominator=newden,check=check,
                        G_minus_previous_interval=[str(new_g[0]-previous_g[1]),str(new_g[1]-previous_g[0])],
                        strictly_below_previous=new_g[1]<previous_g[0])
            results['points'].append(record);rec['status']='COMPLETE_VALIDATED';seen.add(key)
            save('results.json',results);save('run_receipt.json',receipt)
            print(json.dumps({'index':index,'G_display':check['G_display_only'],'gradient_display':check['gradient_midpoint_display_only'],
                              'rejects_g0':check['rejects_g0_strictly'],'native_pid':child.pid}),flush=True)
            if check['rejects_g0_strictly']:
                results['status']='STRICT_G0_COUNTEREXAMPLE_IMMEDIATE_STOP';results['counterexample_index']=index
                break
            if index==3:results['status']='THREE_NEW_POINTS_COMPLETED_G0_STILL_UNRESOLVED'
            nums,den,source_native=newnums,newden,out
        save('results.json',results);receipt['status']='COMPLETED';code=0
    except BaseException as exc:
        receipt['status']='STOPPED_FAILURE_NO_RETRY';receipt['error']=repr(exc);receipt['traceback']=traceback.format_exc()
    finally:
        receipt['exit_code']=code;receipt['elapsed_seconds']=time.monotonic()-start
        receipt['all_child_handles_closed']=True
        if(HERE/'results.json').exists():receipt['results_sha256']=sha(HERE/'results.json')
        save('run_receipt.json',receipt)
    return code

if __name__=='__main__':
    assert sys.argv[1:]==['--run']
    raise SystemExit(run())
