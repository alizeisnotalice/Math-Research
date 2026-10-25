#!/usr/bin/env python3
"""One authorized restricted-face batch; never recompute the saved seed."""
import sys
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
import hashlib,importlib.util,json,os,subprocess,time,traceback
from fractions import Fraction as Q
from pathlib import Path
HERE=Path(__file__).resolve().parent
RULE=HERE.parents[1]/'round030/green_descent_diagnostic/diagnose.py'
spec=importlib.util.spec_from_file_location('frozen_B030_validation',RULE)
rule=importlib.util.module_from_spec(spec);spec.loader.exec_module(rule)
rule.G0=Q(1,65000) # in-memory validation target only; frozen source untouched
DEN=10**14;GAP=Q(1,10**8)
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(n,o):
    p=HERE/(n+'.tmp');p.write_text(json.dumps(o,indent=2)+'\n');os.replace(p,HERE/n)
def check_sources(reg):
    for p,h in reg['source_sha256'].items():assert sha(p)==h,('changed source',p)
    assert sha(__file__)==reg['driver_sha256']

def proposal(out,nums,forms,rays,fixed):
    grad=[sum(rule.dyadic(v))/2 for v in out['gradient']]
    mag=max(abs(grad[4]),abs(grad[5]))
    record={'gradient_exact_midpoints':list(map(str,grad)),'normalization_last_two':str(mag),
            'fixed_first_four_numerators':fixed,'denominator':DEN,'fixed_step':'1/20',
            'repairs':[],'inherits_old_value_certificate':False}
    if not mag:
        record['status']='ZERO_LAST_TWO_GRADIENT_STOP';return record
    direction=[Q(0)]*4+[-grad[j]/mag for j in (4,5)]
    raw=[Q(n,DEN)+d/20 for n,d in zip(nums,direction)]
    clip=raw[:4]
    clip.append(max(2*GAP,min(Q(fixed[3],DEN)-GAP,raw[4])))
    clip.append(max(GAP,min(clip[4]-GAP,raw[5])))
    quant=[rule.prep.nearest(x*DEN) for x in clip]
    record.update(direction=list(map(str,direction)),raw_proposal=list(map(str,raw)),
                  clipped_proposal=list(map(str,clip)),quantized_numerators=quant,
                  quantization_errors=[str(Q(n,DEN)-x) for n,x in zip(quant,clip)])
    assert quant[:4]==fixed and nums[:4]==fixed
    pf=rule.prep.preflight(quant,DEN,forms,rays);record['initial_preflight']=pf
    if not(pf['feasible'] and pf['int64_safe']):
        record['status']='QUANTIZED_INPUT_UNSUPPORTED_STOP';return record
    final=quant
    if pf['zero_form_count']:
        final=None
        for k in range(1,17):
            nn=quant[:];nn[5]+=k
            test=rule.prep.preflight(nn,DEN,forms,rays)
            record['repairs'].append({'k':k,'numerators':nn,'preflight':test})
            if test['feasible'] and test['int64_safe'] and test['zero_form_count']==0:
                final=nn;break
        if final is None:
            record['status']='SIXTEEN_COORDINATE6_REPAIRS_EXHAUSTED_STOP';return record
    record.update(status='SUPPORTED_NEW_EXACT_INPUT',final_numerators=final,
                  final_preflight=rule.prep.preflight(final,DEN,forms,rays),
                  actual_step_gradient_dot=str(sum(grad[j]*Q(final[j]-nums[j],DEN)for j in (4,5))))
    return record

def run():
    assert not(HERE/'run_receipt.json').exists(),'one authorized run only'
    reg=read(HERE/'registration.json');check_sources(reg)
    _,rays=rule.prep.fan_data();forms=rule.prep.all_zero_forms(rays)
    assert len(forms)==9308
    assert hashlib.sha256(json.dumps(forms,separators=(',',':')).encode()).hexdigest()==reg['zero_forms_sha256']
    seen=set(map(tuple,reg['old_point_registry']))
    nums=reg['seed_numerators'];fixed=nums[:4]
    path=Path(reg['seed_native_path']);previous=read(path)
    assert rule.point_key(nums,DEN) in seen
    seedcheck=rule.validate(previous,nums,DEN);assert seedcheck['strictly_above_g0']
    start=time.monotonic()
    receipt={'status':'RUNNING','driver_pid':os.getpid(),'start_epoch':time.time(),
             'max_native_calls':4,'max_wall_seconds':600,'per_native_seconds':120,
             'registration_sha256':sha(HERE/'registration.json'),'seed_recomputed':False,
             'native_calls':[],'new_H_K_p':0,'new_C6_bound':False}
    results={'status':'INCOMPLETE','points':[],'seed_check':seedcheck,'new_C6_bound':False}
    save('run_receipt.json',receipt);save('results.json',results)
    (HERE/'native_outputs').mkdir(exist_ok=False)
    code=1
    try:
        for i in range(1,5):
            assert time.monotonic()-start<600
            c=proposal(previous,nums,forms,rays,fixed)
            c.update(index=i,previous_numerators=nums,previous_native_path=str(path),previous_native_sha256=sha(path))
            save(f'candidate_{i:02}.json',c)
            if c['status']!='SUPPORTED_NEW_EXACT_INPUT':
                results['status']=c['status'];break
            nn=c['final_numerators'];key=rule.point_key(nn,DEN)
            if key in seen:
                results['status']='EXACT_OLD_POINT_DUPLICATE_STOP';break
            c['not_in_old_registry']=True;save(f'candidate_{i:02}.json',c)
            check_sources(reg)
            cmd=[reg['binary_path'],'full',reg['fan_path'],str(DEN),*map(str,nn)]
            rec={'index':i,'numerators':nn,'denominator':DEN,'command':cmd,
                 'candidate_path':f'candidate_{i:02}.json','candidate_sha256':sha(HERE/f'candidate_{i:02}.json'),
                 'stdout':f'native_outputs/{i:02}.stdout.json','stderr':f'native_outputs/{i:02}.stderr.txt',
                 'status':'REGISTERED_BEFORE_NATIVE'}
            receipt['native_calls'].append(rec);save('run_receipt.json',receipt)
            tic=time.monotonic()
            with(HERE/rec['stdout']).open('w')as stdout,(HERE/rec['stderr']).open('w')as stderr:
                child=subprocess.Popen(cmd,stdout=stdout,stderr=stderr,cwd=HERE)
                rec['native_pid']=child.pid;rec['status']='RUNNING';save('run_receipt.json',receipt)
                try:rc=child.wait(timeout=min(120,600-(time.monotonic()-start)))
                except subprocess.TimeoutExpired:
                    child.kill();rc=child.wait();rec['timed_out']=True
                except BaseException:
                    if child.poll()is None:child.kill()
                    rec['exit_code']=child.wait();raise
            rec.update(exit_code=rc,elapsed_seconds=time.monotonic()-tic,
                       stdout_sha256=sha(HERE/rec['stdout']),stderr_sha256=sha(HERE/rec['stderr']))
            save('run_receipt.json',receipt)
            assert rc==0 and not rec.get('timed_out'),'native failure: no retry'
            assert (HERE/rec['stderr']).read_bytes()==b''
            out=read(HERE/rec['stdout']);v=rule.validate(out,nn,DEN)
            oldg=rule.dyadic(previous['G']);newg=rule.dyadic(out['G'])
            row={'index':i,'numerators':nn,'denominator':DEN,'check':v,
                 'strict_descent':newg[1]<oldg[0],
                 'new_minus_previous_interval':[str(newg[0]-oldg[1]),str(newg[1]-oldg[0])]}
            results['points'].append(row);seen.add(key);rec['status']='COMPLETE_VALIDATED'
            save('results.json',results);save('run_receipt.json',receipt)
            print(json.dumps({'index':i,'native_pid':child.pid,'G_display':v['G_display_only'],
                              'last2_gradient_display':v['gradient_midpoint_display_only'][4:],
                              'strict_descent':row['strict_descent'],'rejects_g0':v['rejects_g0_strictly'],
                              'repair_k':c['repairs'][-1]['k']if c['repairs']else 0}),flush=True)
            if v['rejects_g0_strictly']:
                results['status']='STRICT_G0_COUNTEREXAMPLE_IMMEDIATE_STOP';break
            if not row['strict_descent']:
                results['status']='NO_CERTIFIED_STRICT_DESCENT_IMMEDIATE_STOP';break
            if i==4:
                results['status']='FOUR_STEPS_COMPLETE_GLOBAL_FLOOR_UNRESOLVED';break
            nums=nn;previous=out;path=HERE/rec['stdout']
        check_sources(reg);save('results.json',results);receipt['status']='COMPLETED';code=0
    except BaseException as exc:
        receipt.update(status='STOPPED_FAILURE_NO_RETRY',error=repr(exc),traceback=traceback.format_exc())
        results['status']='STOPPED_FAILURE_NO_RETRY';save('results.json',results)
    finally:
        receipt.update(exit_code=code,elapsed_seconds=time.monotonic()-start,all_child_handles_closed=True,
                       results_sha256=sha(HERE/'results.json'))
        save('run_receipt.json',receipt)
    return code

if __name__=='__main__':
    assert sys.argv[1:]==['--run']
    raise SystemExit(run())
