#!/usr/bin/env python3
"""B030: at most two calls to the frozen C013 engine, no Green implementation."""
import sys
sys.dont_write_bytecode = True
sys.set_int_max_str_digits(0)
import hashlib, importlib.util, itertools, json, math, os, subprocess, time, traceback
from fractions import Fraction as Q
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[2]
PREP_PATH = PROJECT/'upper/round015/spatial_minimum_search/prepare.py'
spec = importlib.util.spec_from_file_location('frozen_B015_preflight', PREP_PATH)
prep = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prep)
DEN, GAP, G0 = 10**14, Q(1,10**8), Q(1,60000)

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def read(path):
    return json.loads(Path(path).read_text())

def save(name, obj):
    tmp = HERE/(name+'.tmp')
    tmp.write_text(json.dumps(obj,indent=2,ensure_ascii=False)+'\n')
    os.replace(tmp,HERE/name)

def dyadic(row):
    assert set(row) == {'lower_mantissa','upper_mantissa','common_binary_exponent'}
    a,b,e = map(int,(row['lower_mantissa'],row['upper_mantissa'],row['common_binary_exponent']))
    assert a <= b and abs(e) < 10000
    return Q(a)*Q(2)**e,Q(b)*Q(2)**e

def separate_dyadic(row):
    return (Q(int(row['lower_mantissa']))*Q(2)**int(row['lower_binary_exponent']),
            Q(int(row['upper_mantissa']))*Q(2)**int(row['upper_binary_exponent']))

def point_key(nums,den):
    return tuple(str(Q(n,den)) for n in nums)

def validate(out,nums,den):
    assert out['status']=='generic_nonzero_supported' and out['precision_bits']==512
    assert out['full_group_terms']==967680 and out['computed_sign_pairs']==483840
    assert out['zero_node_terms']==0 and 0<=out['confluent_terms']<=967680
    assert out['confluent_terms']%2==0 and out['unique_node_values']>0
    assert len(out['gradient'])==6
    g=dyadic(out['G']);grad=[dyadic(row)for row in out['gradient']]
    assert g[0]>0
    for lo,hi in [g,*grad]:
        assert hi-lo <= max(Q(1),abs((lo+hi)/2))*Q(2)**-256, 'insufficient precision: stop, no retry'
    euler=dyadic(out['Euler_x_grad_plus5G'])
    assert euler[0]<=0<=euler[1]
    lo,hi=5*g[0],5*g[1]
    for n,(a,b)in zip(nums,grad):
        pair=sorted([Q(n,den)*a,Q(n,den)*b]);lo+=pair[0];hi+=pair[1]
    assert lo<=0<=hi
    # Any S2 point has at most two absolute coordinate values.
    x=list(map(lambda n:Q(n,den),nums))
    distance=max(min(abs(x[i]-x[j]) for i,j in itertools.combinations(triple,2))/2
                 for triple in itertools.combinations(range(6),3))
    assert distance>0
    return dict(status='PASS_FULL_GREEN_AND_GRADIENT_INTERVALS',G_interval=list(map(str,g)),
                gradient_intervals=[list(map(str,p))for p in grad],
                independent_Euler_interval=[str(lo),str(hi)],
                rejects_g0_strictly=g[1]<G0,strictly_above_g0=g[0]>G0,
                G_upper_minus_g0=str(g[1]-G0),G_lower_minus_g0=str(g[0]-G0),
                G_display_only=float((g[0]+g[1])/2),
                gradient_midpoint_display_only=[float((a+b)/2)for a,b in grad],
                exact_dist_infinity_to_S2_lower_bound=str(distance),
                necessary_B_ceiling_for_cost_below64=str(Q(16,3)*g[1]),new_C6_bound=False)

def proposal(out,nums,den,forms,rays):
    grad=[sum(dyadic(row))/2 for row in out['gradient']]
    magnitude=max(map(abs,grad[1:]))
    assert magnitude>0, 'zero free gradient: no supported candidate'
    direction=[Q(0)]+[-v/magnitude for v in grad[1:]]
    raw=[Q(n,den)+d/20 for n,d in zip(nums,direction)]
    clipped=[Q(1)]
    for i in range(1,6):
        independently=max((6-i)*GAP,min(1-i*GAP,raw[i]))
        clipped.append(min(independently,clipped[-1]-GAP))
    assert clipped[-1]>=GAP and all(clipped[i]-clipped[i+1]>=GAP for i in range(5))
    quantized=[prep.nearest(v*DEN)for v in clipped]
    pf=prep.preflight(quantized,DEN,forms,rays)
    record=dict(gradient_exact_dyadic_midpoints=list(map(str,grad)),normalization_max_free_abs=str(magnitude),
                negative_normalized_direction=list(map(str,direction)),fixed_step='1/20',
                raw_rational_proposal=list(map(str,raw)),clipped_rational_proposal=list(map(str,clipped)),
                quantized_numerators=quantized,denominator=DEN,
                quantization_errors=[str(Q(n,DEN)-v)for n,v in zip(quantized,clipped)],
                initial_preflight=pf,repairs=[],inherits_any_prior_point_certificate=False)
    assert pf['feasible'] and pf['int64_safe'], 'quantized candidate unsupported'
    final=quantized
    if pf['zero_form_count']:
        final=None
        for repair in prep.gap_repairs(quantized,DEN):
            test=prep.preflight(repair['numerators'],DEN,forms,rays)
            repair['preflight']=test
            record['repairs'].append(repair)
            if test['feasible'] and test['int64_safe'] and test['zero_form_count']==0:
                final=repair['numerators'];break
        assert final is not None, 'six inward generic repairs exhausted; stop'
    record['final_numerators']=final
    record['final_preflight']=prep.preflight(final,DEN,forms,rays)
    record['clipped_step_gradient_dot_displacement']=str(sum(grad[i]*(Q(final[i],DEN)-Q(nums[i],den))for i in range(1,6)))
    return record

def run():
    assert not(HERE/'run_receipt.json').exists(),'batch already started, no automatic rerun'
    reg=read(HERE/'registration.json')
    for path,digest in reg['source_sha256'].items():
        assert sha(path)==digest,('frozen input changed',path)
    assert sha(HERE/'diagnose.py')==reg['driver_sha256']
    _,rays=prep.fan_data();forms=prep.all_zero_forms(rays)
    digest=hashlib.sha256(json.dumps(forms,separators=(',',':')).encode()).hexdigest()
    assert digest==reg['zero_forms_sha256'] and len(forms)==9308
    seed=reg['seed_numerators'];den=reg['seed_denominator']
    assert prep.preflight(seed,den,forms,rays)==reg['seed_preflight']
    assert reg['seed_preflight']['zero_form_count']==0
    seen=set(tuple(p)for p in reg['old_G_exact_point_keys'])
    assert point_key(seed,den) not in seen
    receipt=dict(status='RUNNING',driver_pid=os.getpid(),start_epoch=time.time(),max_wall_seconds=1200,
                 max_native_calls=2,native_calls=[],registration_sha256=sha(HERE/'registration.json'),
                 new_G_engine=False,new_H=0,new_K=0,new_p=0,new_C6_bound=False)
    save('run_receipt.json',receipt)
    start=time.monotonic()
    outdir=HERE/'native_outputs';outdir.mkdir(exist_ok=False)
    results=dict(status='INCOMPLETE',points=[],new_C6_bound=False)
    code=1
    def call(nums,den,index):
        key=point_key(nums,den)
        assert key not in seen,'duplicate exact Green input: stop'
        pf=prep.preflight(nums,den,forms,rays)
        assert pf['feasible'] and pf['int64_safe'] and pf['zero_form_count']==0
        assert len(receipt['native_calls'])<2 and time.monotonic()-start<1200
        for path,digest in reg['source_sha256'].items():
            assert sha(path)==digest,('source changed before launch',path)
        cmd=[reg['binary_path'],'full',reg['fan_path'],str(den),*map(str,nums)]
        rec=dict(index=index,numerators=nums,denominator=den,exact_key=list(key),preflight=pf,
                 command=cmd,status='REGISTERED_BEFORE_CALL',stdout=f'native_outputs/{index:02d}.stdout.json',
                 stderr=f'native_outputs/{index:02d}.stderr.txt')
        receipt['native_calls'].append(rec);save('run_receipt.json',receipt)
        t=time.monotonic()
        with (HERE/rec['stdout']).open('w')as stdout,(HERE/rec['stderr']).open('w')as stderr:
            child=subprocess.Popen(cmd,stdout=stdout,stderr=stderr,cwd=HERE)
            rec['native_pid']=child.pid;rec['status']='RUNNING';save('run_receipt.json',receipt)
            try:
                rc=child.wait(timeout=min(120,1200-(time.monotonic()-start)))
            except subprocess.TimeoutExpired:
                child.kill();rc=child.wait();rec['timed_out']=True
            except BaseException:
                if child.poll()is None:child.kill()
                rec['exit_code']=child.wait();rec['status']='INTERRUPTED_CHILD_REAPED';raise
        rec.update(exit_code=rc,elapsed_seconds=time.monotonic()-t,
                   stdout_sha256=sha(HERE/rec['stdout']),stderr_sha256=sha(HERE/rec['stderr']))
        save('run_receipt.json',receipt)
        assert rc==0 and not rec.get('timed_out'),'native failed or timed out: no retry'
        out=read(HERE/rec['stdout']);check=validate(out,nums,den)
        results['points'].append(dict(index=index,numerators=nums,denominator=den,check=check))
        rec['status']='COMPLETE_VALIDATED';seen.add(key)
        save('results.json',results);save('run_receipt.json',receipt)
        print(json.dumps(dict(index=index,G=check['G_display_only'],grad=check['gradient_midpoint_display_only'],rejects_g0=check['rejects_g0_strictly'],native_pid=child.pid)),flush=True)
        return out,check
    try:
        out,check=call(seed,den,1)
        if check['rejects_g0_strictly']:
            results['status']='STRICT_G0_COUNTEREXAMPLE_FIRST_POINT_STOP'
        else:
            candidate=proposal(out,seed,den,forms,rays)
            save('candidate.json',candidate)
            out2,check2=call(candidate['final_numerators'],DEN,2)
            results['status']=('STRICT_G0_COUNTEREXAMPLE_SECOND_POINT_STOP'if check2['rejects_g0_strictly']
                              else'TWO_POINT_DIAGNOSTIC_COMPLETE_GLOBAL_G_UNRESOLVED')
        oldc=dyadic(read(reg['old_c_native_path'])['G'])
        oldx=separate_dyadic(read(reg['old_x0_constants_path'])['G_from_old_Gamma'])
        results['old_G_intervals']={'old_c':list(map(str,oldc)),'x0':list(map(str,oldx))}
        for row in results['points']:
            gl,gu=map(Q,row['check']['G_interval'])
            row['strictly_below_old_c']=gu<oldc[0];row['strictly_below_x0']=gu<oldx[0]
            row['strictly_above_x0']=gl>oldx[1]
        if len(results['points'])==2:
            g1=list(map(Q,results['points'][0]['check']['G_interval']))
            g2=list(map(Q,results['points'][1]['check']['G_interval']))
            results['second_strictly_below_first']=g2[1]<g1[0]
            results['second_minus_first_interval']=[str(g2[0]-g1[1]),str(g2[1]-g1[0])]
        results['public_and_internal_C6_upper_at_registration']=reg['bounds_snapshot']['upper_exact']
        save('results.json',results);receipt['status']='COMPLETED';code=0
    except BaseException as exc:
        receipt['status']='STOPPED_NO_RETRY';receipt['error']=repr(exc);receipt['traceback']=traceback.format_exc()
    finally:
        receipt['exit_code']=code;receipt['elapsed_seconds']=time.monotonic()-start
        receipt['all_child_handles_closed']=True
        if(HERE/'results.json').exists():receipt['results_sha256']=sha(HERE/'results.json')
        save('run_receipt.json',receipt)
    return code

if __name__=='__main__':
    assert sys.argv[1:]==['--run']
    raise SystemExit(run())
