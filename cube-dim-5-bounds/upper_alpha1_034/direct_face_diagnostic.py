#!/usr/bin/env python3
"""New2x1024 unitface fullGreen diagnostic; point arithmetic is rigorous, quadrature is not."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,importlib.util,json,os,sys,time

os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('VECLIB_MAXIMUM_THREADS','1')
import numpy as np
import scipy
from scipy.stats import qmc
import sympy as sp

HERE=Path(__file__).resolve().parent
PRIME=2147483647
SEEDS=(34341,34342)
N=1024


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path,obj):
    temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(obj,indent=2)+'\n')
    temp.replace(path)


def main():
    began=time.monotonic()
    regpath=HERE/'registration.json'
    reg=json.loads(regpath.read_text())
    for row in reg['sources']:assert sha(row['source'])==row['sha256'],row['source']
    assert sp.isprime(PRIME)
    source=next(Path(r['source']) for r in reg['sources'] if r['source'].endswith('astra_green018/kernel.py'))
    assert sha(source)=='5045f262ef0e00eb38746f1cf2293ee6cde138a6874e8f6850dc90d28b33ea08'
    spec=importlib.util.spec_from_file_location('frozen_full_green018_I034',source)
    module=importlib.util.module_from_spec(spec)
    sys.modules[spec.name]=module
    spec.loader.exec_module(module)
    module.ctx.prec=256
    engine=module.GreenKernel()
    assert module.sha(module.TABLE)=='8dd2e9895d7465ac599d43feb1b5099506f4e43d1248215eba0a77e11f4d13f9'
    assert len(engine.vertices)==280 and len(engine.simplices)==3504
    fingerprint=hashlib.sha256((sha(regpath)+sha(Path(__file__))+sha(source)+module.TABLE_SHA).encode()).hexdigest()
    cache=Path('/tmp')/('alpha1-I034-'+fingerprint[:16])
    cache.mkdir(exist_ok=True)
    stamp=cache/'fingerprint.json'
    if stamp.exists():assert json.loads(stamp.read_text())['fingerprint']==fingerprint
    else:write_json(stamp,{'fingerprint':fingerprint,'registration':sha(regpath),'script':sha(Path(__file__))})

    all_rows=[]
    completed=[]
    failure=None
    paused=False
    for seed in SEEDS:
        design=qmc.Sobol(d=4,scramble=True,seed=seed).random_base2(10)
        nums=1+np.floor((PRIME-2)*design).astype(np.int64)
        path=cache/f'rows-{seed}.jsonl'
        rows=[]
        if path.exists():
            for line in path.read_text().splitlines():
                row=json.loads(line)
                assert row['seed']==seed and row['index']==len(rows)
                assert row['tail_numerators']==list(map(int,nums[row['index']]))
                rows.append(row)
        assert len(rows)<=N
        with path.open('a',buffering=1) as output:
            for i in range(len(rows),N):
                tail=list(map(int,nums[i]))
                point=(F(1),)+tuple(F(x,PRIME) for x in tail)
                try:
                    assert all(0<x<PRIME for x in tail)
                    assert len(set([PRIME]+tail))==5,'coordinate-collision guard; no substitution allowed'
                    record=engine.evaluate(point,distinct=False)
                    value=module.arb(record['G_ball'])
                    assert value.is_finite() and value>0
                    power=value**(module.arb(5)/4)
                    assert power.is_finite() and power>0
                    row={'seed':seed,'index':i,'tail_numerators':tail,
                         'G_power_5_over4_ball':str(power),'kernel_row':record}
                    output.write(json.dumps(row,separators=(',',':'))+'\n')
                    output.flush()
                    rows.append(row)
                    if seed==SEEDS[0] and i==0:
                        write_json(HERE/'first_new_point_qualification.json',
                                   {'status':'pass','first_point_is_part_of_seed_batch':True,
                                    'extra_old_or_anchor_evaluation':False,'kernel_sha256':sha(source),
                                    'table_sha256':module.TABLE_SHA,'row':row})
                    if (i+1)%256==0:
                        print(json.dumps({'seed':seed,'completed':i+1,'seconds':time.monotonic()-began}),flush=True)
                    if (i+1)%64==0:
                        write_json(cache/'checkpoint.json',{'fingerprint':fingerprint,'seed':seed,'completed':i+1,
                                                          'pid':os.getpid(),'seconds':time.monotonic()-began})
                    if time.monotonic()-began>=110 and i+1<N:
                        paused=True
                        break
                except Exception as exc:
                    failure={'seed':seed,'index':i,'tail_numerators':tail,'point':list(map(str,point)),
                             'exception':type(exc).__name__,'reason':str(exc),
                             'point_dropped_or_replaced':False,'completed_prefix':len(rows)}
                    write_json(cache/'failure.json',failure)
                    break
        all_rows.extend(rows)
        if len(rows)==N:completed.append((seed,rows))
        if failure or paused:break

    def stats(rows):
        summation=sum((module.arb(row['G_power_5_over4_ball']) for row in rows),module.arb(0))
        mean=summation/len(rows)
        blocks=[]
        for j in range(0,len(rows),64):
            piece=rows[j:j+64]
            blocks.append({'start':j,'count':len(piece),
                           'mean_ball':str(sum((module.arb(x['G_power_5_over4_ball']) for x in piece),module.arb(0))/len(piece))})
        maxrow=max(rows,key=lambda x:F(str(module.arb(x['G_power_5_over4_ball']).upper().fmpq())))
        minrow=min(rows,key=lambda x:F(x['kernel_row']['G_upper_dyadic']))
        maximum=module.arb(maxrow['G_power_5_over4_ball'])
        return {'count':len(rows),'discrete_mean_ball':str(mean),
                'discrete_mean_lower_dyadic':str(F(str(mean.lower().fmpq()))),
                'discrete_mean_upper_dyadic':str(F(str(mean.upper().fmpq()))),
                'blocks_descriptive_not_confidence':blocks,
                'maximum_contribution_fraction_of_sum_ball':str(maximum/summation),
                'maximum_contribution_row':maxrow,'minimum_G_upper_row':minrow,
                'min_distinct_projection_count':min(x['kernel_row']['distinct_projection_count'] for x in rows),
                'nonzero_log_arguments_every_record':True,
                'all_records_outside_S2_by_distinct_magnitudes':True}

    done=len(completed)==2 and not failure and not paused
    batch_stats=[{'seed':seed,**stats(rows)} for seed,rows in completed]
    result={'status':'complete_new_unitface_direct_Green_diagnostic' if done else
                     ('failed_no_missing_point_integral_estimate' if failure else 'checkpoint_incomplete'),
            'script_sha256':sha(Path(__file__)),'registration_sha256':sha(regpath),
            'kernel_sha256':sha(source),'table_sha256':module.TABLE_SHA,
            'exact_design':{'seeds':list(SEEDS),'points_each':N,'prime_denominator':PRIME,
                            'mapping':'1+floor((p-2)u), no sorting','fingerprint':fingerprint},
            'completed_kernel_calls':len(all_rows),'completed_batches':len(completed),
            'failure':failure,'paused':paused,'batches':batch_stats,
            'quadrature_error_certified':False,'continuous_I_bound':False,'new_C5_bound':False,
            'normalization':'mean estimates I_face directly;128 timesmean would be positivepartmass atlevel1',
            'arithmetic_scope':'Arb balls cover only the finite rational sample means, not the continuum integral',
            'singular_handling':'No invalid/infinite/zero-node points dropped or replaced; any failure halts full batch.',
            'cache_directory':str(cache),'full_sample_rows_archived':False,
            'runtime':{'pid':os.getpid(),'seconds':time.monotonic()-began,'python':sys.executable,
                       'scipy':scipy.__version__,'numpy':np.__version__,'Arb_bits':module.ctx.prec}}
    if done:
        result['combined_2048_design']=stats(all_rows)
        result['combined_2048_design']['blocks_descriptive_not_confidence']=[]
    result['cache_row_files']=[{'path':str(cache/f'rows-{seed}.jsonl'),
                               'sha256':sha(cache/f'rows-{seed}.jsonl')}
                              for seed in SEEDS if (cache/f'rows-{seed}.jsonl').exists()]
    write_json(HERE/'diagnostic_results.json',result)
    print(json.dumps({'status':result['status'],'calls':len(all_rows),
                      'batch_discrete_means':[x['discrete_mean_ball'] for x in batch_stats],
                      'combined_discrete_mean':result.get('combined_2048_design',{}).get('discrete_mean_ball'),
                      'failure':failure,'seconds':result['runtime']['seconds']},indent=2),flush=True)
    return 0 if done else (2 if failure else 3)


if __name__=='__main__':sys.exit(main())
