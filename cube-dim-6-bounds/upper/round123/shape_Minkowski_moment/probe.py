#!/usr/bin/env python3
"""UNRUN: one fixed determinant-moment MC probe, needs primary_release.json.

One CPU, fresh8192 shapes; no point kernel,109 observable or sample reuse.
Saved MC moments and derived budgets are diagnostics, never certificates.
"""
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import resource
import sys
import time
import traceback
from fractions import Fraction as Q

HERE=Path(__file__).resolve().parent
SOURCE=Path(__file__).resolve()

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def save(path,value):
    temp=path.with_name(path.name+'.tmp')
    temp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n'); temp.replace(path)

def cpu():
    usage=resource.getrusage(resource.RUSAGE_SELF)
    return usage.ru_utime+usage.ru_stime

def same_Y_mp(Y,V,dps):
    import mpmath as mp
    with mp.workdps(dps):
        ys=[mp.mpf(float(y)) for y in Y]; R=sum(ys); A=mp.matrix(6)
        assert all(y>0 and mp.isfinite(y) for y in ys)
        for y,v in zip(ys,V):
            for i in range(6):
                for j in range(6): A[i,j]+=(R/y)*int(v[i])*int(v[j])
        eig=mp.eigsy(A,eigvals_only=True); determinant=mp.det(A)
        assert eig[0]>0 and determinant>0 and mp.isfinite(determinant)
        z=determinant**(-mp.mpf(1)/12)
        assert 0<z<=mp.mpf(1)/32
        zf=float(z); condition=float(eig[5]/eig[0])
        assert math.isfinite(zf) and zf>0 and math.isfinite(condition)
        return zf,condition

def main():
    if (HERE/'state.json').exists():
        print('Existing probe state: refusing overwrite or automatic second run.',file=sys.stderr); return 2
    spec=json.loads((HERE/'inputs.json').read_text())
    release=json.loads((HERE/'primary_release.json').read_text())
    assert release['authorized'] is True
    assert release['source_sha256']==sha(SOURCE) and release['input_sha256']==sha(HERE/'inputs.json')
    started={'pid':os.getpid(),'unix_time':time.time(),'source_sha256':sha(SOURCE),
             'input_sha256':sha(HERE/'inputs.json'),'release_sha256':sha(HERE/'primary_release.json')}
    with (HERE/'started.json').open('x') as handle:
        json.dump(started,handle,indent=2,allow_nan=False)
        handle.write('\n')
        handle.flush()
        os.fsync(handle.fileno())
    assert spec['samples']==8192 and spec['chunk_samples']==1024 and spec['chunks']==list(range(8))
    assert spec['workers']==1 and spec['numerical_threads']==1
    keys=('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS')
    for key in keys: os.environ[key]='1'
    import numpy as np
    start=time.monotonic(); c0=cpu(); code=0; active_chunk=None; active_Y=None; active_index=None
    (HERE/'chunks').mkdir(exist_ok=True)
    state={'status':'RUNNING','source_sha256':sha(SOURCE),'input_sha256':sha(HERE/'inputs.json'),
           'pid':os.getpid(),'completed_chunks':[],'fixed_shapes':8192,'workers':1}
    save(HERE/'state.json',state)
    print(json.dumps({'pid':os.getpid(),'fixed_shapes':8192,'workers':1,'observable':'detA^(-1/12)'}),flush=True)
    all_z=[]; all_cond=[]; means=[]; fallback_records=[]
    try:
        V=np.asarray([(1,)+v for v in itertools.product((-1,1),repeat=5)],dtype=np.float64)
        assert np.array_equal(V.T@V,32*np.eye(6))
        outer=np.einsum('vi,vj->vij',V,V).reshape(32,36)
        assert abs(math.exp(-0.5*math.log(1024))-1/32)<1e-15
        Drat=Q(2*math.factorial(32),4**16*math.factorial(16)*math.factorial(15))
        D=float(Drat)*math.sqrt(math.pi)
        Dgamma=2*math.exp(math.lgamma(16.5)-math.lgamma(16))
        assert abs(D/Dgamma-1)<1e-13
        threshold6=Q(8,9)**6*24576**5*Q(1,65000)**6
        threshold=float(threshold6)**(1/6)
        for chunk in range(8):
            active_chunk=chunk; active_index=None
            if (HERE/'STOP.json').exists(): raise InterruptedError('registered STOP marker')
            seed=[20261027123,123,chunk]
            Y=np.random.default_rng(np.random.SeedSequence(seed)).gamma(0.5,1.0,size=(1024,32))
            active_Y=Y
            valid=np.all(np.isfinite(Y)&(Y>0),axis=1)
            if not np.all(valid):
                active_index=int(np.flatnonzero(~valid)[0]); raise RuntimeError('invalid Gamma draw; no replacement allowed')
            R=Y.sum(axis=1); smallest=Y.min(axis=1); ratios=smallest[:,None]/Y
            rsum=ratios.sum(axis=1); B=((ratios/rsum[:,None])@outer).reshape(1024,6,6)
            logTA=np.log(rsum)+np.log(R)-np.log(smallest)
            eig=np.linalg.eigvalsh(B); sign,ld=np.linalg.slogdet(B)
            with np.errstate(divide='ignore',invalid='ignore',under='ignore',over='ignore'):
                conditions=eig[:,-1]/eig[:,0]
                z=np.exp(-0.5*logTA-ld/12)
            bad=(sign<=0)|(eig[:,0]<=0)|(~np.isfinite(conditions))|(conditions>spec['fallback_condition_threshold'])|(~np.isfinite(z))|(z<=0)|(z>1/32)
            indices=np.flatnonzero(bad)
            if len(indices)>spec['max_fallbacks_per_chunk']:
                active_index=int(indices[0]); raise RuntimeError('registered fallback count limit exceeded')
            local=[]
            for index in indices:
                active_index=int(index)
                z[index],conditions[index]=same_Y_mp(Y[index],V,spec['fallback_dps'])
                local.append({'chunk':chunk,'sample_index':int(index),'condition_after_mp':float(conditions[index]),
                              'original_Y_hex':[float(y).hex() for y in Y[index]]})
            assert np.all(np.isfinite(z)&(z>0)&(z<=1/32)) and np.all(np.isfinite(conditions)&(conditions>0))
            fallback_records.extend(local)
            path=HERE/'chunks'/('chunk_%02d.npz'%chunk); tmp=path.with_name(path.name+'.tmp')
            with tmp.open('wb') as handle:
                np.savez_compressed(handle,Z=z,condition=conditions,
                                    max_condition_Y=Y[int(np.argmax(conditions))],max_Z_Y=Y[int(np.argmax(z))])
            tmp.replace(path)
            meta={'chunk':chunk,'seed':seed,'n':1024,'mean_Z':float(z.mean()),'min_Z':float(z.min()),'max_Z':float(z.max()),
                  'max_condition':float(conditions.max()),'fallbacks':local,'data_sha256':sha(path),
                  'source_sha256':sha(SOURCE),'input_sha256':sha(HERE/'inputs.json')}
            save(path.with_suffix('.json'),meta)
            all_z.append(z); all_cond.append(conditions); means.append(meta['mean_Z'])
            state['completed_chunks'].append(chunk); save(HERE/'state.json',state)
            if chunk in (3,7): print(json.dumps({'completed_chunks':chunk+1,'total_chunks':8}),flush=True)
        values=np.concatenate(all_z); conditions=np.concatenate(all_cond)
        assert len(values)==8192
        mean=float(values.mean()); M=D*mean; g0=1/65000
        IM=(3/16)**(6/5)*M**(6/5)/384
        UM=(9/8)**(6/5)*(M/g0)**(6/5)/384
        result={'status':'COMPLETE_MC_DIAGNOSTIC_NOT_EXPECTATION_CERTIFICATE','n':8192,
                'observable':'Z=detA^(-1/12)','mean_Z':mean,'empirical_variance_Z':float(values.var(ddof=1)),
                'chunk_means_Z':means,'Z_quantiles':np.quantile(values,[0,.5,.9,.99,1]).tolist(),
                'condition_quantiles':np.quantile(conditions,[0,.5,.9,.99,1]).tolist(),
                'analytic_Z_bound':'0<Z<=1/32','all_sampled_bounds_pass':True,
                'D_rational_times_sqrt_pi':str(Drat),'D_display':D,
                'D_gamma_vs_factorial_relative_error':abs(D/Dgamma-1),
                'moment_M_empirical':M,'threshold_M_power6_exact':str(threshold6),
                'threshold_M_display':threshold,'threshold_mean_Z_display':threshold/D,
                'empirical_power6_threshold_ratio':M**6/float(threshold6),
                'I_M_empirical_budget_not_certificate':IM,'U_M_at_unproved_g0_empirical_not_certificate':UM,
                'S_M_at_unproved_g0_empirical_not_certificate':64*(g0/6)**(6/5)-IM,
                'fallback_count':len(fallback_records),'fallback_records':fallback_records,
                'dependencies':{'numpy':np.__version__,'python':sys.version,'numpy_blas':np.__config__.CONFIG['Build Dependencies']['blas']},
                'thread_environment':{key:os.environ[key] for key in keys},
                'not_claimed':['rigorous E Z upper bound','certified global g0','actual model obstruction','new C6 bound','deterministic MC error']}
        state['status']='COMPLETE'; save(HERE/'state.json',state)
    except BaseException as exc:
        code=1
        result={'status':'FAILED_NO_PARTIAL_MEAN','exception':repr(exc),'traceback':traceback.format_exc(),
                'chunk':active_chunk,'sample_index':active_index,'completed_chunks':state['completed_chunks']}
        if active_Y is not None:
            if active_index is not None: result['original_Y_hex']=[float(y).hex() for y in active_Y[active_index]]
            else:
                # Batch eigensolver failure has no identified culprit: retain this one finite input batch.
                with (HERE/'failure_batch.npz').open('wb') as handle: np.savez_compressed(handle,Y=active_Y)
                result['failure_batch_sha256']=sha(HERE/'failure_batch.npz')
        state['status']='FAILED'; save(HERE/'state.json',state)
    result.update({'source_sha256':sha(SOURCE),'input_sha256':sha(HERE/'inputs.json'),'pid':os.getpid(),
                   'wall_seconds':time.monotonic()-start,'cpu_seconds':cpu()-c0,'exit_code':code,
                   'child_processes':0,'time_or_CPU_ceiling':None,'partial_samples_renormalized':False})
    save(HERE/'results.json',result)
    print(json.dumps({k:v for k,v in result.items() if k not in ('fallback_records','dependencies','thread_environment','traceback')},indent=2),flush=True)
    return code

if __name__=='__main__': sys.exit(main())
