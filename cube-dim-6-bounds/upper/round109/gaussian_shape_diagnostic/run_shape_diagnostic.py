#!/usr/bin/env python3
"""UNRUN, needs primary_mc_release.json binding this source and frozen inputs.

One finite 12-process iid diagnostic. Child numerical libraries import only
after the initializer limits their own BLAS environment. No time ceiling.
"""
import argparse
import concurrent.futures as futures
import hashlib
import itertools
import json
import math
import multiprocessing as multiprocessing
import os
from pathlib import Path
import platform
import resource
import sys
import time
import traceback

HERE=Path(__file__).resolve().parent
SOURCE=Path(__file__).resolve()
G={}

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def dump(path,obj):
    tmp=path.with_name(path.name+'.tmp.'+str(os.getpid()))
    tmp.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
    tmp.replace(path)

def cpu_seconds():
    r=resource.getrusage(resource.RUSAGE_SELF)
    return r.ru_utime+r.ru_stime

def init_worker(spec,point_rows,K):
    for key in spec['child_thread_environment']: os.environ[key]='1'
    import numpy as np
    X=np.asarray([row['numerators'] for row in point_rows],dtype=np.float64)/4096
    V=np.asarray([(1,)+v for v in itertools.product((-1,1),repeat=5)],dtype=np.float64)
    signs=np.asarray(list(itertools.product((-1,1),repeat=6)),dtype=np.float64)
    pairs=[(i,j) for i in range(6) for j in range(i,6)]
    def phi(Z):
        return np.asarray([Z[:,i]*Z[:,j] for i,j in pairs]).T
    Z=(X[:,None,:]*signs[None,:,:]).reshape(16*64,6)
    PHI=phi(Z); CPHI=phi(signs)
    assert np.array_equal(V.T@V,32*np.eye(6))
    # A=1024I at P_i=1/32, so its sample integrand is K/(32 |x|^5).
    diag=np.asarray([1.0 if i==j else 0.0 for i,j in pairs])
    q=PHI@diag
    sanity=(K/32)*q**(-2.5)
    reference=K/(32*np.sum(X*X,axis=1)**2.5)
    assert np.max(np.abs(sanity.reshape(16,64)/reference[:,None]-1))<1e-13
    G.update(np=np,X=X,V=V,signs=signs,pairs=pairs,PHI=PHI,CPHI=CPHI,
             OUTER=np.einsum('vi,vj->vij',V,V).reshape(32,36),spec=spec,K=K,
             worker_init_wall=time.time(),identity_sign_index=63)
    dump(HERE/'workers'/('worker_'+str(os.getpid())+'.json'),
         {'pid':os.getpid(),'numpy':np.__version__,'python':sys.version,
          'environment':{key:os.environ[key] for key in spec['child_thread_environment']},
          'numpy_blas':np.__config__.CONFIG['Build Dependencies']['blas'],
          'uniform_shape_sanity_pass':True})

def fallback_same_shape(Y):
    """80 dps recomputation of the SAME float gamma inputs, no resampling."""
    import mpmath as mp
    np=G['np']; spec=G['spec']
    with mp.workdps(spec['high_precision_dps']):
        ys=[mp.mpf(float(y)) for y in Y]; R=sum(ys)
        assert all(y>0 and mp.isfinite(y) for y in ys)
        A=mp.matrix(6)
        for y,v in zip(ys,G['V']):
            w=R/y
            for i in range(6):
                for j in range(6): A[i,j]+=w*int(v[i])*int(v[j])
        det=mp.det(A); assert det>0 and mp.isfinite(det)
        inv=A**-1; pref=mp.mpf(G['K'])/mp.sqrt(det)
        vals=[]
        for x in G['X']:
            reflected=[]
            for signs in G['signs']:
                z=[mp.mpf(float(x[i]))*int(signs[i]) for i in range(6)]
                q=sum(z[i]*inv[i,j]*z[j] for i in range(6) for j in range(6))
                assert q>0 and mp.isfinite(q)
                reflected.append(pref*q**(-mp.mpf('2.5')))
            vals.append(reflected)
        corners=[]
        for signs in G['signs']:
            q=sum(int(signs[i])*inv[i,j]*int(signs[j]) for i in range(6) for j in range(6))
            assert q>0 and mp.isfinite(q)
            corners.append(q)
        floor=pref*max(corners)**(-mp.mpf('2.5'))
        eigen=mp.eigsy(A,eigvals_only=True)
        assert eigen[0]>0
        cond=float(eigen[5]/eigen[0])
        B=A/(sum(A[i,i] for i in range(6))/6)
        residual=float(mp.norm(B*(B**-1)-mp.eye(6),p=mp.inf))
        v=np.asarray([[float(t) for t in row] for row in vals])
        floor=float(floor)
        assert np.all(np.isfinite(v)) and np.all(v>0) and math.isfinite(floor) and floor>0
        assert math.isfinite(cond)
        return v[:,63],np.mean(v,axis=1),floor,cond,residual

def process_chunk(chunk_id):
    np=G['np']; spec=G['spec']; t0=time.monotonic(); c0=cpu_seconds()
    if (HERE/'STOP.json').exists(): raise InterruptedError('cooperative STOP before chunk')
    Y=None; active_index=None
    try:
        seed=[20261027,109,chunk_id]
        rng=np.random.default_rng(np.random.SeedSequence(seed))
        Y=rng.gamma(shape=0.5,scale=1.0,size=(2048,32))
        if not np.all(np.isfinite(Y)) or not np.all(Y>0):
            bad=np.flatnonzero(~np.all(np.isfinite(Y)&(Y>0),axis=1)); active_index=int(bad[0])
            raise RuntimeError('zero or nonfinite gamma: no replacement allowed')
        # Normalize with minY/Y, avoiding formation of an overflowing reciprocal.
        R=Y.sum(axis=1); minY=Y.min(axis=1); Z=minY[:,None]/Y
        Zsum=Z.sum(axis=1); weights=Z/Zsum[:,None]
        logTA=np.log(Zsum)+np.log(R)-np.log(minY)
        B=(weights@G['OUTER']).reshape(2048,6,6)
        eigen=np.linalg.eigvalsh(B)
        sign,logdet=np.linalg.slogdet(B)
        with np.errstate(divide='ignore',invalid='ignore'):
            cond=eigen[:,-1]/eigen[:,0]
        bad=(~np.isfinite(cond))|(eigen[:,0]<=0)|(sign<=0)|(~np.isfinite(logdet))|(cond>spec['high_precision_trigger_condition'])
        plain=np.empty((2048,16)); sym=np.empty((2048,16)); floor=np.empty(2048)
        residual=np.full(2048,np.nan)
        normal=np.flatnonzero(~bad)
        if len(normal):
            inverse=np.linalg.solve(B[normal],np.broadcast_to(np.eye(6),(len(normal),6,6)))
            residual[normal]=np.max(np.abs(B[normal]@inverse-np.eye(6)),axis=(1,2))
            features=np.asarray([inverse[:,i,j] if i==j else inverse[:,i,j]+inverse[:,j,i]
                                 for i,j in G['pairs']]).T
            q=features@G['PHI'].T; cq=features@G['CPHI'].T
            suspicious=(~np.all(np.isfinite(q)&(q>0),axis=1))|(~np.all(np.isfinite(cq)&(cq>0),axis=1))|(residual[normal]>spec['high_precision_trigger_residual'])
            good_local=np.flatnonzero(~suspicious); good=normal[good_local]
            with np.errstate(over='ignore',under='ignore',invalid='ignore'):
                pref=G['K']*np.exp(-0.5*logTA[good]-0.5*logdet[good])
                vals=pref[:,None]*q[good_local]**(-2.5)
                cf=pref*np.max(cq[good_local],axis=1)**(-2.5)
            finite=np.all(np.isfinite(vals)&(vals>0),axis=1)&np.isfinite(cf)&(cf>0)
            valid=good[finite]; vv=vals[finite].reshape(len(valid),16,64)
            plain[valid]=vv[:,:,63]; sym[valid]=vv.mean(axis=2); floor[valid]=cf[finite]
            bad[normal[suspicious]]=True; bad[good[~finite]]=True
        fallback_indices=np.flatnonzero(bad)
        if len(fallback_indices)>spec['maximum_high_precision_shapes_per_chunk']:
            active_index=int(fallback_indices[0]); raise RuntimeError('registered fallback count limit exceeded')
        fallback_records=[]
        for index in fallback_indices:
            active_index=int(index)
            p,s,f,c,r=fallback_same_shape(Y[index])
            plain[index]=p; sym[index]=s; floor[index]=f; cond[index]=c; residual[index]=r
            fallback_records.append({'sample_index':int(index),'condition_after_mp':c,'original_Y_hex':[float(y).hex() for y in Y[index]]})
        assert np.all(np.isfinite(plain)&(plain>0)) and np.all(np.isfinite(sym)&(sym>0))
        assert np.all(np.isfinite(floor)&(floor>0)) and np.all(np.isfinite(residual))
        gridmin=sym.min(axis=1)  # UPPER envelope of the true orbit minimum.
        assert np.all(gridmin>=floor*(1-1e-10))
        bins=np.searchsorted(np.asarray(spec['condition_bin_upper_edges']),cond,side='left')
        bc=np.asarray([np.count_nonzero(bins==b) for b in range(6)])
        bs=np.asarray([sym[bins==b].sum(axis=0) for b in range(6)])
        bp=np.asarray([plain[bins==b].sum(axis=0) for b in range(6)])
        arrays={
            'plain_sum':plain.sum(axis=0),'sym_sum':sym.sum(axis=0),
            'plain_sumsq':np.einsum('ij,ij->j',plain,plain),'sym_sumsq':np.einsum('ij,ij->j',sym,sym),
            'top_plain':np.sort(np.partition(plain,-16,axis=0)[-16:],axis=0),
            'top_sym':np.sort(np.partition(sym,-16,axis=0)[-16:],axis=0),
            'gridmin_sum':np.asarray([gridmin.sum(),gridmin@gridmin]),
            'cornerfloor_sum':np.asarray([floor.sum(),floor@floor]),
            'condition_counts':bc,'condition_sym_sums':bs,'condition_plain_sums':bp,
            'conditions':cond,'solve_residuals':residual,
            'worst_condition_Y':Y[int(np.argmax(cond))],
            'worst_plain_Y':Y[int(np.argmax(plain.max(axis=1)))],
            'worst_sym_Y':Y[int(np.argmax(sym.max(axis=1)))]}
        if (HERE/'STOP.json').exists(): raise InterruptedError('cooperative STOP before checkpoint')
        path=HERE/'chunks'/('chunk_%02d.npz'%chunk_id)
        tmp=path.with_name(path.name+'.tmp.'+str(os.getpid()))
        with tmp.open('wb') as handle: np.savez_compressed(handle,**arrays)
        tmp.replace(path)
        meta={'chunk_id':chunk_id,'pid':os.getpid(),'seed_sequence':seed,'n':2048,'status':'COMPLETE',
              'cpu_seconds':cpu_seconds()-c0,'wall_seconds':time.monotonic()-t0,
              'max_condition':float(cond.max()),'max_solve_residual':float(residual.max()),
              'fallback_count':len(fallback_indices),'fallback_records':fallback_records,
              'checkpoint_sha256':sha(path),'source_sha256':sha(SOURCE),
              'input_sha256':sha(HERE/'mc_inputs.json'),'points_sha256':sha(HERE/'points.json')}
        dump(path.with_suffix('.json'),meta)
        return meta
    except BaseException as exc:
        failure={'status':'SCIENTIFIC_OR_RUNTIME_FAILURE_NO_PARTIAL_MEAN','chunk_id':chunk_id,'pid':os.getpid(),
                 'sample_index':active_index,'exception':repr(exc),'traceback':traceback.format_exc(),
                 'source_sha256':sha(SOURCE),'input_sha256':sha(HERE/'mc_inputs.json'),
                 'scientific_failure':not isinstance(exc,(InterruptedError,KeyboardInterrupt))}
        if Y is not None and active_index is not None:
            failure['original_Y_hex']=[float(y).hex() for y in Y[active_index]]
        dump(HERE/'failures'/('chunk_%02d.json'%chunk_id),failure)
        raise

def aggregate(spec,metas,K):
    import numpy as np
    acc={}; conditions=[]; residuals=[]; ptop=[]; stop=[]
    for meta in sorted(metas,key=lambda t:t['chunk_id']):
        path=HERE/'chunks'/('chunk_%02d.npz'%meta['chunk_id'])
        assert sha(path)==meta['checkpoint_sha256'] and meta['n']==2048
        with np.load(path) as data:
            for key in ('plain_sum','sym_sum','plain_sumsq','sym_sumsq','gridmin_sum','cornerfloor_sum','condition_counts','condition_sym_sums','condition_plain_sums'):
                acc[key]=acc.get(key,0)+data[key]
            conditions.append(data['conditions']); residuals.append(data['solve_residuals'])
            ptop.append(data['top_plain']); stop.append(data['top_sym'])
    n=sum(m['n'] for m in metas); assert n==65536
    pm=acc['plain_sum']/n; sm=acc['sym_sum']/n
    pv=(acc['plain_sumsq']-n*pm*pm)/(n-1); sv=(acc['sym_sumsq']-n*sm*sm)/(n-1)
    assert np.all(pv>=0) and np.all(sv>=0)
    pt=np.sort(np.concatenate(ptop),axis=0)[-16:]; st=np.sort(np.concatenate(stop),axis=0)[-16:]
    chunk_plain=[]; chunk_sym=[]
    for meta in sorted(metas,key=lambda t:t['chunk_id']):
        with np.load(HERE/'chunks'/('chunk_%02d.npz'%meta['chunk_id'])) as data:
            chunk_plain.append(data['plain_sum']/2048); chunk_sym.append(data['sym_sum']/2048)
    chunk_plain=np.asarray(chunk_plain); chunk_sym=np.asarray(chunk_sym)
    perpoint=[]
    for j in range(16):
        perpoint.append({'point_id':j,'plain_mean':float(pm[j]),'sym_mean':float(sm[j]),
                         'plain_empirical_variance':float(pv[j]),'sym_empirical_variance':float(sv[j]),
                         'empirical_variance_ratio_plain_over_sym':float(pv[j]/sv[j]),
                         'plain_max_share':float(pt[-1,j]/acc['plain_sum'][j]),
                         'sym_max_share':float(st[-1,j]/acc['sym_sum'][j]),
                         'plain_top16_share':float(pt[:,j].sum()/acc['plain_sum'][j]),
                         'sym_top16_share':float(st[:,j].sum()/acc['sym_sum'][j]),
                         'chunk_sym_means':chunk_sym[:,j].tolist(),
                         'chunk_plain_mean_range':[float(chunk_plain[:,j].min()),float(chunk_plain[:,j].max())],
                         'chunk_sym_mean_range':[float(chunk_sym[:,j].min()),float(chunk_sym[:,j].max())]})
    empirical_min=float(sm.min()); gridmean=float(acc['gridmin_sum'][0]/n)
    assert empirical_min>=gridmean*(1-1e-12)
    return {'status':'COMPLETE_FINITE_MC_DIAGNOSTIC_ONLY','n':n,'K':K,'per_point':perpoint,
            'finite_grid_min_point_id':int(np.argmin(sm)),
            'min_of_empirical_point_means':empirical_min,
            'empirical_mean_of_per_shape_finite_grid_min':gridmean,
            'finite_grid_order_of_min_loss_ratio':empirical_min/gridmean,
            'empirical_original_corner_floor_mean':float(acc['cornerfloor_sum'][0]/n),
            'known_original_floor_ceiling_display':4.286013605856154e-6,
            'condition_quantiles':np.quantile(np.concatenate(conditions),[0,.5,.9,.99,.999,1]).tolist(),
            'max_solve_residual':float(np.concatenate(residuals).max()),
            'condition_bin_counts':acc['condition_counts'].tolist(),
            'condition_bin_sym_contributions':(acc['condition_sym_sums']/n).tolist(),
            'condition_bin_plain_contributions':(acc['condition_plain_sums']/n).tolist(),
            'fallback_total':sum(m['fallback_count'] for m in metas),
            'worker_pids':sorted({m['pid'] for m in metas}),
            'chunk_cpu_seconds_sum':sum(m['cpu_seconds'] for m in metas),
            'chunk_wall_seconds_sum':sum(m['wall_seconds'] for m in metas),
            'not_claimed':['global g or g0 certificate','lower bound from finite-grid per-shape minimum','deterministic integration error','finite variance theorem','confidence guarantee','I or positive S','new C6 bound']}

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--resume',action='store_true'); args=parser.parse_args()
    if (HERE/'state.json').exists() and not args.resume:
        print('Refusing to overwrite an existing run; preserve state and results.',file=sys.stderr)
        return 2
    t0=time.monotonic(); c0=cpu_seconds(); pool=None; metas=[]; code=0; result=None
    try:
        spec=json.loads((HERE/'mc_inputs.json').read_text())
        points=json.loads((HERE/'points.json').read_text())
        release=json.loads((HERE/'primary_mc_release.json').read_text())
        sh=sha(SOURCE); ih=sha(HERE/'mc_inputs.json'); ph=sha(HERE/'points.json')
        assert release['authorized'] is True and release['source_sha256']==sh and release['input_sha256']==ih
        assert release['points_sha256']==ph==spec['points_sha256']
        assert spec['total_shapes']==65536 and spec['chunk_shapes']==2048 and spec['workers']==12
        assert spec['chunk_ids']==list(range(32)) and len(points['points'])==16 and spec['reflections']==64
        Krat=3*math.factorial(32)/(8*4**16*math.factorial(16)*math.factorial(15))
        K=Krat/math.pi**2
        Cgamma=math.gamma(2.5)/(4*math.pi**3)
        Kgamma=2*Cgamma*math.exp(math.lgamma(16.5)-math.lgamma(16))
        assert abs(K/Kgamma-1)<1e-13
        for name in ('chunks','workers','failures'): (HERE/name).mkdir(exist_ok=True)
        if (HERE/'STOP.json').exists(): raise RuntimeError('STOP marker exists; explicit recovery decision required')
        if (HERE/'state.json').exists():
            state=json.loads((HERE/'state.json').read_text())
            assert args.resume and state['source_sha256']==sh and state['input_sha256']==ih
            assert state['status']=='INTERRUPTED'
            assert all(json.loads(path.read_text())['scientific_failure'] is False for path in (HERE/'failures').glob('*.json'))
            for cid in state['completed_chunks']:
                meta=json.loads((HERE/'chunks'/('chunk_%02d.json'%cid)).read_text())
                assert meta['source_sha256']==sh and meta['input_sha256']==ih and meta['points_sha256']==ph
                assert sha(HERE/'chunks'/('chunk_%02d.npz'%cid))==meta['checkpoint_sha256']; metas.append(meta)
        else:
            assert not args.resume
        state={'status':'RUNNING','source_sha256':sh,'input_sha256':ih,'points_sha256':ph,'completed_chunks':[m['chunk_id'] for m in metas],
               'fixed_total_chunks':32,'workers':12,'constant_independent_gamma_vs_rational_relative_error':abs(K/Kgamma-1),
               'K_rational_over_pi_squared':'3*32!/(8*4^16*16!*15!*pi^2)'}
        dump(HERE/'state.json',state)
        todo=[cid for cid in range(32) if cid not in state['completed_chunks']]
        pool=futures.ProcessPoolExecutor(max_workers=12,mp_context=multiprocessing.get_context('spawn'),
                                         initializer=init_worker,initargs=(spec,points['points'],K))
        submitted={pool.submit(process_chunk,cid):cid for cid in todo}
        for completed in futures.as_completed(submitted):
            meta=completed.result(); metas.append(meta)
            state['completed_chunks']=sorted(m['chunk_id'] for m in metas); dump(HERE/'state.json',state)
            if len(metas)%4==0 or len(metas)==1:
                print(json.dumps({'completed_chunks':len(metas),'total_chunks':32,'latest_pid':meta['pid'],'max_condition_latest':meta['max_condition']}),flush=True)
        pool.shutdown(wait=True,cancel_futures=False); pool=None
        result=aggregate(spec,metas,K)
        state['status']='COMPLETE'; dump(HERE/'state.json',state)
    except (KeyboardInterrupt,InterruptedError):
        code=130; result={'status':'INTERRUPTED_NO_PARTIAL_MEAN','complete_chunks':len(metas)}
        dump(HERE/'STOP.json',{'reason':'KeyboardInterrupt','own_process_tree_only':True})
        if 'state' in locals(): state['status']='INTERRUPTED'; dump(HERE/'state.json',state)
    except BaseException as exc:
        code=1; result={'status':'FAILED_NO_PARTIAL_MEAN','exception':repr(exc),'traceback':traceback.format_exc(),'complete_chunks':len(metas)}
        dump(HERE/'STOP.json',{'reason':repr(exc),'own_process_tree_only':True})
        if 'state' in locals(): state['status']='FAILED'; dump(HERE/'state.json',state)
    finally:
        if pool is not None: pool.shutdown(wait=True,cancel_futures=True)
    result.update({'source_sha256':sha(SOURCE),'input_sha256':sha(HERE/'mc_inputs.json'),'points_sha256':sha(HERE/'points.json'),
                   'wall_seconds_including_worker_cleanup':time.monotonic()-t0,'main_cpu_seconds':cpu_seconds()-c0,
                   'children_cpu_seconds_after_cleanup':resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime+resource.getrusage(resource.RUSAGE_CHILDREN).ru_stime,
                   'platform':platform.platform(),'python':sys.version,'exit_code':code,'no_time_ceiling':True,
                   'pool_shutdown_complete':True,'partial_samples_renormalized':False})
    dump(HERE/'results.json',result)
    print(json.dumps({k:v for k,v in result.items() if k not in ('per_point','condition_bin_sym_contributions','condition_bin_plain_contributions')},indent=2),flush=True)
    return code

if __name__=='__main__':
    sys.exit(main())
