#!/usr/bin/env python3
"""One registered 8192-shape run; no automatic extension or resampling."""
import os
for _k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_k]='1'
import hashlib,json,math,platform,resource,sys,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import scipy
from scipy.linalg import solve_triangular
from scipy.special import gammaln

ROOT=Path(__file__).resolve().parent
N=7; M=64; BLOCKS=8; SIZE=1024; SEED=202610277
V=np.array([[1]+[1 if (b>>k)&1 else -1 for k in range(6)] for b in range(64)],dtype=float)
POINTS=np.array([[1]+[(j+k)/(j+7) for k in range(1,7)] for j in range(8)],dtype=float)
C=1/(2*math.pi**3.5)
D=2*math.exp(gammaln(32.5)-gammaln(32))
K=1/(105*2**(25/6)*math.pi**(13/12))

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def moments(a):
    a=np.asarray(a)
    return dict(mean=float(a.mean()),sd=float(a.std(ddof=1)),sem=float(a.std(ddof=1)/math.sqrt(len(a))),min=float(a.min()),max=float(a.max()))

def highprec(y):
    # Same exact binary floating-point input, not a fresh Gamma or shape.
    import mpmath as mp
    with mp.workdps(80):
        yy=[mp.mpf(float(z)) for z in y]
        rr=sum(yy); weights=[rr/z for z in yy]
        aa=mp.matrix([[sum(weights[k]*int(V[k,i])*int(V[k,j]) for k in range(M)) for j in range(N)] for i in range(N)])
        ll=mp.cholesky(aa); det=mp.fprod([ll[i,i]**2 for i in range(N)])
        zz=det**(-mp.mpf(1)/14)
        tt=[]
        for j in range(8):
            xx=mp.matrix([1]+[mp.mpf(j+k)/(j+7) for k in range(1,7)])
            sol=mp.lu_solve(aa,xx); qq=(xx.T*sol)[0]
            tt.append(det**(-mp.mpf(1)/2)*qq**(-3))
        if not (0<zz<=mp.mpf(1)/64) or not all(mp.isfinite(a) and a>0 for a in tt):
            raise ArithmeticError('80-digit same-input check failed')
        return float(zz),np.array([float(a) for a in tt]),str(zz)

def main():
    if (ROOT/'receipt.json').exists() or (ROOT/'results.json').exists():
        raise RuntimeError('Refusing repeat science run: existing results/receipt')
    reg=json.loads((ROOT/'registration.json').read_text())
    for name,digest in reg['frozen_sources_sha256'].items():
        if sha(ROOT/name)!=digest:
            raise RuntimeError('Frozen source mismatch: '+name)
    start_utc=datetime.now(timezone.utc).isoformat()
    wall0=time.perf_counter(); cpu0=time.process_time()
    zall=[]; tall=[]; kall=[]; mineall=[]; maxeall=[]; resall=[]; methods=[]
    blocks=[]; escalations=[]; diagnostics=dict(gamma_zero=0,gamma_nonfinite=0,spd_eigen_failures=0,cholesky_failures=0,output_nonfinite=0,envelope_violations=0,high_precision_count=0,clipped=0,dropped=0,resampled=0)
    stopped=None
    try:
        for b in range(BLOCKS):
            bw=time.perf_counter(); bc=time.process_time()
            rng=np.random.Generator(np.random.PCG64(np.random.SeedSequence([SEED,b])))
            ys=rng.gamma(shape=0.5,scale=1.0,size=(SIZE,M))
            bz=[]; bt=[]
            for j,y in enumerate(ys):
                idx=b*SIZE+j
                zero=int(np.count_nonzero(y<=0)); bad=int(np.count_nonzero(~np.isfinite(y)))
                diagnostics['gamma_zero']+=zero; diagnostics['gamma_nonfinite']+=bad
                if zero or bad:
                    raise ArithmeticError(f'Invalid Gamma at {idx}; no redraw')
                p=y/y.sum(); w=1/p; wmax=w.max(); ws=w/wmax
                aa=(V.T*ws)@V
                ev=np.linalg.eigvalsh(aa); mine,maxe=float(ev[0]),float(ev[-1])
                eigbad=not(np.all(np.isfinite(ev)) and mine>0)
                diagnostics['spd_eigen_failures']+=int(eigbad)
                cond=maxe/mine if mine>0 else math.inf
                reason=[]; resid=math.nan; zz=math.nan; tt=np.full(8,math.nan)
                if eigbad: reason.append('eigen_not_SPD')
                try:
                    ll=np.linalg.cholesky(aa)
                    resid=float(np.linalg.norm(ll@ll.T-aa,ord='fro')/np.linalg.norm(aa,ord='fro'))
                    logdet=2*np.log(np.diag(ll)).sum()+N*math.log(wmax)
                    xx=solve_triangular(ll,POINTS.T,lower=True,check_finite=True)
                    qt=(xx*xx).sum(axis=0)
                    zz=math.exp(-logdet/14)
                    tt=np.exp(-0.5*logdet-3*(np.log(qt)-math.log(wmax)))
                except np.linalg.LinAlgError:
                    diagnostics['cholesky_failures']+=1; reason.append('cholesky_failure')
                nonfinite=not(math.isfinite(zz) and np.all(np.isfinite(tt)) and zz>0 and np.all(tt>0))
                if nonfinite: diagnostics['output_nonfinite']+=1; reason.append('nonfinite_or_nonpositive')
                violate=math.isfinite(zz) and not(0<zz<=1/M)
                if violate: diagnostics['envelope_violations']+=1; reason.append('Z_envelope')
                if cond>1e10: reason.append('condition_above_1e10')
                method=0
                if reason:
                    zlow=zz; tlow=tt.copy()
                    zz,tt,zstr=highprec(y)
                    method=1; diagnostics['high_precision_count']+=1
                    escalations.append(dict(index=idx,reasons=reason,condition=cond,Z_80digits=zstr,Z_float=zlow,max_relative_T_discrepancy=float(np.max(np.abs(tlow/tt-1))) if np.all(np.isfinite(tlow)) else None))
                zall.append(zz); tall.append(tt); kall.append(cond); mineall.append(mine); maxeall.append(maxe); resall.append(resid); methods.append(method)
                bz.append(zz); bt.append(tt)
            blocks.append(dict(block=b,samples=SIZE,seed_sequence=[SEED,b],Z=moments(bz),M_mean=D*float(np.mean(bz)),G_means=(C*D*np.mean(bt,axis=0)).tolist(),wall_seconds=time.perf_counter()-bw,cpu_seconds=time.process_time()-bc))
            print(json.dumps({'block_complete':b,'M_mean':blocks[-1]['M_mean'],'wall_seconds':blocks[-1]['wall_seconds']}),flush=True)
    except BaseException as e:
        stopped=repr(e)
    elapsed=time.perf_counter()-wall0; cpu=time.process_time()-cpu0
    za=np.array(zall); ta=np.array(tall); ga=C*D*ta
    payload=dict(status='completed' if stopped is None else 'stopped',error=stopped,shape_count=len(za),dimension=7,model='equal weight body64 common-time Cauchy',seed=SEED,blocks=blocks,diagnostics=diagnostics,high_precision_events=escalations,constants=dict(C=C,D=D,I_coefficient=K),point_rationals=reg['point_rationals'],points_status='eight fixed face diagnostics; no global g floor')
    if len(za):
        mh=D*float(za.mean()); gh=ga.mean(axis=0); local=float(gh.min()); im=K*mh**(7/6)
        payload.update(Z=moments(za),M_estimate=mh,G_estimates=[dict(point=j,**moments(ga[:,j])) for j in range(8)],sampled_point_min=local,empirical_I_moment_cost=im,empirical_U_using_sampled_point_min=7**(7/6)*im/local**(7/6),empirical_U_ratio_to_128=7**(7/6)*im/local**(7/6)/128,condition=moments(kall),minimum_scaled_eigenvalue=float(min(mineall)),maximum_cholesky_relative_residual=float(max(resall)))
        np.savez_compressed(ROOT/'observables.npz',Z=za,T=ta,condition=np.array(kall),scaled_min_eigen=np.array(mineall),scaled_max_eigen=np.array(maxeall),cholesky_relative_residual=np.array(resall),method=np.array(methods,dtype=np.uint8),points=POINTS)
    (ROOT/'results.json').write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n')
    receipt=dict(start_utc=start_utc,finish_utc=datetime.now(timezone.utc).isoformat(),wall_seconds=elapsed,cpu_seconds=cpu,thread_limits={k:os.environ[k] for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')},processes=1,python=sys.version,python_executable=sys.executable,numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform(),max_rss_native_units=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,status=payload['status'],shape_count=len(za),frozen_sources_sha256=reg['frozen_sources_sha256'],registration_sha256=sha(ROOT/'registration.json'),results_sha256=sha(ROOT/'results.json'),observables_sha256=sha(ROOT/'observables.npz') if len(za) else None)
    (ROOT/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    if stopped: raise RuntimeError(stopped)

if __name__=='__main__': main()
