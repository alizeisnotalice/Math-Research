#!/usr/bin/env python3
"""Authorized deterministic engineering recovery of the original shapes only.

No new independent science sample, new point, or estimator is produced.
"""
import os
for _k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_k]='1'
import hashlib,json,math,runpy,sys,time
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
from scipy.linalg import solve_triangular,cho_solve

ROOT=Path(__file__).resolve().parent

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    reg=json.loads((ROOT/'cache_rebuild-registration.json').read_text())
    cache=Path(reg['cache_directory']); cache.mkdir(parents=True,exist_ok=True)
    dest=cache/'original-8192-shapes.npz'
    if dest.exists() or (ROOT/'cache_rebuild-receipt.json').exists():
        raise RuntimeError('Existing cache recovery; refusing repeat')
    for name,digest in reg['locked_files_sha256'].items():
        if sha(ROOT/name)!=digest: raise RuntimeError('Locked hash mismatch: '+name)
    frozen=runpy.run_path(str(ROOT/'probe.py'))
    V=frozen['V']; POINTS=frozen['POINTS']; D=frozen['D']; C=frozen['C']
    old=np.load(ROOT/'observables.npz')
    assert np.all(old['method']==0), 'Unexpected prior high-precision inputs; stop'
    before=datetime.now(timezone.utc).isoformat(); w0=time.perf_counter(); c0=time.process_time()
    ps=[]; invs=[]; ds=[]; lds=[]; ls=[]; wmaxs=[]; block_receipts=[]
    checked=0
    try:
        for b in range(8):
            rng=np.random.Generator(np.random.PCG64(np.random.SeedSequence([202610277,b])))
            ys=rng.gamma(shape=0.5,scale=1.0,size=(1024,64))
            for j,y in enumerate(ys):
                idx=b*1024+j
                if np.any(y<=0) or not np.all(np.isfinite(y)): raise RuntimeError(f'Gamma failure at {idx}')
                p=y/y.sum(); w=1/p; wmax=w.max(); ws=w/wmax
                aa=(V.T*ws)@V
                ev=np.linalg.eigvalsh(aa); mine,maxe=float(ev[0]),float(ev[-1]); cond=maxe/mine
                ll=np.linalg.cholesky(aa)
                resid=float(np.linalg.norm(ll@ll.T-aa,ord='fro')/np.linalg.norm(aa,ord='fro'))
                logdet=2*np.log(np.diag(ll)).sum()+7*math.log(wmax)
                xx=solve_triangular(ll,POINTS.T,lower=True,check_finite=True)
                qt=(xx*xx).sum(axis=0)
                zz=math.exp(-logdet/14)
                tt=np.exp(-0.5*logdet-3*(np.log(qt)-math.log(wmax)))
                # Exact binary equality is required: not just approximate agreement.
                tests={'Z':zz==old['Z'][idx],'T':np.array_equal(tt,old['T'][idx]),'condition':cond==old['condition'][idx],'min_eigen':mine==old['scaled_min_eigen'][idx],'max_eigen':maxe==old['scaled_max_eigen'][idx],'residual':resid==old['cholesky_relative_residual'][idx]}
                if not all(tests.values()): raise RuntimeError(f'Row {idx} disagrees: {tests}')
                ai=cho_solve((ll,True),np.eye(7),check_finite=True)/wmax
                if not(np.all(np.isfinite(ai)) and np.all(np.diag(ai)>0)): raise RuntimeError(f'Inverse failure at {idx}')
                ps.append(p); invs.append(ai); ds.append(math.exp(-0.5*logdet)); lds.append(logdet); ls.append(ll); wmaxs.append(wmax)
                checked+=1
            block_receipts.append({'block':b,'rows_checked':1024,'all_original_observables_and_diagnostics_bitwise_equal':True})
            print(json.dumps({'cache_block_verified':b,'original_rows_checked':checked}),flush=True)
        np.savez_compressed(dest,P=np.array(ps),A_inverse=np.array(invs),det_A_minus_half=np.array(ds),logdet_A=np.array(lds),scaled_A_cholesky=np.array(ls),A_scale_wmax=np.array(wmaxs),radial_D=np.array(D),potential_C=np.array(C),original_Z=old['Z'],original_T=old['T'],original_points=POINTS)
        receipt={'status':'completed','kind':'engineering cache recovery of exact original sample; not independent science','authorized_exception':'Parent explicitly allowed one deterministic seed replay solely to persist original shape cache, after original-only observable retention was disclosed. No new points or scientific sample count.','start_utc':before,'finish_utc':datetime.now(timezone.utc).isoformat(),'wall_seconds':time.perf_counter()-w0,'cpu_seconds':time.process_time()-c0,'original_rows_checked':checked,'original_observables_and_shape_diagnostics_bitwise_equal':True,'blocks':block_receipts,'new_independent_science_samples':0,'new_points':0,'cache_file':str(dest),'cache_file_bytes':dest.stat().st_size,'cache_sha256':sha(dest),'cache_rebuild_registration_sha256':sha(ROOT/'cache_rebuild-registration.json'),'locked_files_sha256':reg['locked_files_sha256'],'cache_schema':{'P':'original Dirichlet shape, 8192x64','A_inverse':'inverse of A=sum vvT/P; 8192x7x7; Sigma=A/(4R), with R already integrated out','det_A_minus_half':'positive shape determinant weight','logdet_A':'log determinant before exponentiation','scaled_A_cholesky':'Cholesky of A/wmax, reusable stable evaluation','A_scale_wmax':'A scaling factor for stored Cholesky','radial_D':'2 Gamma(65/2)/Gamma(32), floating value; exact formula in proof','potential_C':'1/(2 pi^(7/2)), floating value; exact formula in proof'},'future_point_formula':'Ghat(x)=C*D*mean(det_A_minus_half*(x^T A_inverse x)^(-3)); better numerically use stored scaled Cholesky and wmax; all floating point candidate values only','python_executable':sys.executable}
    except BaseException as e:
        receipt={'status':'stopped','error':repr(e),'original_rows_checked':checked,'new_independent_science_samples':0,'new_points':0}
        (ROOT/'cache_rebuild-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
        raise
    (ROOT/'cache_rebuild-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    (cache/'README.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':'completed','cache':str(dest),'sha256':receipt['cache_sha256']}))

if __name__=='__main__': main()
