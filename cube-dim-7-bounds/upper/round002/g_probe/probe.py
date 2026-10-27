#!/usr/bin/env python3
"""Finite common-shape candidate search. All outputs are numerical, never floors."""
import os
for _k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_k]='1'
import argparse,hashlib,itertools,json,math,platform,sys,time
from pathlib import Path
from datetime import datetime,timezone
from fractions import Fraction
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
import scipy
from scipy.optimize import minimize

ROOT=Path(__file__).resolve().parent
CACHE=Path('/Users/zhengzhihao/Documents/Codex/2026-10-27/n-n-5-1-1-verify/work/cache/gaussian-round001/original-8192-shapes.npz')
OBS=ROOT.parents[1]/'round001/gaussian/observables.npz'
EXPECTED='6dfa48bf8954a7d02270af5bda9c1d33d4fcaf0ba4f5e486589186b7f7c09552'
DATA=None

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def utc(): return datetime.now(timezone.utc).isoformat()
def write(p,x):
    tmp=Path(str(p)+'.tmp'); tmp.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n'); tmp.replace(p)
def xrecord(x): return {'floathex':[float(v).hex() for v in x],'decimal':[float(v) for v in x]}
def load():
    global DATA
    if DATA is None:
        if sha(CACHE)!=EXPECTED: raise RuntimeError('Cache SHA mismatch')
        with np.load(CACHE) as a:
            DATA={k:a[k] for k in a.files}
        with np.load(OBS) as a: DATA['condition']=a['condition']
        DATA['CD']=float(DATA['potential_C']*DATA['radial_D'])
        DATA['high']=np.flatnonzero(DATA['condition']>1e8)
        if DATA['A_inverse'].shape!=(8192,7,7): raise RuntimeError('Shape mismatch')
    return DATA

def mp_shape_diagnostic(x,indices):
    import mpmath as mp
    d=load(); out=[]
    with mp.workdps(80):
        xx=mp.matrix([mp.mpf(float(v)) for v in x])
        vv=[[1]+[1 if (b>>k)&1 else -1 for k in range(6)] for b in range(64)]
        for i in indices:
            pp=[mp.mpf(float(p)) for p in d['P'][i]]
            aa=mp.matrix([[sum(mp.mpf(vv[k][a]*vv[k][b])/pp[k] for k in range(64)) for b in range(7)] for a in range(7)])
            ll=mp.cholesky(aa); det=mp.fprod(ll[a,a]**2 for a in range(7))
            sol=mp.lu_solve(aa,xx); qq=(xx.T*sol)[0]
            val=det**(-mp.mpf(1)/2)*qq**(-3)
            ax=d['A_inverse'][i]@np.asarray(x); qlow=float(np.asarray(x)@ax)
            vlo=float(d['det_A_minus_half'][i]*qlow**(-3)) if qlow>0 else None
            out.append({'shape_index':int(i),'condition':float(d['condition'][i]),'Q_80digit':str(qq),'T_80digit':str(val),'T_float':vlo,'relative_T_difference':float(mp.mpf(vlo)/val-1) if vlo is not None else None})
    return out

def evaluate(x):
    """One full point invocation computes all 8192 values and analytic gradient."""
    d=load(); x=np.asarray(x,dtype=float)
    ax=np.einsum('nij,j->ni',d['A_inverse'],x,optimize=False)
    q=np.einsum('ni,i->n',ax,x,optimize=False)
    bad=np.flatnonzero((q<=0)|~np.isfinite(q))
    if len(bad):
        detail=mp_shape_diagnostic(x,bad.tolist())
        raise ArithmeticError('Nonpositive/nonfinite Q; same-point 80-digit diagnostic: '+json.dumps(detail))
    contributions=d['CD']*d['det_A_minus_half']*q**(-3)
    grad=-6*np.mean((contributions/q)[:,None]*ax,axis=0)
    value=float(np.mean(contributions))
    if not(value>0 and math.isfinite(value) and np.all(np.isfinite(grad))): raise ArithmeticError('Nonfinite value/gradient; no dropping')
    return value,grad,contributions,{'min_Q':float(q.min()),'max_Q':float(q.max()),'max_original_condition':float(d['condition'].max()),'all_shapes_retained':8192}

class LocalLimit(Exception): pass

def run_start(idx):
    reg=json.loads((ROOT/'registration.json').read_text()); init=reg['starts'][idx]
    out=ROOT/'starts'/f'{idx:03d}.json'; out.parent.mkdir(exist_ok=True)
    # Exclusive persistent handle: unknown jobs are never silently re-dispatched.
    fd=os.open(out,os.O_CREAT|os.O_EXCL|os.O_WRONLY); os.close(fd)
    start=time.time(); wall=time.perf_counter(); cpu=time.process_time(); trace=[]; best=None; calls=0
    state={'handle':f'gprobe-r002-start-{idx:03d}','start_id':idx,'pid':os.getpid(),'status':'running','started_utc':utc(),'started_epoch':start,'start_point':init,'full_point_calls':0}
    write(out,state)
    def fun(y):
        nonlocal best,calls
        if calls>=256: raise LocalLimit('256 full-point per-start cap')
        calls+=1
        x=np.r_[1.0,y]
        try: g,gr,con,diag=evaluate(x)
        except BaseException as e:
            trace.append({'call':calls,'x':xrecord(x),'error':repr(e)}); raise
        item={'call':calls,'x':xrecord(x),'Ghat':g,'gradient_free':gr[1:].tolist(),'min_Q':diag['min_Q']}
        trace.append(item)
        if best is None or g<best['Ghat']: best=item
        if len(trace)%16==0:
            state.update(full_point_calls=calls,best=best); write(out,state)
        return math.log(g),gr[1:]/g
    try:
        result=minimize(fun,np.array(init),jac=True,bounds=[(0,1)]*6,method='L-BFGS-B',options={'maxiter':128,'maxfun':256,'maxls':24,'ftol':1e-13,'gtol':1e-8,'maxcor':12})
        state.update(status='finished',solver_success=bool(result.success),solver_message=str(result.message),solver_nit=int(result.nit),solver_nfev=int(result.nfev),solver_njev=int(result.njev),solver_point=xrecord(np.r_[1.0,result.x]))
    except LocalLimit as e: state.update(status='local_budget',solver_success=False,solver_message=str(e))
    except BaseException as e: state.update(status='error',solver_success=False,solver_message=repr(e))
    state.update(full_point_calls=calls,best=best,trace=trace,finished_utc=utc(),finished_epoch=time.time(),wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu)
    write(out,state)
    return {k:v for k,v in state.items() if k!='trace'}

def verify_registration():
    reg=json.loads((ROOT/'registration.json').read_text())
    if sha(CACHE)!=reg['cache_sha256']: raise RuntimeError('Cache lock mismatch')
    for k,h in reg['sources_sha256'].items():
        if sha(ROOT/k)!=h: raise RuntimeError('Source lock mismatch '+k)
    return reg

def run_stage(stage):
    reg=verify_registration(); first=stage=='first8'; ids=list(range(8)) if first else list(range(8,64)); workers=8 if first else 15
    dest=ROOT/f'{stage}-receipt.json'
    if dest.exists(): raise RuntimeError('Stage receipt exists; no re-dispatch')
    if not first:
        prior=json.loads((ROOT/'first8-receipt.json').read_text())
        if prior['status']!='completed': raise RuntimeError('First stage not complete')
    if any((ROOT/'starts'/f'{i:03d}.json').exists() for i in ids): raise RuntimeError('Existing start handles; inspect status, never re-send')
    wall=time.perf_counter(); results=[]; before=utc()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        fut={pool.submit(run_start,i):i for i in ids}
        for f in as_completed(fut):
            item=f.result(); results.append(item)
            print(json.dumps({'start':item['start_id'],'status':item['status'],'calls':item['full_point_calls'],'best_Ghat':item['best']['Ghat'] if item['best'] else None}),flush=True)
    elapsed=time.perf_counter()-wall; cp=sum(x['cpu_seconds'] for x in results)
    events=sorted([(a['started_epoch'],1) for a in results]+[(a['finished_epoch'],-1) for a in results]); active=0; peak=0
    for _,delta in events: active+=delta; peak=max(peak,active)
    result={'stage':stage,'status':'completed' if all(x['status'] in ('finished','local_budget') for x in results) else 'errors','started_utc':before,'finished_utc':utc(),'workers':workers,'unique_worker_pids':sorted({a['pid'] for a in results}),'peak_concurrent_start_handles':peak,'starts':sorted(results,key=lambda a:a['start_id']),'full_point_calls':sum(a['full_point_calls'] for a in results),'value_and_gradient_joint_calls':True,'wall_seconds':elapsed,'sum_worker_cpu_seconds':cp,'configured_lane_cpu_utilization':cp/(elapsed*workers),'frozen_registration_sha256':sha(ROOT/'registration.json'),'no_new_shape_draws':True}
    write(dest,result)

def det3(a):
    return a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])-a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])+a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0])

def classify(x):
    """Exact rational S3 membership, no floating rank threshold."""
    if len({abs(v) for v in x})>=5: return {'classification':'off_S3','reason':'at least five distinct exact absolute values; every 3-body span has at most four','distinct_abs':len({abs(v) for v in x})}
    vv=[[1]+[1 if (b>>k)&1 else -1 for k in range(6)] for b in range(64)]
    checked=0
    for js in itertools.combinations(range(64),3):
        rows=[[vv[j][r] for j in js] for r in range(7)]
        for a,b in itertools.combinations(range(1,7),2):
            mat=[rows[0],rows[a],rows[b]]; det=det3(mat)
            if det: break
        else: raise ArithmeticError('Three distinct hypercube columns unexpectedly dependent')
        # Cramer's rule with exact Fraction inputs.
        xp=[x[0],x[a],x[b]]; coeff=[]
        for col in range(3):
            mt=[row.copy() for row in mat]
            for rr in range(3): mt[rr][col]=xp[rr]
            coeff.append(det3(mt)/det)
        checked+=1
        if all(sum(rows[r][k]*coeff[k] for k in range(3))==x[r] for r in range(7)):
            return {'classification':'in_S3','body_indices':list(js),'coefficients':[str(t) for t in coeff],'triples_checked':checked,'true_G':'infinity, so finite MC value is not a counterexample'}
    return {'classification':'off_S3','reason':'exhaustive exact Cramer check of all 41664 triples','triples_checked':checked}

def summarize(x,counted):
    g,gr,con,diag=counted(x)
    blocks=[float(a.mean()) for a in con.reshape(8,1024)]
    return {'point':xrecord(x),'Ghat':g,'ratio_to_target':g*179000,'gradient_free':gr[1:].tolist(),'blocks':blocks,'block_max_min_ratio':max(blocks)/min(blocks),'sample_sd':float(con.std(ddof=1)),'empirical_SEM':float(con.std(ddof=1)/math.sqrt(8192)),'max_single_sample':float(con.max()),'largest_single_sample_share_of_total':float(con.max()/con.sum()),'conditioning':diag,'same_point_high_condition_80digit':counted(x,highprec=True)}

def diagnose():
    verify_registration()
    if (ROOT/'diagnostics.json').exists(): raise RuntimeError('Diagnostics already exist')
    receipts=[json.loads((ROOT/f'{a}-receipt.json').read_text()) for a in ('first8','rest56')]
    if not all(a['status']=='completed' for a in receipts): raise RuntimeError('Unfinished/error stage')
    total=sum(a['full_point_calls'] for a in receipts); calls=0; fullcalls=0; mpcalls=0; wall=time.perf_counter(); cpu=time.process_time()
    starts=[a for r in receipts for a in r['starts']]
    selected=[]; seen=set()
    for a in sorted(starts,key=lambda a:a['best']['Ghat']):
        orig=[float.fromhex(h) for h in a['best']['x']['floathex']]
        rational=[Fraction(1)]+[Fraction(round(v*1000000),1000000) for v in orig[1:]]
        key=tuple(rational)
        if key in seen: continue
        seen.add(key); selected.append((a,orig,rational))
        if len(selected)==8: break
    def counted(x,highprec=False):
        nonlocal calls,fullcalls,mpcalls
        if total+calls>=32768 or calls>=8192: raise LocalLimit('Global/diagnostic budget')
        calls+=1
        if highprec:
            mpcalls+=1
            return mp_shape_diagnostic(x,load()['high'].tolist())
        fullcalls+=1
        return evaluate(x)
    candidates=[]
    for item,orig,rat in selected:
        base=summarize(orig,counted); xr=np.array([float(z) for z in rat]); rr=summarize(xr,counted)
        rr['rational']=[str(z) for z in rat]; rr['S3']=classify(rat)
        # One fixed rational near-point perturbation bundle, without new starts.
        eps=Fraction(1,1000003)
        near=[Fraction(1)]+[(1-7*eps)*rat[k]+eps*k for k in range(1,7)]
        nr=summarize([float(z) for z in near],counted); nr['rational']=[str(z) for z in near]; nr['S3']=classify(near)
        f0=rr['Ghat']; finite=[]
        for h in (1e-5,3e-6):
            fd=[]; modes=[]
            for k in range(1,7):
                if xr[k]>=h and xr[k]<=1-h:
                    xp=xr.copy(); xm=xr.copy(); xp[k]+=h; xm[k]-=h
                    fp=counted(xp)[0]; fm=counted(xm)[0]; fd.append((fp-fm)/(2*h)); modes.append('central')
                else:
                    sign=1 if xr[k]<h else -1
                    x1=xr.copy(); x2=xr.copy(); x1[k]+=sign*h; x2[k]+=sign*2*h
                    f1=counted(x1)[0]; f2=counted(x2)[0]; fd.append(sign*(-3*f0+4*f1-f2)/(2*h)); modes.append('second_order_one_sided')
            analytic=np.array(rr['gradient_free']); fda=np.array(fd)
            finite.append({'h':h,'gradient':fd,'mode':modes,'max_abs_error':float(np.max(np.abs(fda-analytic))),'error_over_G':float(np.max(np.abs(fda-analytic))/f0)})
        candidates.append({'source_start_handle':item['handle'],'original_optimized_point':base,'rounded_rational_point':rr,'fixed_nearby_offS3_probe':nr,'finite_difference_checks':finite})
        print(json.dumps({'candidate_complete':len(candidates)-1,'source':item['handle'],'rounded_Ghat':rr['Ghat'],'S3':rr['S3']['classification'],'diagnostic_calls_so_far':calls}),flush=True)
    final={'status':'completed','candidate_bundles':candidates,'optimization_calls':total,'diagnostic_full_point_calls':fullcalls,'high_precision_same_point_subset_verification_calls':mpcalls,'diagnostic_budget_charged_calls':calls,'total_full_point_calls':total+fullcalls,'total_budget_charged_calls':total+calls,'total_full_point_budget':32768,'wall_seconds':time.perf_counter()-wall,'cpu_seconds':time.process_time()-cpu,'strict_global_g':False,'strict_global_minimum':False,'new_C7_bound':False,'new_shapes':0,'independent_science_sample_count':8192,'all_results':'numerical candidates from the single inherited common shape sample; no proof and no MC counterexample to a true floor'}
    write(ROOT/'diagnostics.json',final)

def main():
    p=argparse.ArgumentParser(); p.add_argument('action',choices=['first8','rest56','diagnose']); args=p.parse_args()
    if args.action=='diagnose': diagnose()
    else: run_stage(args.action)
if __name__=='__main__': main()
