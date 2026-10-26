#!/usr/bin/env python3
"""Eight fixed one-h optimizations, at most20000 scalar objective calls."""
import sys
sys.set_int_max_str_digits(0)
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
import hashlib,json,resource,time
from pathlib import Path
from datetime import datetime,timezone
from fractions import Fraction as Q
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def req(x,s):
    if not x:raise ArithmeticError(s)
def utc():return datetime.now(timezone.utc).isoformat()
def atomic(n,d):
    p=HERE/n;t=p.with_name(p.name+'.own-tmp')
    with t.open('w') as f:
        json.dump(d,f,ensure_ascii=False,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
    os.replace(t,p)
class ObjectiveBudgetExhausted(Exception):pass

def optimize_one(core,caps,job_index):
    import numpy as np
    from scipy.optimize import minimize
    import scipy
    from math import floor
    beta=list(map(Q,caps['beta_upper']));p=sum(map(Q,caps['p_interval']))/2
    probs=np.array([float(v) for v in __import__('scalar_core').pmfs(p)],dtype=float)
    feat=np.array([[float(Q(v)/b) for v,b in zip(row,beta)] for row in core.table],dtype=float)
    xstart=np.array([1.,0.,0.,0.,0.,0.,0.]);floor0=2.**-24
    calls=0;best=None;best_value=float('inf')
    def objective(x):
        nonlocal calls,best,best_value
        if calls>=2500:raise ObjectiveBudgetExhausted()
        calls+=1;w=feat@x
        req(bool(np.all(np.isfinite(w))) and bool(np.all(w>0)),'finite positive float objective')
        value=float(np.dot(probs,w**(-.2)))
        grad=-.2*(feat.T@(probs*w**(-1.2)))
        req(np.isfinite(value) and bool(np.all(np.isfinite(grad))),'finite scalar objective/gradient')
        feasible=bool(np.all(x>=-1e-13)) and x[0]>=floor0/2 and abs(float(np.sum(x))-1)<1e-6
        if feasible and value<best_value:best=np.array(x,copy=True);best_value=value
        if calls%100==0:
            atomic('optimization_checkpoint.json',{'job':job_index,'h':caps['h'],'objective_calls':calls,
                'best_float_x':None if best is None else best.tolist(),'best_float_objective':best_value,
                'not_a_certificate':True})
        return value,grad
    try:
        result=minimize(objective,xstart,jac=True,method='SLSQP',
                        bounds=[(floor0,1.)]+[(0.,1.)]*6,
                        constraints=[{'type':'eq','fun':lambda x:float(np.sum(x)-1.),'jac':lambda x:np.ones(7)}],
                        options={'ftol':1e-12,'maxiter':500,'disp':False})
        outcome={'success':bool(result.success),'status':int(result.status),'message':str(result.message),'iterations':int(result.nit)}
    except ObjectiveBudgetExhausted:
        outcome={'success':False,'status':'REGISTERED_OBJECTIVE_CAP','message':'best candidate retained; no restart'}
    req(best is not None and calls<=2500,'one bounded optimization with a candidate')
    clean=np.maximum(best,0.);clean[0]=max(clean[0],floor0);clean/=float(np.sum(clean))
    den=1<<44;integers=[int(floor(float(v)*den)) for v in clean]
    integers[0]+=den-sum(integers)
    req(min(integers)>=0 and integers[0]>0 and sum(integers)==den,'rationalized exact positive simplex')
    return integers,{'method':'one SLSQP with analytic objective/gradient','objective_calls':calls,
            'initial_x':xstart.tolist(),'best_float_x':best.tolist(),'best_float_objective_not_certificate':best_value,
            'rationalization':'nonnegative weight construction; floor to2^-44 and leftover assigned x0',
            'outcome':outcome,'numpy_version':np.__version__,'scipy_version':scipy.__version__,
            'random_samples':0,'optimizer_restarts':0}

def main():
    start,cpu=time.monotonic(),time.process_time()
    reg=json.loads((HERE/'registration.json').read_text());rel=json.loads((HERE/'parent_release.json').read_text())
    req(rel.get('authorized') is True and str(rel.get('round'))=='124','primary one-batch release')
    req(rel['script_sha256']==sha(__file__),'released producer source')
    req(rel['registration_sha256']==sha(HERE/'registration.json'),'released registration')
    for n,h in reg['prepared_pins'].items():req(sha(HERE/n)==h,'prepared pin '+n)
    sources=json.loads((HERE/'source_manifest.json').read_text())['sources']
    for s in sources:req(sha(s['path'])==s['sha256'],'source '+s['role'])
    with (HERE/'started.json').open('x') as f:
        json.dump({'pid':os.getpid(),'utc':utc(),'release_sha256':sha(HERE/'parent_release.json')},f);f.write('\n')
    done=[];status='FAILED';failure=None;calls=0
    try:
        from scalar_core import Core
        inp=json.loads((HERE/'upstream.json').read_text());core=Core(inp)
        atomic('model_constants.json',core.coefficients())
        for j,hs in enumerate(inp['h_values']):
            h=Q(hs);caps=core.caps(h);weights,opt=optimize_one(core,caps,j)
            calls+=opt['objective_calls'];req(calls<=20000,'total registered objective calls')
            cert=core.certify(caps,weights)
            row={'index':j,'h':hs,'caps':caps,'optimization':opt,'certificate':cert}
            atomic('h_'+str(j)+'.json',row);done.append({'index':j,'h':hs,'path':'h_'+str(j)+'.json','sha256':sha(HERE/('h_'+str(j)+'.json')),'status':cert['status'],'M_upper_expression_display':cert['M_upper_expression_display']})
            atomic('checkpoint.json',{'completed_h':done,'objective_calls':calls})
            print(json.dumps({'h':hs,'status':cert['status'],'objective_calls':opt['objective_calls'],'M_upper':cert['M_upper_expression_display']}),flush=True)
        status='COMPLETE_EIGHT_FIXED_SINGLE_H_CERTIFICATES'
        atomic('results.json',{'status':status,'completed_h':done,'objective_calls_total':calls,
            'threshold_M_power6_exact':str(core.threshold),
            'any_strict_moment_threshold_pass':any(x['status']=='STRICT_MOMENT_THRESHOLD_PASS' for x in done),
            'all_fixed_single_h_tables_obstructed':all(x['status']=='FIXED_SINGLE_H_TABLE_RIGOROUSLY_OBSTRUCTED' for x in done),
            'MC_or_native_calls':0,'actual_M_lower_bound':False,'global_g0_proved':False,'new_C6_bound':False,
            'release_sha256':sha(HERE/'parent_release.json')})
    except BaseException as e:failure=type(e).__name__+': '+str(e);raise
    finally:
        atomic('run_receipt.json',{'status':status,'failure':failure,'completed_h':len(done),'completed_h_objective_calls':calls,
            'pid':os.getpid(),'utc':utc(),'CPU_sample_seconds':time.process_time()-cpu,'wall_sample_before_receipt':time.monotonic()-start,
            'scope':'sample excludes final receipt write/exit; observer covers full child','max_rss_bytes_macos':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'science_CPU_cores':1,'children':0,'wall_limit':None,'CPU_limit':None,'old_program_or_peak_runs':0,'Monte_Carlo_samples':0})
if __name__=='__main__':main()
