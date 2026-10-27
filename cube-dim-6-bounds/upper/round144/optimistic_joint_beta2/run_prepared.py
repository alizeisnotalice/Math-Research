#!/usr/bin/env python3
"""One gated optimistic joint optimization; never a strict M certificate."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations_with_replacement
from math import comb,factorial,log,floor
import hashlib,json,os,sys,time,traceback
sys.dont_write_bytecode=True;sys.set_int_max_str_digits(0)
P=Path(__file__).resolve().parent
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def req(ok,msg):
    if not ok:raise ArithmeticError(msg)
def write_new(n,x):
    with(P/n).open('x')as f:json.dump(x,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
def release_requirements():
    return {'round':'144','authorized':True,'scope':'ONE_NEW_OPTIMISTIC_BETA2_COMMON_PARAMETER_FLOAT_DIAGNOSTIC','max_objective_gradient_calls':150,'optimizer_invocations':1,'workers':1,'nice_minimum':10,'time_or_CPU_limit':None,'kernel_erfc_native_calls_permitted':0,'strict_M_or_global_exclusion_claim_permitted':False,'run_sha256':sha(P/'run_prepared.py'),'registration_sha256':sha(P/'registration.json'),'inputs_sha256':sha(P/'optimistic_inputs.json'),'source_manifest_sha256':sha(P/'source_manifest.json'),'proof_sha256':sha(P/'PROOF.md')}
def preflight():
    reg=read(P/'registration.json');inp=read(P/'optimistic_inputs.json');manifest=read(P/'source_manifest.json')
    req(sha(P/'run_prepared.py')==reg['run_sha256']and sha(P/'PROOF.md')==reg['proof_sha256'],'FIXED_CODE_PROOF')
    req(sha(P/'optimistic_inputs.json')==reg['inputs_sha256']and sha(P/'source_manifest.json')==reg['source_manifest_sha256'],'FIXED_PREPARED_INPUTS')
    for row in manifest['sources']:req(sha(row['path'])==row['sha256'],'SOURCE_HASH:'+row['path'])
    source={Path(r['path']).name:Path(r['path'])for r in manifest['sources']}
    old=read(source['upstream.json']);true130=read(source['results.json'])
    req(inp['h_values']==old['h_values']and inp['p_intervals']==old['p_intervals']and inp['Phi_table']==old['Phi_table']and inp['beta0']==old['beta0']and F(inp['beta0'])>0,'UNCHANGED_LAW_AND_PHI')
    req(inp['threshold_M_power6_exact']==old['threshold_M_power6_exact']and inp['prefix_Kmax']==old['prefix_Kmax']and inp['tail_K_values']==old['tail_K_values'],'UNCHANGED_COMPLETE_TARGET_AND_TAIL')
    for j,row in enumerate(inp['beta_by_h']):
        for k,v in enumerate(row):
            req(F(v)>0,'POSITIVE_BETA')
            if(j,k)==(5,1):req(F(v)==F(true130['moment_interval'][0])<F(old['beta_by_h'][j][k]),'SOLE_TRUE_LOWER_REPLACEMENT')
            else:req(v==old['beta_by_h'][j][k],'NO_OTHER_CAP_CHANGE')
    req(all(isinstance(v,int)and v>=0 for row in inp['Phi_table']for v in row),'WHOLE_PARAMETER_DOMAIN_NONNEGATIVE_FEATURES')
    req(all(inp['Phi_table'][n+1][k]>=inp['Phi_table'][n][k]for n in range(32)for k in range(7)),'MONOTONE_COUNT_FEATURES')
    req(inp['warm_start_49']==read(source['candidate.json']),'CURRENT128_EXACT_WARM_START')
    return inp,reg,source
class ObjectiveCap(Exception):pass
def run():
    release=read(P/'parent_release.json')
    for k,v in release_requirements().items():req(release.get(k)==v,'RELEASE:'+k)
    req(release.get('main_approved')is True,'MAIN_RELEASE_REQUIRED');req(not(P/'started.json').exists(),'NO_AUTOMATIC_RESTART')
    os.nice(max(0,10-os.getpriority(os.PRIO_PROCESS,0)))
    req(all(os.environ.get(k)=='1'for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']),'ONE_THREAD_ENV')
    inp,reg,source=preflight()
    write_new('started.json',{'round':'144','pid':os.getpid(),'nice':os.getpriority(os.PRIO_PROCESS,0),'release_sha256':sha(P/'parent_release.json'),'time_limit':None})
    start=time.monotonic();cpu=time.process_time();calls=0;result=None
    try:
        import numpy as np
        import scipy
        from scipy.optimize import minimize
        req(np.__version__==reg['numpy_version']and scipy.__version__==reg['scipy_version'],'PINNED_NUMPY_SCIPY')
        n=735471;nt=n+16;beta0=F(inp['beta0']);beta=[[F(v)for v in row]for row in inp['beta_by_h']]
        phi=np.asarray(inp['Phi_table'],dtype=np.int64);counts=np.empty((n,8),dtype=np.uint8);countK=[0]*17
        for i,state in enumerate(combinations_with_replacement(range(17),8)):counts[i]=state;countK[state[-1]]+=1
        req(i+1==n and countK==[comb(k+7,7)for k in range(17)],'COMPLETE_PREFIX_STATES')
        state_hash=hashlib.sha256(counts.tobytes()).hexdigest()
        req(state_hash==read(source['candidate_certificate.json'])['state_bytes_sha256'],'REUSED_SYMBOLIC_STATE_ORDER')
        cats=np.empty_like(counts);cats[:,0]=counts[:,0];cats[:,1:]=counts[:,1:]-counts[:,:-1]
        mids=[float((F(a)+F(b))/2)for a,b in inp['p_intervals']];q=[mids[0]]+[mids[j]-mids[j-1]for j in range(1,8)]+[1-mids[-1]]
        req(min(q)>0,'POSITIVE_WHOLE_MIDPOINT_LAW')
        lf=np.array([log(factorial(j))for j in range(33)]);last=32-counts[:,-1].astype(np.int64)
        logprob=np.full(n,lf[32])-lf[last]+last*log(q[-1])
        for j in range(8):logprob+=cats[:,j]*log(q[j])-lf[cats[:,j]]
        probs=np.empty(nt);probs[:n]=np.exp(logprob)
        for row,k in enumerate(range(17,33),n):probs[row]=comb(32,k)*mids[-1]**k*(1-mids[-1])**(32-k)
        features=np.zeros((nt,49),dtype=np.float64);features[:,0]=float(1/beta0)
        for j in range(8):
            for k in range(1,7):features[:n,1+6*j+k-1]=phi[counts[:,j],k]/float(beta[j][k-1])
        for row,k in enumerate(range(17,33),n):
            for rank in range(1,7):features[row,1+7*6+rank-1]=float(F(int(phi[k,rank]))/beta[-1][rank-1])
        mass=float(np.sum(probs));req(np.all(features>=0)and np.all(np.isfinite(features))and np.all(probs>0)and abs(mass-1)<1e-10,'HONEST_NONNEGATIVE_DIAGNOSTIC_MODEL')
        warm=inp['warm_start_49'];grid=int(warm['denominator']);xstart=np.array(warm['numerators'],dtype=float)/grid
        floor0=2.**-24;best=None;best_value=float('inf');final_mode=False
        def objective(x):
            nonlocal calls,best,best_value
            if calls>=(150 if final_mode else 149):raise ObjectiveCap()
            calls+=1;w=features@x;req(np.all(w>0)and np.all(np.isfinite(w)),'POSITIVE_FINITE_W')
            f=float(np.dot(probs,w**(-.2)));g=-.2*(features.T@(probs*w**(-1.2)))
            req(np.isfinite(f)and np.all(np.isfinite(g)),'FINITE_DIAGNOSTIC_OBJECTIVE')
            if not final_mode and np.all(x>=-1e-13)and x[0]>=floor0/2 and abs(float(np.sum(x))-1)<1e-6 and f<best_value:best=np.array(x,copy=True);best_value=f
            return f,g
        try:
            res=minimize(objective,xstart,jac=True,method='SLSQP',bounds=[(floor0,1.)]+[(0.,1.)]*48,constraints=[{'type':'eq','fun':lambda x:float(np.sum(x)-1.),'jac':lambda x:np.ones(49)}],options={'maxiter':120,'ftol':1e-12,'disp':False})
            outcome={'success':bool(res.success),'status':int(res.status),'message':str(res.message),'iterations':int(res.nit)}
        except ObjectiveCap:outcome={'success':False,'status':'FINITE_149_CALL_SELECTOR_CAP','message':'best finite candidate kept; no restart'}
        req(best is not None,'ONE_OPTIMIZER_PRODUCED_CANDIDATE')
        clean=np.maximum(best,0.);clean[0]=max(clean[0],floor0);clean/=float(np.sum(clean));grid=1<<44
        nums=[int(floor(float(v)*grid))for v in clean];nums[0]+=grid-sum(nums)
        req(min(nums)>=0 and nums[0]>0 and sum(nums)==grid,'EXACT_NEW_SIMPLEX')
        xexact=np.array(nums,dtype=float)/grid;final_mode=True;f,g=objective(xexact);req(calls<=150,'TOTAL_CALL_CAP')
        theta=float(F(inp['threshold_M_power6_exact']));old6=float(F(read(source['candidate_certificate.json'])['M_power6_upper']))
        tangent=float(f-np.dot(g,xexact)+np.min(g))
        result={'status':'COMPLETE_ONE_OPTIMISTIC_NUMERICAL_DIAGNOSTIC','qualification':'lower beta2 endpoint is optimistic, never a strict M upper; float values/tangent are not certified; local failure excludes no family',
          'outcome':outcome,'objective_gradient_calls':calls,'optimizer_invocations':1,'restarts':0,'candidate':{'denominator':str(grid),'numerators':nums},
          'optimistic_beta2_lower':inp['optimistic_beta2_lower'],'H_float_at_exact_candidate':f,'H_power5_float':f**5,
          'M_expression_float_diagnostic_only':f**(5/6),'H_power5_below_g0_float_diagnostic_only':f**5<theta,'H_power5_below_saved128_float_diagnostic_only':f**5<old6,
          'float_tangent_global_lower_hint_not_certificate':tangent,'float_tangent_power5_if_positive_not_certificate':max(0.,tangent)**5,'Euler_homogeneity_residual_float':float(np.dot(g,xexact)+f/5),
          'prefix_states':n,'tail_rows':16,'features':49,'state_bytes_sha256':state_hash,'probability_mass_unrenormalized_float':mass,'feature_matrix_bytes':int(features.nbytes),
          'new_count_probability_float_exp_calls':n,'new_probability_logfactorial_calls':33,'new_category_probability_float_log_calls':9,
          'new_kernel_erfc_heatkernel_native_MC_calls':0,'old_functional_replays':0,'certified_global_exclusion':False,'new_M_upper':False,'global_g0_proved':False,'new_C6_bound':False,
          'run_sha256':sha(Path(__file__)),'release_sha256':sha(P/'parent_release.json'),'wall_before_receipt_seconds':time.monotonic()-start,'CPU_before_receipt_seconds':time.process_time()-cpu}
        write_new('result.json',result);write_new('exit_receipt.json',{'status':result['status'],'result_sha256':sha(P/'result.json'),'objective_gradient_calls':calls,'optimizer_invocations':1,'new_kernel_erfc_native_calls':0,'new_M_upper':False,'new_C6_bound':False})
        print(json.dumps({k:result[k]for k in ['status','objective_gradient_calls','H_power5_float','H_power5_below_g0_float_diagnostic_only','float_tangent_global_lower_hint_not_certificate','new_M_upper']}))
    except BaseException:
        write_new('exit_receipt.json',{'status':'FAILED_STOP_NO_RETRY','failure':traceback.format_exc(),'objective_gradient_calls':calls,'optimizer_invocations_at_most':1,'new_kernel_erfc_native_calls':0,'new_M_upper':False,'new_C6_bound':False});raise
if __name__=='__main__':
    if sys.argv[1:]==['--saved-only-preflight']:
        inp,reg,_=preflight();print(json.dumps({'status':'SAVED_ONLY_PREFLIGHT_PASS','nonnegative_Phi_whole_domain':True,'positive_beta2_optimistic_direction_checked':True,'new_objective_calls':0,'new_kernel_erfc_calls':0}))
    elif not sys.argv[1:]:run()
    else:raise SystemExit('usage: run_prepared.py [--saved-only-preflight]')
