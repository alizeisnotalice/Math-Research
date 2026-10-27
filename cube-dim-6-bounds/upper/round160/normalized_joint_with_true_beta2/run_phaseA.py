#!/usr/bin/env python3
"""One separately gated new-cap SLSQP selector. Never a strict M certificate."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations_with_replacement
from math import comb,factorial,log,floor
import hashlib,json,os,sys,time,traceback
sys.dont_write_bytecode=True
from common import P,read,sha,req,new,gated,canon
D=P/'phaseA'
class ObjectiveCap(Exception):pass
def run():
    inp,s,reg,release=gated('phaseA');m=inp['model_input']
    new(D,'started.json',{'round':'160','phase':'phaseA','pid':os.getpid(),'release_sha256':sha(D/'parent_release.json'),'nice':os.getpriority(os.PRIO_PROCESS,0),'time_limit':None})
    start=time.monotonic();cpu=time.process_time();calls=0
    try:
        import numpy as np
        import scipy
        from scipy.optimize import minimize
        req(np.__version__==reg['numpy_version']and scipy.__version__==reg['scipy_version'],'PINNED_NUMPY_SCIPY')
        n=735471;nt=n+16;beta0=F(m['beta0']);beta=[[F(v)for v in row]for row in m['beta_by_h']]
        phi=np.asarray(m['Phi_table'],dtype=np.int64);counts=np.empty((n,8),dtype=np.uint8);countK=[0]*17
        for i,state in enumerate(combinations_with_replacement(range(17),8)):counts[i]=state;countK[state[-1]]+=1
        req(i+1==n and countK==[comb(k+7,7)for k in range(17)],'COMPLETE_PREFIX_STATES')
        state_hash=hashlib.sha256(counts.tobytes()).hexdigest();req(state_hash==inp['expected_state_hash'],'UNCHANGED_COMPLETE_STATE_ORDER')
        cats=np.empty_like(counts);cats[:,0]=counts[:,0];cats[:,1:]=counts[:,1:]-counts[:,:-1]
        mids=[float((F(a)+F(b))/2)for a,b in m['p_intervals']];q=[mids[0]]+[mids[j]-mids[j-1]for j in range(1,8)]+[1-mids[-1]]
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
        mass=float(np.sum(probs));req(np.all(features>=0)and np.all(np.isfinite(features))and np.all(probs>0)and abs(mass-1)<1e-10,'HONEST_COMPLETE_DIAGNOSTIC_LAW')
        warm=inp['warm_start_49'];xstart=np.array(warm['numerators'],dtype=float)/int(warm['denominator'])
        floor0=2.**-24;best=None;best_value=float('inf');final_mode=False
        def objective(x):
            nonlocal calls,best,best_value
            if calls>=(150 if final_mode else 149):raise ObjectiveCap()
            calls+=1;w=features@x;req(np.all(w>0)and np.all(np.isfinite(w)),'POSITIVE_FINITE_W')
            f=float(np.dot(probs,w**(-.2)));g=-.2*(features.T@(probs*w**(-1.2)))
            req(np.isfinite(f)and np.all(np.isfinite(g)),'FINITE_NEW_OBJECTIVE')
            if not final_mode and np.all(x>=-1e-13)and x[0]>=floor0/2 and abs(float(np.sum(x))-1)<1e-6 and f<best_value:best=np.array(x,copy=True);best_value=f
            return f,g
        try:
            res=minimize(objective,xstart,jac=True,method='SLSQP',bounds=[(floor0,1.)]+[(0.,1.)]*48,constraints=[{'type':'eq','fun':lambda x:float(np.sum(x)-1.),'jac':lambda x:np.ones(49)}],options={'maxiter':120,'ftol':1e-12,'disp':False})
            outcome={'success':bool(res.success),'status':int(res.status),'message':str(res.message),'iterations':int(res.nit)}
        except ObjectiveCap:outcome={'success':False,'status':'FINITE_149_CALL_SELECTOR_CAP','message':'one best candidate retained; no restart or optimum claim'}
        req(best is not None,'ONE_OPTIMIZER_PRODUCED_CANDIDATE')
        clean=np.maximum(best,0.);clean[0]=max(clean[0],floor0);clean/=float(np.sum(clean));grid=1<<44
        nums=[int(floor(float(v)*grid))for v in clean];nums[0]+=grid-sum(nums)
        req(min(nums)>=0 and nums[0]>0 and sum(nums)==grid,'EXACT_NEW_SIMPLEX')
        candidate={'denominator':str(grid),'numerators':nums};xexact=np.array(nums,dtype=float)/grid;final_mode=True;f,g=objective(xexact);req(calls<=150,'TOTAL_CALL_CAP')
        theta0=float(F(m['threshold_M_power6_exact']));thetas=float(F(inp['threshold_gstar_power6_exact']));old6=float(F(read(s['128_certificate'])['M_power6_upper']))
        result={'status':'COMPLETE_ONE_GENUINE_CAP_NUMERICAL_SELECTOR','qualification':'true upper cap normalized model; float H values are selectors, not strict M upper certificates or global exclusions','outcome':outcome,'objective_gradient_calls':calls,'optimizer_invocations':1,'restarts':0,'candidate':candidate,'candidate_canonical_sha256':hashlib.sha256(canon(candidate)).hexdigest(),'new_beta2_cap_exact':inp['new_beta2_cap_exact'],'raw_a_quantized':False,'cost_cap_exact':'1',
          'H_float_at_exact_candidate':f,'H_power5_float':f**5,'M_expression_float_diagnostic_only':f**(5/6),'H_power5_below_g0_float_diagnostic_only':f**5<theta0,'H_power5_below_gstar_float_diagnostic_only':f**5<thetas,'H_power5_below_saved128_float_diagnostic_only':f**5<old6,'Euler_homogeneity_residual_float':float(np.dot(g,xexact)+f/5),
          'prefix_states':n,'tail_rows':16,'features':49,'state_bytes_sha256':state_hash,'probability_mass_unrenormalized_float':mass,'feature_matrix_bytes':int(features.nbytes),'new_count_probability_float_exp_calls':n,'new_probability_logfactorial_calls':33,'new_category_probability_float_log_calls':9,'new_kernel_erfc_heatkernel_native_MC_density_calls':0,'old_functional_replays':0,'certified_global_exclusion':False,'new_M_upper':False,'global_g_proved':False,'new_C6_bound':False,'run_sha256':sha(Path(__file__)),'release_sha256':sha(D/'parent_release.json'),'wall_before_receipt_seconds':time.monotonic()-start,'CPU_before_receipt_seconds':time.process_time()-cpu}
        new(D,'result.json',result);new(D,'exit_receipt.json',{'status':result['status'],'result_sha256':sha(D/'result.json'),'objective_gradient_calls':calls,'optimizer_invocations':1,'new_kernel_erfc_native_density_calls':0,'new_M_upper':False,'new_C6_bound':False})
        print(json.dumps({k:result[k]for k in ['status','objective_gradient_calls','H_power5_float','H_power5_below_saved128_float_diagnostic_only','H_power5_below_g0_float_diagnostic_only','H_power5_below_gstar_float_diagnostic_only','new_M_upper']}))
    except BaseException:
        new(D,'exit_receipt.json',{'status':'FAILED_STOP_NO_RETRY','failure':traceback.format_exc(),'objective_gradient_calls':calls,'new_M_upper':False,'new_C6_bound':False});raise
if __name__=='__main__':
    req(not sys.argv[1:],'usage: run_phaseA.py (after phaseA dual release only)');run()
