"""One new deterministic joint objective. No saved single-h objective is run."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from math import comb,lgamma,log,floor
import hashlib

class ObjectiveCap(Exception):pass

def optimize(inp,checkpoint):
    import numpy as np
    import scipy
    from scipy.optimize import minimize
    n=735471;nt=n+16
    beta0=Q(inp['beta0']);beta=[[Q(v) for v in row] for row in inp['beta_by_h']]
    phi=np.asarray(inp['Phi_table'],dtype=np.int64)
    counts=np.empty((n,8),dtype=np.uint8)
    countK=[0]*17
    for i,state in enumerate(combinations_with_replacement(range(17),8)):
        counts[i]=state;countK[state[-1]]+=1
    if i+1!=n or countK!=[comb(k+7,7) for k in range(17)]:raise ArithmeticError('complete new joint state coverage')
    state_hash=hashlib.sha256(counts.tobytes()).hexdigest()
    cats=np.empty_like(counts);cats[:,0]=counts[:,0];cats[:,1:]=counts[:,1:]-counts[:,:-1]
    mids=[float((Q(a)+Q(b))/2) for a,b in inp['p_intervals']]
    q=[mids[0]]+[mids[j]-mids[j-1] for j in range(1,8)]+[1-mids[-1]]
    if min(q)<=0:raise ArithmeticError('valid midpoint multinomial categories')
    lf=np.array([lgamma(j+1) for j in range(33)]);last=32-counts[:,-1].astype(np.int64)
    logprob=np.full(n,lgamma(33))-lf[last]+last*log(q[-1])
    for j in range(8):logprob+=cats[:,j]*log(q[j])-lf[cats[:,j]]
    probs=np.empty(nt);probs[:n]=np.exp(logprob)
    for row,k in enumerate(range(17,33),n):probs[row]=comb(32,k)*mids[-1]**k*(1-mids[-1])**(32-k)
    features=np.zeros((nt,49),dtype=np.float64);features[:,0]=float(1/beta0)
    for j in range(8):
        for k in range(1,7):features[:n,1+6*j+k-1]=phi[counts[:,j],k]/float(beta[j][k-1])
    for row,k in enumerate(range(17,33),n):
        for rank in range(1,7):features[row,1+7*6+rank-1]=float(Q(int(phi[k,rank]))/beta[-1][rank-1])
    mass=float(np.sum(probs))
    if not np.all(np.isfinite(probs)) or not np.all(probs>0) or abs(mass-1)>1e-10:raise ArithmeticError('finite unrenormalized complete float probability')
    warm=inp['warm_start_49'];den=int(warm['denominator']);xstart=np.array(warm['numerators'],dtype=float)/den
    floor0=2.**-24;calls=0;best=None;best_value=float('inf')
    def objective(x):
        nonlocal calls,best,best_value
        if calls>=1000:raise ObjectiveCap()
        calls+=1;w=features@x
        if not np.all(np.isfinite(w)) or not np.all(w>0):raise ArithmeticError('positive finite joint float W')
        f=float(np.dot(probs,w**(-.2)));g=-.2*(features.T@(probs*w**(-1.2)))
        if not np.isfinite(f) or not np.all(np.isfinite(g)):raise ArithmeticError('finite new joint objective/gradient')
        if np.all(x>=-1e-13) and x[0]>=floor0/2 and abs(float(np.sum(x))-1)<1e-6 and f<best_value:
            best=np.array(x,copy=True);best_value=f
        if calls%10==0:checkpoint({'phase':'new_joint_optimization','objective_calls':calls,'best_x_float':None if best is None else best.tolist(),'best_objective_float_not_certificate':best_value})
        return f,g
    try:
        res=minimize(objective,xstart,jac=True,method='SLSQP',bounds=[(floor0,1.)]+[(0.,1.)]*48,
                     constraints=[{'type':'eq','fun':lambda x:float(np.sum(x)-1.),'jac':lambda x:np.ones(49)}],
                     options={'maxiter':500,'ftol':1e-12,'disp':False})
        outcome={'success':bool(res.success),'status':int(res.status),'message':str(res.message),'iterations':int(res.nit)}
    except ObjectiveCap:outcome={'success':False,'status':'REGISTERED_1000_OBJECTIVE_CAP','message':'best candidate kept; no restart'}
    if best is None or calls>1000:raise ArithmeticError('one finite new optimization produced no candidate')
    clean=np.maximum(best,0.);clean[0]=max(clean[0],floor0);clean/=float(np.sum(clean))
    grid=1<<44;nums=[int(floor(float(v)*grid)) for v in clean];nums[0]+=grid-sum(nums)
    if min(nums)<0 or nums[0]<=0 or sum(nums)!=grid:raise ArithmeticError('exact49 weight simplex')
    return {'denominator':str(grid),'numerators':nums}, {
        'objective_calls':calls,'optimizer_invocations':1,'restarts':0,'outcome':outcome,
        'warm_start_source':inp['warm_start_49']['source'],'best_float_x':best.tolist(),
        'best_float_objective_not_certificate':best_value,'prefix_states':n,'tail_rows':16,'features':49,
        'state_bytes_sha256':state_hash,'prefix_count_by_K':countK,'float_probability_mass_unrenormalized':mass,
        'feature_matrix_bytes':int(features.nbytes),'float_feature_matrix_saved':False,
        'numpy_version':np.__version__,'scipy_version':scipy.__version__,'Monte_Carlo_samples':0}
