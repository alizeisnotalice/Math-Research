"""Fixed conditional upper formulas; no source evaluation or random sampling."""
import math
import numpy as np
C_LOG=81/20 # pays the exact 1000001/249001 coefficient, independently of data.

def per_scale_cantelli(logmu,logv,h,factor=1.):
    with np.errstate(over='ignore',under='ignore',divide='ignore',invalid='ignore'):
        t=np.maximum(0,-np.expm1(logmu-h))
        rv=np.exp(logv+math.log(factor)-2*h)
        q=np.where(logmu<h,rv/(rv+t*t),1.)
        q=np.minimum(q,np.minimum(1.,np.exp(logmu-h)))
    return np.nan_to_num(q,nan=1.)

def conditional_uppers(logmu,logv,max_logproduct,table0_logproduct,h,R):
    # logmu/logv: (model,N,J), max/table0 logproduct: (N,J).
    # The 0th model is fixed0; the others are joint and actual prior.
    old=[];variance=[];sharp=[]
    for model in range(3):
        m,v=logmu[model],logv[model]
        old_j=per_scale_cantelli(m,v,h)
        tight_j=per_scale_cantelli(m,v,h,1/3)
        S=C_LOG*(table0_logproduct if model==0 else max_logproduct)
        terms=1 if model==0 else R
        gap=np.maximum(0,h-m)
        exponent=np.divide(-2*gap*gap,S,out=np.full_like(gap,-np.inf),where=S>0)
        with np.errstate(under='ignore'):
            exp_j=np.where(h>m,np.minimum(1.,terms*np.exp(exponent)),1.)
        new_j=np.minimum(tight_j,exp_j)
        # Keeping the old bound is also a safe FP diagnostic and retains its
        # Markov/source cap in the real-arithmetic formula.
        new_j=np.minimum(new_j,old_j)
        old.append(np.minimum(1.,old_j.sum(axis=1)))
        variance.append(np.minimum(1.,tight_j.sum(axis=1)))
        sharp.append(np.minimum(1.,new_j.sum(axis=1)))
    return np.array(old),np.array(variance),np.array(sharp)

def micro_checks():
    # Constant, strict-threshold and a nonzero-gap example.
    m=np.full((3,2,2),math.log(.5));v=np.full_like(m,-math.inf);s=np.zeros((2,2))
    a,b,c=conditional_uppers(m,v,s,s,0.,4)
    assert np.array_equal(a,np.zeros((3,2)))and np.array_equal(c,a)
    m.fill(0.);a,b,c=conditional_uppers(m,v,s,s,0.,4)
    assert np.array_equal(c,np.ones((3,2)))
    m.fill(math.log(.5));v.fill(math.log(.02));s.fill(.04)
    a,b,c=conditional_uppers(m,v,s,s,0.,4)
    assert np.all(c<=b)and np.all(b<=a)
    return {'zero_variance_gap':True,'strict_threshold_boundary':True,'upper_chain':True}
