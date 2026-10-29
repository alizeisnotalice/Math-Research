"""Independent continuous-uniform Lipschitz extension; frozen factor2 untouched."""
import math
import numpy as np
C_LOG=81/20

def uniform_upper(logmu,logv,maxlog,firstlog,h,R):
    answer=[]
    floor_count=0
    for model in range(3):
        m,v=logmu[model],logv[model]
        with np.errstate(over='ignore',under='ignore',divide='ignore',invalid='ignore'):
            t=np.maximum(0,-np.expm1(m-h))
            rv=np.exp(v-math.log(3)-2*h)
            cant=np.where(m<h,rv/(rv+t*t),1.)
            cant=np.minimum(cant,np.minimum(1.,np.exp(m-h)))
        cant=np.nan_to_num(cant,nan=1.)
        S=C_LOG*(firstlog if model==0 else maxlog)
        terms=1 if model==0 else R
        gap=np.maximum(0,h-m)
        exponent=np.divide(-6*gap*gap,S,out=np.full_like(gap,-np.inf),where=S>0)
        floor_count+=int(np.sum((S>0)&(h>m)&(exponent<-700)))
        # Flooring the exponent upward avoids converting a positive upper
        # into zero solely through exp underflow. This is not full FP auditing.
        e=np.minimum(1.,terms*np.exp(np.maximum(exponent,-700)))
        e=np.where(S>0,e,0.)
        e=np.where(h>m,e,1.)
        per=np.minimum(cant,e)
        answer.append(np.minimum(1.,per.sum(axis=1)))
    return np.array(answer),floor_count
