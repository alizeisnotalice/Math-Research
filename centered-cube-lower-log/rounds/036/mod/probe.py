"""True mod-cell source and exact residue convolution DP."""
from pathlib import Path
from fractions import Fraction as F
import numpy as np
import itertools
import math
import time
import json
import argparse
import hashlib
import sys

OUT=Path(__file__).resolve().parent
PERM=np.array([11,6,4,5,1,7,2,3,9,0,8,10])
CASES=[(n,r)for n in [3,7,15]for r in [2,3,4,6]]


class Source:
    def __init__(self,n,r):
        self.n,self.r,self.k=n,r,12
        assert self.k%r==0
        self.s=PERM%r
        j=np.arange(self.k)
        self.atoms=np.array([(j+((self.s+b)%r+.5)/r)/self.k
                             for b in range(r)])
        self.delta=1/(64*self.k*r)
        assert self.atoms.min()-self.delta/2>0
        assert self.atoms.max()+self.delta/2<1
        short=np.geomspace(self.delta/16,1/(2*self.k),13)
        long=np.linspace(1/(2*self.k),1.25,41)
        extra=[1/self.k,2/self.k,3/self.k,4/self.k,6/self.k,
               1.,10/9,7/6,11/9]
        self.lengths=np.unique(np.concatenate([short,long,extra]))
        self.cyclic=(np.arange(r)[:,None]-np.arange(r)[None,:])%r
        self.models=["joint"]+[f"fixed{b}"for b in range(r)]+["prior"]
        self.density_max=(64*r)**n


def residue_weights(src,x,length):
    if length==1.:
        return np.full((len(x),src.n,src.r,src.r),1/src.r)
    distance=x[:,:,None,None]-src.atoms[None,None,:,:]
    mass=np.zeros_like(distance)
    for shift in [-1,0,1]:
        d=distance-shift
        mass+=np.clip((d+length/2+src.delta/2)/src.delta,0,1) \
             -np.clip((d-length/2+src.delta/2)/src.delta,0,1)
    return (mass/(src.k*length)).reshape(
        len(x),src.n,src.r,src.k//src.r,src.r).sum(axis=3)


def averages_dp(src,x,length,return_residue=False):
    if length==1. and not return_residue:
        return np.ones((len(x),src.r+2))
    z=residue_weights(src,x,length)
    coeff=np.zeros((len(x),src.r,src.r))
    coeff[:,:,0]=1
    for i in range(src.n):
        coeff=np.sum(coeff[:,:,src.cyclic]*z[:,i,:,None,:],axis=-1)
    fixed=np.prod(z.sum(axis=3),axis=1)
    assert np.allclose(coeff.sum(axis=2),fixed,rtol=1e-10,atol=1e-12)
    joint=np.diagonal(coeff,axis1=1,axis2=2).sum(axis=1)
    result=np.column_stack([joint,fixed,fixed.mean(axis=1)])
    return (result,z)if return_residue else result


def direct_average(src,x,length):
    indices=np.indices((src.k,)*src.n,dtype=np.int16).reshape(src.n,-1).T
    classifier=(indices%src.r).sum(axis=1)%src.r
    assert np.all(np.bincount(classifier,minlength=src.r)==len(indices)//src.r)
    fixed=[];joint=0.
    for b in range(src.r):
        d=x[:,None]-src.atoms[b][None,:]
        mass=np.zeros_like(d)
        for shift in [-1,0,1]:
            mass+=np.clip((d-shift+length/2+src.delta/2)/src.delta,0,1) \
                -np.clip((d-shift-length/2+src.delta/2)/src.delta,0,1)
        w=mass/(src.k*length)
        product=np.ones(len(indices))
        for i in range(src.n):product*=w[i,indices[:,i]]
        fixed.append(product.sum())
        joint+=product[classifier==b].sum()
    return np.array([joint]+fixed+[np.mean(fixed)])


def small_control(r):
    src=Source(3,r)
    centers=np.array([[1/7,2/7,4/7],[.11,.47,.82]])
    cell=np.arange(3)%src.k
    b=int((cell%r).sum()%r)
    centers=np.vstack([centers,src.atoms[b,cell]])
    rows=[]
    for length in [src.delta,1/src.k,3/8,1.,7/6]:
        dp=averages_dp(src,centers,length)
        for i,x in enumerate(centers):
            direct=direct_average(src,x,length)
            error=float(np.max(np.abs(dp[i]-direct)/np.maximum(1,np.abs(direct))))
            assert error<1e-10
            rows.append({"L":length,"center":i,"max_relative_error":error,
                         "joint_DP":float(dp[i,0]),"joint_direct":float(direct[0])})
    return {"n":3,"K":12,"R":r,"cells":12**3,"prior_exact":str(F(1,r)),
            "max_relative_error":max(x["max_relative_error"]for x in rows),
            "rows":rows}


def fraction_control(r):
    k,n=12,2;delta=F(1,64*k*r)
    atoms=[[F(2*r*j+2*((int(PERM[j])+b)%r)+1,2*k*r)
            for j in range(k)]for b in range(r)]
    x=[F(1,7),F(4,7)]
    clip=lambda a:max(F(0),min(F(1),a))
    output=[]
    for length in [F(3,8),F(7,6),F(1)]:
        weights=[]
        for b in range(r):
            weights.append([[sum((clip((xx-a-s+length/2+delta/2)/delta)
                -clip((xx-a-s-length/2+delta/2)/delta)
                for s in [-1,0,1]),F(0))/(k*length)
                for a in atoms[b]]for xx in x])
        direct_joint=F(0);fixed=[F(0)]*r
        for cells in itertools.product(range(k),repeat=n):
            b=sum(j%r for j in cells)%r
            direct_joint+=math.prod(weights[b][i][j]for i,j in enumerate(cells))
            for table in range(r):
                fixed[table]+=math.prod(weights[table][i][j]
                                         for i,j in enumerate(cells))
        dp_joint=F(0);dp_fixed=[]
        for b in range(r):
            coeff=[F(0)]*r;coeff[0]=F(1)
            for i in range(n):
                z=[sum((weights[b][i][j]for j in range(k)if j%r==q),F(0))
                   for q in range(r)]
                nxt=[F(0)]*r
                for q in range(r):
                    for t in range(r):nxt[(q+t)%r]+=coeff[q]*z[t]
                coeff=nxt
            dp_joint+=coeff[b];dp_fixed.append(sum(coeff,F(0)))
        assert direct_joint==dp_joint and fixed==dp_fixed
        output.append({"L":str(length),"joint":str(dp_joint),
                       "prior":str(sum(fixed,F(0))/r),"exact_equality":True})
    return {"n":n,"K":k,"R":r,"x":[str(a)for a in x],"rows":output}


def draw_proposal(src,rng,samples):
    x=rng.random((samples,src.n))
    tilted=rng.random(samples)>=.5;rows=np.flatnonzero(tilted)
    model=rng.integers(src.r+2,size=len(rows))
    li=rng.integers(len(src.lengths),size=len(rows))
    cells=rng.integers(src.k,size=(len(rows),src.n))
    b=np.zeros(len(rows),dtype=int)
    b[model==0]=(cells[model==0]%src.r).sum(axis=1)%src.r
    for table in range(src.r):b[model==table+1]=table
    b[model==src.r+1]=rng.integers(src.r,size=np.count_nonzero(model==src.r+1))
    t=src.atoms[b[:,None],cells]+src.delta*(rng.random(cells.shape)-.5)
    x[rows]=(t-src.lengths[li,None]*(rng.random(cells.shape)-.5))%1
    return x,{"uniform":int((~tilted).sum()),
              "tilted_model_counts":np.bincount(model,minlength=src.r+2).tolist(),
              "tilted_source_b_counts":np.bincount(b,minlength=src.r).tolist()}


def evaluate_queries(src,x):
    maxima=np.zeros((len(x),src.r+2));density=np.zeros(len(x))
    for length in src.lengths:
        values=averages_dp(src,x,float(length))
        assert values.min()>=-1e-12
        maxima=np.maximum(maxima,values);density+=values.sum(axis=1)
    weights=1/(.5+.5*density/((src.r+2)*len(src.lengths)))
    assert weights.max()<=2 and maxima.min()>=1
    cap=2*(src.r+2)*len(src.lengths)
    assert (weights[:,None]*maxima).max()<=cap*(1+1e-8)
    # Structural domination is a pointwise check, not a gain theorem.
    assert np.max(maxima[:,0]-src.r*maxima[:,-1])<=1e-8*max(1,maxima.max())
    return maxima,weights


def fourier_diagnostics():
    output=[]
    for n,r in CASES:
        src=Source(n,r)
        seed=202610291000+10*n+r
        x=np.random.default_rng(seed).random((64,n))
        for length in [10/9,7/6,11/9]:
            values,z=averages_dp(src,x,length,True)
            hats=np.fft.fft(z,axis=3)
            totals=hats[:,:,:,0].real
            normalized=hats/totals[:,:,:,None]
            products=np.prod(normalized,axis=1)
            phase=np.exp(2j*np.pi*np.arange(r)[:,None]*np.arange(r)[None,:]/r)
            nonzero_sum=np.sum(products[:,:,1:]*phase[None,:,1:],axis=2)
            fixed=np.prod(totals,axis=1)
            relative_complex=np.sum(fixed*nonzero_sum,axis=1)/fixed.sum(axis=1)
            product_bound=np.sum(fixed*np.sum(np.abs(products[:,:,1:]),axis=2),
                                 axis=1)/fixed.sum(axis=1)
            actual=values[:,0]/values[:,-1]-1
            error=np.max(np.abs(actual-relative_complex.real))
            assert error<1e-10
            output.append({"n":n,"R":r,"K":12,"seed":seed,"centers":64,
                "L":length,"max_abs_joint_over_prior_minus1":float(np.abs(actual).max()),
                "max_abs_Fourier_relative_residual":float(np.abs(relative_complex).max()),
                "max_weighted_nonzero_character_product_bound":
                    float(product_bound.max()),
                "max_abs_nonzero_normalized_character":
                    float(np.abs(normalized[:,:,:,1:]).max()),
                "max_DP_vs_Fourier_difference":float(error),
                "max_imaginary_residual":float(np.abs(relative_complex.imag).max())})
    return output


def calibrate():
    started=time.time()
    controls=[small_control(r)for r in [2,3,4,6]]
    rational=[fraction_control(r)for r in [2,3,4,6]]
    src=Source(15,6)
    x=(.12345+np.arange(128*15)*((math.sqrt(5)-1)/2))%1
    x=x.reshape(128,15)
    t=time.time();evaluate_queries(src,x);bench=time.time()-t
    chosen=None
    for nmc in [4096,2048]:
        estimate=bench*(nmc/128)*len(CASES)*1.3
        if estimate<=180:
            chosen=nmc;break
    calibration={"direct_controls":controls,"Fraction_controls":rational,
                 "worst_batch_seconds":bench,"batch":128,
                 "selected_N":chosen,"estimated_seconds":estimate,
                 "total_seconds":time.time()-started}
    (OUT/"calibration.json").write_text(json.dumps(calibration,separators=(",",":")))
    print(json.dumps({k:calibration[k]for k in [
        "worst_batch_seconds","selected_N","estimated_seconds","total_seconds"]}),
        flush=True)
    if chosen is None:raise RuntimeError("Frozen cost gate exceeds 180 seconds.")
    configs=[]
    count_models=sum(r+2 for n,r in CASES)
    for i,(n,r)in enumerate(CASES):
        s=Source(n,r)
        configs.append({"n":n,"K":12,"R":r,"seed":202610290361+i,"N":chosen,
                        "delta":s.delta,"s_labels":s.s.tolist(),
                        "lengths":s.lengths.tolist(),"J":len(s.lengths)})
    reg={"status":"frozen after exact checks/cost calibration, before MC",
         "script_SHA256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         "plan_SHA256":hashlib.sha256((OUT/"PLAN.md").read_bytes()).hexdigest(),
         "calibration_SHA256":hashlib.sha256((OUT/"calibration.json").read_bytes()).hexdigest(),
         "permutation_seed":202610290360,"permutation":PERM.tolist(),
         "cases":configs,"batch":128,"alpha_each":.05/(2*count_models),
         "classifier":"sum_i (j_i mod R) mod R, source-cell determined",
         "prior_exact":"1/R, since K is multiple of R",
         "proposal":".5 uniform + .5 all sources/all target lengths",
         "weak_screen":"empirical scans all observations; ESS>=64 primary; "
                       "discovery ESS>=32 then nextafter threshold downward, heldout",
         "CI_scope":"only heldout selected finite-tail coefficient; no fullweak ranking",
         "cost_limit_seconds":180,"RAM_limit_bytes":2_000_000_000,
         "stopping":"one fixed seed/case; no gain-adaptive new R/N/source"}
    if (OUT/"registration.json").exists():raise RuntimeError("Already registered.")
    (OUT/"registration.json").write_text(json.dumps(reg,separators=(",",":")))


def weak_scan(values,weights,ess_min):
    order=np.argsort(values)[::-1];v=values[order];w=weights[order]
    sums=np.cumsum(w);square=np.cumsum(w*w);ess=sums*sums/square
    score=v*sums/len(v)
    ix=int(np.argmax(np.where(ess>=ess_min,score,-np.inf)))
    return {"coefficient":float(score[ix]),"threshold":float(v[ix]),
            "ESS":float(ess[ix]),"min_ESS":ess_min}


def kl(x,mu):
    return (x*math.log(x/mu)if x else 0) \
        +((1-x)*math.log((1-x)/(1-mu))if x<1 else 0)


def confidence(mean,bound,n,alpha,upper):
    x=min(1.,max(0.,mean/bound))
    if x==1 and upper:return bound
    if x==0 and not upper:return 0.
    lo,hi=(x,1.)if upper else(0.,x)
    target=math.log(1/alpha)/n
    for _ in range(80):
        mid=(lo+hi)/2
        if upper:
            if kl(x,mid)>target:hi=mid
            else:lo=mid
        else:
            if kl(x,mid)>target:lo=mid
            else:hi=mid
    return bound*(hi if upper else lo)


def run_case(config,alpha):
    started=time.time();src=Source(config["n"],config["R"])
    assert src.lengths.tolist()==config["lengths"]
    rng=np.random.default_rng(config["seed"])
    mm=[];ww=[];uniform=0;models=np.zeros(src.r+2,dtype=int)
    source_b=np.zeros(src.r,dtype=int)
    for left in range(0,config["N"],128):
        count=min(128,config["N"]-left)
        x,audit=draw_proposal(src,rng,count)
        m,w=evaluate_queries(src,x);mm.append(m);ww.append(w)
        uniform+=audit["uniform"];models+=audit["tilted_model_counts"]
        source_b+=audit["tilted_source_b_counts"]
    m=np.concatenate(mm);w=np.concatenate(ww);cut=len(w)//2
    results={}
    for mi,name in enumerate(src.models):
        discovery=weak_scan(m[:cut,mi],w[:cut],32)
        threshold=float(np.nextafter(discovery["threshold"],0.))
        z=threshold*w[cut:]*(m[cut:,mi]>threshold)
        bound=min(2*threshold,2*(src.r+2)*len(src.lengths))
        assert z.max()<=bound*(1+1e-8)
        squares=float(np.dot(z,z));total=float(z.sum())
        results[name]={"empirical_unrestricted":weak_scan(m[:,mi],w,1),
            "empirical_ESS_ge64":weak_scan(m[:,mi],w,64),
            "discovery":discovery,
            "holdout":{"threshold":threshold,"coefficient":float(z.mean()),
                "SE_descriptive":float(z.std(ddof=1)/math.sqrt(len(z))),
                "ESS":total*total/squares if squares else 0.,
                "contribution_cap":bound,"max_observed":float(z.max()),
                "CI_lower":confidence(float(z.mean()),bound,len(z),alpha,False),
                "CI_upper":confidence(float(z.mean()),bound,len(z),alpha,True)}}
    out={"n":src.n,"K":12,"R":src.r,"N":config["N"],"seed":config["seed"],
        "delta":src.delta,"J":len(src.lengths),"s_labels":src.s.tolist(),
        "source_mass_exact":1,"actual_prior_exact":str(F(1,src.r)),
        "mass_weight_mean":float(w.mean()),
        "mass_SE_descriptive":float(w.std(ddof=1)/math.sqrt(len(w))),
        "weight_max":float(w.max()),"uniform_rows":uniform,
        "tilted_model_counts":models.tolist(),"source_b_counts":source_b.tolist(),
        "max_actual_weak_contribution":float((w[:,None]*m).max()),
        "contribution_cap_all_thresholds":2*(src.r+2)*len(src.lengths),
        "models":results,"runtime_seconds":time.time()-started}
    print(json.dumps({"n":src.n,"R":src.r,"seconds":out["runtime_seconds"],
        "mass":out["mass_weight_mean"],"best_joint":results["joint"]["empirical_ESS_ge64"],
        "best_prior":results["prior"]["empirical_ESS_ge64"]}),flush=True)
    return out


def run():
    reg=json.loads((OUT/"registration.json").read_text())
    assert hashlib.sha256(Path(__file__).read_bytes()).hexdigest()==reg["script_SHA256"]
    if (OUT/"RESULTS.json").exists():raise RuntimeError("Refusing extra run.")
    started=time.time()
    cases=[run_case(c,reg["alpha_each"])for c in reg["cases"]]
    diagnostic=fourier_diagnostics()
    result={"registration_SHA256":hashlib.sha256(
                (OUT/"registration.json").read_bytes()).hexdigest(),
        "python":sys.executable,"python_version":sys.version,"numpy":np.__version__,
        "scope":"finite true mod-cell source screen; no fullweak ranking/asymptotics",
        "alpha_each":reg["alpha_each"],"cases":cases,
        "macroscopic_Fourier_diagnostics":diagnostic,
        "runtime_seconds":time.time()-started}
    (OUT/"RESULTS.json").write_text(json.dumps(result,separators=(",",":")))
    print("COMPLETE",result["runtime_seconds"],flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--small",action="store_true")
    parser.add_argument("--calibrate",action="store_true")
    args=parser.parse_args()
    if args.small:
        result=small_control(2)
        (OUT/"first_small_control.json").write_text(json.dumps(result,indent=2))
        print(json.dumps(result),flush=True)
    elif args.calibrate:calibrate()
    else:run()
