"""Exact quantile-cell source; stable positive PB probability times log mean."""
from pathlib import Path
from fractions import Fraction as F
import numpy as np
import math
import itertools
import json
import time
import argparse
import hashlib
import sys

OUT=Path(__file__).resolve().parent
CASES=[(n,r)for n in [15,31,63]for r in [2,4,8]]


def quantiles(n,r):
    den=2**n
    counts=[math.comb(n,s)for s in range(n+1)]
    cumulative=list(itertools.accumulate(counts))
    upper=[]
    for b in range(1,r):
        upper.append(next(s for s,c in enumerate(cumulative)if r*c>=b*den))
    upper.append(n)
    lower=[0]+[x+1 for x in upper[:-1]]
    numerator=[sum(counts[a:b+1])if a<=b else 0 for a,b in zip(lower,upper)]
    assert sum(numerator)==den
    return lower,upper,numerator,den


class Source:
    def __init__(self,n,k,r):
        assert k%2==0
        self.n,self.k,self.r=n,k,r
        self.lower,self.upper,self.prior_num,self.prior_den=quantiles(n,r)
        self.prior=np.array([float(F(a,self.prior_den))for a in self.prior_num])
        theta=(np.arange(r)+.5)/r
        j=np.arange(k)
        u=np.where(j[None,:]<k//2,theta[:,None],1-theta[:,None])
        self.atoms=(j[None,:]+u)/k
        self.delta=1/(64*k*r)
        assert self.atoms.min()-self.delta/2>0
        assert self.atoms.max()+self.delta/2<1
        short=np.geomspace(self.delta/16,1/(2*k),13)
        linear=np.linspace(1/(2*k),1.25,41)
        self.lengths=np.unique(np.concatenate([short,linear,
            [1/k,2/k,3/k,4/k,6/k,1.,10/9,7/6,11/9]]))
        self.models=["joint"]+[f"fixed{b}"for b in range(r)]+["prior"]
        self.log_density_max=n*math.log(64*r)
        self.max_probability_mass_residual=0.
        self.subnormal_probability_entries=0

    def classify(self,s):
        return np.searchsorted(self.upper,s,side="left")

    def metadata(self):
        return {"n":self.n,"K":self.k,"R":self.r,"delta":self.delta,
            "segments":[{"b":b,"lower":a,"upper":z,"empty":a>z,
                         "prior_exact":str(F(p,self.prior_den))}
                        for b,(a,z,p)in enumerate(zip(
                            self.lower,self.upper,self.prior_num))],
            "active_classes":int(np.count_nonzero(self.prior)),
            "source_mass_exact":1,"J":len(self.lengths),
            "log_density_max":self.log_density_max}


def logsumexp(a,axis=-1):
    maximum=np.max(a,axis=axis,keepdims=True)
    safe=np.where(np.isfinite(maximum),maximum,0)
    with np.errstate(divide="ignore",invalid="ignore",under="ignore"):
        result=safe+np.log(np.sum(np.exp(a-safe),axis=axis,keepdims=True))
    return np.squeeze(result,axis=axis)


def class_weights(src,x,length):
    distance=x[:,:,None,None]-src.atoms[None,None,:,:]
    mass=np.zeros_like(distance)
    for shift in [-1,0,1]:
        d=distance-shift
        mass+=np.clip((d+length/2+src.delta/2)/src.delta,0,1) \
              -np.clip((d-length/2+src.delta/2)/src.delta,0,1)
    w=mass/(src.k*length)
    return w[:,:,:,:src.k//2].sum(axis=3),w[:,:,:,src.k//2:].sum(axis=3)


def log_averages(src,x,length,return_probability=False):
    if length==1. and not return_probability:
        return np.zeros((len(x),src.r+2))
    plus,minus=class_weights(src,x,length)
    total=plus+minus
    with np.errstate(divide="ignore",invalid="ignore"):
        log_fixed=np.log(total).sum(axis=1)
    # Tables with a zero coordinate total have exact product mean zero.
    # Pack only active (query,table) pairs, especially at very short lengths.
    rows,tables=np.where(np.all(total>0,axis=1))
    pp=plus[rows,:,tables]/total[rows,:,tables]
    pm=minus[rows,:,tables]/total[rows,:,tables]
    probability=np.zeros((len(rows),src.n+1));probability[:,0]=1
    for i in range(src.n):
        nxt=np.zeros_like(probability)
        nxt[:,:i+1]+=probability[:,:i+1]*pm[:,i,None]
        nxt[:,1:i+2]+=probability[:,:i+1]*pp[:,i,None]
        probability=nxt
    branch=np.zeros((len(x),src.r))
    for b,(lo,hi)in enumerate(zip(src.lower,src.upper)):
        selected=np.flatnonzero(tables==b)
        if lo<=hi:branch[rows[selected],b]=probability[selected,lo:hi+1].sum(axis=1)
    if len(rows):
        src.max_probability_mass_residual=max(src.max_probability_mass_residual,
            float(np.abs(probability.sum(axis=1)-1).max()))
        src.subnormal_probability_entries+=int(np.count_nonzero(
            (probability>0)&(probability<np.finfo(float).tiny)))
    with np.errstate(divide="ignore"):
        joint=logsumexp(log_fixed+np.log(branch),axis=1)
        prior=logsumexp(log_fixed+np.log(src.prior)[None,:],axis=1)
    result=np.column_stack([joint,log_fixed,prior])
    return (result,probability)if return_probability else result


def direct_average(src,x,length):
    cells=np.indices((src.k,)*src.n,dtype=np.int16).reshape(src.n,-1).T
    b=src.classify((cells<src.k//2).sum(axis=1))
    histogram=np.bincount(b,minlength=src.r)
    assert np.array_equal(histogram,[
        src.k**src.n*p//src.prior_den for p in src.prior_num])
    joint=0.;fixed=[]
    for table in range(src.r):
        d=x[:,None]-src.atoms[table][None,:]
        mass=np.zeros_like(d)
        for shift in [-1,0,1]:
            mass+=np.clip((d-shift+length/2+src.delta/2)/src.delta,0,1) \
                -np.clip((d-shift-length/2+src.delta/2)/src.delta,0,1)
        w=mass/(src.k*length);product=np.ones(len(cells))
        for i in range(src.n):product*=w[i,cells[:,i]]
        fixed.append(product.sum());joint+=product[b==table].sum()
    return np.array([joint]+fixed+[np.dot(src.prior,fixed)])


def small_control(r):
    src=Source(3,4,r)
    cell=np.arange(3)%4;b=int(src.classify(np.count_nonzero(cell<2)))
    centers=np.array([[1/7,2/7,4/7],[.11,.47,.82],src.atoms[b,cell]])
    rows=[]
    for length in [src.delta,.25,3/8,1.,7/6]:
        logs=log_averages(src,centers,length)
        for i,x in enumerate(centers):
            direct=direct_average(src,x,length)
            error=float(np.max(np.abs(np.exp(logs[i])-direct)
                               /np.maximum(1,np.abs(direct))))
            assert error<1e-10
            rows.append({"L":length,"center":i,"relative_error":error})
    return {"metadata":src.metadata(),"cells":64,
            "max_relative_error":max(x["relative_error"]for x in rows),"rows":rows}


def fraction_control(r):
    n,k=3,4;lo,hi,nums,den=quantiles(n,r);delta=F(1,64*k*r)
    atoms=[]
    for b in range(r):
        theta=F(2*b+1,2*r)
        atoms.append([(j+(theta if j<k//2 else 1-theta))/k for j in range(k)])
    x=[F(1,7),F(2,7),F(4,7)];rows=[]
    clip=lambda v:max(F(0),min(F(1),v))
    for length in [F(3,8),F(7,6),F(1)]:
        weights=[[[sum((clip((xx-a-s+length/2+delta/2)/delta)
            -clip((xx-a-s-length/2+delta/2)/delta)for s in [-1,0,1]),F(0))
            /(k*length)for a in atoms[b]]for xx in x]for b in range(r)]
        joint=F(0);fixed=[F(0)]*r
        for cells in itertools.product(range(k),repeat=n):
            s=sum(j<k//2 for j in cells)
            b=next(i for i,cut in enumerate(hi)if s<=cut)
            joint+=math.prod(weights[b][i][j]for i,j in enumerate(cells))
            for bb in range(r):
                fixed[bb]+=math.prod(weights[bb][i][j]for i,j in enumerate(cells))
        dp=F(0);dp_fixed=[]
        for b in range(r):
            coeff=[F(1)]
            for i in range(n):
                plus=sum(weights[b][i][:k//2],F(0))
                minus=sum(weights[b][i][k//2:],F(0))
                nxt=[F(0)]*(len(coeff)+1)
                for q,c in enumerate(coeff):
                    nxt[q]+=c*minus;nxt[q+1]+=c*plus
                coeff=nxt
            dp+=sum(coeff[lo[b]:hi[b]+1],F(0))if lo[b]<=hi[b]else F(0)
            dp_fixed.append(sum(coeff,F(0)))
        assert dp==joint and dp_fixed==fixed
        rows.append({"L":str(length),"joint":str(joint),
            "prior":str(sum((F(p,den)*f for p,f in zip(nums,fixed)),F(0))),
            "exact_equality":True})
    return {"n":n,"K":k,"R":r,"prior_exact":[str(F(x,den))for x in nums],
            "empty_segments":[b for b in range(r)if lo[b]>hi[b]],"rows":rows}


def draw_proposal(src,rng,count):
    x=rng.random((count,src.n))
    tilted=rng.random(count)>=.5;rows=np.flatnonzero(tilted)
    models=rng.integers(src.r+2,size=len(rows))
    li=rng.integers(len(src.lengths),size=len(rows))
    cells=rng.integers(src.k,size=(len(rows),src.n))
    b=np.zeros(len(rows),dtype=int)
    joint=models==0
    b[joint]=src.classify((cells[joint]<src.k//2).sum(axis=1))
    for table in range(src.r):b[models==table+1]=table
    prior=models==src.r+1
    # An independent exact Binomial classifier gives the actual prior,
    # avoiding a floating categorical approximation for dyadic priors.
    auxiliary=rng.integers(2,size=(int(prior.sum()),src.n)).sum(axis=1)
    b[prior]=src.classify(auxiliary)
    assert all(src.prior_num[int(bb)]>0 for bb in b[joint|prior])
    t=src.atoms[b[:,None],cells]+src.delta*(rng.random(cells.shape)-.5)
    x[rows]=(t-src.lengths[li,None]*(rng.random(cells.shape)-.5))%1
    audit=np.zeros((src.r+2,src.r),dtype=int)
    for m,table in zip(models,b):audit[int(m),int(table)]+=1
    return x,{"uniform":int((~tilted).sum()),"source_b_by_model":audit}


def evaluate_queries(src,x):
    maximum=np.full((len(x),src.r+2),-np.inf)
    log_sum=np.full(len(x),-np.inf)
    for length in src.lengths:
        values=log_averages(src,x,float(length))
        maximum=np.maximum(maximum,values)
        log_sum=np.logaddexp(log_sum,logsumexp(values,axis=1))
    log_d=np.logaddexp(math.log(.5),math.log(.5)+log_sum
                      -math.log((src.r+2)*len(src.lengths)))
    assert maximum.min()>=0
    assert maximum.max()<=src.log_density_max+1e-7
    assert np.exp(-log_d).max()<=2
    cap=2*(src.r+2)*len(src.lengths)
    assert np.exp(maximum-log_d[:,None]).max()<=cap*(1+1e-8)
    return maximum,log_d


def calibrate():
    started=time.time()
    controls=[small_control(r)for r in [2,4,8]]
    fractions=[fraction_control(r)for r in [2,4,8]]
    src=Source(63,12,8)
    x=(.12345+np.arange(64*63)*((math.sqrt(5)-1)/2))%1
    x=x.reshape(64,63)
    t=time.time();evaluate_queries(src,x);bench=time.time()-t
    selected=None
    for nmc in [4096,2048,1024]:
        estimate=bench*(nmc/64)*len(CASES)*1.3
        if estimate<=180:
            selected=nmc;break
    calibration={"controls":controls,"Fraction_controls":fractions,
        "worst_n":63,"worst_R":8,"batch":64,"benchmark_seconds":bench,
        "selected_N":selected,"estimated_seconds":estimate,
        "max_PB_mass_residual_benchmark":src.max_probability_mass_residual,
        "positive_subnormal_probability_entries":src.subnormal_probability_entries,
        "total_seconds":time.time()-started}
    (OUT/"calibration.json").write_text(json.dumps(calibration,separators=(",",":")))
    print(json.dumps({k:calibration[k]for k in [
        "benchmark_seconds","selected_N","estimated_seconds","total_seconds"]}),
        flush=True)
    if selected is None:raise RuntimeError("Cost gate exceeds 180 seconds.")
    models_count=sum(r+2 for _,r in CASES)
    configs=[]
    for i,(n,r)in enumerate(CASES):
        s=Source(n,12,r)
        configs.append({"n":n,"K":12,"R":r,"seed":202610290371+i,
            "N":selected,"source":s.metadata(),"lengths":s.lengths.tolist()})
    reg={"status":"frozen after exact controls and cost, before MC",
        "script_SHA256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "plan_SHA256":hashlib.sha256((OUT/"PLAN.md").read_bytes()).hexdigest(),
        "calibration_SHA256":hashlib.sha256(
            (OUT/"calibration.json").read_bytes()).hexdigest(),
        "quantile":"first integer CDF >= b/R; repeated cuts retained",
        "cases":configs,"batch":64,"alpha_each":.05/(2*models_count),
        "prior_sampling":"independent n fair bits then exact classifier",
        "stability":"positive PB probabilities; direct segment sums; "
                    "log product means, logsumexp, no Gaussian approximation",
        "cost_limit_seconds":180,"RAM_limit_bytes":2_000_000_000,
        "stopping":"one fixed seed per case; no gain-adaptive K/R/N/seeds",
        "CI_scope":"heldout selected tail only; not fullweak ordering"}
    if (OUT/"registration.json").exists():raise RuntimeError("Already registered.")
    (OUT/"registration.json").write_text(json.dumps(reg,separators=(",",":")))


def weak_scan(log_values,log_d,ess_min):
    order=np.argsort(log_values)[::-1]
    v=log_values[order];w=np.exp(-log_d[order])
    sums=np.cumsum(w);squares=np.cumsum(w*w)
    ess=sums*sums/squares
    score=v+np.log(sums/len(w))
    ix=int(np.argmax(np.where(ess>=ess_min,score,-np.inf)))
    return {"coefficient":float(np.exp(score[ix])),
            "log_threshold":float(v[ix]),"ESS":float(ess[ix]),"min_ESS":ess_min}


def kl(x,mu):
    return (x*math.log(x/mu)if x else 0) \
        +((1-x)*math.log((1-x)/(1-mu))if x<1 else 0)


def confidence(mean,bound,n,alpha,upper):
    x=min(1.,max(0.,mean/bound))
    if upper and x==1:return bound
    if not upper and x==0:return 0.
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


def macro_diagnostic():
    output=[]
    for n,r in CASES:
        src=Source(n,12,r)
        seed=202610292000+10*n+r
        x=np.random.default_rng(seed).random((64,n))
        for length in [10/9,7/6,11/9]:
            values=log_averages(src,x,length)
            relative=np.expm1(values[:,0]-values[:,-1])
            output.append({"n":n,"R":r,"L":length,"seed":seed,"centers":64,
                "max_abs_joint_over_prior_minus1":float(np.abs(relative).max()),
                "relative_quantiles":np.quantile(relative,[0,.25,.5,.75,1]).tolist()})
    return output


def run_case(config,alpha):
    started=time.time();src=Source(config["n"],12,config["R"])
    assert src.lengths.tolist()==config["lengths"]
    rng=np.random.default_rng(config["seed"])
    maxima=[];densities=[];uniform=0
    audit=np.zeros((src.r+2,src.r),dtype=int)
    for left in range(0,config["N"],64):
        count=min(64,config["N"]-left)
        x,a=draw_proposal(src,rng,count)
        m,d=evaluate_queries(src,x)
        maxima.append(m);densities.append(d)
        uniform+=a["uniform"];audit+=a["source_b_by_model"]
    m=np.concatenate(maxima);d=np.concatenate(densities)
    w=np.exp(-d);cut=len(w)//2;results={}
    cap=2*(src.r+2)*len(src.lengths)
    for mi,name in enumerate(src.models):
        discovery=weak_scan(m[:cut,mi],d[:cut],32)
        log_t=float(np.nextafter(discovery["log_threshold"],-np.inf))
        threshold=math.exp(log_t)
        z=np.exp(log_t-d[cut:])*(m[cut:,mi]>log_t)
        bound=min(2*threshold,cap)
        assert z.max()<=bound*(1+1e-8)
        total=float(z.sum());squares=float(np.dot(z,z))
        results[name]={"empirical_unrestricted":weak_scan(m[:,mi],d,1),
            "empirical_ESS_ge64":weak_scan(m[:,mi],d,64),
            "discovery":discovery,
            "holdout":{"log_threshold":log_t,"coefficient":float(z.mean()),
                "SE_descriptive":float(z.std(ddof=1)/math.sqrt(len(z))),
                "ESS":total*total/squares if squares else 0.,
                "contribution_cap":bound,"max_observed":float(z.max()),
                "CI_lower":confidence(float(z.mean()),bound,len(z),alpha,False),
                "CI_upper":confidence(float(z.mean()),bound,len(z),alpha,True)}}
    # Conservative magnitude for probability-subnormal loss only, NOT a
    # total floating-point error certificate. Other rounding is not certified.
    underflow_log=(math.log(2*src.r*src.n*(src.n+1))
                   +math.log(np.finfo(float).tiny)+src.log_density_max)
    out={"metadata":src.metadata(),"N":config["N"],"seed":config["seed"],
        "mass_weight_mean":float(w.mean()),
        "mass_SE_descriptive":float(w.std(ddof=1)/math.sqrt(len(w))),
        "weight_max":float(w.max()),"uniform_rows":uniform,
        "source_b_counts_by_model":audit.tolist(),
        "max_actual_coefficient_contribution":float(np.exp(m-d[:,None]).max()),
        "contribution_cap_all_thresholds":cap,
        "max_PB_probability_mass_residual":src.max_probability_mass_residual,
        "positive_subnormal_probability_entries":src.subnormal_probability_entries,
        "subnormal_loss_magnitude_budget_only":math.exp(underflow_log),
        "models":results,"runtime_seconds":time.time()-started}
    best_fixed=max(results[name]["empirical_ESS_ge64"]["coefficient"]
                    for name in src.models if name.startswith("fixed"))
    print(json.dumps({"n":src.n,"R":src.r,"seconds":out["runtime_seconds"],
        "active_classes":src.metadata()["active_classes"],
        "mass":out["mass_weight_mean"],"best_joint":
            results["joint"]["empirical_ESS_ge64"],
        "best_prior":results["prior"]["empirical_ESS_ge64"],
        "best_fixed":best_fixed}),flush=True)
    return out


def run():
    reg=json.loads((OUT/"registration.json").read_text())
    assert hashlib.sha256(Path(__file__).read_bytes()).hexdigest()==reg["script_SHA256"]
    if (OUT/"RESULTS.json").exists():raise RuntimeError("Refusing extra run.")
    started=time.time()
    cases=[run_case(c,reg["alpha_each"])for c in reg["cases"]]
    diagnostic=macro_diagnostic()
    result={"registration_SHA256":hashlib.sha256(
                (OUT/"registration.json").read_bytes()).hexdigest(),
        "python":sys.executable,"python_version":sys.version,"numpy":np.__version__,
        "scope":"finite true quantile-cell source screen; no fullweak ordering",
        "alpha_each":reg["alpha_each"],"cases":cases,
        "macroscopic_diagnostics":diagnostic,"runtime_seconds":time.time()-started}
    (OUT/"RESULTS.json").write_text(json.dumps(result,separators=(",",":")))
    print("COMPLETE",result["runtime_seconds"],flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--small",action="store_true")
    parser.add_argument("--calibrate",action="store_true")
    args=parser.parse_args()
    if args.small:
        out={"direct":[small_control(r)for r in [2,4]],
             "Fraction":[fraction_control(r)for r in [2,4]]}
        (OUT/"first_small_control.json").write_text(json.dumps(out,indent=2))
        print(json.dumps(out),flush=True)
    elif args.calibrate:calibrate()
    else:run()
