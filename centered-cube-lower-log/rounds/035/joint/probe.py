"""Joint-cell source: exact finite sums and weighted Poisson-binomial DP."""
from pathlib import Path
from fractions import Fraction as F
import itertools
import math
import json
import time
import sys
import hashlib
import argparse
import numpy as np

OUT=Path(__file__).resolve().parent
MODELS=["joint","fixed0","fixed1","prior_mix"]
CASES=[(n,k)for n in [3,5,7]for k in [4,6]]


class Source:
    def __init__(self,n,k):
        assert n%2==1 and k%2==0
        self.n,self.k=n,k
        self.delta=1/(64*k*k)
        j=np.arange(k)
        u0=(j+.5)/k
        self.atoms=np.array([(j+u0)/k,(j+1-u0)/k])
        self.majority=(n+1)//2
        geometric=np.geomspace(self.delta/16,1/(2*k),17)
        linear=np.linspace(1/(2*k),1.25,49)
        critical=np.array([1.,1/k,2/k,3/k,4/k])
        self.lengths=np.unique(np.concatenate(
            [geometric,linear,critical[critical<=1.25]]))
        assert np.any(self.lengths==1.)
        assert self.atoms.min()-self.delta/2>0
        assert self.atoms.max()+self.delta/2<1
        self.density_max=(64*k)**n
        self.hgrid=np.arange(math.ceil(math.log(self.density_max)/.125)+1)*.125


def coordinate_weights(src,x,length,b):
    # Each entry is (K L)^-1 times exact periodic box mass in centered interval.
    if length==1.:
        return np.full((*x.shape,src.k),1/src.k)
    distance=x[:,:,None]-src.atoms[b][None,None,:]
    mass=np.zeros_like(distance)
    for shift in [-1,0,1]:
        d=distance-shift
        mass+=np.clip((d+length/2+src.delta/2)/src.delta,0,1) \
             -np.clip((d-length/2+src.delta/2)/src.delta,0,1)
    assert mass.min()>=-1e-12
    return mass/(src.k*length)


def averages_dp(src,x,length):
    if length==1.:
        return np.ones((len(x),4))
    result=np.zeros((len(x),4))
    for b in [0,1]:
        w=coordinate_weights(src,x,length,b)
        plus=w[:,:,:src.k//2].sum(axis=2)
        minus=w[:,:,src.k//2:].sum(axis=2)
        coefficient=np.zeros((len(x),src.n+1))
        coefficient[:,0]=1
        for i in range(src.n):
            nxt=np.zeros_like(coefficient)
            nxt[:,:i+1]+=coefficient[:,:i+1]*minus[:,i,None]
            nxt[:,1:i+2]+=coefficient[:,:i+1]*plus[:,i,None]
            coefficient=nxt
        result[:,1+b]=np.prod(plus+minus,axis=1)
        branch=(coefficient[:,:src.majority].sum(axis=1)if b==0 else
                coefficient[:,src.majority:].sum(axis=1))
        result[:,0]+=branch
    result[:,3]=.5*(result[:,1]+result[:,2])
    return result


def cell_indices(src):
    indices=np.indices((src.k,)*src.n,dtype=np.int16).reshape(src.n,-1).T
    gate=(indices<src.k//2).sum(axis=1)>=src.majority
    assert int(gate.sum())*2==len(gate)
    return indices,gate


def averages_direct(src,x,length,indices=None,gate=None):
    if indices is None:indices,gate=cell_indices(src)
    products=[]
    for b in [0,1]:
        w=coordinate_weights(src,x[None,:],length,b)[0]
        product=np.ones(len(indices))
        for i in range(src.n):product*=w[i,indices[:,i]]
        products.append(product)
    f0,f1=products
    return np.array([f0[~gate].sum()+f1[gate].sum(),
                     f0.sum(),f1.sum(),.5*(f0.sum()+f1.sum())])


def fraction_check(k):
    n=3;delta=F(1,64*k*k);majority=2
    atoms=[]
    for b in [0,1]:
        row=[]
        for j in range(k):
            u=F(2*j+1,2*k)
            if b:u=1-u
            row.append((j+u)/k)
        atoms.append(row)
    x=[F(1,7),F(2,7),F(4,7)]
    rows=[]
    clip=lambda z:max(F(0),min(F(1),z))
    for length in [F(3,8),F(5,8),F(1)]:
        weights=[]
        for b in [0,1]:
            a=[]
            for xx in x:
                row=[]
                for atom in atoms[b]:
                    w=sum((clip((xx-atom-shift+length/2+delta/2)/delta)
                           -clip((xx-atom-shift-length/2+delta/2)/delta)
                           for shift in [-1,0,1]),F(0))
                    row.append(w/(k*length))
                a.append(row)
            weights.append(a)
        direct=[F(0)]*3
        for cells in itertools.product(range(k),repeat=n):
            b=int(sum(j<k//2 for j in cells)>=majority)
            direct[0]+=math.prod(weights[b][i][j]for i,j in enumerate(cells))
            for bb in [0,1]:
                direct[bb+1]+=math.prod(weights[bb][i][j]
                                        for i,j in enumerate(cells))
        dp=[F(0)]*3
        for b in [0,1]:
            coeff=[F(1)]
            for i in range(n):
                plus=sum(weights[b][i][:k//2],F(0))
                minus=sum(weights[b][i][k//2:],F(0))
                nxt=[F(0)]*(len(coeff)+1)
                for q,c in enumerate(coeff):
                    nxt[q]+=c*minus;nxt[q+1]+=c*plus
                coeff=nxt
            dp[0]+=sum(coeff[:majority]if b==0 else coeff[majority:],F(0))
            dp[1+b]=sum(coeff,F(0))
        assert direct==dp
        rows.append({"L":str(length),"joint":str(dp[0]),
                     "fixed0":str(dp[1]),"fixed1":str(dp[2]),
                     "prior_mix":str((dp[1]+dp[2])/2),
                     "exact_equality":True})
    return {"K":k,"n":3,"x":[str(a)for a in x],"rows":rows}


def small_control(n,k,seed=202610290350):
    started=time.time();src=Source(n,k);rng=np.random.default_rng(seed)
    indices,gate=cell_indices(src)
    cell=np.arange(n)%k
    b=int(np.count_nonzero(cell<k//2)>=src.majority)
    atom=src.atoms[b,cell]
    centers=np.vstack([rng.random((2,n)),atom])
    errors=[];examples=[]
    for length in [src.delta,1/k,7/8,1.]:
        dp=averages_dp(src,centers,length)
        for row,x in enumerate(centers):
            direct=averages_direct(src,x,length,indices,gate)
            rel=np.max(np.abs(dp[row]-direct)/np.maximum(1,np.abs(direct)))
            assert rel<1e-10
            errors.append(float(rel))
            if length==1/k:
                examples.append({"center":row,"direct":direct.tolist(),
                                 "DP":dp[row].tolist()})
    return {"n":n,"K":k,"cells_enumerated":len(indices),
            "b1_cell_count":int(gate.sum()),"b0_cell_count":int((~gate).sum()),
            "actual_prior":[.5,.5],"queries_checked":len(errors),
            "max_relative_error":max(errors),"examples_at_L_1_over_K":examples,
            "seconds":time.time()-started}


def draw_proposal(src,rng,samples):
    x=rng.random((samples,src.n))
    tilted=rng.random(samples)>=.5
    rows=np.flatnonzero(tilted)
    model=rng.integers(4,size=len(rows))
    li=rng.integers(len(src.lengths),size=len(rows))
    cells=rng.integers(src.k,size=(len(rows),src.n))
    b=np.zeros(len(rows),dtype=int)
    b[model==0]=((cells[model==0]<src.k//2).sum(axis=1)>=src.majority)
    b[model==2]=1
    b[model==3]=rng.integers(2,size=np.count_nonzero(model==3))
    points=src.atoms[b[:,None],cells]+src.delta*(rng.random(cells.shape)-.5)
    x[rows]=(points-src.lengths[li,None]*(rng.random(cells.shape)-.5))%1
    audit={"uniform":int((~tilted).sum()),
           "tilted_model_counts":np.bincount(model,minlength=4).tolist(),
           "tilted_source_b_counts_by_model":[
               np.bincount(b[model==m],minlength=2).tolist()for m in range(4)]}
    return x,audit


def evaluate_queries(src,x):
    maxima=np.zeros((len(x),4))
    density_sum=np.zeros(len(x))
    for length in src.lengths:
        values=averages_dp(src,x,float(length))
        assert values.min()>=0
        maxima=np.maximum(maxima,values)
        density_sum+=values.sum(axis=1)
    weights=1/(.5+.5*density_sum/(4*len(src.lengths)))
    assert weights.max()<=2
    assert maxima.min()>=1
    assert maxima.max()<=src.density_max*(1+1e-8)
    assert (weights[:,None]*maxima).max()<=8*len(src.lengths)*(1+1e-8)
    return maxima,weights


def calibration():
    started=time.time()
    controls=[small_control(n,k)for n,k in CASES]
    fractions=[fraction_check(k)for k in [4,6]]
    src=Source(7,6)
    x=(.12345+np.arange(512*7)*((math.sqrt(5)-1)/2))%1
    x=x.reshape(512,7)
    t=time.time();evaluate_queries(src,x);benchmark=time.time()-t
    estimate_16384=benchmark*(16384/512)*len(CASES)*1.5
    nmc=16384 if estimate_16384<=120 else 8192
    estimated=benchmark*(nmc/512)*len(CASES)*1.5
    result={"controls":controls,"fraction_controls":fractions,
            "benchmark_n":7,"benchmark_K":6,"benchmark_batch":512,
            "benchmark_J":len(src.lengths),"benchmark_seconds":benchmark,
            "estimated_16384_six_case_seconds":estimate_16384,
            "selected_N_per_case":nmc,"estimated_six_case_seconds":estimated,
            "within_180_second_MC_budget":estimated<=180,
            "calibration_total_seconds":time.time()-started}
    (OUT/"calibration.json").write_text(json.dumps(result,indent=2))
    print(json.dumps({k:result[k]for k in [
        "benchmark_seconds","selected_N_per_case","estimated_six_case_seconds",
        "calibration_total_seconds"]}),flush=True)
    if estimated>180:
        raise RuntimeError("Cost estimate exceeds frozen budget; no MC run.")
    configs=[]
    total_h=0
    for i,(n,k)in enumerate(CASES):
        s=Source(n,k);total_h+=4*len(s.hgrid)
        configs.append({"n":n,"K":k,"seed":202610290351+i,"N":nmc,
            "delta":s.delta,"J":len(s.lengths),"lengths":s.lengths.tolist(),
            "hgrid":s.hgrid.tolist(),"density_max":s.density_max,
            "actual_b_prior":[.5,.5]})
    reg={"status":"frozen after exact controls/cost calibration, before MC",
         "script_SHA256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         "plan_SHA256":hashlib.sha256((OUT/"PLAN.md").read_bytes()).hexdigest(),
         "calibration_SHA256":hashlib.sha256((OUT/"calibration.json").read_bytes()).hexdigest(),
         "cases":configs,"batch":512,"model_order":MODELS,
         "b_rule":"majority of j_i<K/2; b is source-cell dependent",
         "u0":"(j+1/2)/K","u1":"1-u0","smoothing":"product box delta=1/(64K²)",
         "proposal":".5 uniform + .5 all four sources/all query lengths",
         "source_sampling":"uniform joint cell, THEN b=f(cell) for joint; "
                           "prior-mixture draws b independently with actual pi=.5",
         "alpha_total":.05,"alpha_each":.05/(2*total_h),
         "full_profile_basis":"M>=1 from L=1; M<=(64K)^n from disjoint box density",
         "contribution_cap":"min(2lambda,8J)",
         "cost":{"memory_limit_bytes":2_000_000_000,
                 "estimated_MC_seconds":estimated},
         "stopping":"one fixed seed/case; no adaptation to gains; no asymptotics"}
    if (OUT/"registration.json").exists():
        raise RuntimeError("Registration already exists.")
    (OUT/"registration.json").write_text(json.dumps(reg,separators=(",",":")))


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
        divergence=kl(x,mid)
        if upper:
            if divergence>target:hi=mid
            else:lo=mid
        else:
            if divergence>target:lo=mid
            else:hi=mid
    return bound*(hi if upper else lo)


def empirical_weak(values,weights,min_ess):
    order=np.argsort(values)[::-1]
    v=values[order];w=weights[order]
    sums=np.cumsum(w);square=np.cumsum(w*w)
    ess=sums*sums/square
    score=v*sums/len(v)
    index=int(np.argmax(np.where(ess>=min_ess,score,-np.inf)))
    return {"coefficient":float(score[index]),"threshold":float(v[index]),
            "tail_ESS":float(ess[index]),"min_ESS_filter":min_ess}


def run_case(config,alpha):
    started=time.time()
    src=Source(config["n"],config["K"])
    assert src.lengths.tolist()==config["lengths"]
    rng=np.random.default_rng(config["seed"])
    all_m=[];all_w=[]
    audit={"uniform":0,"tilted_model_counts":np.zeros(4,dtype=int),
           "tilted_source_b_counts_by_model":np.zeros((4,2),dtype=int)}
    for left in range(0,config["N"],512):
        count=min(512,config["N"]-left)
        x,aa=draw_proposal(src,rng,count)
        values,w=evaluate_queries(src,x)
        all_m.append(values);all_w.append(w)
        audit["uniform"]+=aa["uniform"]
        audit["tilted_model_counts"]+=aa["tilted_model_counts"]
        audit["tilted_source_b_counts_by_model"]+=aa["tilted_source_b_counts_by_model"]
    maxima=np.concatenate(all_m);weights=np.concatenate(all_w)
    n_samples=len(weights);out={}
    for mi,model in enumerate(MODELS):
        curve=[]
        for h in src.hgrid:
            lam=math.exp(float(h));z=lam*weights*(maxima[:,mi]>lam)
            bound=min(2*lam,8*len(src.lengths))
            assert z.max()<=bound*(1+1e-8)
            mean=float(z.mean())
            lo=confidence(mean,bound,n_samples,alpha,False)
            hi=confidence(mean,bound,n_samples,alpha,True)
            if lam>=src.density_max:hi=0.
            else:hi=min(hi,lam)
            total=float(z.sum());squares=float(np.dot(z,z))
            curve.append({"h":float(h),"lambda":lam,"coefficient_MC":mean,
                "SE_descriptive":float(z.std(ddof=1)/math.sqrt(n_samples)),
                "tail_ESS":total*total/squares if squares else 0.,
                "CI_lower":lo,"CI_upper":hi,"contribution_bound":bound,
                "max_observed_contribution":float(z.max())})
        interval_upper=[min(math.exp(curve[i+1]["h"]),src.density_max)
            /curve[i]["lambda"]*curve[i]["CI_upper"]
            for i in range(len(curve)-1)]
        best=int(np.argmax([row["coefficient_MC"]for row in curve]))
        out[model]={"profile":curve,"full_weak_lower_CI":max(
                1.,max(row["CI_lower"]for row in curve)),
            "full_weak_upper_CI":max(1.,max(interval_upper)),
            "best_frozen_grid_MC":curve[best],
            "empirical_weak_unrestricted":empirical_weak(maxima[:,mi],weights,1),
            "empirical_weak_tail_ESS_ge64":empirical_weak(maxima[:,mi],weights,64)}
    result={"n":src.n,"K":src.k,"seed":config["seed"],"N":n_samples,
        "delta":src.delta,"J":len(src.lengths),"density_max":src.density_max,
        "source_mass_exact":1,"actual_b_prior_exact":["1/2","1/2"],
        "mass_weight_mean":float(weights.mean()),
        "mass_weight_SE_descriptive":float(weights.std(ddof=1)/math.sqrt(n_samples)),
        "weight_max":float(weights.max()),
        "max_actual_coefficient_contribution":float(
            (weights[:,None]*maxima).max()),
        "universal_coefficient_contribution_cap":8*len(src.lengths),
        "audit":{"uniform":audit["uniform"],
            "tilted_model_counts":audit["tilted_model_counts"].tolist(),
            "tilted_source_b_counts_by_model":
                audit["tilted_source_b_counts_by_model"].tolist()},
        "models":out,"runtime_seconds":time.time()-started}
    print(json.dumps({"n":src.n,"K":src.k,"seconds":result["runtime_seconds"],
        "mass":result["mass_weight_mean"],"best":{
            m:out[m]["best_frozen_grid_MC"]["coefficient_MC"]for m in MODELS}}),
        flush=True)
    return result


def run():
    reg=json.loads((OUT/"registration.json").read_text())
    assert hashlib.sha256(Path(__file__).read_bytes()).hexdigest()==reg["script_SHA256"]
    if (OUT/"RESULTS.json").exists():raise RuntimeError("Refusing another MC run.")
    started=time.time()
    cases=[run_case(config,reg["alpha_each"])for config in reg["cases"]]
    controls=json.loads((OUT/"calibration.json").read_text())
    result={"registration_SHA256":hashlib.sha256(
                (OUT/"registration.json").read_bytes()).hexdigest(),
        "python":sys.executable,"python_version":sys.version,"numpy":np.__version__,
        "scope":"true nonproduct joint-cell finite n and finite L model, no asymptotics",
        "calibration":controls,"alpha_total":reg["alpha_total"],
        "alpha_each":reg["alpha_each"],"cases":cases,
        "runtime_seconds":time.time()-started}
    (OUT/"RESULTS.json").write_text(json.dumps(result,separators=(",",":")))
    print("COMPLETE",result["runtime_seconds"],flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--small",action="store_true")
    parser.add_argument("--calibrate",action="store_true")
    args=parser.parse_args()
    if args.small:
        result={"control":small_control(3,4),"fraction":fraction_check(4)}
        (OUT/"first_small_control.json").write_text(json.dumps(result,indent=2))
        print(json.dumps(result,ensure_ascii=False),flush=True)
    elif args.calibrate:calibration()
    else:run()
