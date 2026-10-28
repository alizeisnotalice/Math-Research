"""Frozen full weak-profile screen; dedicated G=1 and G=8 proposals.

Registration is deterministic (--register), with no random sampling.
The ordinary run executes exactly the two registered fixed seeds, in batches.
"""
from pathlib import Path
from fractions import Fraction
import argparse
import hashlib
import json
import math
import time
import sys
import block_probe as old
import numpy as np

HERE=Path(__file__).resolve().parent
N=512
BATCH=64
SEEDS={1:202610290343,8:202610290344}
R=4


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def build_cache():
    started=time.time()
    src,grid=old.configure_source(16)
    j=len(grid)
    means=np.empty((j,src.q))
    transition=np.empty((j,src.q),dtype=bool)
    norm=np.empty((j,R))
    logs2=[]
    second_raw=[]
    for li in range(j):
        m,t=src.base_query(li)
        means[li]=m;transition[li]=t
        for r in range(R):
            norm[li,r]=np.roll(m,r*(src.q//R)).mean()
        one_coordinate=math.fsum(
            (m*m*np.where(t,1+src.t1,1+src.t0)).tolist())/src.q
        assert one_coordinate>0
        second_raw.append(one_coordinate)
        logs2.append(src.n*math.log(one_coordinate))
    log_h2=float(old.lse(np.array(logs2)[None,:])[0])
    h2=math.exp(log_h2)
    steps=math.ceil(log_h2/.25)
    hgrid=np.arange(steps+1)*.25
    assert hgrid[0]==0 and hgrid[-1]>=log_h2
    return src,grid,means,transition,norm,{
        "log_H2":log_h2,"H2":h2,
        "log_second_moment_bounds_by_query":logs2,
        "one_coordinate_second_moment_bounds":second_raw,
        "h_grid":hgrid.tolist(),
        "cache_bytes":means.nbytes+transition.nbytes+norm.nbytes,
        "deterministic_preflight_seconds":time.time()-started,
        "max_raw_mean_normalizer_residual":float(np.abs(norm-1).max()),
        "C_Z":math.exp(max(0.,float(src.n*np.log(norm).max())))}


def estimated_peak(src,cache_bytes,g):
    phase=BATCH*g*len(src.lengths)*R*8
    return int(cache_bytes+BATCH*g*src.q*8+4*phase
               +8*src.q*R*8+50_000_000)


def register():
    src,grid,means,transition,norm,det=build_cache()
    prior=["REPORT.md","registration.json","block_probe.py","RESULTS.json"]
    entry={
        "status":"frozen deterministic H2/profile before random sampling",
        "script_SHA256":sha(Path(__file__)),
        "preserved_pilot_files_SHA256":{p:sha(HERE/p)for p in prior},
        "K":255,"v":16,"n":src.n,"B":4,"M":src.m,"p":.5,"w":4,
        "R":4,"models_G":[1,8],"N_per_model":N,"batch":BATCH,
        "seeds":{str(k):v for k,v in SEEDS.items()},
        "J":len(grid),"Q":src.q,
        "query_grid":"natural v16 union original nine; endpoints included",
        "H2_formula":"sum_j [mean_q m_qj0^2(1+t_T)]^n",
        "H2_validity":"Jensen over latent phases then sum over common target L",
        "lambda_outside_range":"lambda<=1 and lambda>=H2 give weak value<=1",
        "profile_interpolation":"lambda_(k+1) P(M>lambda_k), moment clip allowed",
        "CI":"one-sided Bernoulli-KL: lower mean of lower contribution; "
             "upper mean of upper contribution",
        "alpha_total":.05,"alpha_each":.05/(4*len(det["h_grid"])),
        "proposal":"each model separately: .5 uniform + .5 mean_j normalized "
                   "model block-mixture tilt; independent phase per block",
        "lower_and_upper_cap":"min(2lambda,2J C_Z), upper from Markov",
        "deterministic":det,
        "memory_estimate_bytes":{str(g):estimated_peak(src,det["cache_bytes"],g)
                                  for g in [1,8]},
        "cost_gate":{"max_estimated_total_seconds":180,
                     "prediction_seconds":[90,130],
                     "max_major_memory_bytes":1_000_000_000,
                     "no_extra_seed_or_N":True},
        "scope":"finite-grid complete weak-profile screen only; ordinary "
                "floating point is not interval certified; no original C_n/log claim"}
    assert max(entry["memory_estimate_bytes"].values())<1_000_000_000
    output=HERE/"profile_registration.json"
    if output.exists():raise RuntimeError("Registration already exists.")
    output.write_text(json.dumps(entry,indent=2))
    print(json.dumps({"logH2":det["log_H2"],"H2":det["H2"],
                      "h_count":len(det["h_grid"]),
                      "last_h":det["h_grid"][-1],
                      "memory":entry["memory_estimate_bytes"],
                      "preflight_seconds":det["deterministic_preflight_seconds"]}),
          flush=True)


def draw_batch(src,means,g,rng,batch):
    counts=np.empty((batch,g,src.q),dtype=float)
    block_n=src.n//g
    uniform=np.full(src.q,1/src.q)
    rows=[]
    for s in range(batch):
        if rng.random()<.5:
            for h in range(g):counts[s,h]=rng.multinomial(block_n,uniform)
            rows.append({"uniform":True});continue
        li=int(rng.integers(len(src.lengths)))
        phases=rng.integers(R,size=g)
        for h,r in enumerate(phases):
            m=np.roll(means[li],int(r)*(src.q//R))
            counts[s,h]=rng.multinomial(block_n,m/m.sum())
        rows.append({"uniform":False,"j":li,"phase":phases.tolist()})
    assert np.all(counts.sum(axis=2)==block_n)
    return counts,rows


def kl(x,mu):
    total=0.
    if x:total+=x*math.log(x/mu)
    if x<1:total+=(1-x)*math.log((1-x)/(1-mu))
    return total


def ci_lower(mean,bound,n,alpha):
    x=min(1.,max(0.,mean/bound))
    if x==0:return 0.
    target=math.log(1/alpha)/n
    low,high=0.,x
    for _ in range(80):
        mid=(low+high)/2
        if kl(x,mid)>target:low=mid
        else:high=mid
    return low*bound


def ci_upper(mean,bound,n,alpha):
    x=min(1.,max(0.,mean/bound))
    if x==1:return bound
    target=math.log(1/alpha)/n
    low,high=x,1.
    for _ in range(80):
        mid=(low+high)/2
        if kl(x,mid)>target:high=mid
        else:low=mid
    return high*bound


def accumulate(summary,z):
    summary["sum"]+=z.sum(axis=0)
    summary["squares"]+=(z*z).sum(axis=0)
    summary["max"]=np.maximum(summary["max"],z.max(axis=0))
    summary["active"]+=np.count_nonzero(z,axis=0)


def run_model(src,means,transition,norm,det,g,seed,alpha):
    started=time.time();jcount=len(src.lengths)
    hgrid=np.array(det["h_grid"]);lambdas=np.exp(hgrid);hcount=len(hgrid)
    rng=np.random.default_rng(seed)
    summaries=[{"sum":np.zeros(hcount),"squares":np.zeros(hcount),
                "max":np.zeros(hcount),"active":np.zeros(hcount,dtype=int)}
                for _ in [0,1]]
    mass_values=[];latent_rows=[]
    cap=2*jcount*det["C_Z"]
    block_n=src.n//g
    for batch_id,start in enumerate(range(0,N,BATCH)):
        batch=min(BATCH,N-start)
        counts,latent=draw_batch(src,means,g,rng,batch)
        latent_rows.extend(latent)
        flat=counts.reshape(batch*g,src.q)
        z=np.empty((batch,g,jcount,R))
        e=np.empty_like(z)
        for li in range(jcount):
            m=np.column_stack([np.roll(means[li],r*(src.q//R))for r in range(R)])
            t=np.column_stack([np.roll(transition[li],r*(src.q//R))
                                for r in range(R)])
            z[:,:,li,:]=(flat@np.log(m)).reshape(batch,g,R)
            nt=(flat@t.astype(float)).reshape(batch,g,R)
            e[:,:,li,:]=block_n*math.log1p(src.t0) \
                +nt*(math.log1p(src.t1)-math.log1p(src.t0))
        del counts,flat
        block_log=old.lse(z,axis=3)-math.log(R)
        block_weight=np.exp(z-old.lse(z,axis=3)[:,:,:,None])
        block_rv=np.sum(block_weight*np.sqrt(np.expm1(e)),axis=3)**2
        model_log=block_log.sum(axis=1)
        model_rv=np.expm1(np.log1p(block_rv).sum(axis=1))
        normalized=z-block_n*np.log(norm)[None,None,:,:]
        log_target=(old.lse(normalized,axis=3)-math.log(R)).sum(axis=1)
        log_d=np.logaddexp(math.log(.5),math.log(.5)
                           +old.lse(log_target)-math.log(jcount))
        weights=np.exp(-log_d);assert weights.max()<=2
        mass_values.extend(weights.tolist())
        contributions_l=np.empty((batch,hcount))
        contributions_u=np.empty_like(contributions_l)
        for hi,lam in enumerate(lambdas):
            lower=old.lower_factor(model_log,model_rv,float(lam))
            upper=old.upper_factor(model_log,model_rv,float(lam))
            cl=lam*weights*lower;cu=lam*weights*upper
            bound=min(2*lam,cap)
            assert cl.max()<=bound*(1+1e-9)+1e-10
            assert cu.max()<=bound*(1+1e-9)+1e-10
            assert np.max(lower-upper)<=1e-12
            contributions_l[:,hi]=cl;contributions_u[:,hi]=cu
        accumulate(summaries[0],contributions_l)
        accumulate(summaries[1],contributions_u)
        print(json.dumps({"G":g,"batch":batch_id+1,"of":N//BATCH,
                          "seconds":time.time()-started}),flush=True)
    profile=[]
    for hi,lam in enumerate(lambdas):
        entry={"h":float(hgrid[hi]),"lambda":float(lam),
               "contribution_cap":float(min(2*lam,cap))}
        for label,summ in zip(["lower","upper"],summaries):
            mean=float(summ["sum"][hi]/N)
            variance=max(0.,float((summ["squares"][hi]-N*mean*mean)/(N-1)))
            ess=float(summ["sum"][hi]**2/summ["squares"][hi]) \
                if summ["squares"][hi]>0 else 0.
            ci=(ci_lower if label=="lower"else ci_upper)(
                mean,entry["contribution_cap"],N,alpha)
            entry[label]={"mean":mean,"SE_descriptive":math.sqrt(variance/N),
                          "ESS":ess,"active":int(summ["active"][hi]),
                          "max_observed":float(summ["max"][hi]),"CI":ci}
        entry["upper_CI_after_deterministic_clips"]=min(
            entry["upper"]["CI"],float(lam),det["H2"]/float(lam))
        profile.append(entry)
    upper_intervals=[]
    upper_intervals_raw=[]
    for k in range(hcount-1):
        ratio=lambdas[k+1]/lambdas[k]
        upper_intervals_raw.append(ratio*profile[k]["upper"]["CI"])
        upper_intervals.append(min(
            ratio*profile[k]["upper_CI_after_deterministic_clips"],
            det["H2"]/float(lambdas[k])))
    low_index=int(np.argmax([x["lower"]["CI"]for x in profile]))
    empirical_index=int(np.argmax([x["lower"]["mean"]for x in profile]))
    upper_index=int(np.argmax(upper_intervals))
    mass=np.array(mass_values)
    return {"G":g,"seed":seed,"N":N,"batch":BATCH,
            "dedicated_proposal":".5 uniform + .5 all-J normalized G-block model",
            "mass_mean":float(mass.mean()),
            "mass_SE_descriptive":float(mass.std(ddof=1)/math.sqrt(N)),
            "weight_max":float(mass.max()),
            "uniform_rows":sum(x["uniform"]for x in latent_rows),
            "latent_rows":latent_rows,
            "profile":profile,
            "complete_weak_lower_CI":profile[low_index]["lower"]["CI"],
            "lower_CI_argmax_h":profile[low_index]["h"],
            "complete_weak_upper_CI":max(1.,max(upper_intervals)),
            "complete_weak_upper_CI_without_moment_clip":
                max(1.,max(upper_intervals_raw)),
            "upper_CI_argmax_interval":[profile[upper_index]["h"],
                                       profile[upper_index+1]["h"]],
            "empirical_best_lower_integral":profile[empirical_index]["lower"],
            "empirical_best_lower_h":profile[empirical_index]["h"],
            "full_range_outside_cap":1.,
            "runtime_seconds":time.time()-started}


def run():
    reg=json.loads((HERE/"profile_registration.json").read_text())
    assert sha(Path(__file__))==reg["script_SHA256"]
    for name,digest in reg["preserved_pilot_files_SHA256"].items():
        assert sha(HERE/name)==digest
    target=HERE/"profile_results.json"
    if target.exists():raise RuntimeError("Refusing an extra profile run.")
    started=time.time()
    src,grid,means,transition,norm,det=build_cache()
    assert det["h_grid"]==reg["deterministic"]["h_grid"]
    assert det["log_H2"]==reg["deterministic"]["log_H2"]
    results=[]
    for g in [1,8]:
        results.append(run_model(src,means,transition,norm,det,g,SEEDS[g],
                                 reg["alpha_each"]))
    for name,digest in reg["preserved_pilot_files_SHA256"].items():
        assert sha(HERE/name)==digest
    result={"registration_SHA256":sha(HERE/"profile_registration.json"),
            "script_SHA256":sha(Path(__file__)),
            "python":sys.executable,"python_version":sys.version,
            "numpy":np.__version__,"scope":reg["scope"],
            "source":{"K":255,"v":16,"n":src.n,"M":src.m,"B":4,
                      "p":.5,"w":4,"R":4,"J":len(grid),"Q":src.q,
                      "field":src.field,"original_matched_source":False},
            "deterministic":det,"alpha_total":reg["alpha_total"],
            "alpha_each":reg["alpha_each"],"models":results,
            "pilot_hashes_unchanged":True,
            "runtime_seconds":time.time()-started}
    target.write_text(json.dumps(result,separators=(",",":")))
    print(json.dumps({"DONE":True,"seconds":result["runtime_seconds"],
                      "bands":[{"G":x["G"],"lower":x["complete_weak_lower_CI"],
                                "upper":x["complete_weak_upper_CI"],
                                "best_lower_MC":x["empirical_best_lower_integral"]}
                               for x in results]}),flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--register",action="store_true")
    args=parser.parse_args()
    register()if args.register else run()
