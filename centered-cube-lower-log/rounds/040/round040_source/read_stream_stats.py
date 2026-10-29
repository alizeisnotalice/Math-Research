"""Read-only saved-stream payload audit. No source import, cache, or RNG."""
from pathlib import Path
from decimal import Decimal, localcontext
import hashlib
import json
import math
import sys
import time
import numpy as np

HERE=Path(__file__).resolve().parent
INPUT=HERE.parent/"round040_stream"


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def independent_ci(mean,cap,n,alpha,upper):
    with localcontext()as ctx:
        ctx.prec=55
        zero,one=Decimal(0),Decimal(1)
        cap=Decimal(str(cap));x=Decimal(str(mean))/cap;x=max(zero,min(one,x))
        if upper and x==one:return float(cap)
        if not upper and x==zero:return 0.
        target=-Decimal(str(alpha)).ln()/n
        lo,hi=(x,one)if upper else(zero,x)
        for _ in range(120):
            p=(lo+hi)/2
            d=(x*(x/p).ln()if x else zero)+((one-x)*((one-x)/(one-p)).ln()if x<one else zero)
            if upper:
                if d>target:hi=p
                else:lo=p
            else:
                if d>target:lo=p
                else:hi=p
        return float(cap*(lo+hi)/2)


def close(a,b,tolerance=1e-9):
    assert abs(a-b)<=tolerance*max(1,abs(a),abs(b)),(a,b)


def main():
    target=HERE/"STREAM_STATS_AUDIT.json"
    if target.exists():raise RuntimeError("Refusing replacement of readonly audit.")
    started=time.time();plan=json.loads((INPUT/"PLAN.json").read_text())
    receipt=json.loads((INPUT/"RUN_RECEIPT.json").read_text())
    planhash=sha(INPUT/"PLAN.json");assert receipt["plan_SHA256"]==planhash
    for key,path in [("script_SHA256",INPUT/"probe.py"),("helper_SHA256",INPUT/"helpers.py"),
                     ("source_SHA256",HERE/"source_tables.py"),("builder_SHA256",HERE/"build_cache.py")]:
        assert plan[key]==sha(path)
    family=2*len(plan["R_values"])*len(plan["models"])*len(plan["log_thresholds"])
    assert family==504==plan["CI_family_count"]
    close(plan["alpha_each"],plan["alpha_total"]/family,1e-15)
    count_ci=0;max_ci_error=max_summary_error=max_pair_error=0.;cases=[]
    for R,rece in zip(plan["R_values"],receipt["cases"]):
        resultfile=INPUT/f"R{R}_RESULTS.json";d=json.loads(resultfile.read_text())
        assert d["R"]==R==rece["R"] and sha(resultfile)==rece["results_SHA256"]
        assert d["plan_SHA256"]==planhash and d["N"]==plan["N_per_R"]
        assert d["seed"]==plan["seeds"][str(R)]
        payload=INPUT/d["paired_payload"];assert sha(payload)==d["paired_payload_SHA256"]
        with np.load(payload)as data:
            h=data["h"];logd=data["logD"];components=data["component_log_sum_normalized"]
            lower=data["lower"];upper=data["upper"];tags=data["tags"]
        n=d["N"];assert lower.shape==upper.shape==(3,n,21)and tags.shape==(n,2)
        assert h.tolist()==plan["log_thresholds"]
        assert np.isfinite(lower).all()and np.isfinite(upper).all()
        assert np.min(lower)>=0 and np.max(lower-upper)<1e-8
        reconstructed=np.logaddexp(math.log(.5),math.log(.5)-math.log(R*d["J"])+components)
        assert np.array_equal(logd,reconstructed)
        w=np.exp(-logd)
        for a,b in [(float(w.mean()),d["mass_weight_mean"]),
                    (float(w.std(ddof=1)/math.sqrt(n)),d["mass_weight_SE_descriptive"]),
                    (float(w.max()),d["weight_max"])]:close(a,b,1e-13)
        assert w.max()<=2+1e-12
        uniform=(tags[:,0]<0)
        assert np.all(tags[uniform]==-1)and int(uniform.sum())==d["uniform_rows"]
        assert np.all((tags[~uniform,0]>=0)&(tags[~uniform,0]<R))
        assert np.all((tags[~uniform,1]>=0)&(tags[~uniform,1]<d["J"]))
        capall=2*R*d["J"]*d["C_Z"]*d["C_E"];close(capall,d["cap_all"],1e-14)
        models={};ciarrays=np.empty((2,3,21))
        for mi,name in enumerate(plan["models"]):
            rows=d["models"][name]["profile"]
            assert [x["h"]for x in rows]==h.tolist()
            for hi,row in enumerate(rows):
                cap=min(2*math.exp(float(h[hi])),capall);close(cap,row["cap"],1e-14)
                for direction,(key,values,isupper)in enumerate([
                    ("lower",lower[mi,:,hi],False),("upper",upper[mi,:,hi],True)]):
                    mean=float(values.mean());se=float(values.std(ddof=1)/math.sqrt(n))
                    total=float(values.sum());square=float(values@values)
                    ess=total**2/square if square else 0.
                    for field,value in [("mean",mean),("SE_descriptive",se),
                                        ("ESS",ess),("max_contribution",float(values.max()))]:
                        error=abs(value-row[key][field]);max_summary_error=max(max_summary_error,error)
                        close(value,row[key][field],1e-11)
                    assert values.max()<=cap*(1+1e-8)
                    ci=independent_ci(mean,cap,n,plan["alpha_each"],isupper)
                    error=abs(ci-row[key]["CI"]);max_ci_error=max(max_ci_error,error)
                    close(ci,row[key]["CI"],1e-10);ciarrays[direction,mi,hi]=ci;count_ci+=1
            best=int(np.argmax(lower[mi].mean(axis=0)))
            assert rows[best]==d["models"][name]["best_lower_MC"]
            close(max(1,float(ciarrays[0,mi].max())),d["models"][name]["best_lower_CI_with_L1_anchor"],1e-10)
            models[name]={"best_h":float(h[best]),"best_lower_MC":rows[best]["lower"]["mean"],
                          "same_h_upper_MC":rows[best]["upper"]["mean"],
                          "best_lower_CI_with_L1_anchor":d["models"][name]["best_lower_CI_with_L1_anchor"]}
        for hi,pair in enumerate(d["paired_comparisons"]):
            assert pair["h"]==float(h[hi])
            delta=lower[1,:,hi]-upper[2,:,hi]
            calc={"joint_lower_minus_prior_upper_MC":float(delta.mean()),
                  "paired_SE_descriptive":float(delta.std(ddof=1)/math.sqrt(n)),
                  "simultaneous_CI_endpoint_gap":float(ciarrays[0,1,hi]-ciarrays[1,2,hi])}
            for key,value in calc.items():
                max_pair_error=max(max_pair_error,abs(value-pair[key]));close(value,pair[key],1e-10)
        bestpair=max(d["paired_comparisons"],key=lambda x:x["joint_lower_minus_prior_upper_MC"])
        cases.append({"R":R,"N":n,"payload_SHA256_verified":True,"result_SHA256_verified":True,
            "logD_exactly_reconstructed_from_all_component_sum":True,
            "weight_mean_recomputed":float(w.mean()),"weight_SE_recomputed":float(w.std(ddof=1)/math.sqrt(n)),
            "weight_max_recomputed":float(w.max()),"uniform_rows_recomputed":int(uniform.sum()),
            "tilted_b_counts":np.bincount(tags[~uniform,0],minlength=R).tolist(),
            "all_saved_pairs_use_same_h":True,"models":models,"best_descriptive_pair":bestpair})
    assert count_ci==504
    out={"scope":"saved payload audit only; zero source evaluation, cache builds or samples",
        "script_SHA256":sha(Path(__file__)),"plan_SHA256":planhash,
        "python":sys.executable,"numpy":np.__version__,"all_hash_chains_valid":True,
        "one_sided_CIs_independently_recomputed":count_ci,"max_CI_absolute_error":max_ci_error,
        "max_saved_mean_SE_ESS_max_absolute_error":max_summary_error,
        "max_paired_summary_absolute_error":max_pair_error,"cases":cases,
        "seconds":time.time()-started,
        "limitations":["raw counts/kernel statistics not regenerated","FP premise still not interval certified",
                       "paired SE is descriptive; no new paired CI","extra L1 anchor is outside bare natural grid"]}
    target.write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(out,indent=2))


if __name__=="__main__":main()
