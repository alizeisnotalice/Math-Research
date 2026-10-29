"""Independent saved-summary audit. No driver import/cache/RNG/source evaluation."""
from pathlib import Path
from decimal import Decimal, localcontext
import hashlib
import json
import math
import sys
import time

HERE=Path(__file__).resolve().parent
INPUT=HERE.parent/"round038_highdim"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def decimal_ci(mean,cap,n,alpha,upper):
    # Independent high-precision KL inversion, rather than importing driver.
    with localcontext() as ctx:
        ctx.prec=55
        zero,one=Decimal(0),Decimal(1)
        cap=Decimal(str(cap));x=Decimal(str(mean))/cap
        x=max(zero,min(one,x));target=-(Decimal(str(alpha))).ln()/Decimal(n)
        if upper and x==one:return float(cap)
        if not upper and x==zero:return 0.
        lo,hi=(x,one)if upper else(zero,x)
        for _ in range(120):
            p=(lo+hi)/2
            div=(x*(x/p).ln()if x else zero)+((one-x)*((one-x)/(one-p)).ln()if x<one else zero)
            if upper:
                if div>target:hi=p
                else:lo=p
            else:
                if div>target:lo=p
                else:hi=p
        return float(cap*(lo+hi)/2)


def main():
    started=time.time()
    plan=json.loads((INPUT/"HIGH_DIM_PLAN.json").read_text())
    reg=json.loads((INPUT/"HIGH_DIM_REGISTRATION.json").read_text())
    result=json.loads((INPUT/"HIGH_DIM_RESULTS.json").read_text())
    assert plan["script_SHA256"]==sha(INPUT/"probe_highdim.py")
    assert plan["source_SHA256"]==sha(INPUT/"source_tables.py")
    assert plan["builder_SHA256"]==sha(INPUT/"build_cache.py")
    assert reg["plan_SHA256"]==sha(INPUT/"HIGH_DIM_PLAN.json")
    assert result["registration_SHA256"]==sha(INPUT/"HIGH_DIM_REGISTRATION.json")
    assert result["N"]==reg["N"] and result["seed"]==plan["seed"]
    family=2*len(plan["models"])*len(plan["log_thresholds"])
    assert family==336==reg["CI_family_count"]
    assert abs(reg["alpha_each"]-reg["alpha_total"]/family)<1e-20
    n=result["N"]
    cap_all=2*reg["source"]["R"]*reg["J"]*reg["normalizer_C_Z"]*reg["kernel_mean_C_E"]
    assert abs(cap_all-reg["all_threshold_contribution_cap"])<1e-10
    max_ci_error=0.;max_ess_error=0.;checked=0;models={}
    for name in plan["models"]:
        data=result["models"][name];rows=data["profile"]
        assert [r["h"]for r in rows]==plan["log_thresholds"]
        for row in rows:
            h=row["h"];threshold=math.exp(h);cap=min(2*threshold,cap_all)
            assert row["lambda"]==threshold and abs(row["cap"]-cap)<1e-10
            assert row["lower"]["mean"]<=row["upper"]["mean"]+1e-12
            for kind,upper in [("lower",False),("upper",True)]:
                s=row[kind]
                assert 0<=s["mean"]<=s["max_contribution"]+1e-12<=cap+1e-8
                ci=decimal_ci(s["mean"],cap,n,reg["alpha_each"],upper)
                max_ci_error=max(max_ci_error,abs(ci-s["CI"]))
                assert abs(ci-s["CI"])<1e-9
                ess=(n*s["mean"]**2/(s["mean"]**2+(n-1)*s["SE_descriptive"]**2)
                     if s["mean"] else 0.)
                max_ess_error=max(max_ess_error,abs(ess-s["ESS"]))
                assert abs(ess-s["ESS"])<1e-8
                checked+=1
        winner=max(rows,key=lambda r:r["lower"]["mean"])
        ciwinner=max(rows,key=lambda r:r["lower"]["CI"])
        assert winner==data["best_lower_MC"]
        assert max(1.,ciwinner["lower"]["CI"])==data["best_lower_CI"]
        models[name]={"best_lower_MC_h":winner["h"],"best_lower_MC":winner["lower"]["mean"],
                      "lower_MC_SE_descriptive":winner["lower"]["SE_descriptive"],
                      "same_h_lower_CI":winner["lower"]["CI"],
                      "same_h_upper_MC":winner["upper"]["mean"],
                      "same_h_upper_CI":winner["upper"]["CI"],
                      "best_lower_CI_h":ciwinner["h"],
                      "best_lower_CI_with_L1_anchor":data["best_lower_CI"]}
    comparisons=[]
    for r in [2,4]:
        joint=result["models"]["joint"+str(r)];prior=result["models"]["prior"+str(r)]
        h=joint["best_lower_MC"]["h"]
        j=next(x for x in joint["profile"]if x["h"]==h)
        p=next(x for x in prior["profile"]if x["h"]==h)
        comparisons.append({"R":r,"h_selected_by_joint_lower_MC":h,
            "joint_lower_MC_minus_prior_lower_MC_only":j["lower"]["mean"]-p["lower"]["mean"],
            "joint_lower_MC_minus_prior_upper_MC":j["lower"]["mean"]-p["upper"]["mean"],
            "joint_lower_CI_minus_prior_upper_CI":j["lower"]["CI"]-p["upper"]["CI"],
            "paired_difference_SE_recomputable_from_summary":False})
    out={"scope":"saved-summary recomputation only; zero new samples/cache builds",
         "source_input_directory":str(INPUT),"input_hash_chain_valid":True,
         "input_results_SHA256":sha(INPUT/"HIGH_DIM_RESULTS.json"),
         "audit_script_SHA256":sha(Path(__file__)),"python":sys.executable,
         "checked_one_sided_CIs":checked,"max_CI_absolute_recompute_error":max_ci_error,
         "max_ESS_absolute_recompute_error":max_ess_error,
         "family_count":family,"alpha_each":reg["alpha_each"],
         "all_threshold_contribution_cap_recomputed":cap_all,
         "all_stored_max_contributions_below_registered_caps":True,
         "models":models,"joint_prior_same_threshold_comparisons":comparisons,
         "limitations":["saved summaries cannot reconstruct sample RN weights or paired covariances",
                        "mean/BE enclosure validity audited separately",
                        "ordinary FP premise is not interval certification",
                        "L1 anchor augments the natural window; does not belong to its bare grid"],
         "seconds":time.time()-started}
    target=HERE/"HIGH_DIM_SAVED_STATS_AUDIT.json"
    if target.exists():raise RuntimeError("Refusing replacement of audit output.")
    target.write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps({k:v for k,v in out.items()if k not in ["limitations","source_input_directory"]},indent=2))


if __name__=="__main__":main()
