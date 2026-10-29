"""Independent factor6 full-input arithmetic/statistic check, no root imports/CI inversion."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import math
import sys
import time
import numpy as np

HERE=Path(__file__).resolve().parent
POST=HERE.parent/"round041_post"
REPLAY=HERE.parent/"round041_replay"
OLD=HERE.parent/"round040_stream"
C=float(F(81,20))


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def variance_third_upper(mu,lv,h):
    # Different log-odds evaluation from root's rv/(rv+gap^2).
    out=np.ones_like(mu);active=mu<h
    with np.errstate(divide="ignore",invalid="ignore",under="ignore",over="ignore"):
        loggap=h+np.log(-np.expm1(mu[active]-h))
        out[active]=np.exp(-np.logaddexp(0.,2*loggap-lv[active]+math.log(3)))
        markov=np.exp(np.minimum(0.,mu-h))
    return np.minimum(out,markov)


def uniform_six(mu,lv,maxlog,firstlog,h,R):
    result=[];floor_count=0
    for mi in range(3):
        s=C*(firstlog if mi==0 else maxlog)
        active=mu[mi]<h
        positive=active&(s>0)
        exponential=np.ones_like(s)
        exponential[active]=0. # Exact S0, positive gap.
        gap=h-mu[mi][positive]
        raw=(-6*gap*gap)/s[positive]
        floor_count+=int(np.count_nonzero(raw<-700))
        # Exponent is floored BEFORE multiplying the table-union factor.
        terms=1 if mi==0 else R
        exponent=math.log(terms)+np.maximum(raw,-700.)
        exponential[positive]=np.exp(np.minimum(0.,exponent))
        per=np.minimum(variance_third_upper(mu[mi],lv[mi],h),exponential)
        result.append(np.minimum(1.,per.sum(axis=-1)))
    return np.array(result),floor_count


def assert_close(a,b,absolute=2e-9,relative=2e-12):
    assert abs(a-b)<=absolute+relative*max(abs(a),abs(b)),(a,b)


def moments(values):
    n=len(values);total=float(values.sum());squared=float(values@values)
    return {"mean":float(values.mean()),"SE_descriptive":float(values.std(ddof=1)/math.sqrt(n)),
            "ESS":total*total/squared if squared else 0.,"max_contribution":float(values.max())}


def main():
    target=HERE/"uniform_offline_check.json"
    if target.exists():raise RuntimeError("Refusing replacement of offline audit.")
    started=time.time()
    plan=json.loads((POST/"FINAL_PLAN.json").read_text());final=json.loads((POST/"FINAL_RESULTS.json").read_text())
    receipt=json.loads((REPLAY/"RECONSTRUCTION_RECEIPT.json").read_text())
    assert final["plan_SHA256"]==sha(POST/"FINAL_PLAN.json")
    for key,path in [("script_SHA256",POST/"aggregate.py"),("uniform_upper_SHA256",POST/"uniform_upper.py"),
        ("uniform_proof_SHA256",HERE/"uniform_PROOF.md"),
        ("method_uniform_CI_proof_SHA256",HERE.parent/"round041_audit/METHOD_UNIFORM_CI.md"),
        ("reconstruction_receipt_SHA256",REPLAY/"RECONSTRUCTION_RECEIPT.json"),
        ("previous_lower_plan_SHA256",POST/"PLAN.json")]:assert plan[key]==sha(path),(key,path)
    assert receipt["all_acceptance_checks_passed"]and receipt["new_independent_samples"]==0
    assert final["new_independent_samples"]==0
    oldplan=json.loads((OLD/"PLAN.json").read_text())
    models=oldplan["models"];all_rows=0;stat_count=0;pair_count=0;floor_total=0
    max_stat_error=max_pair_error=max_pointwise_leak=0.;cases=[]
    for accepted,case in zip(receipt["cases"],final["cases"]):
        R=case["R"];assert R==accepted["R"]
        inpfile=REPLAY/accepted["upper_inputs_file"]
        assert sha(inpfile)==plan["inputs"][inpfile.name]==accepted["upper_inputs_SHA256"]
        assert accepted["max_checked_source_r"]<=1e-6
        meta=json.loads((OLD/f"R{R}_RESULTS.json").read_text())
        repeated=json.loads((REPLAY/"base_replay"/f"R{R}_RESULTS.json").read_text())
        for key in ["count_stream_SHA256","sampled_anchor_nominal_stream_SHA256","stream_statistics_SHA256"]:
            assert meta[key]==repeated[key]
        with np.load(OLD/f"R{R}_paired.npz",allow_pickle=False)as a,\
             np.load(REPLAY/"base_replay"/f"R{R}_paired.npz",allow_pickle=False)as b:
            assert a.files==b.files
            for key in a.files:assert np.array_equal(a[key],b[key])
            oldlower=a["lower"].copy();oldupper=a["upper"].copy()
            oldh=a["h"].copy();oldd=a["logD"].copy()
        with np.load(inpfile,allow_pickle=False)as data:
            mu=data["logmean_upper"];lv=data["logvariance_old"]
            sm=data["max_table_logproduct"];s0=data["table0_logproduct"]
            h=data["h"];logd=data["logD"];factor2=data["new_upper"];varonly=data["variance_only_upper"]
        assert np.array_equal(h,oldh)and np.array_equal(logd,oldd)
        n=case["N"];j=mu.shape[-1]
        assert mu.shape==lv.shape==(3,n,j)and sm.shape==s0.shape==(n,j)
        assert h.tolist()==oldplan["log_thresholds"]
        multiplier=np.exp(h[None,None,:]-logd[None,:,None])
        B=oldlower/multiplier;assert B.min()>=0 and B.max()<=1+1e-12
        B=np.clip(B,0,1);sharp=multiplier*3*B/(1+2*B)
        upper=np.empty_like(oldlower);floors=0
        max_chain={name:0. for name in ["lower_to_upper","upper_to_variance","variance_to_old","factor2_floor_leak","caps"]}
        for hi,hh in enumerate(h):
            u,fc=uniform_six(mu,lv,sm,s0,float(hh),R);floors+=fc
            upper[:,:,hi]=u*multiplier[:,:,hi]
            violations=[("lower_to_upper",sharp[:,:,hi]-upper[:,:,hi]),
                        ("upper_to_variance",upper[:,:,hi]-varonly[:,:,hi]),
                        ("variance_to_old",varonly[:,:,hi]-oldupper[:,:,hi])]
            for name,difference in violations:
                value=float(difference.max());max_chain[name]=max(max_chain[name],value)
                assert value<=2e-9,(R,hh,name,value)
            leak=upper[:,:,hi]-factor2[:,:,hi]
            max_chain["factor2_floor_leak"]=max(max_chain["factor2_floor_leak"],float(leak.max()))
            max_pointwise_leak=max(max_pointwise_leak,float(leak.max()))
            # Upward flooring may exceed the old factor2 zero by a tiny term.
            # Bound it analytically, independently for fixed and mixture models.
            for mi in range(3):
                terms=1 if mi==0 else R
                slack=multiplier[0,:,hi]*j*terms*math.exp(-700)
                assert np.max(leak[mi]-slack)<=2e-9
            cap=min(2*math.exp(float(hh)),meta["cap_all"])
            max_chain["caps"]=max(max_chain["caps"],float(upper[:,:,hi].max())-cap)
            assert upper[:,:,hi].max()<=cap+2e-8
            all_rows+=3*n*j
        assert floors==case["upper_exp_floor_applications"]
        floor_total+=floors;model_summary={}
        for mi,name in enumerate(models):
            rows=case["models"][name]["profile"]
            assert [x["h"]for x in rows]==h.tolist()
            for hi,row in enumerate(rows):
                assert row["lambda"]==math.exp(float(h[hi]))
                for key,array in [("old_lower",oldlower),("sharp_lower",sharp),
                                  ("old_upper",oldupper),("sharp_upper",upper)]:
                    recomputed=moments(array[mi,:,hi])
                    for field,value in recomputed.items():
                        max_stat_error=max(max_stat_error,abs(value-row[key][field]))
                        assert_close(value,row[key][field])
                    stat_count+=1
                assert_close(float(factor2[mi,:,hi].mean()),row["mcd2_upper_mean_descriptive"])
                assert_close(float(varonly[mi,:,hi].mean()),row["variance_only_upper_mean_descriptive"])
            best=int(np.argmax(sharp[mi].mean(axis=0)))
            assert rows[best]==case["models"][name]["best_sharp_lower"]
            model_summary[name]={"best_h":float(h[best]),"sharp_lower":float(sharp[mi,:,best].mean()),
                "same_h_sharp_upper":float(upper[mi,:,best].mean())}
        positive=0;positive_endpoints=0
        for hi,pair in enumerate(case["pairs"]):
            assert pair["h"]==float(h[hi])
            delta=sharp[1,:,hi]-upper[2,:,hi]
            values={"sharp_joint_lower_minus_sharp_prior_upper":float(delta.mean()),
                "paired_SE_descriptive":float(delta.std(ddof=1)/math.sqrt(n)),
                "simultaneous_endpoint_gap":case["models"]["joint"]["profile"][hi]["sharp_lower"]["CI"]
                  -case["models"]["prior"]["profile"][hi]["sharp_upper"]["CI"]}
            for key,value in values.items():
                max_pair_error=max(max_pair_error,abs(value-pair[key]));assert_close(value,pair[key])
            positive+=values["sharp_joint_lower_minus_sharp_prior_upper"]>0
            positive_endpoints+=values["simultaneous_endpoint_gap"]>0
            pair_count+=1
        assert positive==case["positive_same_h_MC"]and positive_endpoints==case["positive_same_h_endpoint_gaps"]
        cases.append({"R":R,"N":n,"J":j,"same_old40_data_verified_directly":True,
            "all_input_hashes_valid":True,"floor_count_recomputed":floors,
            "max_pointwise_chain_violations":max_chain,"models":model_summary,
            "positive_same_h_MC":int(positive),"positive_same_h_endpoint_gaps":int(positive_endpoints)})
    assert floor_total==final["upper_exp_floor_count"]and stat_count==1008 and pair_count==84
    out={"scope":"factor6 independent full-input arithmetic and saved moments/pairs; no root math imports/source/RNG/CI inversion",
        "script_SHA256":sha(Path(__file__)),"FINAL_PLAN_SHA256":sha(POST/"FINAL_PLAN.json"),
        "FINAL_RESULTS_SHA256":sha(POST/"FINAL_RESULTS.json"),"all_hash_chains_valid":True,
        "rows_recomputed":all_rows,"profile_statistics_rows_recomputed":stat_count,
        "same_h_pairs_recomputed":pair_count,"floor_count_recomputed":floor_total,
        "max_saved_statistics_absolute_error":max_stat_error,"max_pair_absolute_error":max_pair_error,
        "max_new_upper_minus_factor2":max_pointwise_leak,"cases":cases,
        "python":sys.executable,"numpy":np.__version__,"seconds":time.time()-started,
        "qualification":"r0/C/continuous-uniform Lipschitz and mean/BE premises retained; full FP not certified; CI audited elsewhere"}
    target.write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(out,indent=2))


if __name__=="__main__":main()
