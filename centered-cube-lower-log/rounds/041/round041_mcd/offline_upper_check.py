"""Independent offline upper audit; never imports upper_math/source/driver.

Run only after root accepts RECONSTRUCTION_RECEIPT. Reads saved upper inputs,
old40 arrays and replay identity files; no samples or kernel evaluations.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import math
import sys
import time
import numpy as np

HERE=Path(__file__).resolve().parent
REPLAY=HERE.parent/"round041_replay"
OLD=HERE.parent/"round040_stream"
C=float(F(81,20))


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def cantelli_scale(mu,v,h,variance_factor):
    """Independent log-odds formula, not upper_math's linear rv fraction."""
    answer=np.ones_like(mu)
    active=mu<h
    if np.any(active):
        # lambda-u = exp(h)*[-expm1(log(u)-h)], retaining close gaps.
        with np.errstate(divide="ignore",invalid="ignore",over="ignore",under="ignore"):
            log_gap=h+np.log(-np.expm1(mu[active]-h))
            log_odds=2*log_gap-v[active]-math.log(variance_factor)
            answer[active]=np.exp(-np.logaddexp(0.,log_odds))
    # The independent Markov guard includes zero upper-mean cases.
    with np.errstate(under="ignore"):
        markov=np.exp(np.minimum(0.,mu-h))
    return np.minimum(answer,markov)


def mcd_scale(mu,log_product,h,terms):
    answer=np.ones_like(mu)
    active=mu<h
    answer[active]=0. # S=0 and positive gap: deterministic empty tail.
    s=C*log_product
    positive=active&(s>0)
    with np.errstate(under="ignore",over="ignore",invalid="ignore"):
        answer[positive]=np.minimum(1.,np.exp(math.log(terms)
            -2*(h-mu[positive])**2/s[positive]))
    # At equality the frozen implementation conservatively uses 1,
    # including S=0. This is safe for strict > and matches its convention.
    return answer


def independent_uppers(log_mean,log_v,max_cost,fixed_cost,h,R):
    old=[];variance=[];new=[]
    for m in range(3):
        a=cantelli_scale(log_mean[m],log_v[m],h,1.)
        b=cantelli_scale(log_mean[m],log_v[m],h,1/3)
        d=mcd_scale(log_mean[m],fixed_cost if m==0 else max_cost,h,1 if m==0 else R)
        old.append(np.clip(a.sum(axis=-1),0,1))
        variance.append(np.clip(b.sum(axis=-1),0,1))
        new.append(np.clip(np.minimum(np.minimum(a,b),d).sum(axis=-1),0,1))
    return np.array(old),np.array(variance),np.array(new)


def close_arrays(a,b,name,errors):
    assert a.shape==b.shape,(name,a.shape,b.shape)
    difference=float(np.max(np.abs(a-b)))
    errors[name]=max(errors.get(name,0.),difference)
    assert np.allclose(a,b,rtol=2e-12,atol=2e-9,equal_nan=False),(name,difference)


def main():
    target=HERE/"OFFLINE_UPPER_AUDIT.json"
    if target.exists():raise RuntimeError("Refusing repeated audit output.")
    receipt_path=REPLAY/"RECONSTRUCTION_RECEIPT.json"
    if not receipt_path.exists():raise RuntimeError("Wait for root's accepted reconstruction receipt.")
    receipt=json.loads(receipt_path.read_text())
    assert receipt["all_acceptance_checks_passed"] and receipt["same_protocol_replay"]
    assert receipt["new_independent_samples"]==0
    plan=json.loads((REPLAY/"PLAN.json").read_text());oldplan=json.loads((OLD/"PLAN.json").read_text())
    oldreceipt=json.loads((OLD/"RUN_RECEIPT.json").read_text())
    assert receipt["plan_SHA256"]==sha(REPLAY/"PLAN.json")
    assert plan["old_plan_SHA256"]==sha(OLD/"PLAN.json")
    assert plan["base_driver_SHA256"]==sha(OLD/"probe.py")
    assert plan["observer_SHA256"]==sha(REPLAY/"reconstruct.py")
    assert plan["upper_math_SHA256"]==sha(REPLAY/"upper_math.py")
    assert plan["MCD_proof_SHA256"]==sha(HERE/"PROOF.md")
    assert plan["C_log_bound"]==C and F(plan["source_r0"])==F(1,10**6)
    assert F(81,20)>F(1000001,249001)
    started=time.time();cases=[];global_errors={};all_rows=0
    for accepted,oldrece in zip(receipt["cases"],oldreceipt["cases"]):
        R=accepted["R"];assert R==oldrece["R"]
        assert all(accepted[k]for k in ["count_stream_matches","sampled_anchor_stream_matches",
            "source_statistics_stream_matches","all_old_payload_arrays_identical"])
        assert accepted["max_checked_source_r"]<=1e-6
        old_result_path=OLD/f"R{R}_RESULTS.json"
        assert sha(old_result_path)==oldrece["results_SHA256"]
        original=json.loads(old_result_path.read_text())
        replayed=json.loads((REPLAY/"base_replay"/f"R{R}_RESULTS.json").read_text())
        for key in ["count_stream_SHA256","sampled_anchor_nominal_stream_SHA256","stream_statistics_SHA256"]:
            assert original[key]==replayed[key],(R,key)
        with np.load(OLD/original["paired_payload"],allow_pickle=False)as a,\
             np.load(REPLAY/"base_replay"/replayed["paired_payload"],allow_pickle=False)as b:
            assert a.files==b.files
            for key in a.files:assert np.array_equal(a[key],b[key]),(R,key)
            old_lower=a["lower"].copy();old_upper=a["upper"].copy()
            old_h=a["h"].copy();old_d=a["logD"].copy()
        assert sha(OLD/original["paired_payload"])==original["paired_payload_SHA256"]
        input_path=REPLAY/accepted["upper_inputs_file"]
        assert sha(input_path)==accepted["upper_inputs_SHA256"]
        with np.load(input_path,allow_pickle=False)as data:
            mu=data["logmean_upper"];v=data["logvariance_old"]
            maximum=data["max_table_logproduct"];fixed=data["table0_logproduct"]
            h=data["h"];logd=data["logD"]
            saved_new=data["new_upper"];saved_var=data["variance_only_upper"]
        n=original["N"];j=original["J"]
        assert mu.shape==v.shape==(3,n,j)and maximum.shape==fixed.shape==(n,j)
        assert np.array_equal(h,old_h)and np.array_equal(logd,old_d)
        assert np.all(maximum>=fixed-1e-14)and np.min(fixed)>=0
        errors={};violations={k:0. for k in ["old_lower_to_phi","phi_to_new","new_to_variance","variance_to_old","caps"]}
        max_B=0.
        for hi,hh in enumerate(h):
            a,b,c=independent_uppers(mu,v,maximum,fixed,float(hh),R)
            multiplier=np.exp(hh-logd)[None,:]
            computed_old=a*multiplier;computed_var=b*multiplier;computed_new=c*multiplier
            close_arrays(computed_old,old_upper[:,:,hi],"old_upper",errors)
            close_arrays(computed_var,saved_var[:,:,hi],"variance_only_upper",errors)
            close_arrays(computed_new,saved_new[:,:,hi],"new_upper",errors)
            B=old_lower[:,:,hi]/multiplier
            assert B.min()>=0 and B.max()<=1+1e-10
            max_B=max(max_B,float(B.max()))
            sharp=multiplier*(3*B/(1+2*B))
            for name,left,right in [("old_lower_to_phi",old_lower[:,:,hi],sharp),
                ("phi_to_new",sharp,saved_new[:,:,hi]),("new_to_variance",saved_new[:,:,hi],saved_var[:,:,hi]),
                ("variance_to_old",saved_var[:,:,hi],old_upper[:,:,hi])]:
                violation=float(np.max(left-right));violations[name]=max(violations[name],violation)
                assert violation<=2e-9,(R,hh,name,violation)
            cap=min(2*math.exp(float(hh)),original["cap_all"])
            violation=max(float(sharp.max()),float(saved_new[:,:,hi].max()),float(saved_var[:,:,hi].max()))-cap
            violations["caps"]=max(violations["caps"],violation)
            assert violation<=2e-8
            all_rows+=3*n*j
        for key,val in errors.items():global_errors[key]=max(global_errors.get(key,0.),val)
        cases.append({"R":R,"N":n,"J":j,"upper_input_SHA256":sha(input_path),
            "three_stream_hashes_independently_match":True,"all_old_arrays_independently_bit_identical":True,
            "receipt_max_source_r":accepted["max_checked_source_r"],"maximum_recovered_conditional_B":max_B,
            "max_absolute_recompute_errors":errors,"max_chain_or_cap_violations":violations,
            "all_model_sample_scale_threshold_rows_checked":3*n*j*len(h)})
    out={"scope":"independent offline formulas and saved-data identity; no source/RNG/upper_math imports",
        "receipt_SHA256":sha(receipt_path),"script_SHA256":sha(Path(__file__)),
        "accepted_same_sample_receipt_verified_directly":True,"C":C,"r0":"1/1000000",
        "all_rows_checked":all_rows,"max_absolute_recompute_errors":global_errors,"cases":cases,
        "python":sys.executable,"numpy":np.__version__,"seconds":time.time()-started,
        "qualification":"mathematical r0/C and kernel/BE premises retained; complete FP not interval certified; no CI recomputation here"}
    target.write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(out,indent=2))


if __name__=="__main__":main()
