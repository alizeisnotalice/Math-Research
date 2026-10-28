"""One registered, independent finite-grid single-tail confirmation.

Does not modify the prior probe or results. No seed search. Output refuses
overwrite; --output may be used for an explicitly authorized archival rerun.
"""
from pathlib import Path
from fractions import Fraction
import argparse
import hashlib
import json
import math
import sys
import time
import numpy as np
import natural_probe as base

HERE=Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def log_e(d,theta,b,rho):
    assert -b<theta<b
    factors=1+rho*(d-theta)/(b+theta)
    assert np.all(factors>0)
    return math.fsum(np.log(factors).tolist())


def inverse_e(d,b,rho,alpha):
    target=math.log(1/alpha)
    low=float(np.nextafter(-b,b))
    high=float(np.nextafter(b,-b))
    if log_e(d,low,b,rho)<target:
        return -b
    if log_e(d,high,b,rho)>target:
        return b
    for _ in range(100):
        mid=(low+high)/2
        if log_e(d,mid,b,rho)>target:
            low=mid
        else:
            high=mid
    # Conservative endpoint of the monotone-root bracket.
    return low


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",default="confirm_results.json")
    args=parser.parse_args()
    target=HERE/args.output
    if target.exists():
        raise RuntimeError("Refusing to overwrite a prior confirmation result.")
    reg=json.loads((HERE/"confirm_registration.json").read_text())
    assert sha(Path(__file__))==reg["confirm_script_SHA256"]
    for name,digest in reg["preserved_prior_files_SHA256"].items():
        assert sha(HERE/name)==digest
    started=time.time()
    k,v,rmax,n_samples=255,8,4,reg["N"]
    src=base.Source(8,0x11d,v,rmax)
    natural={Fraction(10,9)+Fraction(j,k*v) for j in range(k*v//9+1)}
    old={Fraction(80+j,72) for j in range(9)}
    union=sorted(natural|old)
    assert len(natural)==227 and len(natural&old)==3 and len(union)==233
    assert union[0]==Fraction(10,9) and union[-1]==Fraction(11,9)
    assert src.q==48960 and src.n==6704800
    shifts=[length*src.q for length in union]
    assert all(x.denominator==1 for x in shifts)
    src.lengths=np.array([float(x)for x in union])
    src.length_q=np.array([int(x)for x in shifts],dtype=np.int64)
    jcount=len(union)
    rng=np.random.default_rng(reg["seed"])
    counts,audit=base.draw_counts(src,rng,n_samples)
    log_a=np.empty((n_samples,jcount,rmax))
    rel_a=np.empty_like(log_a)
    norms=[];sum_order=[]
    for li in range(jcount):
        means,transition=src.query(li)
        norm=np.array([np.ascontiguousarray(means[:,r]).mean()
                       for r in range(rmax)])
        norms.append(norm)
        sum_order.append(src.n*np.max(np.abs(np.log(means.mean(axis=0))
                                             -np.log(norm))))
        log_a[:,li]=counts@np.log(means)
        nt=counts@transition.astype(float)
        rel_a[:,li]=np.expm1(src.n*math.log1p(src.t0)
            +nt*(math.log1p(src.t1)-math.log1p(src.t0)))
        if li in [0,jcount//2,jcount-1]:
            print(json.dumps({"completed_L":li+1,"J":jcount,
                              "seconds":time.time()-started}),flush=True)
    norms=np.array(norms)
    normalized=log_a-src.n*np.log(norms)[None,:,:]
    log_d=np.logaddexp(math.log(.5),math.log(.5)
        +base.lse(normalized.reshape(n_samples,-1))-math.log(jcount*rmax))
    weights=np.exp(-log_d)
    assert weights.max()<=2
    log_m=base.lse(log_a,axis=2)-math.log(rmax)
    mix_weight=np.exp(log_a-base.lse(log_a,axis=2)[:,:,None])
    r_m=np.sum(mix_weight*np.sqrt(rel_a),axis=2)**2
    lam=math.exp(2)
    b_shared=base.cantelli(log_m,r_m,lam).max(axis=1)
    u_product=base.product_upper(log_a[:,:,0],rel_a[:,:,0],lam)
    d=lam*weights*(b_shared-u_product)
    b=2*lam
    assert d.min()>=-b and d.max()<=b
    rho,alpha=reg["rho"],reg["alpha"]
    log_e0=log_e(d,0.,b,rho)
    lower=inverse_e(d,b,rho,alpha)
    mean=float(d.mean())
    lower_hoeffding=mean-2*b*math.sqrt(math.log(1/alpha)/(2*n_samples))
    c_z=math.exp(max(0.,float(src.n*np.log(norms).max())))
    contribution=lam*weights*b_shared
    assert contribution.max()<=min(2*lam,2*jcount*c_z)*(1+1e-9)
    signed_se=float(d.std(ddof=1)/math.sqrt(n_samples))
    label_counts=audit.pop("label_counts_L_by_phase")
    audit["nonzero_label_counts_L_phase_count"]=[
        [li,r,count]for li,row in enumerate(label_counts)
        for r,count in enumerate(row)if count]
    result={
        "registration_SHA256":sha(HERE/"confirm_registration.json"),
        "script_SHA256":sha(Path(__file__)),
        "source_script_SHA256":sha(HERE/"natural_probe.py"),
        "python":sys.executable,"python_version":sys.version,
        "numpy":np.__version__,
        "seed":reg["seed"],"N":n_samples,
        "K":k,"v":v,"n":src.n,"B":4,"M":src.m,"p":.5,"w":4,
        "delta":src.delta,"R":rmax,"phases":[0,.25,.5,.75],
        "source_field":src.field,"source_background_lower":src.bg_lower,
        "original_matched_lower_bound_claimed":False,
        "lambda":lam,"b":b,"rho":rho,"alpha":alpha,
        "natural_count":len(natural),"old_count":len(old),
        "intersection_count":len(natural&old),"J_union":jcount,
        "query_lengths_exact":[str(x)for x in union],
        "Q":src.q,"C_Z":c_z,
        "max_raw_mean_mass_residual":float(np.abs(norms-1).max()),
        "max_n_log_sum_order_difference":float(max(sum_order)),
        "mass_mean":float(weights.mean()),
        "mass_SE_descriptive":float(weights.std(ddof=1)/math.sqrt(n_samples)),
        "weight_max":float(weights.max()),"proposal_audit":audit,
        "paired_D_mean":mean,"paired_D_SE_descriptive":signed_se,
        "paired_D_min":float(d.min()),"paired_D_max":float(d.max()),
        "log_e_at_theta_0":log_e0,
        "e_at_theta_0":math.exp(log_e0)if log_e0<700 else None,
        "log_e_rejection_cutoff":math.log(1/alpha),
        "e_inverted_one_sided_lower":lower,
        "log_e_at_inverted_lower":log_e(d,lower,b,rho),
        "Hoeffding_one_sided_lower":lower_hoeffding,
        "e_value_finite_grid_confirmation_positive":bool(lower>0),
        "shared_Cantelli_coefficient_MC":float(contribution.mean()),
        "product_conditional_upper_coefficient_MC":
            float(np.mean(lam*weights*u_product)),
        "max_shared_coefficient_contribution":float(contribution.max()),
        "shared_lower_contribution_bound":min(2*lam,2*jcount*c_z),
        "product_upper_contribution_bound":2*lam,
        "primary_D":d.tolist(),
        "interpretation":"Single frozen finite-grid, single-tail paired integral "
            "difference. Statistical e-value statement assumes exact computed "
            "integrands. Floating-point certification remains open. Not a "
            "weak-norm, C_n, or logarithmic asymptotic improvement.",
        "runtime_seconds":time.time()-started,
        "preserved_prior_hashes_verified":True,
    }
    for name,digest in reg["preserved_prior_files_SHA256"].items():
        assert sha(HERE/name)==digest
    target.write_text(json.dumps(result,separators=(",",":")))
    print(json.dumps({key:result[key]for key in [
        "runtime_seconds","mass_mean","paired_D_mean","paired_D_SE_descriptive",
        "log_e_at_theta_0","e_inverted_one_sided_lower",
        "Hoeffding_one_sided_lower"]}),flush=True)


if __name__=="__main__":
    main()
