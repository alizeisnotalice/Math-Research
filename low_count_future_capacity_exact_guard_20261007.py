#!/usr/bin/env python3
"""Exact RN mask-weight / union envelope / source-only entropy guards.

No actual cube/FIRST/history model. The checker epsilon is 1/4096, not a
numerically certified value of the theory's ln(n+2)^(-4). Counts here are
Bernstein diagnostic labels, not early-continuation or lower-source labels.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import comb, isqrt
from pathlib import Path
import hashlib
import json

CHECKS = 0
EPS = F(1,4096)


def check(truth, message):
    global CHECKS
    CHECKS += 1
    if not truth:
        raise AssertionError(message)


def small(x):
    x = F(x)
    return f"{x.numerator}/{x.denominator}"


def encoded(x):
    x = F(x)
    def frame(value):
        raw = abs(value).to_bytes(max(1,(abs(value).bit_length()+7)//8),"big")
        return (b"-" if value < 0 else b"+")+len(raw).to_bytes(8,"big")+raw
    out = {"numerator_bits":abs(x.numerator).bit_length(),
           "denominator_bits":x.denominator.bit_length(),
           "framed_pair_sha256":hashlib.sha256(frame(x.numerator)+frame(x.denominator)).hexdigest()}
    if max(abs(x.numerator).bit_length(),x.denominator.bit_length()) <= 256:
        out["exact"] = small(x)
    return out


@lru_cache(None)
def weight(n,alpha,sigma,k):
    rho = (1-sigma)/sigma
    term = F(1,comb(n-k+alpha,alpha))
    total = term
    for j in range(k):
        term *= F(k-j,j+1)*rho*F(alpha+j,alpha+n-k+j+1)
        check(term >= 0, "positive beta polynomial term")
        total += term
    return total


def direct_weight(n,alpha,sigma,k):
    rho = (1-sigma)/sigma
    total = F(0)
    for j in range(k+1):
        beta_term = F(alpha,(alpha+j)*comb(n-k+alpha+j,alpha+j))
        total += comb(k,j)*rho**j*beta_term
    return total


def rn_identity(n,alpha,sigma,d):
    gamma = F(3,4)
    initial = 1-sigma+sigma*gamma
    p = sigma*gamma/initial
    a = gamma/initial
    via_mask = sum((F(comb(d,k))*p**k*(1-p)**(d-k)*weight(n,alpha,sigma,k)
                    for k in range(d+1)), F(0))
    # Independent replacement expansion of
    # (1-v)^(n-d)*(1-v+v a)^d, integrated against alpha v^(alpha-1).
    via_future = sum((F(comb(d,j))*a**j*
                      F(alpha,(alpha+j)*comb(n+alpha,alpha+j))
                      for j in range(d+1)), F(0))
    check(via_mask == via_future, "full-dimension exact RN identity by independent expansion")
    check(sum((F(comb(d,k))*p**k*(1-p)**(d-k) for k in range(d+1)),F(0)) == 1,
          "reference coefficient mask normalization")
    return {"supported_soft_coordinates":d,"gamma":small(gamma),
            "current_count_p":small(p),"future_RN_ratio":encoded(via_mask),
            "scope":"finite nonnegative hardpersistent coefficient fixture, not original phi input"}


def union_guard(weights,values,S,eps,label):
    check(sum(weights,F(0)) == 1, "one original source mask mass")
    Z = sum((a*b for a,b in zip(weights,values)),F(0))
    check(S >= Z > 0, "common positive envelope covers averaged response")
    source_high = S >= eps
    new_label_high = sum((a for a,b in zip(weights,values) if b >= eps),F(0))
    old_paid = F(1) if source_high else F(0)
    new_paid = F(0) if source_high else new_label_high
    combined = old_paid+new_paid
    intermediate = F(1) if source_high else Z/eps
    check(combined <= intermediate <= S/eps, "joint source/label marginal uses one envelope fee")
    low_mass = F(0) if source_high else sum((a for a,b in zip(weights,values) if b<eps),F(0))
    check(old_paid+new_paid+low_mass == 1, "direct strict-complement label masses partition original traffic")
    gated = sum((mass*F(2+i,4) for i,mass in enumerate(weights)
                 if source_high or values[i] >= eps),F(0))
    check(gated <= combined <= S/eps, "pair/label dependent bounded gate domination")
    return {"label":label,"epsilon":encoded(eps),"S":encoded(S),"Z":encoded(Z),
            "source_high":source_high,"old_paid_mass":small(old_paid),
            "new_label_paid_mass":small(new_paid),"combined_paid_mass":small(combined)}


def one_softness(n,alpha,sigma,d):
    counts = sorted(set([0,1,2,alpha//2,alpha,2*alpha,4*alpha+64]))
    sampled = []
    last = None
    for k in counts:
        w = weight(n,alpha,sigma,k)
        if last is not None:
            check(w > last, "strict monotonicity on sampled count range")
        last = w
        if k <= 2 or k == alpha//2:
            check(w == direct_weight(n,alpha,sigma,k), "independent positive beta-sum formula")
        sampled.append({"k":k,"w":encoded(w),"w_ge_checker_epsilon":w>=EPS})
    w0,w1 = weight(n,alpha,sigma,0),weight(n,alpha,sigma,1)
    check(w0 == F(1,comb(n+alpha,alpha)), "w0 exact beta endpoint")
    check(w0 < EPS/8, "small exact lowest count weight relative to registered epsilon")
    high_k = counts[-1]
    high_w = weight(n,alpha,sigma,high_k)
    if high_w >= EPS:
        lo,hi = 0,high_k
        while hi-lo>1:
            mid = (lo+hi)//2
            if weight(n,alpha,sigma,mid) >= EPS:
                hi=mid
            else:
                lo=mid
        check(weight(n,alpha,sigma,hi) >= EPS and weight(n,alpha,sigma,hi-1) < EPS,
              "exact rational-epsilon cutoff with inclusive paid boundary")
        cutoff={"exact_within_registered_range":hi,"epsilon":small(EPS)}
        high_mass = EPS/(8*high_w)
        masses=[1-high_mass,high_mass]
        Z=sum((a*b for a,b in zip(masses,[w0,high_w])),F(0))
        check(2*Z < EPS, "low-source with positive high-label scalar guard")
        union = [union_guard(masses,[w0,high_w],2*Z,EPS,"source-low with label-high"),
                 union_guard([F(1)],[w0],EPS,EPS,"source threshold equality paid")]
    else:
        check(high_w < EPS, "registered-range strict lower cutoff bound")
        cutoff={"strictly_greater_than":high_k,"analytic_upper_bound":n,"epsilon":small(EPS)}
        union=[union_guard([F(1)],[w0],2*w0,EPS,"source-low all registered low labels"),
               union_guard([F(1)],[w0],EPS,EPS,"source threshold equality paid")]
    # Exact label equality guard at a separate registered rational epsilon_eq.
    eps_eq=w1
    Z_eq=F(3,4)*w0+F(1,4)*w1
    S_eq=(Z_eq+eps_eq)/2
    eq=union_guard([F(3,4),F(1,4)],[w0,w1],S_eq,eps_eq,"label threshold equality paid")
    check(eq["new_label_paid_mass"] == "1/4" and not eq["source_high"],
          "equality assigned label paid rather than strict complement")
    union.append(eq)
    return {"sigma":small(sigma),"sample_weights":sampled,"cutoff":cutoff,
            "independent_RN_identity":rn_identity(n,alpha,sigma,d),"joint_fee_models":union}


def endpoints(n,alpha):
    w0=F(1,comb(n+alpha,alpha))
    check(0 < w0 < 1, "sigma0 surviving hard coefficient is not total future RN")
    eta_n=F(alpha,n+alpha)
    check(w0+eta_n <= 1 and eta_n>0, "sigma0 future has positive new masks outside original support")
    for k in [0,1,alpha,n]:
        check(weight(n,alpha,F(1),k) == F(1,comb(n-k+alpha,alpha)),
              "sigma1 formal weights and allsoft support")
    check(weight(n,alpha,F(1),n) == 1, "sigma1 genuine support weight")
    tiny=F(1,n*n)
    base=F(1,comb(n-1+alpha,alpha))
    leading=base*F(alpha,n+alpha)
    check(tiny*weight(n,alpha,tiny,1)-leading == tiny*base*F(n,n+alpha),
          "sigma down to zero weight singularity exact identity")
    return {"sigma0_w0":encoded(w0),"sigma0_new_future_allsoft_prior_mass":small(eta_n),
            "sigma0_full_RN_identity":"not applicable; future has new support",
            "sigma1_genuine_K":n,"sigma1_w_n":"1/1"}


def source_entropy(n):
    rootn=isqrt(n)
    b=1<<n
    m=F(1,b)
    opposite=1-m
    p=F(1,rootn*(b-1))
    check(b*m == 1, "exact compressed 2^n-child source normalization")
    check(p*opposite == m/rootn, "source entropy fixture preserves original capacity")
    check(b*p*opposite == F(1,rootn), "one complete mass per source layer")
    inverse=p.denominator  # p numerator is one
    check(inverse >= 2**(n-1), "u lower bound (n-1)ln2 exact certificate")
    check(inverse <= rootn*(1<<n), "u upper bound .5lnn+nln2 exact certificate")
    depths=[]
    for depth in [1,2,8]:
        mass_reuse=F(0)
        for level in range(1,depth+1):
            child_mass=F(1,2**level)
            check((2**level)*child_mass == 1, "binary layer source masses sum to one")
            check(child_mass/(rootn*child_mass) == F(1,rootn), "binary child coin constant across depth")
            mass_reuse += (2**level)*child_mass
        check(mass_reuse == depth, "explicit full-tree entropy depth multiplier")
        depths.append({"depth":depth,"full_mass_layer_reuse":small(mass_reuse),
                       "exact_entropy_sum":"depth * ln(sqrt(n))"})
    return {"occupied_children":"2^n, compressed exact symmetry","mass":"1/1",
            "source_p":encoded(p),"entropy_sum":"ln(sqrt(n)*(2^n-1))",
            "entropy_lower":"(n-1)*ln2","entropy_upper":".5lnn+nln2",
            "binary_depths":depths,"scope":"source-only range certificate, not actual input/output"}


def main():
    here=Path(__file__).parent
    registration=here/"low_count_future_capacity_guard_registration_20261007.json"
    reg=json.loads(registration.read_text())
    check(reg["dimensions"] == [1024,4096,16384], "frozen registered rounds")
    check(reg["status"] == "REGISTERED_BEFORE_EXECUTION", "preregistration state")
    rounds=[]
    for n,d in zip(reg["dimensions"],[8,16,32]):
        before=CHECKS
        alpha=isqrt(n)
        check(alpha*alpha == n, "registered square dimension alpha")
        models=[one_softness(n,alpha,sigma,d) for sigma in [F(2,n),F(1,4*alpha),F(1,2)]]
        rounds.append({"n":n,"alpha":alpha,"checker_epsilon":small(EPS),"models":models,
                       "endpoints":endpoints(n,alpha),"source_entropy":source_entropy(n),
                       "exact_checks":CHECKS-before,"status":"PASS_EXACT"})
    result={"status":"PASS_EXACT","seed":"none; integers and Fraction",
            "scope":reg["scope"],"theory_epsilon_unchanged":reg["theory_epsilon"],
            "checker_epsilon":small(EPS),"cutoffs_are_not_theory_log_epsilon_certificates":True,
            "registration_sha256":hashlib.sha256(registration.read_bytes()).hexdigest(),
            "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "rounds":rounds,"total_exact_checks":CHECKS}
    output=here/"low_count_future_capacity_exact_guard_results_20261007.json"
    output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps({"status":result["status"],"rounds":3,"checks":CHECKS,"output":str(output),
                      "cutoff_table":[{"n":r["n"],"sigma":m["sigma"],**m["cutoff"]}
                                      for r in rounds for m in r["models"]]},ensure_ascii=False))


if __name__ == "__main__":
    main()
