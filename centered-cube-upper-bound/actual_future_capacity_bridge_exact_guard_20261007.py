#!/usr/bin/env python3
"""Finite exact guards for the new unnormalized future/capacity tail bridge.

Not a cube input, actual FIRST sample, path mask, or numerical space experiment.
All continuous-future coefficient fixtures have C_k=gamma**k, 0<gamma<=1;
their full response decreases with softness, so their row ceiling is analytic.
The physical hard column and original source tree are upstream theorem inputs.
"""
from fractions import Fraction as F
from math import comb, isqrt
from pathlib import Path
import hashlib
import json

CHECKS = 0


def check(truth, label):
    global CHECKS
    CHECKS += 1
    if not truth:
        raise AssertionError(label)


def ceil_f(x):
    return -((-x.numerator)//x.denominator)


def ceil_sqrt(x):
    k = isqrt(x.numerator//x.denominator)
    return k + int(k*k*x.denominator < x.numerator)


def small(x):
    x = F(x)
    return f"{x.numerator}/{x.denominator}"


def big_pair(numerator, denominator):
    # The integers need not be reduced. Each reported pair is a genuine exact
    # numerator/denominator certificate; bit framing distinguishes all pairs.
    chunks = []
    for v in (numerator, denominator):
        sign = b"-" if v < 0 else b"+"
        binary = abs(v).to_bytes(max(1,(abs(v).bit_length()+7)//8),"big")
        chunks.append(sign+len(binary).to_bytes(8,"big")+binary)
    out = {"numerator_bits": abs(numerator).bit_length(),
           "denominator_bits": denominator.bit_length(),
           "framed_integer_pair_sha256": hashlib.sha256(b"".join(chunks)).hexdigest()}
    if max(abs(numerator).bit_length(),denominator.bit_length()) <= 256:
        out["unreduced_exact"] = f"{numerator}/{denominator}"
    return out


def certified_cutoff(n, sigma, p):
    if p == 0:
        return n+1, None, None
    uint = 0
    while (2**uint)*p.numerator < p.denominator:
        uint += 1
    # e>=2 and 2**uint>=1/p give uint>=ln(1/p), with no float log.
    v = n*sigma
    root_upper = ceil_sqrt(2*v*uint)
    k = ceil_f(v+root_upper+F(2*uint,3))
    check(F(root_upper**2) >= 2*v*uint, "integer sqrt upper certificate")
    check(k >= v+root_upper+F(2*uint,3), "integer threshold upper certificate")
    return k, uint, root_upper


def exact_tail(n, probability, k):
    """Integer common-denominator CDF; correlated parity acceptance <=1.

    Acceptance is 1/3 on even K and 1/2 on odd K. This only tests positive
    gate domination, not an actual history distribution.
    """
    if k > n:
        return 0, 1, 0, 12
    if probability == 0:
        return (1,1,4,12) if k <= 0 else (0,1,0,12)
    if probability == 1:
        gate_num = 6 if n % 2 else 4
        return (1,1,gate_num,12) if k <= n else (0,1,0,12)
    a, d = probability.numerator, probability.denominator
    b = d-a
    denominator = d**n
    term = b**n
    cdf = weighted_cdf = 0
    for j in range(max(0,k)):
        if j in {0,1,k-1}:
            check(term == comb(n,j)*a**j*b**(n-j), "independent binomial coefficient endpoint")
        cdf += term
        weighted_cdf += (3 if j%2 else 2)*term  # denominator 6D
        if j+1 < k:
            numerator = term*(n-j)*a
            divisor = (j+1)*b
            check(numerator % divisor == 0, "exact binomial recurrence")
            term = numerator//divisor
    tail = denominator-cdf
    # Full parity-weighted total is [5D-(d-2a)^n]/(12D).
    gated_tail = 5*denominator-(d-2*a)**n-2*weighted_cdf
    check(0 <= tail <= denominator, "unnormalized plain tail positivity")
    check(0 <= gated_tail <= 12*tail, "positive correlated gate tail domination")
    return tail, denominator, gated_tail, 12*denominator


def coefficient_fixture(n, sigma, gamma, source_p, arity):
    k, uint, root_upper = certified_cutoff(n,sigma,source_p)
    v = n*sigma
    initial_base = 1-sigma+sigma*gamma
    check(0 < initial_base <= 1, "positive original contact response")
    # A common full-source normalized Bernstein polynomial is
    # [(1-s+s gamma)/(1-sigma+sigma gamma)]**n. All future s>=sigma
    # obey the same ceiling q=1, because gamma<=1.
    check(0 < gamma <= 1, "analytic continuous future monotonicity certificate")
    if k > n:
        z = None
        pref_num, pref_den = 0, 1
    elif k == v or sigma == 0:
        z = F(1)
        pref_num, pref_den = 1, 1
    else:
        check(0 < sigma < 1 and v < k < n, "interior rational optimal tilt")
        z = F(k)*(1-sigma)/(sigma*(n-k))
        check(z >= 1, "future tilt direction")
        B = 1-sigma+sigma*z
        future_s = sigma*z/B
        check(sigma <= future_s <= 1, "same original physical L and future softness")
        future_base = (1-future_s+future_s*gamma)/initial_base
        original_count_p = sigma*gamma/initial_base
        check(1-original_count_p+original_count_p*z == B*future_base,
              "exact Bernstein tilt with original q normalization")
        check(future_base <= 1, "whole-source continuous future ceiling")
        pref_num = B.numerator**n*z.denominator**k
        pref_den = B.denominator**n*z.numerator**k
        check(pref_num*source_p.denominator <= pref_den*source_p.numerator,
              "certified rational tail prefactor <= original source capacity")
    posterior_count_p = sigma*gamma/initial_base
    tail_num, tail_den, gated_num, gated_den = exact_tail(n,posterior_count_p,k)
    check(tail_num*source_p.denominator <= tail_den*source_p.numerator,
          "exact full-source coefficient tail <= p q")
    check(gated_num*source_p.denominator <= gated_den*source_p.numerator,
          "gate-weighted tail <= p q")
    # arity equal full-born source children. Each group's response is 1/arity
    # of the whole coefficient polynomial; source masses are also 1/arity.
    rootn = isqrt(n)
    h_opp = F(arity-1,arity)
    check(source_p*h_opp == F(1,rootn*arity), "source-fixed forward capacity exact equality")
    check(gated_num*source_p.denominator <= arity*gated_den*source_p.numerator,
          "individual unnormalized group tail <= p q, not posterior-conditioned")
    # Divide by q=1, sum all children, and use their normalized hard column
    # Hopp. Multiplication by the real theorem's N_h is left symbolic.
    check(rootn*(arity-1)*gated_num <= arity*gated_den,
          "combined source-once normalized hard-column charge <= 1/sqrt n")
    check(arity*source_p*h_opp == F(1,rootn), "all child source masses paid exactly once")
    return {"sigma": small(sigma), "gamma": small(gamma), "p_C": small(source_p),
            "u_integer_upper": uint, "sqrt_integer_upper": root_upper,
            "certified_integer_cutoff": k, "tail_empty": k>n,
            "rational_tilt_z": None if z is None else small(z),
            "prefactor_certificate": big_pair(pref_num,pref_den),
            "exact_plain_tail": big_pair(tail_num,tail_den),
            "exact_gate_tail": big_pair(gated_num,gated_den),
            "normalized_combined_charge_bound": small(F(1,rootn)),
            "status": "PASS_EXACT"}


def endpoint_checks(n,p):
    result = []
    for sigma,cap in [(F(0),p),(F(1),p),(F(0),F(1)),(F(1),F(1)),(F(1,2),F(0))]:
        k, uint, root_upper = certified_cutoff(n,sigma,cap)
        if sigma in (0,1):
            K = 0 if sigma == 0 else n
            tail = int(K >= k)
        else:
            check(cap == 0 and k == n+1, "zero-capacity mass-null tail convention")
            tail = 0
        check(F(tail) <= cap, "sigma/capacity endpoint tail response")
        equality = sigma in (0,1) and K == k
        if equality:
            check(tail == 1 and not (K < k), "computed threshold equality is paid, not strict low")
        result.append({"sigma":small(sigma), "p":small(cap), "k":k,
                       "tail_mass":tail, "threshold_equality":equality,
                       "status":"PASS_EXACT"})
    return result


def main():
    rounds = []
    for n in [1024,4096,16384]:
        before = CHECKS
        rootn, arity = isqrt(n), 512
        check(rootn*rootn == n, "integer square dimension")
        p = F(1,rootn*(arity-1))
        cases = []
        for sigma in [F(2,n),F(1,4*rootn),F(1,2)]:
            for gamma in [F(1),F(1,2)]:
                cases.append(coefficient_fixture(n,sigma,gamma,p,arity))
        endpoints = endpoint_checks(n,p)
        rounds.append({"n":n, "sqrt_n":rootn, "source_arity":arity,
                       "source_mass": "1/1", "models":cases, "endpoints":endpoints,
                       "exact_checks":CHECKS-before, "status":"PASS_EXACT"})
    result = {"scope":"finite positive Bernstein/tilt/tail and source-capacity coefficient certificates only",
              "not_actual_FIRST_or_geometric_source_fixture":True,
              "not_actual_early_continuation_mask":True,
              "real_hard_column_Nh_is_upstream_analytic_input":True,
              "canonical_threshold":"n sigma + sqrt(2 n sigma ln(1/p)) + (2/3)ln(1/p)",
              "guard_threshold":"certified integer upper cutoff; integer u>=ln(1/p), integer sqrt upper, inclusive tail",
              "seed":"none; exact integers and Fraction",
              "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "rounds":rounds, "total_exact_checks":CHECKS, "status":"PASS_EXACT"}
    output = Path(__file__).with_name(Path(__file__).stem.replace("_exact_guard","_exact_guard_results")+".json")
    output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps({"status":result["status"],"rounds":3,"models":18,"endpoint_models":15,
                      "exact_checks":CHECKS,"output":str(output)},ensure_ascii=False))


if __name__ == "__main__":
    main()
