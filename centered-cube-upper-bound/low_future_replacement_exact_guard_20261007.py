#!/usr/bin/env python3
"""Three exact finite replacement-spectrum guards, never a FIRST sample."""
from fractions import Fraction as F
from hashlib import sha256
from math import comb, isqrt, log10
from pathlib import Path
import json


def receipt(q):
    def digest(k):
        return sha256(k.to_bytes(max(1,(k.bit_length()+7)//8),"big")).hexdigest()
    return {"numerator_bits":q.numerator.bit_length(),
            "denominator_bits":q.denominator.bit_length(),
            "numerator_sha256_binary":digest(q.numerator),
            "denominator_sha256_binary":digest(q.denominator),
            "log10_diagnostic":log10(q.numerator)-log10(q.denominator) if q else None}


def spectrum(n,alpha):
    denominator=comb(n+alpha,alpha)
    weights=[comb(k+alpha-1,alpha-1) for k in range(n+1)]
    assert sum(weights)==denominator
    assert F(sum((n-k)*weights[k] for k in range(n+1)),denominator)==F(n,alpha+1)
    return weights,denominator


def uniform_beta_average(n,alpha,a,weights,denominator):
    p,d=a.numerator,a.denominator
    value=F(sum(weights[k]*p**k*d**(n-k) for k in range(n+1)),
            denominator*d**n)
    if a==F(3,4):
        # Independent incomplete-beta tail representation for (1-v/4)^n.
        N=n+alpha
        tail=4**N-sum(comb(N,j)*3**(N-j) for j in range(alpha))
        independent=F(tail,comb(N,alpha)*4**n)
        assert value==independent
    # Jensen in a radical-free exact form.
    terminal=a**n
    assert value**(alpha+1)>=terminal**alpha
    for v in [F(0),F(1,4),F(1,2),F(3,4),F(1)]:
        likelihood=(1-v+v*a)**n
        assert likelihood**v.denominator>=terminal**v.numerator
    # Original mask count, conditional on y, has p_i=sigma*a_i.
    sigma=F(1,2) if a<=2 else F(1,4)
    p_soft=sigma*a
    assert 0<=p_soft<=1
    for v in [F(1,4),F(1,2),F(3,4)]:
        z=1+v/(sigma*(1-v))
        original=(1-v)**n*(1-p_soft+p_soft*z)**n
        assert original==(1-v+v*a)**n
    return {"a":str(a),"Z_exact":receipt(value),
            "original_K_mean":str(n*p_soft),
            "replacement_K_mean":str(F(n*alpha,alpha+1)),
            "AMGM_Jensen_Beta_status":"PASS_EXACT"}


def mixed_coefficient(n,k,inside,outside):
    # a_i=3/4 for inside, a_i=2 for outside; positive elementary-symmetric
    # numerator with a common 4^k denominator.
    total=sum(comb(outside,j)*comb(inside,k-j)*8**j*3**(k-j)
              for j in range(max(0,k-inside),min(k,outside)+1))
    return F(total,4**k*comb(n,k))


def exact_round(n):
    root=isqrt(n)
    alpha=root if root*root==n else root+1
    weights,denominator=spectrum(n,alpha)
    tail_checks=[]
    for M in sorted({alpha,2*alpha,min(n-1,alpha*n.bit_length())}):
        exact=F(sum(weights[k] for k in range(n-M)),denominator)
        closed=F(comb(n-M+alpha-1,alpha),denominator)
        upper=F(n+alpha-M-1,n+alpha)**alpha
        assert exact==closed and exact<=upper
        tail_checks.append({"M":M,"tail":receipt(exact),
                            "hockey_stick_and_product_status":"PASS_EXACT"})
    uniform=[uniform_beta_average(n,alpha,a,weights,denominator)
             for a in [F(3,4),F(1),F(2)]]
    outside=n//8
    inside=n-outside
    terminal=F(3,4)**inside*2**outside
    # Scalar likelihood row is log-concave and initially decreasing. This
    # verifies only that scalar condition, never the full input geometry.
    initial_slope=F(3,4)*inside+2*outside-n
    assert initial_slope<0 and terminal<1
    mixed=[]
    for k in sorted({0,1,n-alpha,n-1,n}):
        B=mixed_coefficient(n,k,inside,outside)
        assert B**n>=terminal**k
        omission=n-k
        inverse_total=sum(comb(outside,j)*comb(inside,omission-j)
                          *3**j*8**(omission-j)
                          for j in range(max(0,omission-inside),min(omission,outside)+1))
        # b=1/a: inside 4/3, outside 1/2; common 6^m numerator.
        C=F(inverse_total,6**omission*comb(n,omission))
        assert B==terminal*C
        mixed.append({"k":k,"m":omission,"B_k":receipt(B),
                      "AMGM_and_omission_identity_status":"PASS_EXACT"})
    for v in [F(1,4),F(1,2),F(3,4)]:
        likelihood=(1-v+v*F(3,4))**inside*(1-v+2*v)**outside
        assert likelihood**v.denominator>=terminal**v.numerator
        p_inside=F(3,8)
        z=1+v/(F(1,2)*(1-v))
        original=(1-v)**n*(1-p_inside+p_inside*z)**inside*z**outside
        assert likelihood==original
    return {"n":n,"alpha":alpha,"prior_mean_omitted":str(F(n,alpha+1)),
            "normalization_status":"PASS_EXACT","tail_checks":tail_checks,
            "uniform_profiles":uniform,"mixed_inside":inside,"mixed_outside":outside,
            "mixed_initial_slope":str(initial_slope),"mixed_coefficients":mixed,
            "scope":"Positive ratio-coordinate finite models; mixed scalar monotonicity only. No common mu/hardband/Lwinner/FIRST realization or low-residual geometry certified."}


def main():
    result={"name":"low_future_replacement_exact_guard_20261007",
            "seed":None,"randomness":"none",
            "rounds":[exact_round(n) for n in [512,1024,4096]],
            "status":"PASS_EXACT",
            "scope":"Exact AMGM/Jensen, positive Beta identity, replacement tail, omission identity, and original-K distinction. No original kernel quadrature, no actual FIRST sample, no coverage or spatial order inference."}
    path=Path(__file__).with_name("low_future_replacement_exact_guard_20261007_results.json")
    path.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print("PASS_EXACT: 3 rounds; 9 uniform profiles; 15 mixed spectrum identities; 9 tails.")
    print(path)


if __name__=="__main__":
    main()
