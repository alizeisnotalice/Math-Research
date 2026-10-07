#!/usr/bin/env python3
"""Exact finite coefficient guards; no original FIRST/geom samples."""
from fractions import Fraction as F
from hashlib import sha256
from math import comb, isqrt, log10
from pathlib import Path
import json


def receipt(q):
    def digest(k):
        return sha256(k.to_bytes(max(1, (k.bit_length()+7)//8), "big")).hexdigest()
    return {"numerator_bits": q.numerator.bit_length(),
            "denominator_bits": q.denominator.bit_length(),
            "numerator_sha256_binary": digest(q.numerator),
            "denominator_sha256_binary": digest(q.denominator),
            "log10_diagnostic": (log10(q.numerator)-log10(q.denominator)
                                   if q else None)}


def log_upper(q, terms=24):
    """Rational upper certificate using log series after dyadic scaling."""
    assert q >= 1
    k = 0
    while q >= 2:
        q /= 2
        k += 1
    def series_upper(x):
        z = (x-1)/(x+1)
        partial = sum((2*z**(2*j+1)/F(2*j+1) for j in range(terms)), F(0))
        tail = 2*z**(2*terms+1)/(F(2*terms+1)*(1-z*z))
        return partial+tail
    return k*series_upper(F(2))+series_upper(q)


def thinned_beta_coefficient(n, alpha, m, r, counts, denominator):
    # W~Beta(1,alpha), J~Bin(n,W), K|J~Bin(J,r).
    # This positive thinning computation is independent of the density peak.
    if r == 0:
        return F(int(m == 0))
    if r == 1:
        return F(counts[m], denominator)
    a, b = r.numerator, r.denominator
    total = sum(counts[j]*comb(j,m)*a**m*(b-a)**(j-m)*b**(n-j)
                for j in range(m,n+1))
    return F(total, denominator*b**n)


def exact_round(n):
    root = isqrt(n)
    alpha = root if root*root == n else root+1
    assert 2 <= alpha <= n and alpha*alpha >= n
    denominator = comb(n+alpha, alpha)
    counts = [comb(n-m+alpha-1,alpha-1) for m in range(n+1)]
    assert sum(counts) == denominator
    assert F(sum(m*counts[m] for m in range(n+1)), denominator) == F(n,alpha+1)

    # Bin(n,1/alpha) common-denominator integer weights.
    bin_den = alpha**n
    bin_counts = [comb(n,m)*(alpha-1)**(n-m) for m in range(n+1)]
    assert sum(bin_counts) == bin_den
    tail = bin_den-bin_counts[0]
    inner_D, inner_M = F(0), F(0)
    majors = {0:F(1)}
    for m in range(1,n+1):
        inner = F(tail, m*bin_den)
        majors[m] = inner+F(counts[m],denominator)
        inner_D += inner
        inner_M += m*inner
        tail -= bin_counts[m]
    assert tail == 0 and inner_M == F(n,alpha)
    D_major = 1+inner_D+F(denominator-counts[0],denominator)
    M_major = inner_M+F(n,alpha+1)
    D_safe = 3+log_upper(F(n,alpha))
    assert D_major <= D_safe
    assert M_major <= F(2*n,alpha)
    assert M_major*M_major <= 4*n

    coefficient_checks = []
    for m in sorted({0,1,alpha,2*alpha,n}):
        for r in [F(0),F(1,alpha),F(1,2),F(1)]:
            exact = thinned_beta_coefficient(n,alpha,m,r,counts,denominator)
            assert exact <= majors[m]
            coefficient_checks.append({"m":m,"r":str(r),
                                       "coefficient":receipt(exact),
                                       "majorant_status":"PASS_EXACT"})

    # Uniform-in-t hard-persistent tensor ratio <=(1-v/9)^n.
    # Its Beta(alpha,1) average is independently computed from the
    # incomplete-beta/binomial-tail identity, not numerical quadrature.
    N = n+alpha
    tail9 = 9**N-sum(comb(N,j)*8**(N-j) for j in range(alpha))
    tensor_average = F(tail9, comb(N,alpha)*9**n)
    split_upper = F(1,alpha)**alpha+F(9*alpha-1,9*alpha)**n
    assert 0 < tensor_average <= split_upper < F(1,8)
    tiny_full_row_lower = F(1,comb(n+alpha,alpha))
    assert tiny_full_row_lower < tensor_average
    return {"n":n,"alpha":alpha,"normalization_status":"PASS_EXACT",
            "D_major":receipt(D_major),"D_safe":receipt(D_safe),
            "M_major":str(M_major),"M_safe":str(F(2*n,alpha)),
            "D_M_status":"PASS_EXACT",
            "coefficient_checks":coefficient_checks,
            "hardpersistent_tensor_average_upper":receipt(tensor_average),
            "hardpersistent_split_upper":receipt(split_upper),
            "hardpersistent_status":"PASS_EXACT_LESS_THAN_ONE_EIGHTH",
            "universal_tiny_row_lower":receipt(tiny_full_row_lower)}


def main():
    result = {"name":"future_softness_moment_exact_guard_20261007",
              "seed":None,"randomness":"none",
              "scope":"Exact finite Beta/Binomial coefficient, normalization, peak-majorant and tensor-kernel guards. No mu, hardband, FIRST, future-cap or geom input realization is certified. No global coverage inference.",
              "rounds":[exact_round(n) for n in [512,1024,4096]],
              "status":"PASS_EXACT"}
    path = Path(__file__).with_name("future_softness_moment_exact_guard_20261007_results.json")
    path.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print("PASS_EXACT: 3 rounds; 60 coefficient guards; D/M identities; 3 hardpersistent tensors.")
    print(path)


if __name__ == "__main__":
    main()
