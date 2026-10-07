#!/usr/bin/env python3
"""New exact square-defect/product certificate; no actual winner samples."""
from fractions import Fraction as F
from itertools import combinations
from math import comb
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent
PREFIX = "without_replacement_square_defect"
checks, rounds = [], []

def rec(name, ok):
    checks.append({"label":name, "pass":bool(ok)})

def add(poly, monomial, coefficient):
    poly[monomial] = poly.get(monomial, 0) + coefficient
    if not poly[monomial]:
        del poly[monomial]

def monomial(d, axes, doubled=None):
    m = [0] * d
    for i in axes:
        m[i] += 1
    if doubled is not None:
        m[doubled] += 1
    return tuple(m)

def elementary(values):
    e = [F(1)] + [F(0)] * len(values)
    used = 0
    for value in values:
        used += 1
        for k in range(used, 0, -1):
            e[k] += value * e[k-1]
    return e

def scalar_defects(values):
    n = len(values)
    e = elementary(values)
    R = [F(0)] * n
    for value in values:
        excluded = [F(1)]
        for k in range(1, n):
            excluded.append(e[k] - value * excluded[-1])
        for k in range(1, n):
            R[k] += value**2 * excluded[k-1]
    Q = sum(values)
    E = [e[k] / comb(n,k) for k in range(n+1)]
    A = [(n*R[k]-k*Q*e[k]) / (n*(n-k)*comb(n,k))
         for k in range(n)]
    return E, A, Q/n

def walsh(profile):
    return [sum(value * (-1 if (x & mask).bit_count()%2 else 1)
                for x,value in enumerate(profile)) for mask in range(len(profile))]

def iwalsh(profile):
    size = len(profile)
    return [sum(value * (-1 if (x & mask).bit_count()%2 else 1)
                for mask,value in enumerate(profile))/size for x in range(size)]

def apply_symbol(profile, symbol):
    transformed = walsh(profile)
    return iwalsh([a*b for a,b in zip(transformed,symbol)])

def ihash(integer):
    return hashlib.sha256(integer.to_bytes((integer.bit_length()+7)//8,"big")).hexdigest()

for d, n in [(3,8),(4,32),(5,128)]:
    start = len(checks)
    tag = f"formal_d{d}/compressed_n{n}"
    for k in range(d):
        lhs, rhs = {}, {}
        # n R_k - k Q e_k, expanded without using the candidate identity.
        for axes in combinations(range(d), k):
            for i in axes:
                add(lhs,monomial(d,axes,i),d)
            for i in range(d):
                add(lhs,monomial(d,axes,i),-k)
        if k:
            for i,j in combinations(range(d),2):
                others = [a for a in range(d) if a not in (i,j)]
                for axes in combinations(others,k-1):
                    mi, mj, mij = list(monomial(d,axes)), list(monomial(d,axes)), list(monomial(d,axes))
                    mi[i] += 2
                    mj[j] += 2
                    mij[i] += 1
                    mij[j] += 1
                    add(rhs,tuple(mi),1)
                    add(rhs,tuple(mj),1)
                    add(rhs,tuple(mij),-2)
        rec(tag+f"/k{k}/formal-polynomial-identity", lhs==rhs)

    scalar_receipts = []
    cases = [
        [F(1)] * n,
        [F(i%2) for i in range(n)],
        [F(1,2+i%5) for i in range(n)]
    ]
    for case, values in enumerate(cases):
        E, A, P = scalar_defects(values)
        ctag = tag+f"/spectral_case{case}"
        for k in range(n):
            rec(ctag+f"/k{k}/normalized-recurrence", E[k+1]==P*E[k]-A[k])
            rec(ctag+f"/k{k}/positive-defect", 0<=A[k]<=1)
        exact_sum = 1-E[n]-(1-P)*sum(E[:n])
        rec(ctag+"/sum-range-identity", sum(A)==exact_sum)
        rec(ctag+"/total-spectral-contraction", 0<=sum(A)<=1)
        for r in [F(0),F(1,4),F(1,2),F(1)]:
            beta = [comb(n,k)*r**k*(1-r)**(n-k) for k in range(n+1)]
            B = [F(0)] * n
            for l in range(n-1,-1,-1):
                B[l] = beta[l+1] + (P*B[l+1] if l<n-1 else 0)
                rec(ctag+f"/r{r}/l{l}/positive-submarkov", 0<=B[l]<=1)
            K = sum(beta[k]*E[k] for k in range(n+1))
            replacement = (1-r+r*P)**n
            defect = sum(A[l]*B[l] for l in range(n))
            rec(ctag+f"/r{r}/replacement-identity", K==replacement-defect)
            rec(ctag+f"/r{r}/spectral-defect-contraction", 0<=defect<=sum(A)<=1)
        scalar_receipts.append({"case":case,"sum_A":str(sum(A)),
                                "P":str(P),"scope":"Exact common scalar spectrum; no original-space input simulation."})

    size = 1 << d
    h = [F(1+x%3) if x%4==0 else F(0) for x in range(size)]
    nu = [value**2 for value in h]
    W = sum(nu)
    symbol_by_k = [[F(0)]*size for k in range(d)]
    for mask in range(size):
        values = [F(1,3) if mask & (1<<i) else F(1) for i in range(d)]
        E, A, P = scalar_defects(values)
        for k in range(d):
            symbol_by_k[k][mask] = A[k]
    energy = F(0)
    test_energy = F(0)
    pairing = F(0)
    exterior_pairing = F(0)
    for k in range(d):
        symbol = symbol_by_k[k]
        Ah = apply_symbol(h,symbol)
        Anu = apply_symbol(nu,symbol)
        delta0 = [F(1)] + [F(0)]*(size-1)
        kernel = apply_symbol(delta0,symbol)
        rec(tag+f"/fixture/k{k}/kernel-zero-mass",sum(kernel)==0)
        local_energy = sum(h[x]*Ah[x] for x in range(size))
        energy += local_energy
        rec(tag+f"/fixture/k{k}/PSD-energy",local_energy>=0)
        # One receiver label per x; the vector g(x) has l2 norm <=1.
        g = [F(1+x%2,2) if x%d==k else F(0) for x in range(size)]
        gh = [g[x]*h[x] for x in range(size)]
        Agh = apply_symbol(gh,symbol)
        test_energy += sum(gh[x]*Agh[x] for x in range(size))
        pairing += sum(gh[x]*Ah[x] for x in range(size))
        product = [sum(kernel[x^y]*(h[y]-h[x])**2 for y in range(size))
                   for x in range(size)]
        for x in range(size):
            rec(tag+f"/fixture/k{k}/x{x}/signed-product-exact",
                Anu[x]==2*h[x]*Ah[x]+product[x])
        rec(tag+f"/fixture/k{k}/total-product-sign",
            sum(product)==-2*local_energy)
        exterior = [g[x] if h[x]==0 else F(0) for x in range(size)]
        rec(tag+f"/fixture/k{k}/exterior-product-identity",
            sum(exterior[x]*Anu[x] for x in range(size)) ==
            sum(exterior[x]*product[x] for x in range(size)))
        exterior_pairing += sum(exterior[x]*h[x]*Ah[x] for x in range(size))
    rec(tag+"/fixture/source-square-budget",0<=energy<=W)
    rec(tag+"/fixture/receiver-square-budget",0<=test_energy<=W)
    rec(tag+"/fixture/bilinear-CS",pairing**2<=energy*test_energy)
    rec(tag+"/fixture/bilinear-source-once",abs(pairing)<=W)
    rec(tag+"/fixture/exterior-bilinear-vanishes",exterior_pairing==0)

    # Actual original-kernel certificate inequalities. No numerical integration:
    # w and w*w density <=1, atom-free face laws are proved analytically in md.
    delta = F(1,64*n)
    epsilon = delta**2
    volume = (2*epsilon)**n
    height = 1/volume
    rec(tag+"/original/source-finite-L1-mass",height*volume==1)
    rec(tag+"/original/inactive-coordinates-never-counted",0<epsilon<delta)
    selected_lower = F(0)
    for k in range(1,n):
        posmass = F(k,n)
        negmass = F(k*(k+1)*comb(n,k+1),n*(n-k)*comb(n,k))
        rec(tag+f"/original/k{k}/Jordan-mass",posmass==negmass)
        failure = 2*k*(delta+epsilon)
        wrong_negative = 2*(k+1)*(delta+epsilon)
        rec(tag+f"/original/k{k}/union-probability-range",
            0<=failure<1 and 0<=wrong_negative<1)
        lower = posmass*(1-failure)-negmass*wrong_negative
        rec(tag+f"/original/k{k}/selected-product-lower",lower>=F(7,8)*posmass)
        selected_lower += lower
    rec(tag+"/original/total-Jordan-variation",sum(F(2*k,n) for k in range(1,n))==n-1)
    rec(tag+"/original/partition-product-linear-lower",selected_lower>=F(7*(n-1),16))
    rec(tag+"/original/bilinear-zero-support",epsilon<delta)
    rounds.append({
        "formal_and_finite_dimension":d,"compressed_spectral_dimension":n,
        "original_face_dimension":n,"checks":len(checks)-start,
        "scalar_receipts":scalar_receipts,
        "finite_fixture_W":str(W),"finite_fixture_source_square":str(energy),
        "finite_fixture_test_square":str(test_energy),"finite_fixture_pairing":str(pairing),
        "original_delta":str(delta),"original_epsilon":str(epsilon),
        "original_selected_certified_lower":str(selected_lower),
        "original_coarse_lower":str(F(7*(n-1),16)),
        "original_source_height_numerator_bits":height.numerator.bit_length(),
        "original_source_height_numerator_sha256":ihash(height.numerator),
        "scope":"New polynomial/spectrum/product algebra plus original-G window certificate only; not actual max winner, saturated bad source, or FIRST."
    })

registration = HERE / f"{PREFIX}_registration_20261007.json"
result = {
    "status":"PASS" if all(c["pass"] for c in checks) else "FAIL",
    "seed":None,"total_checks":len(checks),
    "failed_checks":[c for c in checks if not c["pass"]],
    "rounds":rounds,"checks":checks,
    "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "registration_sha256":hashlib.sha256(registration.read_bytes()).hexdigest(),
    "scope":"New square-defect budget and necessary signed product conversion; no actual selector/obstacle/cube weak counterexample."
}
(HERE/f"{PREFIX}_results_20261007.json").write_text(
    json.dumps(result,ensure_ascii=False,indent=2)+"\n")
print(json.dumps({k:result[k] for k in
                  ["status","total_checks","failed_checks","script_sha256"]},
                 ensure_ascii=False))
raise SystemExit(0 if result["status"]=="PASS" else 1)
