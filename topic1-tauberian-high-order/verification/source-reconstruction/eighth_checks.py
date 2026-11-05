"""Independent exact finite-series arithmetic for the current top-level proof.

No geometric hypothesis, full-input inequality, or infinite-dimensional
optimization is certified by this calculation.
"""
from fractions import Fraction as F
from math import comb
import json

def zero(N): return [F(0)] * (N + 1)
def constant(c, N):
    out = zero(N); out[0] = F(c); return out
def monomial(k, c, N):
    out = zero(N)
    if k <= N: out[k] = F(c)
    return out
def add(*args):
    return [sum(values, F(0)) for values in zip(*args)]
def scale(a, c): return [x * F(c) for x in a]
def mul(a, b):
    N = len(a) - 1
    out = zero(N)
    for i, x in enumerate(a):
        if not x: continue
        for j, y in enumerate(b[:N-i+1]):
            out[i+j] += x*y
    return out
def power(a, n):
    if n < 0: return power(inverse(a), -n)
    out = constant(1, len(a)-1)
    for _ in range(n): out = mul(out, a)
    return out
def inverse(a):
    N = len(a)-1
    assert a[0]
    out = zero(N); out[0] = 1/a[0]
    for k in range(1, N+1):
        out[k] = -sum((a[j]*out[k-j] for j in range(1,k+1)), F(0))/a[0]
    return out
def divide(a, b): return mul(a, inverse(b))
def shift(a, k): return ([F(0)]*k+a)[:len(a)]

def constants(n):
    K = F(comb(n,2)*n*n,27)
    c8 = -F(13,18)-F(n*(n-2),3)+F(400,243)
    polynomial = F((n-1)*n**3*(2609+1971*n+81*n*n),26244)
    J = polynomial+(F(8,27) if n == 2 else 0)
    Jhat = J+(F(8,27) if n == 2 else 0)
    return K,c8,J,Jhat

def upper_series(n):
    N=n+9
    one=constant(1,N)
    pv=add(one,monomial(1,1,N))
    A=power(pv,n)
    a=monomial(n,1,N)
    K,c8,J,Jhat=constants(n)
    Q=scale(add(shift(power(pv,-4),6),
                scale(shift(power(pv,-5),7),F(10,9)),
                scale(shift(power(pv,-2),8),c8)),K)
    if n == 2: Q=add(Q,monomial(8,F(8,27),N))
    L=divide(A,add(A,scale(a,-1)))
    U=divide(A,add(A,scale(mul(a,add(one,Q)),-1)))
    diff=add(U,scale(L,-1))
    got=[diff[n+j] for j in (6,7,8)]
    want=[K,-(n+F(26,9))*K,Jhat]
    assert got == want
    baseline=[]
    for j in range(9):
        closed=sum(((-1)**(j-(q-1)*n)*comb(n+j-1,j-(q-1)*n)
                    for q in range(1,2+j//n)),0)
        assert L[n+j] == closed
        baseline.append(closed)
    true_v=[F(x) for x in baseline]
    for j, corr in zip((6,7,8),want): true_v[j]+=corr
    return {
        "n":n,"upper_normalized_d6_d7_d8":[str(x) for x in got],
        "baseline_ell_j0_to_j8":baseline,
        "ordinary_t_coefficients_j0_to_j8":[str(x*2**(n+j)) for j,x in enumerate(true_v)]
    }

def macro_series(n, safety, micro=False):
    N=n+9
    one=constant(1,N)
    g=inverse(add(one,monomial(1,2,N)))
    a=monomial(n,2**n,N)
    K,c8,J,Jhat=constants(n)
    B=comb(n,2)
    kappa=add(scale(power(g,2),F(8*n,3)),
              scale(shift(power(g,3),1),F(320*n,81)))
    Lambda0=scale(mul(power(kappa,3),power(g,-1)),F(B,8*n*n))
    Lambda1=scale(mul(power(kappa,4),power(g,-2)),-F(5*B,96*n**3))
    Lambda2=add(scale(mul(power(kappa,5),power(g,-3)),F(13*B,512*n**4)),
                scale(mul(power(kappa,4),power(g,-2)),F(B*(n-2),32*n*n)))
    Lambda=add(Lambda0,shift(Lambda1,1),shift(Lambda2,2),
               monomial(3,safety,N))
    if micro:
        assert n == 2
        Lambda=add(Lambda,scale(mul(a,Lambda0),-1))
    e=shift(kappa,3)
    w=add(shift(g,1),shift(Lambda,6))
    rho=add(one,scale(w,-2))
    M=mul(a,power(rho,n))
    det=mul(power(add(one,scale(e,-1)),n-1),add(one,scale(e,n-1)))
    ratio=divide(det,add(det,scale(M,-1)))
    baseline=inverse(add(one,scale(mul(a,power(g,n)),-1)))
    diff=add(ratio,scale(baseline,-1))
    got=[diff[n+j]/2**(n+j) for j in (6,7,8)]
    want=[K,-(n+F(26,9))*K,Jhat if micro else J]
    assert got == want
    return {"n":n,"safety_constant":safety,"micro_budget":micro,
            "normalized_d6_d7_d8":[str(x) for x in got]}

upper=[upper_series(n) for n in range(2,13)]
macro=[macro_series(n,7) for n in range(2,13)]
micro=[macro_series(2,c,True) for c in (1,100)]
assert upper[1]["baseline_ell_j0_to_j8"][5] == 0
assert upper[1]["ordinary_t_coefficients_j0_to_j8"][5] == "0"

# Independent polynomial-in-n check after dividing out K_n.
# Coefficients are in ascending powers of n; dimension-specific terms are
# retained separately, since the indicator n=2 is not a polynomial in n.
from_upper=[F(40,9)-F(13,18)+F(400,243),F(2,3)+F(26,9)+F(1,2),-F(1,3)+F(1,2)]
from_closed=[F(2609,486),F(1971,486),F(81,486)]
assert from_upper == from_closed
result={
    "arithmetic":"Python fractions.Fraction; rational-series inversion and exact determinant",
    "variable_conventions":"Upper series variable v; lower series variable u, v=2u, a=(2u)^n.",
    "upper_rows":upper,
    "macro_rows":macro,
    "micro_rows":micro,
    "general_n_polynomial_coefficients_after_dividing_by_K":[str(x) for x in from_upper],
    "n3_j5_zero":True,
    "scope":"Finite exact algebra checks plus a general polynomial coefficient identity; no geometric proof is supplied by these computations."
}
print(json.dumps(result,ensure_ascii=False,indent=2))
