# R042 complete-cost comparison: independent audit

2026-10-28. **ACCEPT with the source's explicit validity-domain restrictions.**
For L=log n>=2^240, the new displayed bound U42 is strictly below every
listed *valid* old complete cost at t_old=n/(L-3log L). The comparison
against the all-t unit-core infimum is separately eventual and has no
certified explicit starting point in this appendix.

## 1. Source identity and scope

Source read in full:
`/Users/zhengzhihao/Desktop/T/cube-r030-parameters/cube-dim-7-bounds/upper/round042/comparison/PROOF.md`,
SHA256 `e9dc339e64164f4afb51e20ee81dd03d38f49132353469849de906ec84844c6a`.
The actual R040 parameter definitions and the archived old sources were
also checked; their identities, the exact R034 matrix identity and both
new frozen proofs are in `exact_receipt.json`.

This review does not redo the old operator proofs, parameter certificates,
LP solves, original numerical runs or production computations. It checks
the new comparison algebra, real parameters, uniform domains and infimum
logic. Only `audit/round042/comparison/` is written; there is no Git,
n=7, MATLAB, parameter search or decimal reconstruction.

## 2. The old lower bound is for complete cost expressions

For t_old=n/x_old, x_old=L-3log L in [L/2,L], the accepted integral
upper bound and nonnegative core give

```
B_old>=n log t_old-n^2/[(n-1)x_old]
     >=n(L-log x_old)-2n/x_old>=T-4n/L.
```

The subtracted integral has the correct sign here. For every old mechanism
in the table, on that mechanism's valid domain and with B_old-D>0,

```
F_old(D)=A_old+B_old-D+2sqrt(A_old(B_old-D))
        >=B_old-D>=T-4n/L-2n sqrt(l/L),           (2.1)
l=log L.
```

Thus the estimate is a lower bound for a proved *upper-bound formula*,
not for the actual C_n. It keeps enough of the complete expression to
compare different x scales safely: A and the root cross term are
nonnegative, so discarding them is a lower-bound step.

## 3. Each old deduction is within the common envelope only when valid

The common deduction bound 0<=D<=2n sqrt(l/L) is correct for l>=1:

- The cap formula has D=0, and R033 has coefficient 1/128 on n/sqrt L.
- Every legal R034 q satisfies the proved family bound
  q<q_inf<39/320<1/8. The radical comparison is exact:
  sqrt41>6 and 0<11-sqrt41<5<9. Thus this covers the whole theoretical
  family, not just the twelve reported candidates.
- R035/R036 have coefficient 2^(-16) on n sqrt(l/L).
- Reading the actual R040 source confirms x=L-3l, with a full-cost
  mechanism coefficient gamma<c(1-2epsilon)<sqrt2<2. Its simpler
  T-gamma*n sqrt(l/L) corollary is obtained using a larger intermediate
  mechanism coefficient and absorbing a remainder. The appendix does
  not confuse that display coefficient with the intermediate complete
  cost. Both are strictly below sqrt2, so the uniform deduction bound
  remains valid whichever legitimate mechanism instance is compared.

R034 still requires L>=max(64,Y^2) and all its parameter/capacity tests.
R040 still requires its own parameter-dependent sufficiently-large-n
threshold for fixed c,epsilon,gamma. The source correctly intersects these
domains with L>=2^240; that number does not make all R040 parameter choices
simultaneously valid. It also does not automatically satisfy Y^2 for
arbitrarily large Y. R033 and the explicit R036 formula, as well as the
twelve fixed R034 rows, do have their stated explicit gates covered.

## 4. The explicit infinite-tail gap

U42=T-2^(-32)nL^(-1/3) is below (2.1) if

```
G(L)=4L^(-2/3)+2sqrt(log L)L^(-1/6)<2^(-32).
```

The first term decreases. The logarithmic derivative of the second
is [1/(2log L)-1/6]/L, negative whenever log L>3. At L0=2^240,
the elementary inequalities 2/3<log2<1 give
160<log L0<240<16^2. The lower log bound follows by differentiating
log x-2(x-1)/(x+1), whose derivative is positive for x>1; the upper
bound follows from e>2.

Consequently, at that exact base,

```
G(L0)<2^(-158)+2^(-35)<2^(-32).
```

Monotonicity covers every larger L, so

```
C_n<=F42<U42<T-4n/L-2n sqrt(l/L)<=F_old(D).       (4.1)
```

The new theorem's own explicit L>=2^108 requirement is included in this
domain. This is not a deduction-only comparison of changed x parameters:
the new *full cost* is upper-bounded by U42, and all listed old full costs
are lower-bounded by the common expression in (4.1).

The explicit gap is independent of the chosen valid old deduction. Hence
taking an infimum over these valid fixed-t_old instances preserves an
infimum at least the common lower bound, still strictly above U42. No
pointwise-strict-to-infimum error or optimum-attainment assumption occurs.

## 5. All-t unit-core family: eventual only

For every strictly feasible unit-core pair, the existing payment condition
implies E=B-D>A and 0<D<1. Therefore

```
F(A,B)-F(A,B-D)
 =D[1+2sqrt A/(sqrt B+sqrt E)]<2D<2.
```

Taking the infimum over every legal t and sigma gives only the non-strict
bound V_unit>=U_old-2. An empty feasible set has infimum +infinity.
Cost limits of legal feasible parameters preserve that non-strict bound;
this is not permission to run the finite algorithm at an uncertified
guard boundary or to assert validity for every formal weak-boundary point.

The accepted fine asymptotic gives

```
U_old-T=(3l+1)n/L+O(nl^2/L^2).
```

Its error-to-positive-main ratio tends to zero, and that main term tends
to infinity. Thus U_old-2>T eventually, proving U42<T<V_unit once both
the new explicit theorem and this old eventual statement are valid.
No explicit constant or starting point for the old O term was provided;
the appendix correctly does not assign L>=2^240 to this separate result.
These remain formula-family comparisons, not C_n lower bounds.

## 6. Actual minimal exact checks

`python3 check_exact.py` actually returned PASS, PID 77641. It read the
twelve exact `certificate.q` strings from the frozen original JSON
(SHA256 `87ac7a68bf3979f6b654f7954044feac4d0d26031e14dc38152c8e48684d4a3d`).
It confirmed every q is below 1/8 and the maximum exact fraction is
74275164113/10^12 at Y=256, with required L>=65536. No optimizer value,
rounded q_approx or raw decimal was used, and no earlier certificate
verification or LP was repeated.

The same small Fraction check verified the exact positive gap
2^(-32)-2^(-35)-2^(-158), the radical family's rational upper comparison,
the R036 guard 240*(2/3)>128 and the new/fixed-row explicit domains.
The actual infinite-tail conclusion comes from Section 4, not these
finite arithmetic checks. The frozen R042 full and threshold source
identities remain unchanged.

## 7. Final verdict

The complete-cost direction, actual parameter choices, exact old constants,
validity domains and asymptotic-versus-explicit distinction all pass.
No source correction is required. The appendix does not assert an
exhaustive public comparison, publication status, global parameter
optimality or a world record. Reported C_n bounds should still be the
minimum of every bound valid at the dimension in question.
