# R030 independent uniform-upper-bound audit

2026-10-27. **ACCEPT.** The proposed bound and its strict comparison with the
previous uniform upper expression are correct. No mathematical repair of the
proposal was required. This report supplies explicit proofs of every threshold.

## 1. Scope, sources, and notation

This audit independently reconstructs the scalar improvement from the already
accepted operator estimate. It does not re-audit the full covering theorem.
The audited claim is frozen in `proposal_snapshot.md`; its SHA256 and those of
all read sources are in `check_result.json`. The independent clone has HEAD
`f65a5bf92911ebed4cf803071bd7693978e5e582`.

Sources read: repository `AGENTS.md`; `revision_v3/cover/audit.md`, specifically
the definition and positivity of the limiting scalar budget;
`revision_v3/cover/asymptotic_cross_review.md`, specifically the previous uniform
bound and the definition of U_n; and `revision_v3/asymptotic/proof.md`.
All paths in this paragraph are relative to `cube-dim-7-bounds/` except AGENTS.
Only `audit/round030/uniform/` was written. No commit, push, MATLAB run, lower-bound
research, or change to another worker's files was performed.

Let n be an integer, L=log n>=18, l=log L, and T=nL-nl. For t>1 define

```
ell_t(rho)=1+1/(2t)-1/(2rho),       h=ell_t(1)=(1+1/t)/2,
I=integral_1^t ell_t(rho)^(n-1) d rho/rho,
A=(1+1/t)^n,
B=1-h^n+n log t-(n/2)I,
F_n(t)=(sqrt(A)+sqrt(B))^2,
U_n=inf_{t>1} F_n(t).
```

The inherited estimate is C_n<=F_n(t) for the original all-positive-scale,
Lebesgue, centered-cube maximal operator. Also A>0 and B>0: since 0<ell<=1,

```
B=1-h^n+n integral_1^t [1-(1/2)ell_t(rho)^(n-1)] d rho/rho
 >=1-h^n+(n/2)log t>0.
```

We prove, with no hidden n-dependent constants,

```
C_n <= U_n <= T+(3l+1)n/L+(9l^2+1)n/L^2.             (1)
```

## 2. Exact substitution and a general integral lower estimate

For any n>=1 and 0<x<n take t=n/x, and put

```
u=(n/2)(1/rho-1/t),      a=(n-x)/2,
w(u)=(1-u/n)^(n-1),      0<=u<=a.
```

The endpoints rho=t,1 give u=0,a respectively; ell_t(rho)=1-u/n and
d rho/rho=-du/(u+x/2). Consequently the proposed transformation is exact:

```
I=integral_0^a w(u)/(u+x/2) du.                      (2)
```

Define the unnormalized mass and first moment

```
m=integral_0^a w(u)du=1-(1-a/n)^n=1-h^n,
r=integral_0^a u w(u)du.
```

Here h=(1+x/n)/2=1-a/n. Positivity gives 0<m<1 and

```
0<=r<=integral_0^n u(1-u/n)^(n-1)du
      =n^2 integral_0^1 v(1-v)^(n-1)dv
      =n/(n+1)<1.
```

The last integral follows directly by writing v=1-(1-v) and integrating the
two powers; it does not require a probability normalization. Cauchy-Schwarz
applied to sqrt(w/(u+x/2)) and sqrt(w(u+x/2)) yields

```
m^2 <= I [r+(x/2)m],
I >= m^2/[r+(x/2)m] >= 2m^2/(x+2).                 (3)
```

All denominators are positive; the direction of each inequality is preserved.

## 3. Every threshold for the chosen parameter

Choose x=L-3l. The following facts hold for all real L>=18, hence for every
integer n in the claimed range.

First, e<3: for k>=2, k!>=2^(k-1), with strict inequality for k>=3; summing
the resulting geometric series gives e<1+1+1=3. Also

```
e^3 > sum_{k=0}^5 3^k/k! = 92/5 >18,
```

so log18<3. Therefore at L=18,
L/2-3log L>9-9=0, and the derivative 1/2-3/L is positive for L>=18.
It follows that

```
x=L-3l>L/2>=9>1,       0<y:=3l/L<1/2.             (4)
```

Because L>=18>e, we have l>1, hence x+2=L-3l+2<L. Furthermore e^L>2L:
the first three terms of its series give

```
e^L >=1+L+L^2/2>2L,
```

because 1-L+L^2/2=((L-1)^2+1)/2>0. Thus 0<x<L<n/2, t=n/x>1 is finite,
and h=(1+x/n)/2<3/4. Also n=e^L>1+L>=19, so integer n>=20 and in particular
n>=16. For integer n>=16,

```
16(3/4)^16<1,
[(n+1)(3/4)^(n+1)]/[n(3/4)^n]
    =(3/4)(1+1/n)<=(3/4)(17/16)=51/64<1.
```

Induction gives h^n<=(3/4)^n<=1/n. Hence m^2=(1-h^n)^2>=1-2/n, and (3)
implies

```
I >= 2(1-2/n)/(x+2).                               (5)
```

Finally, the constant absorption required below is valid uniformly:

```
1+2/L <= n/L^3.                                    (6)
```

Indeed e^L>=L^5/120, while p(L)=L^3-120L-240 satisfies p(18)=3432>0 and
p'(L)=3L^2-120>0 for L>=18. Therefore
L^5/120>=L^3+2L^2, proving (6). These are analytic statements for every
L>=18, not conclusions from a finite sample of dimensions.

## 4. Budget and cross-term estimates

Using (5), dropping the nonpositive term -h^n, and then x+2<=L and n-2>0,

```
B <=1+n(L-log x)-(n-2)/(x+2)
  <=1+n(L-log x)-n/L+2/L.                          (7)
```

For 0<=y<=1/2,

```
-log(1-y)<=y+y^2.
```

To verify the whole interval, the derivative of
g(y)=y+y^2+log(1-y) is y(1-2y)/(1-y)>=0 and g(0)=0. With y=3l/L and
log x=l+log(1-y), (7) becomes

```
B <= T+(3l-1)n/L+9nl^2/L^2+1+2/L.                (8)
```

For the cross term, log(1+z)<=z for z>=0 yields

```
A=(1+x/n)^n <=e^x=n/L^3.
```

Since x>=9>e, log x>1 and n log x>1. Thus the original definition also gives

```
0<B<=1+n(L-log x)<=nL,
2sqrt(AB)<=2n/L.                                  (9)
```

Expanding F=A+B+2sqrt(AB), combining (8)-(9), and then applying (6) and L>=2,

```
F_n(n/x)
 <= T+(3l+1)n/L+9nl^2/L^2 +1+2/L+n/L^3
 <= T+(3l+1)n/L+9nl^2/L^2 +2n/L^3
 <= T+(3l+1)n/L+(9l^2+1)n/L^2.
```

The admissible-parameter definition of U_n proves (1).

## 5. Strict improvement of the old upper expression

The old expression is T+(6l+5)n/L and was proved on L>=32. The old expression
minus the new expression is

```
(n/L)[3l+4-(9l^2+1)/L].                            (10)
```

Already on L>=18, (4) gives 9l^2/L=(3l)(3l/L)<=3l/2, so the bracket in (10)
is at least 3l/2+4-1/L>0. In particular the new upper expression is strictly
smaller throughout the common range L>=32. This compares explicit proved
upper expressions; it does not assert strict inequality for the unknown U_n.

## 6. Execution evidence and limits of the conclusion

`python3 check.py` was actually executed in this independent clone and returned
PASS; `check_result.json` contains runtime information and source hashes.
It uses only stdlib `fractions.Fraction` for exact arithmetic. It checks the
threshold bases 92/5>18, 16(3/4)^16<1, the induction ratio 51/64<1, and
p(18)=3432; it separately integrates the finite binomial density exactly for
n=1,2,7,16 and three rational x values per n, checking mass, full first moment,
the transformed integrand coordinate and Jacobian. A default-Python SymPy
import was attempted and was unavailable (`ModuleNotFoundError`); no package
installation was needed. No floating-point calculation is used as proof.
The finite checks diagnose implementation and algebra; Sections 2-5 supply the
general proof.

The new result sharpens the explicit upper remainder and extends the stated
sufficient threshold from L>=32 to L>=18. It retains leading coefficient 1.
It proves no lower bound for the true C_n, no globally optimal t, no uniqueness
of a stationary point, and no public-priority claim. The pre-existing lower
side for U_n, if used elsewhere, is logically separate from this audit.

The source repository has missing unrelated Git objects; sparse checkout
reported the missing `cube-lowdim-bounds/AGENTS.md` blob and one unrelated
`centered-cube-upper-bound` blob, leaving a pre-existing out-of-scope deletion
status. All sources named above were fully read and hashed. No repair,
staging, commit or push of that unrelated status was attempted.
