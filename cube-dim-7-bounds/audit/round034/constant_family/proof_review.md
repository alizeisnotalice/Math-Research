# R034 general proof-constant family and C34A: independent audit

2026-10-28. **ACCEPT.** The general-dimension theorem is valid throughout its
stated parameter domain, including its first strict inequality. The fixed
C34A matrix has **12/12 independently verified exact rational certificates**.
This review certifies the twelve displayed candidates and the theorem that
makes them operator bounds; it does not certify global optimality of any
candidate or any numerical optimization run.

## 1. Scope and locked sources

The mathematical source read in full is
`/Users/zhengzhihao/Desktop/T/cube-r030-parameters/cube-dim-7-bounds/upper/round034/constant_theory/PROOF.md`,
SHA256 `44ca63fa21fe3dce5af17d794ff255ec4cd155c9c4194aafb9c775cead7a3c82`.
Its `verify.py` and `verification.json` are hashed for bundle identity in
`matrix_receipt.json`; they were not imported or executed by this audit.

The frozen candidate matrix read is
`/Users/zhengzhihao/Desktop/T/cube-r030-compute/cube-dim-7-bounds/upper/round034/compute/c34a_results.json`,
SHA256 `87ac7a68bf3979f6b654f7954044feac4d0d26031e14dc38152c8e48684d4a3d`.

Dependencies are the accepted original finite centered-cube covering
interface, the R033 actual-measure variable-P lemma and explicit capture
shrinkage, and the R030 scalar upper expression W_n. Parameter-sensitive
steps are reconstructed below. Only `audit/round034/constant_family/` was
written. R033 is unchanged. No n=7 work, new search, MATLAB, commit or push
was performed.

## 2. Precise domain and all-y guards

Let Y>=8 and positive rational a0,b0,c0 satisfy

```
0<a0<b0<Y/2, 0<c0<Y/2, g=b0+c0, gY>=1,
f1=1-4g/Y>0,
f2=1-2g^2/(1-b0/Y)>0.
```

Define

```
KA=c0/[4(1+c0/Y)], CL=f1*f2,
D=CL-(1/2)exp(-9a0Y/32),
KP=(b0-a0)D-1/Y-3/(8Y^3).
```

For 0<q<min(KA,KP), the claim is, for every integer n with
L=log n>=max(64,Y^2), and y=sqrt L, x=L-3log L, t=n/x,

```
C_n<F(A_n(t),B_n(t)-q n/y)<W_n-q n/y.             (2.1)
```

At L=64, log2<3/4 gives L/4-3log L>5/2, and its derivative is positive
thereafter. The log2 bound follows from the positive exponential partial
sum 1+3/4+(3/4)^2/2=65/32>2. Hence 3L/4<x<L and x>1 throughout this range.
Since n=e^L>=L^2/2>=2L,

```
t>=n/L>=L/2>=32, t<n.
```

For s=1+c0/y, eta=1-b0/y, eta'=1-a0/y, d=g/y, all y>=Y give

```
1/2<eta<eta'<1, 1<s<3/2,
eta*t>L/4>=16, rho in (eta*t,eta'*t) implies rho>4>s,
nd=ng/y>=2y/Y>=2.
```

Thus the needed strong radius separation uses n>=L^2/2, rather than the
insufficient weaker lower bound n>=2L alone. The closed good shell set
uses 1<=||w||_infty<=s, a negative maximal coordinate, and
S=sum(|w_i|-eta)_+<=4nd^2, K=#{|w_i|>eta}<=8nd. Its definition and all
actual-increment radius statements are the same as R033 with these constants.

## 3. Shell mass and local capture constants

On each negative maximal face at radius r in [1,s], the fixed coordinate
contributes at most d to S and 1 to K; each free uniform coordinate contributes
at most d^2/2 and d respectively. Therefore ES<=nd^2 and EK<=2nd, using
nd>=2. Two Markov bounds give at least half good points on each negative
face. Each such face has exact weighted contribution dr/(2r), so

```
lambda_infty(A)>= (n/4)log(1+c0/y)
 >= [c0/(4(1+c0/y))]n/y>=KA*n/y.                 (3.1)
```

Every A point belongs to the limiting corner cap because its negative
maximal coordinate has 1+r>beta r for r<=s<t. The finite-a pointwise
inclusion of that cap still gives lambda_a(A)=alpha^(-n)lambda_infty(A).

For any w in A and eta*t<rho<eta'*t, the normalized local side lengths are

```
ell_i=1-(|w_i|-rho/t)_+/(2rho).
```

Their total deficit is at most 2g^2/(1-b0/Y), so every (n-1)-side product
is at least f2. At least 2n-K>=2n(1-4g/y)>=2n f1 maximal faces are captured.
Including the exact factor 1/2 per face gives local full-shell coefficient
at least n CL. This is valid because f1,f2 are positive, so multiplying
their separate uniform lower bounds preserves direction.

The original corner's captured coefficient is
(n/2)ell_t(rho)^(n-1). Its exponent satisfies

```
(n-1)(1/rho-1/t)/2
 >=(1-1/n)x*a0/(2y*eta')>=9a0y/32>=9a0Y/32.
```

Here n>=4, x>=3L/4 and eta'<1 suffice. Subtracting the entire corner-captured
region from the local captured region is a lower bound for their part in
the corner cap. The resulting shell coefficient is at least nD.
The assumption KP>q>0 forces D>0, because b0-a0>0 and both other KP
subtractions are positive. Therefore one may integrate with the lower
band length without reversing an inequality:

```
log(eta'/eta)=integral_a0^b0 dv/(y-v)>=(b0-a0)/y,
Phi_beta(w)>=(b0-a0)D*n/y.                       (3.2)
```

## 4. Uniform finite-a transfer and exact strictness

Take zeta=n^(-3), kappa=beta-zeta. The relation x*n^2>=2 guarantees
0<zeta<=(beta-1)/2, so kappa>1. The R033 explicit face-availability and
side-product estimates apply unchanged since rho>s and
log(eta'/eta)<log2<1. Their loss is at most 3/(2n). Normalizing it by n/y,

```
3/(2n)=[3y/(2n^2)]n/y
 <=[3/(8y^3)]n/y<=[3/(8Y^3)]n/y.                (4.1)
```

This step is uniform in y: it uses n>=2y^2 at the actual y, not at Y only.
For every real a>=8n^4 the finite guard
alpha beta-kappa-2/(a+1)=n^(-3)-(beta+2)/(a+1)>0 holds. An actual item
with increment w in A has radius r_j<alpha+s<5/2, whereas actual items in
P_kappa(w) have radius r_i>=rho>eta*t>16. Strict radius priority places all
such i before the same actual selected j, and the guard makes all their
bQ_i capture its center. Its one eligibility condition pays their total
actual mass. There is no measurable-selector or history-convergence step.

The original scalar estimates at this x give B>=T-4n/L>8n, A<=n/L^3 and
B<=nL. Since q<KA<c0/4<Y/8, the removed budget delta=q n/y<n/8 is safely
positive and less than B. For tau=sqrt(A(B-delta)),

```
tau<=n/L<=(1/Y)n/y.
```

Thus (3.2), (4.1), the finite-a cap inclusion, and the actual variable-P
deficit lemma give

```
lambda_a(A)>=KA*n/y>delta,
lambda_a(P_kappa(w))-tau>=KP*n/y>delta,
actual uncaptured mass<=alpha^(-n)J_a-delta.
```

The old finite covering interface, including tau and the terminal a^(-n)
self-mass, now gives C_n<=F(A,B-delta) after taking only a scalar a-to-infinity
limit. The full original general-L1 and all-positive-common-radius scope is
preserved by the already accepted finite-cover/inner-regularity interface.

For the *first strict inequality* in (2.1), choose a fixed number
q<qstar<min(KA,KP) and apply that non-strict operator bound to delta_star=
qstar*n/y. Positivity is still assured, and F(A,B-delta_star)<F(A,B-delta)
by strict monotonicity in its positive second budget. Merely proving a
non-strict bound at q would not justify this strictness. Finally,

```
F(A,B)-F(A,B-delta)>delta,
F(A,B)<=W_n,
```

proves the second strict inequality in (2.1).

All replacements of y by Y are now explicit: KA's normalized shell capacity
increases with y; the face factor and side-product factor increase; the
subtracted exponential decreases; normalized tau<=1/y and normalized
shrinkage<=3/(8y^3) decrease. No conclusion is inferred solely from a
numerical check at y=Y.

## 5. Independent check of the hand-picked q=1/20 example

For Y=64, a0=1/10, b0=29/100, c0=21/100, exact rational arithmetic gives

```
KA=336/6421>1/20, CL=98301/203872,
sum_{k=0}^6 (9/5)^k/k!=7542629/1250000>6.
```

Therefore exp(-9/5)<1/6, and the coarser rational lower bound

```
KP>(19/100)(98301/203872-1/12)-1/64-3/(8*64^3)
   =12055437337/200414330880>1/20
```

has exact margin 2034720793/200414330880 above 1/20. This example alone
rigorously improves the R033 coefficient on L>=4096, independently of any
optimizer or matrix candidate.

## 6. Fixed 12-row matrix verification and resources

`python3 check_matrix.py` was actually run. It imported no producer helper,
ran no search, and used only standard-library Fraction for certificate
arithmetic. For each row it reconstructed r=9a0Y/32 and

```
Eupper=1/[sum_{k=0}^128 r^k/k!] > exp(-r),
KP_lower=(b0-a0)(CL-Eupper/2)-1/Y-3/(8Y^3).
```

The strict exponential upper direction follows from omitted positive series
terms. Every row passed the complete parameter domain and both positive
factor conditions, 0<q<KA, and q<KP_lower. Its independently computed KA,
KP_lower and the two strict rational margins exactly equal the matrix's
reported fractions. For every row the rational midpoint between q and
min(KA,KP_lower) supplies an explicit valid qstar for the strict theorem.

| Y | L minimum | Certified q (decimal spelling of exact rational) | Producer PID | Objective calls | Wall s | CPU s | Peak RSS bytes |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 14 | 196 | 0.005610906253 | 59046 | 2772 | 0.0396243 | 0.037527 | 76791808 |
| 16 | 256 | 0.011980898994 | 59045 | 2052 | 0.0291787 | 0.028133 | 76824576 |
| 20 | 400 | 0.021616142707 | 59047 | 2520 | 0.0351660 | 0.034767 | 76824576 |
| 24 | 576 | 0.028627506156 | 59049 | 2232 | 0.0318360 | 0.031045 | 76972032 |
| 28 | 784 | 0.034006773928 | 59048 | 2448 | 0.0315147 | 0.030980 | 77086720 |
| 32 | 1024 | 0.038292347062 | 59050 | 2304 | 0.0303180 | 0.030117 | 76857344 |
| 40 | 1600 | 0.04474464185 | 59051 | 2268 | 0.0281335 | 0.027432 | 77119488 |
| 48 | 2304 | 0.049415666219 | 59053 | 2160 | 0.0275735 | 0.027014 | 76693504 |
| 64 | 4096 | 0.055808586202 | 59052 | 2340 | 0.0299300 | 0.029017 | 76595200 |
| 96 | 9216 | 0.063065593838 | 59056 | 2376 | 0.0288713 | 0.028486 | 76562432 |
| 128 | 16384 | 0.067166220283 | 59054 | 2268 | 0.0281196 | 0.027421 | 77135872 |
| 256 | 65536 | 0.074275164113 | 59055 | 2016 | 0.0258157 | 0.025748 | 76644352 |

These producer PIDs and resources are recorded from the locked matrix;
they are not independently observed process timings. It reports twelve
distinct PIDs, workers_max=12 and total batch wall time 0.44844050001120195 s.
The independent audit process itself was PID 59969, wall 0.027897250023670495 s,
CPU 0.047111 s and peak RSS 19316736 bytes. Full unrounded resources and exact
rational inputs are retained in `matrix_receipt.json`.

The Y=256 row has the largest certified q among these twelve rows and applies
only once L>=65536. This comparison is finite and does not assert that its
parameters are globally optimal, even at that Y. Raw optimizer objectives
are not substituted for the certified rational q values.

## 7. Optional analytic ceiling and interpretation

The source's optional family ceiling is also algebraically consistent:
q<min(c0/4,b0(1-2g^2)), 0<g<1/sqrt2. For fixed g the upper envelope equals
Q(g)=g(1-2g^2)/(5-8g^2), whose derivative sign is
5-22g^2+16g^4. The interior maximum occurs at
g^2=(11-sqrt41)/16, yielding
(19-sqrt41)sqrt(11-sqrt41)/320. This is an analytic ceiling for this proof
family, not proof that any displayed finite row attains it and not an upper
limit on improvements obtainable from all covering arguments. The stated
approach by Y tending to infinity, a0 tending to zero with a0Y tending to
infinity, and rational b0,c0 tending to the balancing values is compatible
with all strict domain constraints for large enough Y.

All results remain operator *upper bounds*. Neither the ceiling nor the
matrix gives a lower bound for C_n, a globally optimal covering theorem, or
a leading n log n coefficient smaller than 1.
