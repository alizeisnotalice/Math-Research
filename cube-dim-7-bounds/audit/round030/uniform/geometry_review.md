# R030 independent audit: genuine geometric deficit and uniform consequence

2026-10-27. **ACCEPT.** The genuine finite-history deficit is valid. Replacing
the single negative-coordinate slab by the union over all negative coordinate
slabs gives the stated D. The same argument closes on every integer n with
log n>=18 and proves C_n<=W_n-D, with W_n the newly audited uniform expression.
This is a new budget for the original operator, obtained from a condition on
the shared selection history. It is not merely a reoptimization of the old
two-parameter scalar budget.

## 1. Sources and scope

The independently read geometry source was
`/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round030/geometry_theory/REPORT.md`,
SHA256 `c4c2ba570855414474cdbdceee71e9a630c1064dbb0960f6a1c663180a7c1ea9`.
That source version states the general two-region lemma and the uniform corner
envelope; its numerical application still uses one slab. The larger union-A
and the general-dimension extension were supplied in the parent task message
and are reconstructed below. This audit does not certify a later changed
source file merely because its path is unchanged.

The inherited finite selection/covering interface was read in this clone's
`cube-dim-7-bounds/revision_v3/cover/audit.md`, SHA256
`d7cf4cda032cf8324983f4db1dbf77c9b025a89f70f485f8214b33ce7460f358`.
The independently audited scalar uniform result is `proof.md` in this audit
directory, SHA256 `8b981f3ba43094a429bd93a51cebc92c3a631a72d884f8d41a36afd6a51fcb4f`.
The original covering theorem remains an explicit dependency; its full proof
is not repeated here. Only this `geometry_review.md` was added in this follow-up.
No lower-bound calculation, MATLAB, commit or push was performed.

## 2. Actual measure, its domination, and corner set inclusion

Fix a finite original selection history, a load point at 0, and the last index
k among the smallest-radius positive contributors. Normalize r_k=1. Every
other positive contributor has j<k and r_j>=1. Write

```
beta=1+1/t, b=a beta, alpha=a/(a+1),
c=(a beta-1)/(a+1)=alpha beta-1/(a+1), a>2t,
v_j=y_j/(a+1), u=y_k/(a+1), m(z)=max(1,||z||_infty),
I_j=Q_j minus union_{i<j}Q_i,
w_j=|I_j|/(2a r_j)^n.
```

Let U index the contributors not captured by the last center, i.e.
||y_k-y_j||_infty>=b r_j. Their radii satisfy r_j<t: otherwise the common
load point gives ||y_k-y_j||_infty<a(1+r_j)<=b r_j, a contradiction.
For every z in I_j/(a+1),

```
||z||_infty<r_j<t,    m(z)<=r_j,
||u-z||_infty>c r_j>=c m(z),    ||u||_infty<alpha.
```

The exact actual increment measure is

```
d mu(z)=alpha^(-n) 2^(-n)
        sum_{j in U} r_j^(-n) 1_{I_j/(a+1)}(z) dz.
```

Its j-th integral is
alpha^(-n)2^(-n)r_j^(-n)|I_j|/(a+1)^n=w_j. Thus mu(total)=sum_{j in U}w_j
exactly. The sets I_j are disjoint and use the same scaling map, so their
indicators never add at a point. With r_j>=m(z), this proves

```
mu <= alpha^(-n)2^(-n)m(z)^(-n)1_{H_u}(z)dz,
H_u={||z||_infty<t, ||u-z||_infty>=c m(z)}.        (2.1)
```

Simultaneously reflect the coordinates of the entire history so that
0<=u_i<alpha. Since c>1 and z_i-u_i<=m(z)<c m(z), only the negative coordinate
alternative can realize the distance condition. Hence

```
H_u={||z||_infty<t, some z_i<=u_i-c m(z)}
    subset H_{alpha e}.
```

This is a pointwise set inclusion, not just monotonicity of an integral. Define

```
H=H_{alpha e},
d lambda=alpha^(-n)2^(-n)m(z)^(-n)1_H(z) dz.
```

Then mu<=lambda and lambda(H)=alpha^(-n)J_a(t), with the accepted finite cap
integral J_a. The corner is used as an envelope; no actual terminal center is
assumed to attain it.

## 3. The common-prefix argument

Take R>1, 0<sigma<1, a measurable A subset H with ||w||_infty<sigma on A,
and P=H intersect {m(z)>R}. Assume the strict guard

```
1+sigma/R+2/(a+1)<alpha beta.                     (3.1)
```

If mu(A)=0, the nonnegative deficit measure lambda-mu has total mass at least
lambda(A). Otherwise some actual uncaptured item j has a positive-measure
increment in A; take an actual w in I_j/(a+1) intersect A. Uncapturedness and
the center-to-increment error give

```
c r_j < ||u-w||_infty < alpha+sigma.
```

Multiplying (3.1) by R and subtracting R/(a+1) gives
Rc>R+sigma+R/(a+1)>alpha+sigma. Thus r_j<R.

For any actual item i whose increment has positive measure in P, choose
z in I_i/(a+1) intersect P. Then r_i>=m(z)>R>r_j. The radius priority rule
therefore forces i<j; there is no tie-order issue. Moreover

```
||v_i-v_j||_infty
 < ||z-w||_infty+(r_i+r_j)/(a+1)
 < m(z)+sigma+2r_i/(a+1)
 < [1+sigma/R+2/(a+1)]r_i
 < alpha beta r_i.
```

In the original coordinates this is ||y_i-y_j||_infty<b r_i. Thus *every*
such i is included in the eligibility test of the *same* selected item j.
Its actual selection condition gives sum_i w_i<=tau over these indices, and
mu(P)<=sum_i w_i<=tau. Consequently

```
sum_{j in U}w_j
 <= alpha^(-n)J_a(t)-min{lambda(A),lambda(P)-tau}. (3.2)
```

It is essential that this step uses actual increments and one actual j; it
does not combine hypothetical extremizers from different histories. A need
not lie in every H_u. A part of A absent from the actual H_u has zero actual
mass, which is itself a deficit in the common envelope lambda.

## 4. Union over all negative faces and an explicit sufficient criterion

For t>7/4 choose

```
R=7/4, sigma=5/(3t),
A=(-sigma,sigma)^n intersect {min_i z_i<-1/t},
D=(5^n-4^n)/(3^n t^n)>0.
```

Here sigma<1 and m=1 on A. Also alpha-c=-1/t+beta/(a+1)>-1/t. Therefore
every point in A belongs to H. The complement of A inside the small cube is
the cube [-1/t,sigma)^n up to null boundaries, so

```
2^(-n)|A|=sigma^n-[(sigma+1/t)/2]^n
         =(5^n-4^n)/(3^n t^n)=D,
lambda(A)=alpha^(-n)D.                            (4.1)
```

This counts the union once; there is no summation of overlapping face slabs.
The guard is exact:

```
alpha beta-1-sigma/R-2/(a+1)
 =1/(21t)-(beta+2)/(a+1)>0
```

whenever a>63t+20. The full cube's m<=1 weighted mass is alpha^(-n), and its
1<m<=R weighted shell mass is alpha^(-n)n log R. Removing these from the cap
gives the safe lower bound

```
lambda(P)>=alpha^(-n)[J_a(t)-1-n log R].           (4.2)
```

Let A_n=(1+1/t)^n and B=B_n(t), the limiting scalar quantities. Suppose

```
B>D,
B-D-1-n log R > sqrt(A_n(B-D)).                   (4.3)
```

Set the fixed positive threshold tau=sqrt(A_n(B-D)). Since J_a(t)->B and
alpha->1, the strict margin (4.3) implies, for every sufficiently large finite
a, lambda(P)-tau>alpha^(-n)D. Thus (3.2) yields

```
sum_{j in U}w_j<=alpha^(-n)[J_a(t)-D]
```

uniformly for every original finite instance. The previously paid/captured
contributors have mass at most tau, and the terminal self-item at most a^(-n).
The inherited finite covering interface therefore gives

```
C_n <=(1+A_{n,a}/tau)
      [tau+a^(-n)+alpha^(-n)(J_a(t)-D)],
A_{n,a}=(1+1/t+1/a)^n.
```

Only now take a->infinity in these already proved scalar operator bounds:

```
C_n <=(1+A_n/tau)(tau+B-D)
     =(sqrt(A_n)+sqrt(B-D))^2.                    (4.4)
```

There is no convergence assumption on selected histories. The sufficient
finite-a lower threshold can depend on n,t and the strict scalar margin, but
is independent of the input function, load point and selection instance.

## 5. Seven dimensions and the new uniform bound

At t0=3391749/524288>6, the enlarged D is still less than 1/1000. The supplied
baseline inequalities B>54/5, A_n B<121/4 give tau<11/2. Since
exp(3/5)>1+3/5+(3/5)^2/2=89/50>7/4, log R<3/5, and hence

```
B-1-7 log R-tau >54/5-1-21/5-11/2=1/10>D.
```

Thus the criterion (4.3) also closes for the union-A upgrade at this t0. This
audit independently checks D and the algebraic implication; it relies on the
previous baseline B and AB certificates for the two displayed inequalities
and does not replace their outward interval computation or certify a final
rounded root value.

For the uniform statement, put L=log n>=18, l=log L, x=L-3l, t=n/x, and

```
W_n=nL-nl+(3l+1)n/L+(9l^2+1)n/L^2.
```

The preceding independent uniform audit proves x>L/2, t>2,
A_n<=n/L^3, B<=nL, and F_n(t)<=W_n. In particular 0<D<sigma^n<1.
For completeness the needed lower bound for B is independently obtained as
follows. In the exact u integral (Section 2 of `proof.md`), use
(1-u/n)^(n-1)<=exp[-(n-1)u/n] and u+x/2>=x/2, and enlarge the positive
integration interval to [0,infinity). This gives I<=2t/(n-1). Since n>=2,

```
B>=n log t-nt/(n-1)
 >=nL-n log x-2n/x
 >=T-4n/L.
```

The function L-log L-4/L increases for L>=18, and log18<3 gives its value
at 18 greater than 18-3-2/9=133/9>8. Consequently B>=8n. Now B>D, and

```
tau_new=sqrt(A_n(B-D))<=sqrt(A_n B)<=n/L.
```

Using only log(7/4)<1, D<1, integer n>=20 and L>=18,

```
B-D-1-n log R >7n-2>n/L>=tau_new.
```

Thus (4.4) applies for every dimension in the claimed uniform range. Finally,

```
(sqrt(A_n)+sqrt(B-D))^2
 =A_n+B-D+2sqrt(A_n(B-D))
 <=A_n+B+2sqrt(A_n B)-D
 =F_n(t)-D<=W_n-D.                               (5.1)
```

This is a strict positive improvement of an operator upper bound. It does not
prove an equality or a lower bound for C_n, does not identify a global optimal
parameter, and does not lower the asymptotic leading coefficient 1. The
geometric deduction uses intermediate eligibility, which was absent from the
old scalar cap envelope; its validity is not contradicted by cap-saturation
examples that omit this common-history constraint.

## 6. Actual local exact checks

An inline Python 3 stdlib Fraction calculation was actually run on this
computer and returned `FINITE_EXACT_CHECKS_PASS`. It checked the union-volume
identity for n=1,2,7,16 at the exact t0; the finite guard identity and
(alpha+sigma)/c<R at rational t=2,t0,1000 with a=63t+21; and D7<1/1000.
The exact union D7 was

```
672299955714101173526150772398233636592156672
/
11293140269216990566386715906741754400397983486063.
```

These finite exact checks support algebra and transcription. The general
set inclusion, density inequality, common-prefix argument, all-L thresholds,
and scalar limiting argument are proved analytically above. No numerical
sample or unproved compactness statement is being used as a general theorem.

## 7. General sigma and exact-tail corollary: ACCEPT

The following is also valid, and supplies the precise legal parameter domain
for a finite numerical candidate search. Fix an integer n>=1, a finite t>1,
and any real sigma satisfying

```
1/t<sigma<1.
```

Use the union A=(-sigma,sigma)^n intersect {min_i z_i<-1/t}. Its exact limiting
capacity and the tail capacity at the limiting admissible radius are

```
D=sigma^n-[(sigma+1/t)/2]^n>0,
P=n log(1/sigma)
  -(n/2) integral_{t sigma}^t ell_t(rho)^(n-1) d rho/rho.
```

For the same A_n=(1+1/t)^n and B=B_n(t), the sufficient conditions are

```
B>D,
P-D>sqrt(A_n(B-D)).                               (7.1)
```

Under these strict conditions the original operator satisfies

```
C_n<=(sqrt(A_n)+sqrt(B-D))^2.                      (7.2)
```

To prove this, sigma>1/t ensures 1<t sigma<t and

```
t sigma-(1+sigma)/beta=(t sigma-1)/beta>0.
```

For any R in (t sigma,t), the limiting corner-cap tail is exactly

```
P(R)=n log(t/R)
     -(n/2) integral_R^t ell_t(rho)^(n-1) d rho/rho.
```

This follows from the same n negative-face shell contributions and the
n positive-face captured subtraction used in the accepted cap formula; here
R>1 so the core does not enter. In particular P(t sigma)=P and P(R) is
continuous. Fix tau=sqrt(A_n(B-D)). The strict second condition in (7.1)
allows one to choose a *fixed* R>t sigma sufficiently near t sigma so that
P(R)-D>tau. No claim is made at the inadmissible equality R=t sigma.

For this fixed R, the limiting guard gap is

```
beta-1-sigma/R=(R-t sigma)/(tR)>0.
```

Consequently all sufficiently large finite a satisfy a>2t and the strict
guard (3.1). Since R>t sigma>1, the argument in Section 3 also gives the
necessary true-radius inequality r_j<R; no assumption r_j=m(w) is used.
The finite-a corner tail has exact capacity

```
lambda(P_R)=alpha^(-n) [n log(t/R)
  -(n/2) integral_R^t ell_a(rho)^(n-1) d rho/rho],
ell_a=alpha ell_t.
```

It converges to P(R), while lambda(A)=alpha^(-n)D. Therefore the strict
margin P(R)-D>tau gives lambda(P_R)-tau>alpha^(-n)D for all sufficiently
large finite a. The true-history cut and the finite covering bound now give
(7.2) by exactly the scalar limiting argument of Section 4. One first fixes
R and then chooses large a; no limit of a family of selected cubes or
interchange of instance-dependent extrema is required.

An additional inline stdlib Fraction run on this computer returned
`GENERAL_SIGMA_RATIONAL_GUARDS_PASS`. For t=2,t0,1000 and
sigma=(1+1/t)/2, it checked the two displayed rational guard identities,
the strict radius comparison, and a finite rational a making (3.1) hold.
This finite check is supplementary; the identities and limiting argument
above prove the corollary for every real parameter in the stated domain.

For numerical candidates, membership in the open parameter domain and both
strict inequalities (7.1) require rigorous verification for each selected
candidate before claiming an operator bound. A floating-point maximum, an
uncertified integral value, or equality in the second condition does not
meet this sufficient criterion. The numerical batch may explore candidates
without treating those exploratory evaluations as proofs.

## 8. Final geometry-source version check

At the parent's final archive request, the geometry author's REPORT.md was
read again at the path in Section 1. Its final checked SHA256 is
`b899300ca0150993c738e7312eb0d8dd9367b32c2f5faf65ed995782deb0dff1`.
Sections 3.1.1-3.1.3 now contain the full-face union, the general-dimension
criterion and uniform corollary, and the general-sigma exact-tail corollary.
Their statements and proofs agree with the independently reconstructed
arguments in Sections 2-7 of this review. In particular Section 3.1.3 retains
the strict payment condition P-D>tau, fixes R>t sigma before taking large a,
and uses the exact finite-a tail. **ACCEPT applies to this locked final
version of those sections as well.**

The previous source hash in Section 1 intentionally records the first
single-slab version that was reviewed; it is not being silently replaced.
This final version check required no new mathematical or numerical run.
The final rounded numerical root bound is owned by the parent's independent
certificate and is outside this report's numerical scope.

## 9. Barrier for this unit-core-cut budget family: ACCEPT

Define the feasible pairs (t,sigma) by the open domain in Section 7 and its
strict inequalities B>D and P-D>sqrt(A E), where E=B-D. Define

```
V_n=inf_{feasible (t,sigma)} (sqrt A+sqrt E)^2,
```

with the extended-real convention inf(empty set)=+infinity. We do not assume
a feasible pair, a minimizing pair, or uniqueness for every small n. Because
P is the mass of a tail subset of the limiting corner cap, P<=B. Thus every
feasible pair satisfies E=B-D>=P-D>sqrt(A E)>0, and therefore E>A. Also
0<D<sigma^n<1. Rationalizing the root difference gives the two-line barrier

```
F_old(t)-F_cut(t,sigma)
 =D[1+2sqrt A/(sqrt B+sqrt E)]<2D<2,
V_n>=U_n-2.                                     (9.1)
```

For the first strict inequality, sqrt B+sqrt E>2sqrt A follows from B>E>A.
For the infimum statement, each feasible F_cut>F_old-2>=U_n-2; taking an
infimum preserves only the non-strict bound displayed. It also holds with
the stated empty-set convention.

For all sufficiently large n, Section 5 supplies the concrete feasible pair
t=n/(L-3l), sigma=5/(3t). Indeed its conservative tail at R=7/4 pays the
full D; the exact Section 7 tail starts at t sigma=5/3<7/4 and is at least
that mass. Thus in this large-n range V_n is finite and

```
U_n-2<=V_n<=W_n.
```

Combining this with the independently accepted fine expansion of U_n and
the definition of W_n yields

```
V_n=T+(3l+1)n/L+O(nl^2/L^2),  n tending to infinity.
```

The fixed 2 is absorbed since nl^2/L^2 tends to infinity. Consequently this
specified unit-core-cut budget family cannot lower the leading coefficient
1, nor change the displayed fine term (3l+1)n/L beyond the proved remainder.
This is a statement about a family of *upper budget expressions*, not a
lower bound for the actual C_n or an obstruction to all covering arguments.
It uses no assumption that any infimum is attained. This addendum is purely
analytic and required no further numerical run.
