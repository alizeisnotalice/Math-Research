# R033 typical outer shell: independent adversarial audit

2026-10-28. **ACCEPT for the precise proposal reconstructed below**, conditional
only on the already accepted original finite covering interface. The selected
outer-shell mechanism gives a deficit growing like n/sqrt(log n). It is not
subject to the earlier deficit<2 barrier, whose scope required the removed
region to have unit-core capacity D<1.

## 1. Scope and version identity

This report independently reconstructs the parent-provided R033 proposal and
its explicit finite-a sharpening. At initial review the geometry author's
`upper/round033/typical_shell/` manuscript was not yet present in its clone;
no unobserved author version is certified by this report. The entire precise
mathematical proposal audited is stated here. A final author-source hash can
be appended after that file is available.

The finite actual-increment measure and covering interface are inherited from
`cube-dim-7-bounds/revision_v3/cover/audit.md`, SHA256
`d7cf4cda032cf8324983f4db1dbf77c9b025a89f70f485f8214b33ce7460f358`.
The previously independently audited uniform scalar expression is
`audit/round030/uniform/proof.md`, SHA256
`8b981f3ba43094a429bd93a51cebc92c3a631a72d884f8d41a36afd6a51fcb4f`.
Only `audit/round033/typical_shell/` is written in this task. No n=7 work,
large-dimensional grid, optimization, MATLAB, commit or push is performed.

## 2. Definitions and elementary guards

Fix an integer n with L=log n>=4096. Define

```
y=sqrt L, l=log L, x=L-3l, t=n/x, beta=1+1/t,
eta=1-1/(2y), eta'=1-1/(4y),
s=1+1/(16y), d=s-eta=9/(16y),
delta=n/(128y).
```

We have y>=64, 3L/4<=x<L, t<n, and t>=n/L>=L/2. Indeed at L=4096,
L/4-3log L>1024-36>0, since log4096=12log2<12; its derivative
1/4-3/L is positive thereafter. Also e^L>=L^2/2. It follows that
eta t>4>s, beta<2, and n>=4. The inequalities nd>=2 and nd>=1 used below
follow already from n>=L=y^2 and nd>=9y/16>=36.

For z in R^n put m(z)=max(1,||z||_infty), and write the limiting corner cap
and its weighted volume as

```
D_infty={||z||_infty<t, ||e-z||_infty>=beta m(z)},
d lambda_infty=2^(-n)m(z)^(-n)1_{D_infty}(z) dz.
```

Boundaries may be closed without changing this absolutely continuous measure.
Let A be the *closed* set of w satisfying

```
r=||w||_infty in [1,s],
some coordinate w_i=-r,
S(w)=sum_i (|w_i|-eta)_+ <=4nd^2,
K(w)=#{i:|w_i|>eta} <=8nd.
```

The strict unique-maximum-negative-face version differs only by ties, a
finite union of zero-volume hyperplanes. Using the closed version removes
any compactness ambiguity: the negative-face condition is closed, S is
continuous, and the sublevel set of K is closed because K is lower
semicontinuous. Thus A is compact. No compactness theorem is needed later.
Every A point lies in D_infty: a negative maximal coordinate gives
||e-w||_infty>=1+r>beta r since r<=s<t.

## 3. Correct shell-face normalization and Markov estimate

On a fixed negative maximal face at radius r, the other n-1 coordinates are
independent uniform variables on [-r,r]. The fixed coordinate contributes
r-eta<=d to S and 1 to K. For each free coordinate V,

```
E(|V|-eta)_+=(r-eta)^2/(2r)<=d^2/2,
Pr(|V|>eta)=(r-eta)/r<=d.
```

Therefore ES<=d+nd^2/2<=nd^2 and EK<=1+nd<=2nd. Markov and a union bound
give at least 1/2 good face proportion for every r in [1,s]. The weighted
shell contribution from one face is exactly

```
2^(-n)r^(-n)(2r)^(n-1)dr=dr/(2r).
```

There are n negative faces; overlaps of maximal faces have zero volume.
Consequently

```
lambda_infty(A)>= (n/4)log s>=n/(128y)=delta.      (3.1)
```

For the last step log(1+v)>=v/(1+v)>=v/2 for v=1/(16y)<=1. Thus no Euclidean
surface factor, extra power of r, or missing factor 2 occurs here.

## 4. Uniform lower mass of the w-dependent preceding region

For each fixed w in A define

```
P_beta(w)=D_infty intersect {eta t<rho=||z||_infty<eta' t}
          intersect {||z-w||_infty<beta rho}.
```

For rho in this band, rho>s>=||w||. The normalized length in coordinate i
cut out from [-rho,rho] by the beta-rho cube centered at w is exactly

```
ell_i=min(1,(1+beta)/2-|w_i|/(2rho))
     =1-(|w_i|-rho/t)_+/(2rho).
```

The sum of its side deficits is bounded uniformly by

```
sum_i(1-ell_i)<=S(w)/(2rho)
 <=2nd^2/(eta t)=2xd^2/eta<=81/(128eta)<=81/127<2/3.
```

Thus each product of n-1 side lengths is at least 1/3. Both maximal faces
in coordinate i are captured if |w_i|<=eta<rho/t. Every remaining coordinate
still has at least one captured face; the far endpoint constraint is
automatic from rho>s. Hence at least 2n-K faces are captured. Since
8d=9/(2y)<=9/128<1/8,

```
2n-K>=15n/8.
```

Each captured face has coefficient (1/2) times a product of n-1 normalized
side lengths in the weighted d rho/rho integral. Thus the full-shell
capture coefficient is at least 5n/16.

For the original corner e, its captured coefficient on 1<rho<t is exactly
(n/2)ell_t(rho)^(n-1), with ell_t(rho)=1+1/(2t)-1/(2rho). In this band,

```
1-ell_t(rho)=(1-rho/t)/(2rho)>=1/(8yt),
(n-1)(1-ell_t(rho))>=9y/128.
```

The second inequality uses (n-1)/n>=3/4 and x>=3L/4. Moreover
exp(9y/128)>=exp(9/2)>50; the first six positive series terms alone sum to
80993/1280>50. Therefore the corner-captured coefficient is at most n/100.
Subtracting it from the full-shell captured mass gives a safe coefficient
at least 5n/16-n/100>n/4 for P_beta(w). Since

```
log(eta'/eta)=integral_eta^eta' dv/v>=eta'-eta=1/(4y),
```

uniformly for every w in A,

```
Phi_beta(w):=lambda_infty(P_beta(w))>=n/(16y)=8delta. (4.1)
```

The subtraction direction is correct: the full local capture minus all
corner-captured mass is a lower bound for local capture inside D_infty.

## 5. Explicit shrinkage: availability signs and face coefficients

Put zeta=1/n^3 and kappa=beta-zeta. Because t<n and n>=4,
0<zeta<=(beta-1)/2, so 1<kappa<beta. Define P_kappa(w) by replacing beta
by kappa only in the local distance condition in Section 4. Let Phi_kappa
be its lambda_infty mass. We prove the explicit bound

```
0<=Phi_beta(w)-Phi_kappa(w)<=3/(2n)<=delta          (5.1)
```

uniformly for every w in A.

For a local capture parameter v in [kappa,beta], the normalized side is
ell_i(v)=min(1,(1+v)/2-|w_i|/(2rho)), which lies in (0,1] and is Lipschitz
with constant 1/2 in v. At rho>s, one face per coordinate is always captured.
The other is captured exactly when |w_i|<(v-1)rho. To check the signs,
the positive face requires w_i>-(v-1)rho and the negative face requires
w_i<(v-1)rho; their opposite far endpoint inequalities hold automatically.
If w_i=0 both faces stay captured. Otherwise at most one face is lost on
shrinking beta to kappa, on a radial interval of logarithmic length at most

```
log[(beta-1)/(kappa-1)]
 =-log[1-zeta/(beta-1)]<=2t zeta.
```

Each face has normalized coefficient at most 1/2. Summing over n coordinates
therefore bounds the availability loss by nt zeta. For surviving faces,
the product inequality for factors in [0,1] bounds each (n-1)-factor product
loss by (n-1)zeta/2. With at most 2n faces and the same coefficient 1/2,
the total product loss is at most

```
[n(n-1)zeta/2] log(eta'/eta)<=n^2 zeta/2,
```

since log(eta'/eta)<1. Adding both losses and using t<n yields
at most 3n^2 zeta/2=3/(2n). Filtering both nested local captures by the fixed
set D_infty can only decrease their difference measure. Finally n>=2L=2y^2
implies n^2>=4y^4>=192y, proving (5.1). The exceptional equality radii
rho=|w_i|/(v-1) are finite for fixed w,v and have zero d rho measure; no
entire face at such a radius was incorrectly assigned zero face area.

Combining (4.1) and (5.1) gives Phi_kappa(w)>=7delta uniformly. This explicit
estimate eliminates all need for joint-continuity, compactness, total-variation
limits, or measurable-selection arguments in the finite-a transfer.

## 6. Finite-a domination is pointwise, and the variable-P cut is valid

Use the original normalization
alpha=a/(a+1), c=alpha beta-1/(a+1), with a>2t. Reflect the complete history
so that 0<=u_i<=alpha. The actual uncaptured increment measure mu satisfies

```
mu<=lambda_a,
d lambda_a=alpha^(-n)2^(-n)m(z)^(-n)1_{D_a}(z)dz,
D_a={||z||<t, some z_i<=alpha-c m(z)},
lambda_a(total)=alpha^(-n)J_a(t).
```

This is the actual increment measure: its j-th integral is precisely w_j,
and its density bound follows from disjoint actual I_j and r_j>=m(z).
The limiting cap uses threshold 1-beta m(z). Their difference is exactly

```
alpha-c m(z)-[1-beta m(z)]
 =[(beta+1)m(z)-1]/(a+1)>0.
```

Hence D_infty subset D_a pointwise. In particular A subset D_a and

```
lambda_a(A)=alpha^(-n)lambda_infty(A)>=delta,
lambda_a(P_kappa(w))=alpha^(-n)Phi_kappa(w)>=7delta. (6.1)
```

Choose the fixed tau=sqrt(A_n(B_n-delta)), where A_n=(1+1/t)^n and
B_n=B_n(t). The previously proved elementary estimates at this x give
B_n>=T-4n/L>=8n>delta, A_n<=n/L^3, and B_n<=nL. Thus tau is real positive
and tau<=n/L<=n/(64y)=2delta.

For every real a>=8n^4 we have a>2t, c>1 and

```
alpha beta-kappa-2/(a+1)
 =1/n^3-(beta+2)/(a+1)>0.                         (6.2)
```

Set Delta=lambda_a(total)-mu(total)>=0. Suppose Delta<delta. Then (6.1)
gives mu(A)>0. A finite real history therefore supplies an actual
uncaptured item j and an actual increment point w in I_j/(a+1) intersect A.
Its true radius is bounded by

```
c r_j<||u-w||<=alpha+s<3,
r_j<3.
```

For any actual item i having positive increment measure in P_kappa(w), take
a corresponding actual point z. Its radius satisfies r_i>=m(z)=rho>4,
so strict radius priority forces i<j. The centers obey

```
||v_i-v_j||
 <||z-w||+(r_i+r_j)/(a+1)
 <kappa rho+2r_i/(a+1)
 <=[kappa+2/(a+1)]r_i<alpha beta r_i.
```

Thus all these items are captured in this single selected j's eligibility
test. It follows that mu(P_kappa(w))<=tau<=2delta. But domination and the
total deficit give

```
mu(P_kappa(w))>=lambda_a(P_kappa(w))-Delta>6delta,
```

a contradiction. Therefore Delta>=delta for *every* original history and
every a>=8n^4. P is allowed to depend on the one actual chosen w; no global
selector or common P for all w is being assumed.

## 7. Operator upper bound and its large-L consequence

The finite uncaptured budget is at most alpha^(-n)J_a(t)-delta. Inserting
this into the inherited covering interface, including the captured mass
tau and terminal self-mass a^(-n), gives

```
C_n <=[1+(1+1/t+1/a)^n/tau]
       [tau+a^(-n)+alpha^(-n)J_a(t)-delta].
```

Only now take a to infinity in scalar operator bounds that already hold
for every a>=8n^4. Since J_a->B_n, the same fixed tau yields

```
C_n<=(sqrt(A_n)+sqrt(B_n-delta))^2
    <F_n(t)-delta<=W_n-delta,                    (7.1)
W_n=T+(3l+1)n/L+(9l^2+1)n/L^2.
```

Strictness follows because A_n>0 and delta>0: rationalizing the cross-term
difference gives F_n(t)-F_cut=delta+2sqrt(A_n)(sqrt(B_n)-sqrt(B_n-delta))
strictly greater than delta.

For L>=2^30, every summand of

```
256[(3log L+1)/sqrt L+(9(log L)^2+1)/L^(3/2)]
```

is decreasing: the only nontrivial factors log L/sqrt L and
(log L)^2/L^(3/2) have negative derivatives for log L>2 and >4/3
respectively. At L=2^30, use sqrt L=32768 and log L=30log2<30 to bound
the expression above by the exact rational number

```
256[91/32768+8101/32768^3]
 =97710514085/137438953472<1.
```

Thus W_n-T<=n/(256sqrt L). Since delta=n/(128sqrt L), (7.1) gives

```
C_n<T-n/(256sqrt L),    L>=2^30.                 (7.2)
```

This improves the actual operator upper bound beyond the unit-core deficit
family. It does not lower the n log n leading coefficient, prove a lower
bound for C_n, or establish optimality of this new mechanism.

## 8. Actual bounded local computation

`python3 check_constants.py` was actually run and returned PASS; the
`constant_receipt.json` records its runtime and exact outputs. It uses only
stdlib Fraction. It verifies the small rational constants, the exponential
series partial sum, the final large-L base comparison, and representative
rational side-length and finite-a guard identities. It neither evaluates
the exponentially large dimension n=exp(4096) nor samples a high-dimensional
shell. General-dimensional validity comes from Sections 2-7, not these
finite checks.

## 9. Author-source version lock

After the independent reconstruction was completed, the geometry author's
`/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round033/typical_shell/PROOF.md`
became available and was read in full. Its checked SHA256 is
`762a36853cd9996fe213bdeba336240dee449ee5a224430ab0498c455ddfc543`.
Its closed outer-shell definition correctly uses 1<=||w||_infty<=s rather
than the vacuous lower bound 1<=m(w). Its variable-P lemma, per-negative-face
Markov normalization, shell subtraction, explicit kappa shrinkage and
finite-a guard all agree with this report's independent proof. **ACCEPT
also applies to that locked author version.**

The author's final L>=2^30 estimate uses the sharper log2<7/10 instead of
this report's sufficient log2<1 bound; its positive exponential partial-sum
argument and derivative signs are correct. No numerical rerun was needed
to compare these two sufficient threshold arguments.
