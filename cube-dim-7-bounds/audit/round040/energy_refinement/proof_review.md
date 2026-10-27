# R040 integrated-energy refinement: independent scope audit

2026-10-28. **ACCEPT for fixed parameters and their common sufficiently-
large-n threshold.** For each fixed 0<c<sqrt2, 0<epsilon<1/2 and
0<gamma<c(1-2epsilon), the complete operator inequality at the stated
old x=L-3log L scale is valid. The simpler T-gamma correction is valid
only after using a strictly larger intermediate mechanism coefficient,
as the source correctly specifies. No explicit 128/256-type threshold
is established.

## 1. Source and review scope

Source read in full and checked against its actual identity:
`/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round040/energy_refinement/PROOF.md`,
SHA256 `5e52e7ce1f2cb9d094e6d6e8e32622c0dabcd73fb0dcbab16add999dcf13b2af`.
The previous comparison audit checked its actual x and coefficient domain;
this report completes an independent analytical check of the entire chain.

Dependencies are the previously accepted finite covering/actual-increment
framework, ambient energy-deficit lemma and same-parameter W_n bound.
No production or old numerical check was rerun, and no new search,
n=7 work, MATLAB or Git operation was performed. Only
`audit/round040/energy_refinement/` was written.

## 2. Fixed-parameter domain and exact good-face centering

Use L=log n, l=log L, x=L-3l, t=n/x,

```
d=c sqrt(l/L), h=1/(L d l),
I=[epsilon*d,(1-epsilon)*d],
D=(1-2epsilon)*d, p=L^(-c^2/4), k=n h p.
```

For fixed c>0 and epsilon>0, d tends to zero, h/d=1/(c^2l^2) tends
to zero, t tends to infinity and zeta*t=o(h). The target/source annuli
are disjoint, with r<2 and rho>4; the actual kappa perturbation has
1<kappa<beta and the same finite-a capture guard as the original energy
construction. These are uniform radial statements for fixed parameters.

On any target face, the fixed coordinate contributes d and the other
coordinates are independent uniform variables. Therefore m_w is the
exact conditional mean, with variance at most nd^3/(3r). The tolerance
t/l yields bad probability O(L^2d^3l^2/n)=o(1). Likewise K_w has mean
(2+o(1))nd and variance O(nd), so its 4nd cutoff has bad probability
O(1/(nd))=o(1). The exponentially small target-corner capture does not
change this conclusion. Thus the good target's radial density is

```
n[1-o(1)] du/(1+u), uniformly on I,
a_star=nD[1+o(1)].                              (2.1)
```

The target uses all 2n faces. No negative-face-only factor one half
may be inserted in this density or later column estimates.

## 3. Whole-source deletion is small relative to k precisely for c^2<2

The fixed source-face center m_z is also exact. Variance O(Ld^3) at
tolerance 1/l gives source bad-S probability O(Ld^3l^2). The K_z mean
is (2+o(1))Ld with variance O(Ld); its 8Ld cutoff loses O(1/(Ld)).
Multiplying by the full annulus mass O(nd) gives

```
b_stat=O_c(nL d^4 l^2+n/L)=O_c(nl^4/L+n/L).
```

The corner-captured source contribution is
O(nd exp(-a*epsilon*Ld)). Since

```
k=n/[c l^(3/2)L^(1/2+c^2/4)],
```

the statistical ratios to k are respectively
O_c(l^(11/2)L^(-1/2+c^2/4)) and
O_c(l^(3/2)L^(-1/2+c^2/4)). Both vanish because c^2<2 is fixed.
The corner ratio vanishes as well: its exponential in epsilon*Ld
dominates all remaining powers and logarithms for fixed epsilon>0.
Consequently the total deleted-source mass b_n=o(k), uniformly over
the whole radial source band.

This is the ambient full-source mass. For a row its lost capture is
at most b_n because K<=1; integrated over good targets the loss is
at most a_star*b_n=o(nDk). No conditional independence of source
statistics given target capture, or actual-mu concentration, is assumed.
The estimate need not be uniform as c approaches sqrt2 or epsilon
approaches zero, and the theorem does not claim such uniformity.

## 4. Relative products and full face coefficients

For a matched radial pair, threshold movement b=e+zeta*rho has size
at most h+zeta*t=o(d), and affects only coordinates counted by K_w
or K_z. The row error is bounded by

```
[t/l+4nd(h+zeta*t)+2d]/(2rho)=O(1/l)+o(1),
```

while the column counterpart is
[1/l+8Ld(h+zeta*t)+2d]/(2r)=O(1/l)+o(1).
Here d*h=1/(L*l); both shifts really vanish after their appropriate
normalization. The common principal loss is

```
nd^2/(4r*rho)=c^2 l/4+O_c(l*d+l^2/L)=c^2 l/4+o(1).
```

Maximum row and column losses are O(d/t) and O(d), while total losses
are O_c(l). Their squared-loss sums and any removed coordinate therefore
have vanishing logarithmic effect. This proves, uniformly on every
good point, matching radius and deleted index, the *relative* product
formula p[1+o(1)]. It is stronger than a fixed e^(+/-2) enclosure and
is justified by the smaller tolerances and bandwidth.

Each row's missing second faces are counted by K_w<=4nd=o(n), and
each column's missing second faces by K_z<=8Ld=o(n). Combining these
counts with the same deleted-coordinate products and coefficient 1/2
per face gives full-shell coefficients np[1+o(1)] in both directions.
Restricting targets to the good cap only lowers a column upper bound.

Hence the same kernel has ambient maximum degree at most

```
M_n<=2n h p[1+o(1)]=2k[1+o(1)].                 (4.1)
```

The source and target logarithmic Jacobians are 1+O(d), so replacing
ordinary matching length by its logarithmic length introduces only the
uniform relative o(1) already displayed, not a fixed missing factor.

## 5. Integrated strip energy and endpoint accounting

The interval I has length D and its endpoint sum is d. For h<D, the
antidiagonal strip in I times I has exact area

```
2hD-h^2=D^2-(D-h)^2.
```

It is legitimate to integrate the full-source coefficient over good
targets: that coefficient estimate is uniform at every good target,
not an unconditioned angular average substituted on a conditional set.
Using (2.1) and the uniform radial Jacobian gives raw pre-source-deletion
energy

```
E0=n^2 p(2hD-h^2)[1+o(1)].
```

Since h/D=1/[(1-2epsilon)c^2l^2] tends to zero, and the deletion loss
is o(nDk), the retained energy satisfies

```
E_lambda,lambda>=2k nD[1-o(1)].                 (5.1)
```

This uses the same kernel as (4.1). It does not insert the shortest
endpoint row length h into a global energy estimate; that would lose
a permanent factor one half. The exact integrated area restores it,
and the leftover h^2 term has vanishing relative size.

## 6. Deficit, actual finite histories and general input

The target mass is nD[1+o(1)]. For any tau<=n/L,

```
tau/k<=c l^(3/2)L^(-1/2+c^2/4)->0.
```

The ambient bipartite energy lemma, (4.1) and (5.1) then yield

```
Delta>=(E_lambda,lambda-tau*a_star)/M_n
       >=nD[1-o(1)].                            (6.1)
```

Only the known ambient degrees are used. The actual measure may be
arbitrarily concentrated under its domination constraint.

This remains valid for every original finite history and every real
a>=8n^4. On both good supports, lambda_a=c_a*lambda exactly,
c_a=alpha^(-n)>=1. Energy scales by c_a^2, degrees and target mass by
c_a; hence the finite-a deficit is no smaller than the same quotient
in (6.1). The o(1) bounds do not depend on a in that range.

An actual target item represented in the target good annulus has
radius below 3; all represented source items have radius above 4.
Strict radius priority therefore puts each source before the same
actual target. The center-displacement guard
kappa+2/(a+1)<alpha*beta puts all such sources into that target's
one original eligibility test. Radial matching only removes sources,
so the actual row integral is at most tau. There is no selector limit
or independently chosen per-row history.

For fixed gamma<c(1-2epsilon), nD=c(1-2epsilon)n sqrt(l/L) and the
fixed positive gap in (6.1) gives deficit at least
delta=gamma*n sqrt(l/L) eventually. B0>8n, B0<=nL and A0<=n/L^3
make tau=sqrt(A0(B0-delta)) positive with tau<=n/L for large n.
The original finite-cover/witness integration then gives

```
C_n<=[1+(beta+1/a)^n/tau]
       [tau+a^(-n)+alpha^(-n)J_a(t)-delta].
```

For each fixed n only its scalar constants are sent to the limit
a to infinity. The same compact-strict-level-set and inner-regularity
bridge retains arbitrary L1 input and every positive common cube radius.
Thus the full F(A0,B0-delta) bound, followed by its strict W_n-delta
comparison at the *same old x*, is sound.

## 7. T-gamma corollary and its precise threshold scope

The same-gamma mechanism only gives C_n<W_n-gamma*n sqrt(l/L);
the positive W_n-T remainder may not be dropped at that coefficient.
The source avoids this error. For any desired fixed gamma<sqrt2,
choose c and epsilon with c(1-2epsilon)>gamma, then choose an intermediate
gamma_plus strictly between them. Apply the mechanism at gamma_plus
and absorb W_n-T=o(n sqrt(l/L)) into gamma_plus-gamma. This proves
the stated T-gamma conclusion after a common sufficiently-large-n
threshold depending on the fixed choices.

The specialization c=1, epsilon=1/8 has deficit
(3/4-o(1))n sqrt(l/L); its coefficient 1/2 is therefore valid eventually.
No explicit log log n>=128 or >=256 threshold follows from this report
or the reviewed source. No uniform starting point as the parameters
approach their endpoints is asserted. The endpoint gamma=sqrt2 and
global optimality over covering methods remain unproved.

## 8. Final verdict

The relative Gaussian products, bad-source/k ratio, integrated area,
ambient maximum degree, common-history finite-a transfer and same-scale
scalar comparison all pass. This is a valid older complete-cost
interface, even though the later R042 negative correction eventually
has larger asymptotic order. It gives no C_n lower bound and retains
leading coefficient 1. No source correction is required.
