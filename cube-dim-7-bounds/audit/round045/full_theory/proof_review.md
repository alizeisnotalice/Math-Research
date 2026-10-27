# R045 soft inverse-column mechanism: independent full-theory audit

2026-10-28. **ACCEPT for each fixed 0<d<=2^(-24), at a common sufficiently-
large-n threshold depending on d.** The locked proof establishes

```
C_n<=F(A0,B0-nd/65536),
x=(1-32d^2)log n, t=n/x.
```

For the explicit fixed choice d=2^(-30), it consequently proves

```
C_n<=n log n-n log log n-2^(-47)n
```

for all sufficiently large n. This is a genuine negative absolute-linear
correction, not an explicit-dimension numerical bound. No explicit
threshold or global optimality is certified. The main coordinator retains
responsibility for integration and final reporting.

## 1. Locked sources, division of audit work and status

Main source read independently in full:
`/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round045/soft_inverse_column/PROOF.md`,
SHA256 `add435fdf2f55993eba9064bd50c493e4c42769c1c48463329c4015be2e66ab9`.

The assigned one-dimensional auditor's formal local proof was also read:
`/Users/zhengzhihao/Desktop/T/cube-r030-roots-audit/cube-dim-7-bounds/audit/round045/one_coordinate/LEMMA.md`,
SHA256 `3fa91ce60a9ef05df270d1fa85a314de6a8573f8223fe1ea69f8433d6b6caa39`.
Its scope is the exact J formula, baseline moments, pinned coordinate,
piecewise q derivative/absolute continuity and uniform matching-product
bridge. This audit does not duplicate its tool execution; it independently
checks the matching of those hypotheses to the full construction and all
remaining operator links. The local lemma alone does not certify the
operator theorem; Sections 4-9 below close the rest.

The earlier interim interface note is superseded by this complete review.
Only `audit/round045/full_theory/` is written. No author source or frozen
earlier proof was changed, and no production batch, parameter search,
n=7 work, shared MATLAB or Git operation was performed. No numerical
diagnostic is used as evidence for the general-dimensional theorem.

## 2. A bounded source factor belongs in the kernel

The ambient cap measure lambda and actual increment measure mu remain
the original ones. The source factor

```
b(z)=min(1,theta/G(z)), theta=exp(-12L d^2)
```

is multiplied into the post-selection kernel, not into an unbounded
actual algorithmic weight. The normalized side lengths in G are positive,
so the inverse and b are Borel and well defined. For row weights
0<=a(w)<=1 as in the general preliminary lemma, actual eligibility
gives row_mu<=a(w)tau. With a positive degree bound M, the same deficit
identity therefore gives the displayed weighted-energy inequality. The
application has a=1 and M=32k>0, so no zero-degree division occurs.

Most importantly, b<=1 preserves the original eligibility bound tau.
There is no theta^(-1) factor in actual payment. All inverse factors
are used only to analyze the ambient weighted row.

## 3. Pointwise column cap and source geometry

At source radius rho=tq, R=q+d-h>1 and
g=d-h+zeta*t*q. Every matched target radius belongs to [R,R+2h].
For each coordinate the complete normalized target overlap decreases
with target radius, and the side length at R is at least 1/2. The
complete-column envelope is consequently

```
C_raw(z)<=(4nh/R)G(z).
```

This includes the full 2n faces, their factor 1/2, the deleted-product
factor at most two and matched log-width at most 2h/R. Thus

```
b(z)C_raw(z)<=(4nh/R)min(G,theta)<=4nh*theta=4k.   (3.1)
```

This is a pointwise bound, not a source average. Keeping b in the
kernel is essential: changing only the ambient source measure would
not multiply a column integral over the target variable by b.

At fixed d>0 and sufficiently large n, h=epsilon0/(Ld),
epsilon0=2^(-20), has h/d tending to zero; zeta*t=o(h).
Targets lie below 2 and sources above 4. All source q satisfy q>1/2,
R>1, 0<g<=d and rho>=4R. The limiting cap contains both supports;
their radial separation makes them disjoint.

## 4. The local inverse lemma covers the actual matching extension

For each target radius r, the exact baseline q0=r-d+h has R=r.
Actual matching q lies in [q0-2h,q0]. Even when q0 falls just outside
the source band, the entire extended interval has q,R>=1/2,
rho>=4R, 0<g<=d and r-R<=2h eventually. The source coefficient
zeta*t is constant while q varies. These are precisely the hypotheses
of the formal local lemma, whose d<=1/16 domain includes the author's
much smaller fixed d.

At the baseline, Fubini cancellation gives E_W J=1 for
W uniform on [-r,r]. Its log mean and variance bounds are
-6d^3/rho0^2<=E log J<=0 and Var log J<=12d^3/rho0^2.
The target face's pinned coordinate contributes a deterministic
logarithm in [-d/rho0,0]; it is not given another uniform law.
The full target-face log sum consequently has mean and variance
O_d(L^2/n), tending uniformly to zero over its fixed radial interval.
Chebyshev gives the baseline event with probability 1-o(1).

The K_w count has mean at most 1+2nd and variance at most 2nd;
its 4nd cutoff has probability 1-o(1). No independence between the
two good events is needed. The target corner-captured part is
exponentially small, hence

```
a_star=lambda(A)>=nd/16                          (4.1)
```

eventually. This lower bound remains valid with d fixed: its annular
log-length is bounded below by a positive constant times d, rather
than incorrectly replacing it by an exact asymptotic d/2.

For a given good target, the local derivative bound is a.e.; its
absolute continuity explicitly handles the capture-onset kink. Summing
and integrating over the whole extended q interval gives

```
|sum log J(q)-sum log J(q0)|<=1280Ldh=5/4096.
```

Together with the 1/8 baseline condition this implies the required
product between 1/2 and 2 for *every* matched q. No finite sample of
q, no interchange of uncountably many random events, and no assumption
that q0 remains inside the source band are used.

## 5. Full inverse angular rows retain the pinned source coordinate

At a fixed source face, the inverse integrals of free coordinates
factor into the J_i. The pinned coordinate has denominator
ell_*=1-g/(2R), independent of the face sign. Its capture count in
coordinate i is N_i=1+1_{|w_i|<mathfrak a}. Hence the exact full inverse
angular coefficient is

```
U(w,rho)=(product J_i)/(2ell_*) * sum_i N_i/J_i.
```

This is a shell coefficient of order n; integrating a matched logarithmic
radial interval of width h gives an order-nh row. On the actual domain
J_i lie between 1/2 and 2, ell_* between 3/4 and 1, and the product
between 1/2 and 2. The source's bounds

```
n/8<=U(w,rho)<=8n                                (5.1)
```

therefore have safe constants. The fixed source face has not been
treated as an independent uniform free coordinate.

## 6. Clipping loss is controlled under the correct inverse probability

Normalize each captured source face by density capture/G. Its free
coordinates are independent with densities proportional to
capture/ell_R, with J_i>=1/2 as normalizers. The pinned source
log-loss is at most d. For any free coordinate,

```
E_inv[-log ell_R]
 <=2 E_full[(-log ell_R)/ell_R]
 <=8g^2/(9Rrho)<=d^2/rho.
```

The final step uses actual R>1, not merely the weaker R>=1/2 allowed
in the local derivative lemma. The factor 8/9 follows from both
ell>=3/4 and -log ell<=(4/3)(1-ell), together with the exact full
excess mean. Thus, on each face and their inverse-row mixture,

```
E_inv[-log G]<=d+nd^2/rho<=d+2Ld^2<=3Ld^2
```

eventually, including the simple necessary eventual guard Ld>=1.
Markov therefore bounds the inverse-row probability of G<theta
by 1/4. An ordinary raw-source small probability would not imply
this statement; the source uses the correct inverse law.

On G>=theta, b=theta/G; everywhere b<=theta/G. Therefore (5.1)
and that tail bound give full soft angular coefficients between
3n theta/32 and 8n theta. Integrating matching radii yields full
rows between 3k/32 and 32k.

Only after clipping is the corner complement subtracted by its
unweighted mass, using b<=1. Its relative cost is bounded by a
constant times

```
(d/h)exp[-L(d/32-12d^2)] ->0.
```

The exponent coefficient is strictly positive for every fixed
d<=2^(-24), while d/h is only polynomial in L. This leaves

```
3k/64<=R_lambda(w)<=32k, C_lambda(z)<=4k.         (6.1)
```

No old unweighted energy is held fixed while its column bound shrinks.
Both ambient sides have actually been recomputed at the new scale k.

## 7. Energy deficit and original finite histories

Choose tau<=3k/128. Equations (4.1) and (6.1) imply

```
Delta>=[a_star(3k/64-tau)]/(32k)
      >=3a_star/4096>=3nd/65536.                 (7.1)
```

The weaker proposed deduction delta=nd/65536 therefore has fixed margin.
The energy identity uses only mu<=lambda and disjoint supports; actual
mu may be concentrated. Its ambient maximum degree is the correctly
weighted one in (6.1).

For every real a>=8n^4, the finite actual envelope is c_a*lambda on
both supports, c_a=alpha^(-n)>=1. The kernel and b are fixed independently
of a, so energy scales by c_a^2 and degrees and target mass by c_a.
The finite-envelope deficit is no smaller than (7.1).

Actual target items represented in A have radius below 3 and source
items represented in P above 4. Strict radius priority makes every
included source earlier than the same target. The explicit guard
kappa+2/(a+1)<alpha*beta transfers increment capture to its actual
center. Because b<=1, each partial source increment contributes at
most its original frozen weight, and that target's single eligibility
test pays the row by tau. There is no unbounded inverse factor in
this real-history payment and no modified selection weight.

Thus the proof holds for all original histories and all sufficiently
large n, uniformly in a>=8n^4. No selector convergence or distribution
hypothesis on the actual increments is introduced.

## 8. New scalar budget and general-input/all-scale bridge

Set tau=sqrt(A0(B0-delta)). Since x lies in [L/2,L], the accepted
integral bounds give B0>8n eventually, whereas delta<n. Also
A0<=n exp(-32d^2L), B0<=nL, so

```
tau/k<=sqrt L/h * exp(-4d^2L)->0.
```

This ensures tau<=3k/128 at a threshold depending on the fixed d.
Captured and terminal-self weights retain their original bounds tau
and a^(-n), while the uncaptured bound is reduced by delta. The same
finite witness-cover argument consequently gives

```
C_n<=[1+(beta+1/a)^n/tau]
       [tau+a^(-n)+alpha^(-n)J_a(t)-delta].
```

This concerns arbitrary original L1 inputs and all positive common
cube radii; the soft weight is used only in the post-selection proof.
The compact strict-level-set and inner-regularity bridge is unchanged.
For fixed n, taking only the scalar a-to-infinity limit yields the
proposed F(A0,B0-delta) operator bound.

Finally x=(1-32d^2)L gives
B0<=T+64nd^2+1. The remaining A0 and cross term are o(n) at fixed
d. With d=2^(-30), delta=2^(-46)n and 64nd^2=2^(-54)n=delta/256.
The other positive terms are eventually below the positive remaining
margin delta/2-2^(-54)n. Hence the full cost is at most
T-delta/2=T-2^(-47)n. The changed-x cost has been paid explicitly;
the old-x W_n formula is not reused at the new parameter.

## 9. Final interpretation

The one-dimensional exact cancellation, its actual matched-radial
application, inverse-weighted clipping loss, recomputed ambient energy,
finite-a history and complete scalar/operator bridges are all closed.
The conjectural status of the interim note is therefore resolved for
the two locked sources above. No numerical batch, including independent
cross-face diagnostics, is needed for this proof.

The result improves the asymptotic negative correction to an absolute
multiple of n. It retains leading coefficient 1 and the same -n log log n
second term. It gives no explicit numerical dimension threshold, C_n
lower bound, globally optimal d, or immediate dominance over every older
bound at its original valid dimensions. Until thresholds are quantified,
applicable finite-dimensional bounds should continue to be combined
by taking their minimum.
