# R035 antidiagonal energy cut: independent full-scope audit

2026-10-28. **ACCEPT, for a common sufficiently-large-n threshold.** The
proposed energy mechanism proves

```
C_n<=F(A0,B0-nd/65536)<W_n-nd/65536,
C_n<=T-(nd)/131072,
d=sqrt(log L/L), L=log n, T=nL-n log L.
```

No explicit numerical threshold is certified here. In particular neither
this review nor the reviewed theorem establishes a claim at log L>=512.
The proof applies to general L1 inputs, all positive common cube radii,
and every real sufficiently large finite covering parameter a before its
scalar limit. It does not require any nonconcentration hypothesis on the
actual increment measure mu. No substantive defect was found.

## 1. Reviewed source and scope

The source read in full, after direct communication with the author, is
`/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round035/energy_theory/PROOF.md`,
SHA256 `f7d98eb21b89d83264053fda2c594387387facfdcbdbec7cf83a41e5dd55f0bd`.
This report independently reconstructs the whole new mechanism rather than
using the source's status label as evidence. Its dependencies are the
previously accepted original finite-cover/load interface, actual increment
domination, and the R030 W_n expression. Frozen R033/034 baselines were not
edited or rerun. Production calculations and the separate small-algebra
audit remain with their assigned owners; this review performs analytical
verification only and does not repeat that numerical work.

Only `audit/round035/energy_theory/` is written. No parameter search, n=7
calculation, shared MATLAB, or Git operation was performed.

## 2. Parameters and simultaneous asymptotic domain

Use the source's exact parameters

```
l=log L, x=L-3l, t=n/x, beta=1+1/t,
d=sqrt(l/L), h=1/(100Ld), I=[d/4,3d/4],
zeta=n^(-3), kappa=beta-zeta,
k=n h L^(-1/4)=n/[100 L^(3/4)sqrt l].
```

For all sufficiently large n: d tends to zero, x/L tends to one,
h/d=1/(100l) tends to zero, t tends to infinity, and
zeta*t/d tends to zero. Thus target radii r=1+u, u in I, lie below 2,
source radii rho=tq=t(1-v), v in I, exceed 4, and the target/source
annuli are disjoint. Also 1<kappa<beta<2 and a>=8n^4 implies a>2t.
Every condition invoked below is an eventual inequality or convergence
with absolute constants. There are finitely many such conditions, so their
thresholds may be replaced by a single common sufficiently-large-n threshold.

The matching condition |u+v-d|<h has, for each fixed u or v in I, partner
interval length between h and 2h once h<|I|. The antidiagonal center d-u
or d-v stays in I, including its endpoints. Source logarithmic radial
length is between h and 4h, while target logarithmic radial length is at
most 2h. These are genuine uniform endpoint bounds, not interior-only
matching approximations.

## 3. Correct ambient face law and good-set mass

On a fixed maximal cube face at radius R, the remaining n-1 coordinates
are uniform on [-R,R] and the weighted face mass is dR/(2R). The fixed
coordinate contributes d, not a random cube-coordinate average, to the
exceedance statistic. This makes both centers
m_R=d+(n-1)d^2/(2R) exact. The free-coordinate first and second moments
are d^2/(2R) and d^3/(3R); their use in the source is correct.

For targets the variance bound is O(nd^3), and the tolerance is t/100.
The resulting probability O(L^2d^3/n) tends uniformly to zero. The K_w
Markov estimate loses at most one quarter of each face. Therefore at
least half of every target face is statistically good for large n.
The target annulus has full mass at least nd/4. Its corner-captured
part is O(nd(3/4)^(n-1)), which is negligible. Thus the cap-and-statistics
target set A has

```
a_star=lambda(A)>=nd/16.                         (3.1)
```

This uses all target faces. It does not assume that only negative faces
carry target mass, and no factor 1/2 can silently be transferred from a
negative-face-only model.

For sources, the S_z variance is at most (2/3)Ld^3, and its fixed 1/100
tolerance gives bad probability at most (20000/3)Ld^3. On each face,
K_z=1+Bin(n-1,2d/rho), with mean at most 5Ld and variance at most 4Ld.
The cutoff 16Ld then gives bad probability at most 4/(121Ld). The source
annulus has full mass at most nd. Thus its total statistical bad mass is
at most

```
(20000/3)n l^2/L+(4/121)n/L=o(k).                 (3.2)
```

Indeed the two ratios to k are bounded by absolute multiples of
l^(5/2)L^(-1/4) and sqrt(l)L^(-1/4), both tending to zero. This is a bound
on the whole bad source set under the full source-shell measure, not a
conditional probability given a particular target or capture event.

The corner-captured source mass is O(nd exp(-Ld/32))=o(k). The exponent
bound is uniform because v>=d/4, q<=1, x>=3L/4 and (n-1)/n>=3/4
eventually. Consequently, for *every* row, bad-source and cap-complement
removal costs o(k), simply because the capture indicator is at most one.
This fully pays for the loss even if all removed mass is concentrated
in that row. No conditional independence or spreading hypothesis is used.

## 4. Deterministic row and column product estimates

Let e=r-q-d and b=e+zeta*rho. Matching gives |b|<=h+zeta*t=o(d).
For a fixed good target w, the exact full-source side length is

```
L_i_row=1-(|w_i|-(q-zeta*rho))_+/(2rho)
       =1-(|w_i|-(r-d-b))_+/(2rho).
```

Changing the threshold affects only the coordinates counted by K_w,
because the lower of the two thresholds is strictly above r-2d.
The total exceedance change is bounded by 12nd(h+zeta*t). Hence, with
T_w=sum_i(1-L_i_row), its distance from nd^2/(4r*rho) has uniform limsup
at most 13/200. The fixed-face mean correction is O(d/rho), so it does
not alter this constant. The shared principal center satisfies

```
nd^2/(4r*rho)=l/4+o(1)
```

uniformly: x/L=1-O(l/L), rq=1+O(d), and l*d+l^2/L tends to zero.
The maximum row loss is O(d/t), so the sum of squared losses is at most
O((d/t)l)=o(1). Deleting any coordinate changes the sum by o(1).

For a fixed good source z, the exact target-side length is

```
L_i_col=1-(|z_i|-(kappa*rho-r))_+/(2r)
       =1-(|z_i|-(rho-d-b))_+/(2r).
```

Again threshold changes affect only the K_z coordinates. Its sum-distance
from the same principal center has uniform limsup at most 17/200.
The maximum column loss is O(d); the squared-loss sum is O(dl)=o(1),
and deletion changes the sum by o(1). All side lengths remain positive
because kappa*rho>rho>4>r.

Applying the elementary logarithmic product bounds now proves both
deleted-coordinate products lie between

```
e^(-2)L^(-1/4) and e^2 L^(-1/4).                 (4.1)
```

The factor e^2 leaves fixed slack beyond the displayed limiting errors.
These are deterministic consequences of the two good-statistics conditions,
uniform over all good points, matched radii and deleted indices. The
unconditional face distribution was used only for good-set masses, not
substituted for these deterministic bounds.

## 5. Ambient degrees and their normalization

For rows, each source coordinate has at least one captured face and at
most two; its exact coefficient is 1/2 times its deleted-coordinate
product. The second face indicator is |w_i|<(kappa-1)rho. Thus the full
matched-row shell coefficient lies between
(n/2)e^(-2)L^(-1/4) and n e^2 L^(-1/4). Integrating the matching interval
and then subtracting (3.2) and the source cap complement gives uniformly

```
r_star=e^(-2)k/4<=R_lambda(w)<=4e^2 k.             (5.1)
```

The raw lower bound before removal is e^(-2)k/2; the whole-source o(k)
loss can be made less than e^(-2)k/4 at a common threshold.

For columns, each target coordinate has at most two captured faces, with
second-face indicator |z_i|<kappa*rho-r. The exact coefficient is again
1/2 times its product. The matched target log-length is at most 2h,
so restricting to good targets and the cap can only decrease mass:

```
C_lambda(z)<=2e^2 k.                             (5.2)
```

Taking M=4e^2k bounds both ambient degrees. They are degrees of the
original cap envelope restricted to A and P, not degrees of the unknown
actual measure mu. This distinction is indispensable and is correctly
maintained in the proposed proof.

## 6. Energy-to-total-deficit argument

For disjoint A,P, a kernel 0<=K<=1 supported in A times P, and mu<=lambda,
let d0=lambda-mu and Delta=d0(total). Direct expansion gives

```
E_lambda,lambda-E_mu,mu
 =E_d0|A,lambda|P+E_mu|A,d0|P
 <=M[d0(A)+d0(P)]<=M*Delta.                      (6.1)
```

The second degree is bounded by the lambda-column degree since mu<=lambda.
Disjointness is what allows the final sum of deficits to be bounded by
one total deficit. No assumption on the concentration of mu is needed.
If actual history gives row integral at most tau for mu-almost every
target, E_mu,mu<=tau*mu(A)<=tau*a_star. Thus

```
Delta>=(E_lambda,lambda-tau*a_star)/M.            (6.2)
```

Since k/(n/L)=L^(1/4)/(100sqrt l) tends to infinity, tau<=n/L is eventually
at most r_star/2. Together with (3.1) and (5.1), this yields

```
Delta>=e^(-4)nd/512>nd/65536.                    (6.3)
```

The last strict comparison follows from e<3, hence e^4<81<128. The
counterexample to free addition in the source has energy discrepancy
N(1-tau) and column degree N, giving exactly its true deficit 1-tau.
It correctly illustrates why an ambient maximum degree is required,
and is not asserted to be an actual cube history.

## 7. Finite-a transfer to actual histories

Set alpha=a/(a+1), c=alpha*beta-1/(a+1), a>=8n^4. The previous actual
increment domination gives mu_a<=lambda_a on the original finite corner
envelope H_a. Pointwise on the entire limiting cap,

```
alpha-c*m-(1-beta*m)=[(beta+1)m-1]/(a+1)>0.
```

Thus A,P lie in H_a and lambda_a is exactly c_a*lambda on them, with
c_a=alpha^(-n)>=1. No total-variation convergence is needed. Their
ambient upper degrees scale by c_a, their energy by c_a^2 and their
target mass by c_a. Applying (6.2) to lambda_a gives

```
Delta_a>=[c_a E_lambda,lambda-tau*a_star]/M
        >=[E_lambda,lambda-tau*a_star]/M>nd/65536. (7.1)
```

The *actual* row bound required here follows pointwise in a real history.
An actual uncaptured item j represented at w in A has
c*r_j<||u-w||<=alpha+(1+3d/4)<3, hence r_j<3. Every represented source
item i has r_i>=rho>4, so strict radius priority forces i<j. Its normalized
center displacement is less than
[kappa+2/(a+1)]r_i<alpha*beta*r_i. All relevant source items are therefore
paid in this same j's single eligibility condition. Matching further
restricts them, so its row integral is at most tau. The disjoint actual
increments identify the item j for mu_a-almost every w; no measurable
selector theorem, common hypothetical history, or extra density premise
is invoked.

## 8. General L1 and all-positive-radius conclusion

Let delta=nd/65536 and choose the fixed dimension-dependent threshold
tau=sqrt(A0(B0-delta)). The accepted estimates B0>8n, B0<=nL and
A0<=n/L^3 give positivity and tau<=n/L once n is large enough. This choice
is independent of a, input f, its level, load point and selection history.
At every load point with a positive contribution, the normalized argument
above proves the same finite bound on uncaptured weight. Captured and
self weights are bounded by tau and a^(-n), exactly as in the original
finite covering proof. Consequently its constant is at most

```
[1+(beta+1/a)^n/tau]
[tau+a^(-n)+alpha^(-n)J_a(t)-delta].
```

This holds for every real a>=8n^4 and every finite true-witness family.
The existing compact-subset-of-strict-level-set argument and Lebesgue
inner regularity therefore apply to arbitrary L1 inputs and all positive
common cube radii; no truncation or compact-support assumption is added.

Only after these operator inequalities are proved is a taken to infinity
in their scalar constants. No selected histories need converge. The
limit yields F(A0,B0-delta); its improvement over F(A0,B0) is strictly
greater than delta, and the accepted W_n upper expression proves the
first claimed comparison.

Finally W_n-T=O(nl/L)=O(nd^2)=o(nd). It is therefore at most delta/2
at a common sufficiently-large-n threshold, yielding the second displayed
theorem. This is an existence threshold argument, not a quantified
log L>=512 certificate.

## 9. Final interpretation

The full theorem, uniform good-set product estimates, ambient degree
payment, finite-a histories and general-input/all-scale bridge are accepted
for the locked source version. The mechanism improves the negative
correction's order beyond the previously certified unit-core and typical
single-shell results. It still retains leading n log n coefficient 1.
No lower bound for the actual C_n, global parameter optimum, or impossibility
of a stronger future covering argument is asserted. All production
calculations and numerical threshold certification remain outside this
analytical audit's scope.
