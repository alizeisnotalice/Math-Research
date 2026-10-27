# R042 second-moment cut: independent full-theory audit

2026-10-28. **ACCEPT for all sufficiently large integer n.** The complete
proof establishes

```
C_n<=F(A0,B0-2^(-31)n L^(-1/3)),
C_n<=T-2^(-32)n L^(-1/3),
L=log n, T=nL-n log L.
```

The first bound uses x=L-L^(1/3), not the older x=L-3log L. Its scalar
cost is essential to the second conclusion and is checked below. This
does improve the negative correction's asymptotic order, while retaining
the n log n leading coefficient 1. No explicit numerical dimension
threshold, global parameter optimum or C_n lower bound is certified.

## 1. Locked sources and scope

Main source read independently in full:
`/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round042/second_moment_cut/PROOF.md`,
SHA256 `722764654cd15a8bf7fb18a799c298b54ced2f0188f5554ccae173240ddf90a2`.

The separate one-dimensional moment audit was also read:
`/Users/zhengzhihao/Desktop/T/cube-r030-parameters/cube-dim-7-bounds/upper/round042/second_moment_audit/PROOF.md`,
SHA256 `425a383952f341f42ff0dee9651924b00fff61e958c2957f355e8cc9fdc72d00`.
Its stronger constant 48 and L>=64 moment interface are consistent, but
the main argument needs only its separately proved constant 144. The
operator theorem is not claimed at L>=64 merely because that moment
sublemma is valid there.

The full chain reviewed is: target good mass, raw row energy, the pointwise
full-column envelope Y, its fixed-face second moment, high-column energy
deletion, ambient-degree deficit, actual finite histories, general L1 and
all positive common radii, and the new scalar scale. The accepted original
finite covering framework and increment domination remain dependencies;
they are not replaced by numerical data. Only `audit/round042/full_theory/`
is written. No production run, optimization, n=7 calculation, MATLAB or
Git operation was performed.

## 2. Geometry and the raw row lower bound

Use s=L^(1/3), d=1/s, h=d^2/100, p=exp(-s/4), M0=n h p,
x=L-s, t=n/x, beta=1+1/t, zeta=n^(-3), kappa=beta-zeta.
Targets have radius r=1+u and sources radius rho=tq=t(1-v), with
u,v in I=[d/4,3d/4] and |u+v-d|<h. At a common sufficiently-large-n
threshold, targets lie below 2, sources above 4, r and rho are separated,
q>1/2, 1<kappa<beta<2 and h+zeta*t<d/2. Matching retains a v-interval
of length between h and 2h even at I's endpoints; its source logarithmic
length is between h and 4h.

On each target face, the fixed coordinate contributes d and the free
coordinates have exact excess mean d^2/(2r) and second moment d^3/(3r).
Thus the stated m_w is exact, not a full-cube substitute. Chebyshev with
t/100 loses probability O(L^2d^3/n)=O(L/n), uniformly tending to zero;
the K_w Markov bound loses at most one quarter. At least half each
target face remains eventually. The annulus mass is at least nd/4,
and its captured-corner part is exponentially negligible. Consequently

```
nd/16<=a_star=lambda(A)<=nd.                      (2.1)
```

For a good target, threshold perturbation b=e+zeta*rho is supported only
on the coordinates counted by K_w. The exceedance-sum error bound is
the source's exact inequality, with uniform limiting error at most
13/200 after normalization. The principal loss center is

```
nd^2/(4r*rho)=(s-d)/(4rq).
```

Since |rq-1|<=d and rq>=1-d, its distance from s/4 is at most
(1+d)/[4(1-d)], tending to 1/4 rather than zero. This bounded
compensation was retained. Each individual loss is O(d/t), so the
squared-loss sum is O(1/t) and deleting one coordinate changes the
logarithm by o(1). It follows uniformly on all good targets, matched
radii and deleted indices that the side product is between e^(-2)p
and e^2p. The available two-sided source-face count is between n and
2n, with exact normalized coefficient 1/2 per face.

Before cap restriction the row bounds therefore have the stated scale
M0. The whole source corner complement has mass
O(nd exp(-cLd))=o(M0): here Ld=s^2, whereas p has exponent -s/4.
Subtracting this whole mass is uniform in the target, since the kernel
is bounded by one. No source-statistics good set or global bad-source
concentration argument is used. Hence

```
e^(-2)M0/4<=R_raw(w)<=4e^2M0<36M0,
E_raw>=e^(-2)nd M0/64>c_E nd M0,
c_E=2^(-10).                                    (2.2)
```

This distinguishes the valid cap-complement removal from the older
global bad-source mass deletion, which is not available at this scale.

## 3. Y bounds the complete column, pointwise

For each fixed source z put r_-=q+d-h>1 and a_i=kappa*rho-|z_i|>0.
Every matched target radius belongs to (r_-,r_-+2h), possibly cut down
further by the target annulus. The normalized overlap with its coordinate
interval [-r,r] is exactly

```
ell_i(r)=min(1,(1+a_i/r)/2).
```

Because a_i>0, this function is nonincreasing in r. Every side length
at r_- is at least 1/2. Thus a deleted-coordinate product at any actual
radius r is at most twice the full product at r_-. The complete target
shell has at most 2n faces, each contributing the product times 1/2,
so its coefficient is at most 2n times that full product. The logarithmic
radial interval has length at most 2h/r_-. Therefore

```
C_raw(z)<=Y(z)=(4nh/r_-) product_i ell_i(r_-).     (3.1)
```

The direction does not follow from unnormalized lengths, which can
increase; it follows from the normalized decreasing functions above.
If r_- lies just outside the target annulus, using the larger interval
still only enlarges the upper bound. The proof is pointwise in rho and
never differentiates r_-(rho). Target statistics and cap restriction
can only lower the column. Y and the later cutoff are Borel functions.

## 4. Full-source conditional L2 bound, including its pinned coordinate

Let g=d-h+zeta*rho. Eventually zeta*t<=h/2, so 0<d-h<=g<=d and both
g<rho and g<r_- hold. At a fixed maximal source face, the pinned side is
ell_pin=1-g/(2r_-). Each other absolute coordinate is independently
uniform on [0,rho], yielding the exact moment

```
E ell^2=1-g^2/(2r_-rho)+g^3/(12r_-^2rho),
E_face Y^2=(4nh/r_-)^2 ell_pin^2 (E ell^2)^(n-1). (4.1)
```

The fixed face is not treated as another uniform coordinate. Discarding
its squared factor is valid only upward. Nor is (E ell)^2 substituted
for E ell^2; the cubic term is a constant-order compensation at Ld^3=1.

For clarity, the main source's coarser uniform estimates are sufficient.
Since r_-q<=1+3d/4 and g>=d(1-d/100), reciprocal and product bounds give

```
(n-1)g^2/(2r_-rho)
 >=s/2-[d/2+s/(2n)+77/200]>=s/2-1
```

eventually. The positive cubic term is at most xd^3/(12q)<=1/6.
Thus log u<=u-1 gives (E ell^2)^(n-1)<=exp(-s/2+7/6)<9p^2.
With r_->1 and ell_pin<=1 this proves the conditional bound 144M0^2.
The whole source annulus has mass at most nd. Integrating over all
faces and radii, then restricting to the cap, gives

```
integral Y^2 d lambda<=144 nd M0^2.               (4.2)
```

All constants are uniform over the displayed source radii and faces.
The stronger independent C=48 proof uses the same pinned-factor
structure and keeps an explicit positive O(1) exponent; it is compatible
with this argument and is not needed to improve the operator constants.

## 5. Column cutoff removes energy in the correct direction

Take H=2^20 and P=P_raw intersect {Y<=H M0}. The removed raw energy is
the integral of *actual raw columns*, not the deleted source mass:

```
E_raw-E_good
 =integral_{Y>HM0} C_raw(z) d lambda(z)
 <=integral_{Y>HM0}Y(z) d lambda(z)
 <=[1/(HM0)]integral Y^2 d lambda
 <=(144/H)nd M0<=(c_E/4)nd M0.                  (5.1)
```

The last rational comparison is H>=4*144/c_E. On the retained set,
columns are at most HM0; rows remain at most their raw bound 36M0,
which is also at most HM0. Thus

```
E_good>=(3c_E/4)nd M0,
max(R_lambda,C_lambda)<=HM0.                    (5.2)
```

This needs no assertion that the high-column source mass is o(M0), and
no spreading assumption on the unknown actual measure mu. Deleting
columns by this energy estimate is precisely what repairs that unavailable
old argument at the new scalar scale.

## 6. Deficit constant and the target-mass normalization

For mu<=lambda, disjoint target/source sets, and an actual row bound
row_mu<=tau, the accepted energy comparison is

```
Delta>=(E_good-tau*lambda(A))/max(R_lambda,C_lambda).
```

The term lambda(A) is bounded *above* by nd from (2.1). Therefore
tau/M0<=c_E/4 implies tau*lambda(A)<=(c_E/4)nd M0, and (5.2) gives

```
Delta>=c_E/(2H) nd=2^(-31)nd.                    (6.1)
```

The upper target mass, rather than its lower bound, is the one required
for this subtraction. The lower bound was used earlier for raw energy.
The dimensions of every energy, degree and mass term match.

## 7. Finite-a actual histories and full operator scope

For every a>=8n^4, the accepted actual-increment measure mu_a is dominated
by its finite corner envelope lambda_a. On the limiting cap, and hence
on A and retained P, lambda_a=c_a lambda exactly, c_a=alpha^(-n)>=1.
The kernel uses the same cutoff Y, matching and kappa as in the limiting
estimates; finite-a degrees therefore scale by c_a, energy by c_a^2,
and target mass by c_a. This gives a finite-envelope deficit at least

```
[c_a E_good-tau*a_star]/(HM0)
 >=[E_good-tau*a_star]/(HM0)>=2^(-31)nd.          (7.1)
```

The actual eligibility hypothesis is preserved. A target increment
point w can belong only to a true uncaptured item j with radius below
3, by c*r_j<=||u-w||<alpha+r<3 and c>1. Every represented source item
i has r_i>=rho>4, so radius priority forces i<j. The explicit guard
kappa+2/(a+1)<alpha*beta converts the increment capture into
bQ_i capture of the same actual center y_j. Integrating only a portion
of an item's increment contributes at most its full frozen weight;
that j's single eligibility condition thus bounds the actual row by
tau. Radial matching and the Y restriction only reduce that row.

No actual measure density beyond mu_a<=lambda_a, conditional source
law, measurable selector, or history convergence is needed. The
dimensionless construction applies at each normalized load point.
All its constants and tau depend only on n and the chosen dimension
parameters, not on f or the load point.

Adding the captured budget tau and terminal self-weight a^(-n), then
using the already proved finite-cover inequality, preserves arbitrary
L1 inputs and every positive common radius. The same compact-strict-
level-set and inner-regularity bridge remains in place. Fixed n is
held first while a tends to infinity in these uniform scalar bounds;
there is no exchange with n-to-infinity or convergence of histories.

## 8. The new scalar choice closes both tau and the final cost

Set delta=2^(-31)nd and tau=sqrt(A0(B0-delta)). The old integral upper
bound with x in [L/2,L] gives B0>=T-4n/L>8n eventually, so positivity
is safe. Also A0<=ne^(-s), B0<=nL. Consequently

```
tau/M0<=100L^(7/6)exp(-s/4) ->0.
```

This proves the required tau/M0<=c_E/4 at a common sufficiently-large-n
threshold. Inserting (7.1) and taking the finite-a scalar limit gives
C_n<=F(A0,B0-delta).

To compare with T, the exact identity x=L(1-d^2) gives

```
B0<=1+n log t<=T+2nd^2+1
```

eventually. Meanwhile
A0+2sqrt(A0B0)<=ne^(-s)+2n sqrt L e^(-s/2)=o(nd).
Since B0-delta remains positive,

```
F(A0,B0-delta)
 <=A0+B0-delta+2sqrt(A0B0)
 <=T-2^(-31)nd+o(nd)
 <=T-2^(-32)nd
```

eventually. This explicitly pays for the changed x; the older W_n
bound at x=L-3log L has not been imported at an inapplicable parameter.

The ratio of this correction's order to the R035 one is
L^(1/6)/sqrt(log L), which tends to infinity. The result is therefore
an order improvement in the negative correction, subject to the stated
asymptotic threshold. It does not assert dominance over every finite-
dimensional certificate; their minimum remains a valid reported bound.

## 9. Final status

All requested links are closed: pointwise Y monotonicity and face factors,
fixed-face L2 moment, energy rather than mass deletion, correct target
upper-mass subtraction, finite-a scaling and actual eligibility, original
general-L1/all-radius bridge, and the new scalar cost. No additional
unproved concentration statement or production evidence is needed.
The theorem is accepted for the locked source versions and for a common
sufficiently-large-n threshold only.
