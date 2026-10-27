# Fresh second review: cube-cap geometry and its exact formulas

Date: 2026-10-27. Reviewer: fresh Sol geometry audit. Decision: **ACCEPT within the scope stated below**.

## Scope and evidence boundary

The sole mathematical source read for this review was
`upper/round027/cube_cap_theory/proof.md` at repository HEAD
`79f8acb20c418e60c5196a3c830103253e34c900`, SHA256
`6ecbf72fa9213677dc99e2e2afaa510df240b3d1b08047366a8d0103d010c80f`.
No previous review or PASS decision was read. This review independently checks
the integral H, its coordinate monotonicity, the corner evaluation J_a, the
finite-to-limit formulas, positivity of B_n for every t>1, and the n=7 sixth-power
expansion. It is not a separate approval of the operator selection/covering
argument, the numerical upper certificate, or a claimed lower bound. No new
bound, optimization, parameter search, or research continuation was performed.

The proofs below are the basis of the decision. The accompanying exact Python
Fraction checks are finite audit evidence, not substitutes for those proofs.

## 1. The object, including equality conventions

For an integer n>=1, t>1, c>1, and |u_i|<=alpha<1, put

    m(z)=max(1,||z||_infinity),
    H(u;c,t)=2^(-n) integral_{||z||_infinity<t}
                  1_{||u-z||_infinity>=c m(z)} m(z)^(-n) dz.

Call the complementary condition `||u-z||_infinity<c m(z)` captured. The
inequality is strict because H includes equality. Coordinate reflections
preserve the cube, the weight, and Lebesgue measure, and send u_i to -u_i.
Thus H depends only on |u_i|; below 0<=u_i<=alpha<1.

On the unit cube, the captured interval in coordinate i is

    (-1,1) intersect (u_i-c,u_i+c).

Its length is `L_i(1)=min(2,1+c-u_i)`: its upper endpoint is 1 since
u_i+c>1; its lower endpoint is max(-1,u_i-c). Consequently the ball contribution
is exactly `1-2^(-n) product_i L_i(1)`, and is nondecreasing in each u_i.
All endpoints have zero n-dimensional measure here.

## 2. Shell decomposition and the dangerous equality case

For 1<rho<t, assign each z with ||z||_infinity=rho to its unique coordinate
with largest absolute value. Ties belong to a finite union of sets
`z_i=+/-z_j`, hence have zero n-dimensional measure. For each signed face use
coordinates `(rho,z_1,...,omit z_i,...,z_n)` with `z_i=+/-rho` and other
coordinates in `(-rho,rho)`. The absolute determinant is 1. This establishes
the volume decomposition directly; it does not differentiate a moving volume.

The other-coordinate captured lengths on an i-face are

    L_j(rho)=min(2rho,(1+c)rho-u_j).

They are positive and nonincreasing in u_j. On the positive face,
`|rho-u_i|=rho-u_i<c rho`, so that face is captured in coordinate i.
On the negative face the strict condition is

    rho+u_i<c rho  iff  u_i<(c-1)rho.

Thus the captured cross-sectional face measure is exactly

    S(rho,u)=sum_i [1+1_{u_i<(c-1)rho}] product_{j!=i} L_j(rho).  (A)

This formula also works when n=1, with each empty product equal to 1.
When u_i=(c-1)rho, the entire negative i-face is uncaptured in coordinate i.
It would be wrong to declare its face measure zero merely because it is a cap
equality. Formula (A) correctly excludes it. For any fixed u and c>1 this can
occur at at most n radii; those exceptional radii have zero d rho measure, so
their treatment cannot change the full n-dimensional integral. No positive
measure equality region has been discarded.

It follows, with the ball-shell boundary rho=1 also harmless, that

    H=1-2^(-n) product_i L_i(1)
      +2^(-n) integral_1^t rho^(-n)
           [2n(2rho)^(n-1)-S(rho,u)] d rho.                 (B)

The tempting replacement of S by the derivative of the captured intersection
volume is false. Already for n=1 in the active truncation regime that volume
is `(1+c)rho-u`, with derivative 1+c, whereas the captured face count is 1.
The source explicitly uses (A), avoiding this error.

## 3. Coordinate monotonicity

For fixed rho, each term of S is nonnegative. As u_k increases, its k-th term
has the nonincreasing factor `1+1_{u_k<(c-1)rho}`; every other term has the
nonincreasing factor L_k. Therefore S is coordinatewise nonincreasing, even
at its downward threshold jumps. The shell integrand in (B) is consequently
coordinatewise nondecreasing. The same was shown for the ball term. Hence

    H(u;c,t)<=H(alpha*1;c,t).

This is a maximum over the closed comparison box. Actual selected centers may
lie only in the open box; no realization of the corner is needed. The proof
also covers u_i=0, a completely empty cap (e.g. c>=1+alpha), and n=1.

## 4. Independent dimension-one boundary check

For 0<=u<1, positive z never contributes to H. Direct integration along the
negative half-line gives the exact formula

    H_1(u;c,t)=(u-c+1)_+/2
              +(1/2) log(max(1,min(t,u/(c-1)))).            (C)

The first term is the surviving segment in (-1,0); on z=-rho<-1 the condition
is `rho<=u/(c-1)`. Formula (C) is nondecreasing in u. In particular it is zero
at u=c-1, and remains well-defined when the cutoff is at rho=1 or rho=t.
Evaluating (B) for n=1 gives exactly (C). Six fixed rational cases covering
zero u, those two endpoint thresholds, an interior threshold, the full active
shell, and an empty cap were checked in exact arithmetic by the attached
script. Three fixed face calculations in n=1,2,7 additionally include equality
faces and mixed clipped/unclipped coordinate intervals; all matched (A).

## 5. Finite a: the corner and every normalization factor

Use the source parameters t>1 and a>2t,

    alpha=a/(a+1),  b=a(1+1/t),  c=(b-1)/(a+1).

Then `c-1=(a/t-2)/(a+1)>0`, and

    alpha/(c-1)=a t/(a-2t)>t.

At u=alpha*1, for every rho in [1,t] each negative face is uncaptured; exactly
one face per coordinate is captured. Also `L(rho)=(1+c)rho-alpha<2rho`.
Put `ell_a(rho)=(1+c-alpha/rho)/2`. Formula (B) becomes

    J_a(t)=H(alpha*1;c,t)
          =1-ell_a(1)^n+n log t
             -(n/2) integral_1^t ell_a(rho)^(n-1) d rho/rho.

Indeed the full shell contributes
`2^(-n)rho^(-n) 2n(2rho)^(n-1)=n/rho`, while the n captured faces contribute
`2^(-n)rho^(-n)n L(rho)^(n-1)=(n/2)ell_a(rho)^(n-1)/rho`.
Thus neither a factor 2 nor a shell Jacobian is missing.

The affine change of variables from the source unpaid increments is
`y=(a+1)z`, so their normalized `2^(-n)a^(-n)` density picks up exactly
`P=((a+1)/a)^n`. The finite correction is therefore

    E_a(t)=a^(-n)+P J_a(t).

The self-term a^(-n) is separate from H and must remain at finite a. Nothing
in the shell computation pays for it. This review does not remove it.

One useful check of all finite parameters is the exact identity

    ell_a(rho)=alpha [1+1/(2t)-1/(2rho)].                   (D)

In particular `0<ell_a(rho)<=alpha<1` on [1,t], including the upper endpoint.
The ball term is positive and the shell integrand is at least n/(2rho),
including n=1. Therefore J_a(t)>0.

## 6. Limit and positivity for every t>1

For fixed t and n, (D) gives uniform convergence of ell_a and its integer
powers on [1,t]. Also alpha->1, P->1 and a^(-n)->0. Hence

    J_a(t)->B_n(t)
      =1-[(1+1/t)/2]^n+n log t
        -(n/2) integral_1^t [1+1/(2t)-1/(2rho)]^(n-1) d rho/rho.

Write `ell(rho)=1+1/(2t)-1/(2rho)`. For 1<=rho<=t,
`0<ell(rho)<=1`, with ell(t)=1 and ell(1)=(1+1/t)/2<1. Therefore

    B_n(t) >= 1-[(1+1/t)/2]^n+(n/2)log t > 0.

This proves positivity for all t>1 and all integers n>=1, without relying on
the sign of a rational-plus-log coefficient or a finite sample. There is no
uniform-in-t limiting claim here. The limit passage is for each fixed t.
For n=1 it reduces to `(1-1/t)/2+(1/2)log t`, consistent with (C).

## 7. The n=7 sixth-power formula

Hold t fixed, and put `q=1+1/(2t)`. Expanding in the integration variable rho,

    (q-1/(2rho))^6
      =q^6+sum_{j=1}^6 binom(6,j)q^(6-j)(-1/2)^j rho^(-j).

Integration against d rho/rho from 1 to t yields the coefficient list

    -3, 15/8, -5/6, 15/64, -3/80, 1/384

in front of `q^(6-j)(1-t^(-j))`, respectively. This reproduces every term,
sign, and denominator in source formula (14). It follows exactly that

    B_7(t)=1-[(1+1/t)/2]^7+7 log t-(7/2)I_6(t).

The log coefficient is `7-(7/2)q^6`; it is not necessarily positive. The
source's warning about this is correct. The finite-a formula follows from
the identical binomial calculation with `q_a=(1+c)/2`, `p_a=alpha/2` and
`(-p_a)^j`; its distinct self-term a^(-7) remains in E_a.

The Python run independently computed the six coefficients as exact
Fractions and compared them with the printed list; all agreed. SymPy was
unavailable both in the default and bundled Python, so no claim of a SymPy
run is made. Standard-library exact arithmetic was sufficient for these
finite coefficient identities. No scalar operator bound was evaluated.

## Verdict and recommendations

**ACCEPT** the source cube-cap geometric lemma, coordinate monotonicity,
the finite J_a/P/self-term bookkeeping, B_n formula and all-t positivity,
and the n=7 rational-plus-log expansion. No mathematical correction is
required in this audited part. Publication wording should preserve the
strict captured condition, explicitly distinguish face measure from full
volume measure at cap equality, and keep the finite a^(-n) term until taking
the limit. Identity (D) is an optional simplification for presentation, not
a missing assumption or repair.

Whole-chain validity and the final decimal certificate need their separate
reviews. The present ACCEPT must not be cited as approval of those objects.
Work stops here; no paper was generated and no source file was modified.
