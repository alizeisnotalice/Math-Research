# A074 independent proof and certificate review

This review concerns only the original body32 equal-weight Cauchy model,
q=6/5, rho=1/64, and the union D of the 31 middle tubes with mean_B in
(1/4,3/4]. It accepts the previously sealed 043/049/051 geometry and
coefficients as inputs, checks their exact source bindings, and reviews
the genuinely new 074 time estimates. The author manifest and report
are pinned in registration.json. No 072 claim is used here.

The target is the local long norm of H=integral_(1/1000)^infinity p_t dt
on this D. The author does not claim that an upper expression larger
than a sufficient allowance is a lower bound for that norm.

## 1. Actual events and time denominators

For each A, take precisely the 051 event partition E, E_v, E_w. E requires
both |a_ones|>=1/2 and |a_w|>=1/16; E_v requires |a_ones|<1/2; E_w requires
|a_ones|>=1/2 and |a_w|<1/16. They partition parameter space, and their
nonnegative densities sum to p_t. Enlarging E_w for its pointwise upper
bound is harmless. Only full p_t is used for the subsequent mass bound.

The frozen two-column coarea calculation gives coefficient R2/(16pi^6)
with complete factors t^-2/[(1/4+t^2)(1/256+t^2)]. The frozen
D_big=64R2/pi^6 therefore requires exactly division by 1024. No k-dependent
Gram determinant is lost: that cancellation already occurs in the genuine
two-column coarea calculation of 051.

For the complementary terms, the first two denominators retain
h_b0^2+992t^2 and h_b1^2+960t^2. They hold for the retained small parameter,
before extending its nonnegative integral to all real a. The actual
third denominator is |z-a b'|^2+928t^2, not a false uniform gap on this
entire shifted fiber.

Write lambda=ell^2, L=|b'|>0, gamma=ell*t/L+u with u>=0, and
d=|z|^2+ell^2*t^2. The exact 043 average has numerator 1+t/gamma and
denominator d+L^2*t^2+2L^2*t*gamma. The difference between
(1+L/ell)/d and this exact expression has positive denominator and numerator

    L*u*d+(ell+L)*gamma*(L^2*t^2+2L^2*t*gamma) >= 0.

Thus keeping lambda*t^2 in the averaged denominator is valid. For L=0,
the integral is exactly 1/(|z|^2+lambda*t^2). The original (unshifted) pair
gap g_k can now be substituted. The projected norm bound L<=sqrt(6),
the omitted-parameter peak loss (10/9)^3=1000/729, and the full ordered-triple
coefficient Crem all remain necessary and are present. In particular
074 does not replace this quotient density by the original one with no loss.

It follows from the accepted 051 input that the stated F_k=P_k+B_k bounds
the complete p_t on each tube for every positive t. Every denominator
increases with t; each F_k is decreasing. Dropping only positive spatial
parts gives F_k(t)<=C_k*t^-6, with

    C_k=D_big,k/1024+2K/(992*960*928).

## 2. True marginal, all-time primitive, and one whole-union mass

The 049 complete marginal is 2048t^2/(pi^2 Dplus Dminus), where
Dplus=(z+1)^2+(32t)^2 and Dminus=(z-1)^2+(32t)^2. Reflecting only the
other four coordinates divides its mass by 16, giving coefficient 128.
This reflection statement is for p_t, not a restricted event density.

An independent partial fraction decomposition is

    1/(Dplus Dminus)
      =[(z+2)/Dplus+(2-z)/Dminus]/[4(1+a^2)], a=32t>0.

Indeed (z+2)Dminus+(2-z)Dplus=4(1+a^2). Splitting z+2=(z+1)+1 and
2-z=-(z-1)+1 and integrating gives exactly the author's log ratio and
the sum of the two real arctangents. There is no branch jump: a>0 and
both quadratic denominators are strictly positive. The complete time
denominators have not been replaced by their zero-time limits.

For one fixed j in B, its true residual coordinate norm is
sqrt((m-1)/m); m=1 gives zero and the exact interval [1/4,3/4]. The frozen
intervals [L_m,U_m] are therefore legitimate tube mass bounds. No factor m
is needed, because a single coordinate is chosen for each tube.

Every point of the union has a nonempty B excluding coordinate 1. Any
coordinate in that B is within rho of its mean and hence lies in
[15/64,49/64]. A union over the five possible physical coordinates
therefore covers D, even when several choices or tubes cover a point.
Positivity gives the 5*L_star bound. The full positive face gives the
second bound t/[pi(1+1024t^2)]. Their minimum is a valid upper bound on
the mass of full p_t over D. There is no extra factor 31 here, nor an
extra reflection factor applied to a branch.

## 3. Volume, multiplicity, overlap, and q powers

Use zero-sum residuals e on A and z on B and mean u on B. Since x1=1,
x_A=(1-e1)1_A+e. The true Jacobian for the five-dimensional face is
sqrt(km), while |e|^2+|z|^2<rho^2. The box upper constraints on A imply
e1>=e_i. Permuting the A coordinates partitions the radial residual ball
into k congruent maximal-coordinate cones (ties are null); their fraction
is 1/k. For k=1 the residual space on A is zero-dimensional and the
fraction is 1. Dropping the remaining box constraints is a positive
enlargement. Multiplying the four-ball volume (pi^2/2)rho^4 by u-length
1/2 and by sqrt(km)/k proves

    V_k=(pi^2/4)rho^4*sqrt(m/k).

There are binomial(5,k-1) choices containing coordinate 1. The five weights
are (1,5,10,10,5), totaling 31. k and 6-k cannot be identified when the
coordinate 1 constraint is fixed. No disjointness of tubes is used.

On each tube, p_t^q<=F_k^(q-1)*p_t and p_t^q<=F_k^q. Summing nonnegative
tube integrals gives the first two q-power envelopes. On the union,
p_t<=max_k F_k, giving a third q-power bound Fmax^(1/5)*L_D. After taking
the power 5/6, its pointwise factor is Fmax^(1/6). This verifies all
three norm bounds with the same q and the same full p_t.

The author's implementation selects the candidate having the smallest
certified upper endpoint. It need not compute an interval containing the
mathematical minimum: each candidate separately bounds the true norm,
so the selected candidate's upper endpoint suffices. Its saved lower
endpoint is only a lower endpoint of that selected upper-expression
evaluation; it is not a lower bound on the true norm or on the minimum
of all possible certificates.

## 4. Continuous time coverage and the complete tail

The rational endpoints start at epsilon=1/1000 and multiply by 16/15
until the final endpoint is clipped to T=1/4. For t in [a,b], put r=b/a.
Monotonicity gives F_k(t)<=F_k(a). The joint marginal prefactor t^2
is at most r^2*a^2, while each denominator is at least its value at a;
hence L_m(t)<=r^2*L_m(a), including L_star. The full-face bound similarly
is at most r*a/[pi(1+1024a^2)]. The minimum of these valid slab bounds
is valid throughout the whole slab. Thus a saved norm upper endpoint
times the exact width b-a pays a continuous interval, not just a sample.

The finite checker verifies every saved row, the exact recurrence, absence
of time gaps, and accurate rational width aggregation. It does not evaluate
any of the 86 F/joint/log/atan expressions again. Consequently the numerical
dependency on the author's outward Arb evaluation of each rectangle remains
explicit; that code and its positive formulas are reviewed statically and
bound to the saved source. This is a saved-certificate audit, not a second
independent semigroup calculation.

The remaining half-line follows from the volume union and C_k*t^-6:

    integral_T^infinity ||p_t||q dt
      <=(sum_k n_k V_k C_k^(6/5))^(5/6)/(5*T^5).

The audit uses the complete saved positive dyadic V_k and C_k endpoints.
It independently propagates exact Fraction bounds through each root,
checking lower^n<=input^m<upper^n separately at the lower and upper input
endpoints. It also checks each saved C_k contains the exact rational
interval formed from frozen D_big, Crem and F. No old M or short-time
integral is recomputed. A broad independent tail enclosure only needs to
overlap the author's narrower enclosure; the final strict rational budget
is verified using the independent upper tail and exact widths, not by an
invalid claim that one independently widened root interval must contain
another.

Minkowski on bounded positive-time truncations, followed by monotone
convergence, gives the claimed local norm estimate on [epsilon,infinity).

## 5. Exact conclusion and limitations

Subject to passing the fixed six saved checks, the reviewed certificate
has strict long upper 49091761943697/10^17. The selected finite budget
expression itself exceeds the frozen sufficient D+F allowance. This
means this fixed expression does not close the ledger; it does not
exclude the model or lower-bound the actual Green norm.

The saved 051 short upper is only read and added on the same D, giving
495607457962773/10^18 before any q power is placed in the spatial ledger.
The priority-subtracted middle region of 059 is a subset of this D, so
restriction is allowed. Far, global g, and the complete face I still
require their own certificates. This audit supplies no C6/Cn value and
no new semiconcavity theorem. The parent retains final proof responsibility.
