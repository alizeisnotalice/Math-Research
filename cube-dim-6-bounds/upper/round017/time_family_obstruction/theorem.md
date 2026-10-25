# A continuous time interval obstruction and an all-time acceptance theorem

Fix the body32 model only and an off-S2 point x. Put

    K(T)=integral_0^T p_t(x)dt, q(T)=K(T)/T, c=1/4096.

Here K is finite, p_t is positive and smooth for t>0, and K(0+)=0.
The B015 point used here is strictly positive with six distinct absolute
coordinates, so it is off S2. No assertion is made for the augmented
beta=1 frame.

## 1. What a point obstruction means for arbitrary probability mixtures

Let mu be any probability measure on (0,infinity), and define the
normalized positive kernel

    Q_mu(y)=integral K_T(y)/T dmu(T).

Tonelli gives mass one. Its spatial dilations are probability mixtures
of the same positive semigroup's ergodic averages, so they are bounded
by the ergodic maximal function and have weak (1,1) factor at most one.
The cube-domination method therefore has cost 1/(64 essinf_cube Q_mu).

If q(T)<=c0 for every T in a closed bounded interval J contained in
(0,infinity), then every mu supported in J has Q_mu(x)<=c0. This is
also an essential-minimum obstruction, not merely a zero-measure point:
off S2, K_T(y)/T is jointly continuous and locally uniformly bounded for
y near x and T in J. Dominated convergence makes Q_mu continuous near
x, including from the cube interior. Hence essinf_cube Q_mu<=c0.

For a genuinely all-T bound q_x(T)<=c0, the same conclusion for arbitrary
probability mu requires care at unbounded T supports. A sufficient
uniform-neighborhood domination is available here: for a compact
neighborhood E disjoint from S2, p_t(y)<=C_E for t<=1 by the independent
large-triple rank argument, and p_t(y)<=C t^-6 for t>=1. Thus
K_T(y)/T<=C_E for T<=1 and <=G_E for T>=1. Joint continuity and this
uniform bound imply continuity of Q_mu near x by dominated convergence,
for every probability measure mu. Therefore a true uniform all-T
point bound would obstruct every continuous probability mixture as well.

An unnormalized finite positive time measure leads to the same statement
after division by total kernel mass. The relevant probability weighting
for a mixture of K_T itself is proportional to T times the mixing measure.

## 2. A true time comparison from the 32 Cauchy factors

For c_t(u)=t/[pi(t^2+u^2)] and s>=t,

    c_s(u)/c_t(u)=(s/t)(t^2+u^2)/(s^2+u^2)<=s/t.

Choose any fixed six-column basis and solve those variables in the
positive density formula; all 32 Cauchy factors remain present. The
parameter-space pointwise comparison, integrated over the same remaining
26 variables with the same determinant factor, proves

    p_s(x)<=(s/t)^32 p_t(x)       (s>=t),
    p_s(x)>=(s/t)^32 p_t(x)       (0<s<=t).              (1)

These bounds do not assume that p or q is unimodal, nor that p decreases
after the recorded time. They use the true positive Cauchy integrand.

Suppose at t0 the saved exact intervals are
K(t0) in [L,U] and p_t0(x) in [l,u]. Set N=32. Integrating (1) gives

    K(T)<=U-l t0/(N+1) [1-(T/t0)^(N+1)]  (T<=t0),   (2)
    K(T)<=U+u t0/(N+1) [(T/t0)^(N+1)-1]  (T>=t0).   (3)

After division by T, each right side is A/T+B T^N, with B>=0.
If U-u t0/(N+1)>0, both coefficients A are positive and both envelopes
are convex in T>0. Their maximum on any compact interval is therefore
at an endpoint. In particular, checking only

    U/t0<=c0,
    [U-l t0/33*(1-(a/t0)^33)]/a<=c0,
    [U+u t0/33*((b/t0)^33-1)]/b<=c0                  (4)

proves q(T)<=c0 simultaneously for all a<=T<=b. These are an exact
continuous-envelope certificate, not sampling between time nodes.

The scalar checker applies (4) to the existing B015 intervals, with
t0=6/125, a=223/5000, b=251/5000. Its chosen common rational upper bound
and the signs of all checks are in `scalar_certificate.json`.

## 3. Why the recorded negative q derivative is not an all-time proof

Differentiation gives q'(T)=[p_T(x)-q(T)]/T. Thus the existing interval
p_t0<q(t0) only proves a strict local negative derivative. Neither that
fact nor positivity proves a global maximum has already been passed.
The above power envelopes grow without bound toward at least one end
of the unbounded time range, so they do not close T>0. No unimodality
property is assumed or inferred from finite data.

Likewise an exact Green upper bound G(x)<=G_U would prove q(T)<=c for
T>=4096 G_U, but says nothing adequate near T=0. It supplies a tail,
not the missing compact-time and small-time certificate.

## 4. Exact rational density and a sufficient all-T root certificate

The old exact fan gives, without using a spatial derivative formula,

    p_T(x)=(2pi)^-6 R(T^2),
    R(z)=sum_R alpha_R [z^3-e2(b_R)z^2+e4(b_R)z-e6(b_R)]
                              /product_(j=1)^6(z+b_(R,j)^2).       (5)

Here e_k is the kth elementary symmetric polynomial of the six exact
rational nodes b_(R,j)=x.r_(R,j)/q_(R,j). All 967680 terms, including
any repeated or zero nodes, belong to (5). The formula is simply the
real part of product_j(T-i b_j)^-1; repeated factors are retained with
their multiplicity. It is a rational function of z with rational
coefficients, and R(z)>0 for z>0 by the positive density interpretation.

Choose certified rational numbers 0<ell<pi<u_pi. Define

    Pplus(T)=R(T^2)/(64 ell^6)=(pi/ell)^6 p_T(x),
    Fplus(T)=c T-(pi/ell)^6 K(T).

Then Fplus'=c-Pplus is a rational function of T^2. Write R=N/D with
D(z)>0 for z>0, keeping every noncancelled repeated factor. The positive
critical roots are the positive roots of the rational polynomial

    N(z)-64 ell^6 c D(z).                             (6)

A complete exact root-isolation certificate for (6), together with
strict rational brackets [a_j,b_j] for every positive critical time,
has the following finite acceptance rule:

    c a_j-(u_pi/ell)^6 K_upper(b_j)>0  for every bracket.           (7)

Since K is increasing, (7) proves Fplus>0 throughout every root bracket.
Between brackets its derivative has fixed sign, checked by exact
polynomial signs, so its minimum is at an endpoint. Also Fplus(0+)=0
and Fplus(T) tends to infinity as T tends to infinity because G(x)<infinity.
Therefore these finite certificates imply Fplus(T)>=0 for all T>0,
and hence q(T)<=c for all T>0. Multiple roots and stationary inflections
must be included; the proof does not assume a root count, unimodality,
or simple roots. An identically zero polynomial would require separate
handling, and cannot be silently treated as an empty root set.

The tiny dyadic uncertainty in pi can be made negligible by ordinary
outward Arb bounds. The role of ell is to turn the root question into
exact rational arithmetic; it does not justify omitting any cone terms.

This is a rigorous acceptance theorem, **not an executed certificate**.
The existing native record reports 16931 distinct cached node values.
Naive global denominator expansion could therefore be very expensive.
No common denominator, coefficients, root count, root isolation, new G,
or new K value was computed in this round.

## 5. A concrete bounded next-stage proposal (not released or run)

If a separate batch is authorized, use only this same rational x and:

1. One full exact fmpq G pass, preserving principal parts, with a
   1200-second wall cap. Retain G_U; this pays the large-time tail
   T>=4096 G_U without evaluating a time grid.
2. An exact rational-density aggregation pass in the variable z=T^2,
   with no new spatial point. Aggregate all full-group terms, cancel
   rational factors exactly, and stop after 1200 seconds or 1 GiB.
   Do not expand a denominator of degree above 512: if exact compression
   does not bring it below that cap, return INCONCLUSIVE and stop this
   root-certificate route. The old node count makes this a material risk.
3. Only after an exact degree<=512 polynomial is available, isolate all
   positive roots of (6) with a complete rational Sturm certificate,
   capped at 1200 seconds and at most four positive root intervals.
   If either cap fails, stop without any K calls. Tighten the accepted
   root brackets to rational time endpoints within 2^-30 relative width.
4. At most four new complete Arb K passes, one at each retained upper
   root-bracket endpoint, followed by (7). Existing T0 values are reused
   if an endpoint equals T0 exactly. Every principal/zero/repeated-node
   rule stays in force. The first failed sign, unresolved enclosure,
   unsupported input, or budget limit stops with no all-T claim.

This is deliberately bounded and may fail at step 2. It is preferable
to claiming an unproved single peak or extending a finite time scan.
If a lower-degree exact positive decomposition is later found, it can
replace step 2 only after a new explicit registration. No such change
is authorized here. An all-T certificate for body32 would still say
nothing adverse about the independent beta=1 frame.
