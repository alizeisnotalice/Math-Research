# 104 prepared proof: restore the retained parameter after two deletions

Status: analytic submission and unexecuted finite script. No104 mathematical
program, LP, kernel, peak, or grid has been run. The104 script requires a
separately bound parent release. Acceptance of the098 time-first identity
and review of the new source argument remain explicit dependencies.

The target and model remain unchanged: Lebesgue weak(1,1), dimension6,
centered axis-parallel cubes, a common radius, general L1, every positive
radius; original32 unoriented body columns, weights1, alpha1, common
scalar Cauchy time t, Fourier normalization(2pi)^-6. The actual upper64
and strict lower2073015624830/10^12 are unchanged. No lower-bound work
is proposed or assigned.

## 1. Keep the actual exhaustive source split

For nonempty proper A containing1 put k=|A|, B=A complement, m=6-k,
v=ones, w=1_A-1_B, U_A=span(v,w). Use exactly043/049/051:

 E={|a_v|>=1/2, |a_w|>=1/16},
 E_v={|a_v|<1/2},
 E_w={|a_v|>=1/2, |a_w|<1/16}.

Let P_t be the density on E and B_t the actual disjoint complement.
For an upper bound B_t is at most the sum of the two enlarged densities
f_(b,t) on |a_b|<theta_b, with(b,theta_b)=(v,1/2),(w,1/16).
No disjointness or reflection symmetry is assigned to these enlarged
events. All other parameters remain unrestricted.

Fix one such b. In the31-column source with b omitted, perform only the
first two positive deletions from043. At stage0/1 the counts are31/30,
the frame operator is32I-(P b)(P b)^T<=32I, and the denominators are

 |x-a b|²+992t²,
 dist(x-a b,span i)²+960t².

Let h_b0(x),h_b1(x)>0 be valid lower bounds for these distances for every
|a|<theta_b and every relevant body line i!=b. This is precisely the
first-two-gap interface already proved in043/051, or its100 u/rho
refinement. Each deleted ordered pair(i,j) has i,j!=b, i!=j, with true
coarea J_ij=sqrt(det Gram(i,j)). The projected remaining density is the
actual four-dimensional q^omit_b_(ij),t.

Pull just these two denominators outside the restricted a integral. Now
positivity permits extending that integral to all real a. The resulting
identity is the ordinary source convolution

 int_R c_t(a) q^omit_b_(ij),t(P_(ij)perp x-a P_(ij)perp b) da
                  =q_full_(ij),t(P_(ij)perp x).          (1)

This restores the original body32 four-dimensional quotient **before**
taking a peak or integrating time. No third distance has been introduced.
Its possible zeros along the old a fiber therefore cause no obstruction.

The resulting new positive bound is

 f_(b,t)(x) <= [1024t²/pi²]/
                [(h_b0²+992t²)(h_b1²+960t²)]
        * sum_(ordered i!=j, i,j!=b)
                q_full_(ij),t(P_(ij)perp x)/J_ij.         (2)

All930 ordered pairs are present. There is no division by a pair count,
probability-independence assumption, or unique source label assumption.
All configurations with three or more cancelling normal sources remain
inside each genuine restored quotient.

## 2. Pointwise, coarea, and every-pair domain details

Every three distinct body lines are independent, so the first two
Gram-Schmidt coarea factors are nonzero and multiply to1/J_ij. Before
restoration, the actual projected omitted-column frame has Gram
32I-(P b)(P b)^T>=26I. Its continuous densities make the positive
test-function deletion inequalities pointwise, as in043/051. The
restricted branch density is also continuous: it is an integral of a
continuous bounded31-column full density against a restricted scalar
Cauchy subprobability. Finally(1) holds pointwise by positive convolution
and its continuous Fourier density. Thus(2) is not merely an almost-
everywhere claim inferred by conditioning on a zero-probability fiber.

On the051 full middle tube

 D_A={x=(1,y),y in[0,1]^5: r=dist(x,U_A)<1/64,
                                      u=mean_B x in[1/4,3/4]},

every other actual body pair plane satisfies

 dist(x,span(i,j))>=g_k=1/(4sqrt2)-sqrt(k)/64>1/64.       (3)

The only body lines in U_A are v and w. Any pair avoiding b=v, or
avoiding b=w, cannot span U_A. Consequently **each of the930 pairs in
each branch** satisfies(3), including all mixed signed pairs and all
pair orbits. No pair with rho_ij=0 is hidden inside a finite budget.

Outside a qualified domain, formula(2) is used only where its two first
gaps are valid; a time-integrated expression with rho_ij=0 can be infinity
and is not claimed to pay a finite spatial budget. The present104 local
and full-middle integration domains all satisfy the strict pair gaps(3).

## 3. Time first: an actual weighted line trace, then a coarse paid version

Write z_ij=P_(ij)perp x=r_ij e_ij. The098 identity applies to the actual
restored quotient:

 int_0^infinity t²q_full_(ij),t(z_ij)dt
                         =Q3_ij(e_ij)/(2r_ij), r_ij>0.  (4)

Q3_ij(e) is the time-one peak of the actual full-frame projection onto
span(i,j,e)^perp. It is not an omitted-column peak, not necessarily an033
body-triple orbit peak, and not a pointwise isotropic heat-kernel bound.
The proof of(4) is scaling t^-4 and the continuous true line marginal,
not an oscillatory Fourier interchange.

Keeping the two time denominators from(2) gives the stronger positive
expression

 int f_(b,t)(x)dt <=(1024/pi²)sum_(i,j) 1/(J_ij r_ij)
  *int_0^infinity q_full_(ij),1(s e_ij)
    /[(h_b0²+992r_ij²/s²)(h_b1²+960r_ij²/s²)] ds.         (5)

The prepared finite batch uses the simpler consequence

 int f_(b,t)(x)dt
    <=[512/(pi² h_b0²h_b1²)]sum_(i,j)
                           Q3_ij(e_ij)/(J_ij r_ij).      (6)

The full-body coercivity already proved in098, c=32/sqrt3, implies
Q3_ij(e)<=1/(pi²c³). This is the same full-frame bound, not100's worse
omitted-frame coercivity replacement. The scientific change occurs before
this estimate: two deletions, restoration, then actual time integration
eliminate the entire26970-triple Crem sum, the1000/729 omission loss,
and the old Cauchy-average factor.

For one fixed retained body line the pair orbit counts are

      (180,450,300),  J_ij²=(20,32,36).

Indeed the full ordered counts(192,480,320) lose twice the fixed-line
partner counts(6,15,10). Their sum is930. Define

 Sigma=180/sqrt20+450/sqrt32+300/6,
 d=512 Sigma/(pi^4 c³), D=sqrt2 d=3sqrt6 Sigma/(64pi^4).

Then the complete complement obeys

 int B_t(x)dt <=d[1/(h_v0²h_v1²)+1/(h_w0²h_w1²)]/g(x),  (7)

whenever g(x) is a simultaneous positive lower bound on all relevant
r_ij. This is a lawful new replacement for the old Crem-based full-time
complement budget. It is an upper estimate, not a claim about the actual
loss or sharpness of930 positive terms.

## 4. Prepared joint spatial core budget at original radius

First use the same100 region R with k=m=3, u in[3/8,5/8], r<=delta=1/64.
Its actual center gaps from100 are

 h_v0>=sqrt3/2, h_v1>=sqrt(3/2)mu,
 h_w0>=sqrt822/16, h_w1>=sqrt(3/2)(7/8-u),
 g>=(1-u)/sqrt2, mu=min(u,1-u).

After subtracting sqrt3 r, all five gaps are at least(1-7r) times these
center values. Their total degree in(7) is5. Consequently

 int B_tdt <=D(1-7r)^-5 Z_core(u),
 Z_core=8/[9mu²(1-u)]+256/[1233(7/8-u)²(1-u)].           (8)

Both positive summands are reciprocal products of positive affine
functions on each half interval. Their logarithmic second derivatives
are nonnegative, hence they are convex there. Let M_core be the maximum
at u=3/8,1/2,5/8, and J_core the exact integral over that interval.
It is evaluated by rational/log partial fractions, with no quadrature.

Reuse the exact100 radial identity
J7=int_0^delta r³(1-7r)^-6 dr. The face Jacobian3 and maximal-A-coordinate
cone fraction1/3 cancel. The complete new core budget prepared for testing is

 B_core=D[4v4 J7 M_core^(1/5)J_core]^(5/6), v4=pi²/2,
 P_core=(sqrt3/(128pi^4))[(5v4/14)delta^(14/5)]^(5/6),
 N_core=P_core+B_core, I_core=N_core^(6/5).              (9)

The mixed terms of the two enlarged small events are retained inside
Z_core before raising to q. The only spatial Minkowski step adds the
already full-time large branch and complement. The ten symmetric core
regions are disjoint by(3). The script will compare10I_core with the
saved100 same-region allocation, and compute

 S_core_plan=64(g0/6)^(6/5)-I_paid_CEH-10I_core,
 g0=1/65000, still unproved globally.                    (10)

The remaining middle strips, other pair types, far region, and global g
remain unpaid in(10). A positive value there is remaining allocation,
not a complete positive S or C6 breakthrough.

## 5. A conservative complete31-middle budget, including all u endpoints

The prepared second spatial test covers all051 middle domains, u in
[1/4,3/4], and all0<r<=1/64. For every k,m use center lower gaps

 h_v0>=sqrt(k)/2, h_v1>=sqrt(km/6)mu,
 h_w0>=15sqrt(k)/16, h_w1>=sqrt(km/6)(7/8-u),
 g>=mu/sqrt2.

The actual perturbation is at most sqrt(k)r. Relative losses are at most
2r,4sqrt(6/m)r,(16/15)r,8sqrt(6/m)r, and4sqrt(2k)r, respectively.
All are bounded by20r for k,m in1,...,5. Hence(7) gives

 int B_tdt <=[D/(k²m)](1-20r)^-5 Z_all(u),
 Z_all=24/mu³+512/[75mu(7/8-u)²].                       (11)

This is intentionally conservative. It is not assumed competitive before
execution. M_all is the convex-piece endpoint maximum at1/4,1/2,3/4,
J_all is the exact rational/log u integral, and J20 is the new exact
radial integral int_0^delta r³(1-20r)^-6 dr. It is finite because
20delta=5/16<1. For either lambda7 or20,

 Jlambda=[F(1)-F(1-lambda delta)]/lambda^4,
 F(w)=-1/(5w^5)+3/(4w^4)-1/w³+1/(2w²).

The true face Jacobian times cone fraction is sqrt(m/k). Thus

 B_k=[D/(k²m)]
       [4v4 sqrt(m/k)J20 M_all^(1/5)J_all]^(5/6),
 P_k=[256/(sqrt(km)pi^4c³)]
       [(5v4/7)sqrt(m/k)delta^(14/5)]^(5/6),
 I_middle,U=sum_(k=1)^5 binom(5,k-1)(P_k+B_k)^(6/5).     (12)

Positivity pays the union whether or not overlaps are invoked. No extra
051 short term is added: (12) covers the full G at all positive times.
The script compares the identical full-middle allocation with only the
saved051 short upper plus074 long upper, without evaluating their old
functions, source programs, or time table. It then records

 I_plan=I_paid_CEH+I_middle,U,
 S_plan=64(g0/6)^(6/5)-I_plan.                            (13)

The far allocation is set to zero only in this explicitly optimistic
certificate-expression diagnostic. Its actual payment and global g are
still necessary. Either sign of(13) is retained without a rerun or scan.
An oversized constructed upper allocation is never a lower bound on true I.

## 6. A strict actual-source counterexample to a tempting pair-only cover

A cancellation in four-dimensional normal space need not be a pair of
antiparallel large projections. At the k3 middle x=(1,1,1,1/2,1/2,1/2),
take a_v=0, a_w=5/8 and the following three actual body representatives:

 d1=(1,-1,-1,-1,-1,-1), a_d1=-3/8,
 d2=(1,-1,1,1,1,1),   a_d2=3/8,
 d3=(1,1,-1,1,1,1),   a_d3=3/8.

All other parameters are zero. The resulting physical sum is exactly x,
and this lies in the actual v-small complementary event. Its three
weighted nonzero normal projections, each above the fixed parameter
threshold1/4 in absolute value, sum to zero, have equal lengths,
and mutual cosine-1/2. No pair is parallel or cancels. Their source
densities are positive; nearby parameters preserve three-large labels
and the absence of any antiparallel large pair. Adjusting a_w to keep
x1=1 and choosing the first signed coefficient no larger than the other
two keeps x2,x3<=1 and yields nearby nonzero normal distances in the tube.

Thus a cover whose only active large-label edges are antiparallel pairs
fails on a genuine source pattern; it cannot be a pointwise fractional
cover. A unique-largest-label partition itself is valid with tie breaking,
but its conditional source constraints cannot be discarded or replaced by
independence. The restored quotient in(1)-(6) automatically includes this
triangle cancellation and all larger networks. The fourth prepared exact
configuration checks this witness; it does not fit an invalid cover by
assuming away overlaps.

## 7. Run contract: source only until parent release

The four configurations are fixed in registration.json. The only script,
run_prepared.py, is single process. Its first executable operations set
SIGALRM to SIG_DFL and arm177 seconds. The timer is never cancelled,
including during final result/receipt writing, print, and interpreter exit.
RLIMIT_CPU is60 seconds for this process. No child exists to clean up.
Checkpoints preserve completed configurations; the default alarm can kill
the process even inside a finally block. A timeout may therefore leave
no final receipt, which must not be described as successful completion.

No parent_release.json is supplied. Before any mathematical imports or
configuration it requires an explicit release bound to the exact script
and registration hashes. Existing run artifacts refuse a second attempt.
After a mathematical failure it writes a first-failure receipt and exits;
there is no retry, parameter expansion, or fallback scientific batch.

The177-second fuse starts at its top-level arming point after importing
signal; interpreter initialization before that point is not independently
timed by this script. The intended headroom to180 seconds is not a proof
about launch startup. A parent launcher may add an independent180-second
startup-to-exit cap if that complete outer timing contract is required.

The script and all104 mathematical configurations remain UNEXECUTED.
