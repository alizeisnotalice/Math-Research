# 100: a full-radius time-first budget and a precise remaining bottleneck

The normalized model remains original body32, unit weights, alpha1,
psi(xi)=sum_(v=(1,+/-1,...,+/-1))|v.xi|, scalar Cauchy common time t,
Fourier factor(2pi)^-6, and Lebesgue measure. Put q=6/5 and
G=int_0^infinity p_tdt. No generator, time convention, or spectral frame
is changed. New proof statements are author submissions pending independent
review, including acceptance of the098 time-first interface.

## 1. Deduplication and the reused three-resolvent formula

The full formula below, including its coincident-root extension, was already
proved and actually used in062/theorem.md section5, evaluate.py line92,
and check_saved.py line94. It is not a new100 discovery. It was not used
to integrate the043/074 middle complement. No old source program was run.

For a,b,c>0 let

 J(a,b,c)=int_0^infinity dt/[(t²+a²)(t²+b²)(t²+c²)]
        =pi(a+b+c)/[2abc(a+b)(a+c)(b+c)].                 (1)

For completeness, the two-factor integral is pi/[2ab(a+b)]: subtract
1/(t²+a²) and1/(t²+b²), divide by b²-a², and use the elementary
arctangent integral pi/(2a). Subtracting the two-factor formulas for
(a,b) and(a,c), then dividing by c²-b², gives(1) when the roots are
distinct. On any compact positive parameter set, the integrand is bounded
by(t²+m²)^-3 for a fixed m>0. Dominated convergence extends the integral
to all equalities among positive roots; the right side is already
continuous there. At a=b=c it is3pi/(16a^5). This provides a self-contained
reuse proof without a new executable configuration for the old identity.

Keep the043 k3 pair
v=(1,1,1,1,1,1), w=(1,1,1,-1,-1,-1), U=span(v,w),
and its actual split p_t=P_t+B_t, with P on
{|a_v|>=1/2, |a_w|>=1/16} and B the actual disjoint complement. The old
043/074 positive bound on B is

 B_t <= Krem/(g²+928t²)
   [1/((hv0²+992t²)(hv1²+960t²))
    +1/((hw0²+992t²)(hw1²+960t²))],                      (2)

where, for the outer tube rho<1/64 and mean_B u in[1/4,3/4],

 eta=sqrt3/64,
 hv0=sqrt3/2-eta, hv1=sqrt(3/32)-eta,
 hw0=15sqrt3/16-eta, hw1=sqrt(3/128)-eta,
 g=1/(4sqrt2)-eta,
 Crem=44743079/(2970pi^6),
 Krem=(1000/729)Crem(1+sqrt(6/928)).

The missing-column peak loss and complete ordered-triple sum have not been
removed. Apply(1) directly with normalized roots h/sqrt(A) and divide
by928*960*992. This newly pays the full middle complement by

 B_middle=Krem/(928*960*992)
   [J(g/sqrt928,hv0/sqrt992,hv1/sqrt960)
    +J(g/sqrt928,hw0/sqrt992,hw1/sqrt960)].               (3)

Its outward value is about179.519136791234, versus098's saved coarse
B0 about343.378053668633, ratio about0.522803175314. The new payment
keeps all three resolvents; no074 endpoint function or time grid is replayed.
It is still a uniform worst-middle bound and remains expensive.

## 2. A fixed joint mean/normal region at the original radius

Freeze

 R={x=(1,y), y in[0,1]^5: 3/8<=u=mean_(4,5,6)x<=5/8,
                                    0<r=dist(x,U)<=1/64}.

This is not098's exponentially thin region. The radius is the original
1/64 and the u interval has length1/4, half the original middle interval.
Write x(u)=(1,1,1,u,u,u). The true034/043 face coordinates give
|x-x(u)|<=sqrt3 r. At the center, retain the restricted-parameter geometry
instead of using every worst endpoint at once:

  h_v0(u,0)>=sqrt3/2,
  h_v1(u,0)>=sqrt(3/2)min(u,1-u),
  h_w0(u,0)>=sqrt(3[(15/16)²+(u+1/16)²]),
  h_w1(u,0)>=sqrt(3/2)(7/8-u),
  g(u,0)>=(1-u)/sqrt2.                                  (4)

For the v-small event, write the center residual as(p1_A,q1_B),
p=1-a_v, q=u-a_v. Here p>=1/2. Its nearest body-line distance is
sqrt(3/2)||p|-|q||, whose infimum over a_v<1/2 is
sqrt(3/2)min(u,1-u). For the w-small event p=1-a_w,q=u+a_w are positive;
the squared norm decreases up to a_w=1/16, giving the h_w0 bound, while
p-q>=7/8-u gives h_w1. The other-plane gap is the043 exact center
minimum min((1-u)/sqrt2,sqrt2 u); u>=3/8>1/3 selects the first term.
All lines/pairs used in the complementary branch are included as in043.

Distances are1-Lipschitz, so subtract sqrt3 r from each bound in(4).
Uniformly on3/8<=u<=5/8, each center gap h is at least sqrt3/7. For
example the smallest g is(3/8)/sqrt2, so sqrt3/g<=8sqrt6/3<7;
the smallest w-line gap is sqrt(3/2)/4, giving4sqrt2<7. The other three
inequalities are stronger. Therefore every needed gap satisfies

       h(u,r)>=(1-7r)h(u,0)>0.                           (5)

The use of a common factor is a proved lower bound on the actual gaps,
not a spatially changed kernel or probability time.

## 3. Positive rational time integrals retain the u dependence

For the three normalized roots define the following rational lower bounds:

 a(u)=7(1-u)/310,
 bv(u)=6min(u,1-u)/155,
 bw(u)=6(7/8-u)/155,
 cv=3/112, cw=7/128.                                    (6)

Indeed1/sqrt2>7/10, sqrt(3/2)>6/5,
sqrt928<31, sqrt960<31, sqrt992<32, sqrt3/2>6/7,
and sqrt(3[(15/16)²+(u+1/16)²])>=sqrt822/16>7/4.
All these are elementary square comparisons between positive quantities.
Thus(5)-(6) bound each physical normalized root below by(1-7r) times
its displayed rational counterpart.

Put

 T(a,b,c)=(a+b+c)/[abc(a+b)(a+c)(b+c)],
 Z(u)=T(a(u),bv(u),cv)+T(a(u),bw(u),cw),
 C=pi Krem/[2*928*960*992].

The integral J is decreasing separately in its positive roots and is
homogeneous of degree-5. Integrate the real B_t first, using(1), then
apply the root lower bounds. This gives on R

       int_0^infinity B_t(x)dt <= C(1-7r)^-5 Z(u).        (7)

No per-time Lq estimate, instantaneous mass relaxation, shifted-residual
uniform third gap, or source-event reflection symmetry is used. The
actual split remains exhaustive. Both small-parameter branches and all
other source parameters are paid, even though their enlarged upper events
can overlap.

Z is rational on each of[3/8,1/2] and[1/2,5/8]. Its mixed terms are
retained before raising to q. A useful exact convexity certificate is

 T=1/[bc(a+b)(a+c)(b+c)]
   +1/[ac(a+b)(a+c)(b+c)]
   +1/[ab(a+b)(a+c)(b+c)].                              (8)

Each denominator is a product of positive affine functions of u on the
given piece. A reciprocal product f has
f''=f[(sum l_i'/l_i)^2+sum(l_i'/l_i)^2]>=0. Consequently each piece
of Z is convex. Its maximum is at the two piece endpoints, so

       M=max{Z(3/8),Z(1/2),Z(5/8)}                       (9)

is a rigorous exact rational upper bound. These three fixed analytic
endpoints are not a quadrature grid or an adaptive scan.

The exact positive mean integral

       J_u=int_(3/8)^(5/8) Z(u)du                       (10)

is stored as a rational part plus rational coefficients times logarithms
of positive rational numbers. The minimal code uses exact partial
fractions of the three rational pieces. Each term n/(a u+b)^k is
integrated by its elementary rational primitive for k>1 or
(n/a)log|(a r+b)/(a l+b)| for k=1. No numerical quadrature is used.
All logs and signed sums are enclosed by Arb384. The value is about
83399927.400866. Since0<=Z<=M, Z^q<=M^(1/5)Z pays the complete
combined branch q power without dropping its mixed terms.

## 4. Analytic radial integration and the full local G budget

The original face Jacobian is3 and the maximal-A-coordinate cone has
radial fraction1/3, as explicitly derived in098. They cancel for every
radial integrand that depends on r and u. Other face clipping shrinks
the positive integrals. Let delta=1/64 and v4=pi²/2. Then

 |R|<=v4 delta^4/4,
 int_R r^-q dy <=(5v4/14)delta^(14/5).                   (11)

From098's true normal-marginal identity and the actual two-large-parameter
coarea, the full-time large branch is bounded by K/r with
K=sqrt3/(128pi^4). That identity holds in every normal direction; it does
not assert an isotropic pointwise density bound.

The radial factor in the q power of(7) is now exactly integrable:

 J_r=int_0^delta r³(1-7r)^-6 dr
     =[F(1)-F(1-7delta)]/7^4,
 F(w)=-1/(5w^5)+3/(4w^4)-1/w³+1/(2w²).                 (12)

This follows from w=1-7r; the derivative is(1-w)^3/w^6. The stored
J_r is an exact positive rational. Combining(7)-(12), and only now using
spatial Minkowski, proves

 P_R=K[(5v4/14)delta^(14/5)]^(5/6),
 B_R=C[4v4 J_r M^(1/5) J_u]^(5/6),
 ||G||_Lq(R)<=N_R:=P_R+B_R,
 int_R G^q <=I_R,U:=N_R^(6/5).                          (13)

All t>0 and all source events are included. Zero normal distance is a
null spatial set and can retain an infinite pointwise majorant.

The outward numerical values are

 P_R about1.359641231748e-8,
 B_R about2.796713872380e-5,
 N_R about2.798073513612e-5,
 I_R,U about3.437406901520e-6.

The corresponding volume upper is about7.353428305007e-8. On exactly
this same new R, transfer of the old098 uniform formulas gives norm
about0.000581685625 using512 or about0.000390117911 using its saved
B0 upper endpoint. The new full-radius joint norm is about0.07172379
times the latter. Merely inserting(3) uniformly on R gives the larger
norm about0.000203961371. These old expressions on R are newly compared
analytic allocations, not historically run R certificates.

## 5. Ten disjoint k3 central regions and the conditional S diagnostic

There are binom(5,2)=10 choices of A containing coordinate1 with |A|=3.
Coordinate permutations fixing1 preserve the model, face, and estimates.
For a point in any one such R_A, its distance to every other body pair
plane is at least

       (3/8)/sqrt2-sqrt3/64 >1/64.

Hence these ten R_A are pairwise disjoint. This is a geometric statement
for these central regions, not a disjointness assumption for all31 tubes.
Their constructed total integral budget is

       I_ten,U=10 I_R,U about3.437406901520e-5.           (14)

With the unproved candidate g0=1/65000, the whole-face threshold is
B_gate=4/48750 and I_target=64(g0/6)^(6/5), about1.249978219301e-5.
Even setting all other spatial allocations to zero, the present constructed
budget has

 S_plan=I_target-I_ten,U about-2.187428682219e-5,
 U_plan=I_ten,U/(g0/6)^(6/5) about175.9983000507.          (15)

This rejects competition of this particular upper-certificate expression.
I_ten,U is an upper allocation, not a lower bound on true I. Equations
(14)-(15) neither exclude the true model nor lower-bound actual C6.

B_R is proportional to the complementary coefficient Crem, whereas P_R
is independent of it. If a lawful future estimate reduced the whole
complementary coefficient by a common factor gamma, even these ten regions
alone would require

 gamma < [B_gate/10^(5/6)-P_R]/B_R about0.430143565158.    (16)

The other face regions and a matching global g_L would impose further
requirements. The complement accounts for about99.9514% of the current
norm budget, so within this expression it is the principal loss. This
does not prove that the true complementary source potential is large.

For context,098's thin-Q true integral upper endpoint is at most about
3.016816337966e-21 of I_target. Even eliminating its entire actual
contribution could have only that relative impact. This round instead
pays the full original normal radius on a fixed central u interval, but
it still does not yield a competitive complete face budget.

## 6. A lawful omitted-column coercivity replacement is worse, not better

The coefficient Crem comes from26970 ordered triples of non-retained
directions after three positive deletions. The original exact033 peaks
are used after the correct coarea1/J3. In the middle bound, the retained
small column b is still omitted from that quotient when its peak is
taken; psi-|b.xi|>=(9/10)psi pays the factor1000/729. The subsequent
true Cauchy average pays shifted third-distance zeros and contributes
1+sqrt(6/928). Restoring b before all distances are pulled out, as in062,
cannot be assumed here because a shifted third distance can vanish.

It is legitimate to derive a different coercivity for the actual omitted
law. Set n=|xi|, s=|b.xi|, and x=s²/n² in[0,6]. The true31-form sums
after omitting b satisfy

 sum_(v!=b)(v.xi)²=32n²-s²,
 sum_(v!=b)(v.xi)^4<=96n^4-s^4.

The finite-sum Holder inequality used in098 gives, for n>0,

 psi_omit(xi)>=n sqrt((32-x)^3/(96-x²)) >=c31 n.           (17)

Here c31 is obtained at x*=-32+4sqrt82 in(0,6). The derivative of the
squared expression has the sign of x²+64x-288, so this is its minimum.
The resulting c31 is about16.558397480103. The zero case is immediate.
This proof concerns the actual deleted body32 law; restriction to any
three-dimensional quotient preserves it. No peak of another model is used.

Therefore a lawful universal omitted-quotient peak bound is
1/(pi²c31³). Reuse the exact five triple classes from033/043, with
unordered counts N=(435,1740,435,580,1305) and Gram determinants
D=(64,96,128,160,192). The complete replacement deletion coefficient is

 C_moment=(196608/(pi^5 c31³)) sum_j N_j/sqrt(D_j).        (18)

The factor6 converts to26970 ordered triples. The saved exact R_j also
verify Crem*pi^6=24576 sum N_j R_j=44743079/2970. No tuple enumeration
is repeated. The new lawful coefficient(18) is about58.083412007944,
whereas the existing exact-peak coefficient after omission loss is
(1000/729)Crem about21.495257546860. The ratio is about2.702150085028.
The common Cauchy-average factor does not improve this comparison.

Thus the new uniform31-form moment coercivity cannot simply replace
the old exact projected peaks to improve the budget. This is a comparison
of these two valid bounds, not a claim that all sharper quotient coercivity
or direct joint-source estimates are impossible. A useful next mechanism
must recover structure lost by the positive ordered-triple/global-peak
sum, or retain the omitted parameter and its third denominator jointly.

## 7. Global scope remains unpaid

The target remains the fixed-Lebesgue weak(1,1) constant for centered,
axis-parallel cubes, all coordinates sharing a radius, general L1 inputs,
and all positive radii. The final criterion remains

       S=64(g_L/6)^(6/5)-I_U>0.

No global g_L, complete competitive I_U, or positive S is supplied. The
actual upper64 and strict lower2073015624830/10^12 are unchanged. No
model-cost certificate failure here is an actual C6 lower bound. The four
finite arithmetic configurations and all source bindings are adjacent;
their successful execution does not replace independent proof review.
