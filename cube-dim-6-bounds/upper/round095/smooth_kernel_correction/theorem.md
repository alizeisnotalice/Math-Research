# C095: an explicit small smooth subtraction with whole-space positive source

This is an author proof for the original six-dimensional body32 Cauchy generator.
It awaits independent review. No competitive C6 upper bound is asserted.
The operation of subtracting spatial bumps is an old mechanism (see prior_lookup.md).
The new assertion here is a fully explicit positive admissible amplitude for the
particular off-S3, outside-cube bump below, including all of its nonlocal tails.

## 1. Fixed operator and inherited whole-space inputs

Let V={(1,s2,...,s6):sj=±1}, m=32, Q=[-1,1]^6, and

    A=-psi(D), psi(xi)=sum_(v in V)|v.xi|,
    A f(x)=(1/(2pi)) sum_v integral_R
      [f(x+sv)+f(x-sv)-2f(x)]/s^2 ds.                         (1)

All32 independent scalar Cauchy variables have the same scale t, with density
t/[pi(t^2+z^2)]. The columns are not unit-normalized. The Fourier inverse factor
is (2pi)^-6. Write p_t for their true probability density and G=integral_0^infty p_t dt.
S_k is the union of spans of at most k body directions; no artificial DD stratum
is used. The inherited six-dimensional proof gives:

* p_t is positive/continuous and p_t<=C t^-6, C=1/(160 pi^6), using a body basis
  with determinant160; G is lower semicontinuous everywhere, continuous finite
  off S2, locally L1, and homogeneous of degree-5.
* For every c>0, K0=(G-c)+ is globally L1, and
  A K0=nu0-delta0 on the WHOLE space, nu0 a positive probability measure with
  no atom at0. Equivalently nu0=-A h, h=min(G,c). No singular plane or interface
  was removed in obtaining this identity.
* The obstacle decomposition for this same A holds in dimension6. The original
  finite-source comparison lemma is for ANY K>=0 in L1 with
  AK=a(mu-delta0), a>0 and mu a positive probability measure. It is not restricted
  to K=(G-c)+. Its H1/Fourier proof, one common obstacle set for rational then
  real r>0, L1 dilation continuity and general-input monotone truncation apply
  unchanged to the K constructed below. No new generator is introduced.

These are actual read inputs: six-dimensional theory001 and the original
paper_finite_source.md sections A-C; file hashes are in input_manifest.json.

## 2. A useful explicit uniform tail away from S2

Every three distinct unoriented body directions are independent. Indeed if a
third ±1 vector belongs to span(u,v), coordinates of u/v include both signs;
|a+b|=|a-b|=1 implies ab=0, so the third vector is ±u or ±v. Any independent
triple extends to a six-column body basis, whose nonzero integer determinant
has absolute value at least1.

Fix 0<delta<=1, R>=1 and x with |x|>=R and dist(x,S2)>=delta. Set

    d=delta/(2m sqrt(6)), ell=R/(m sqrt(6)), N=m binom(m-1,2)=14880.

For any coefficient vector z mapping to x, at least three coefficients have
absolute value>d: otherwise the sum of the other at most30 terms lies within
m sqrt(6) d=delta/2 of a two-direction span. A coefficient of maximal magnitude
is at least ell, by |x|<=sqrt(6) sum|z_i|. As ell>=d, select this maximal index
i and two other large indices j,k. Thus the restricted law on the open region
|x|>R, dist(x,S2)>delta is dominated by the sum of N unnormalized events

    |Z_i|>=ell, |Z_j|>d, |Z_k|>d, with j<k and distinct from i.

For each event choose a basis containing i,j,k. In its positive coarea density,
the three restricted Cauchy factors are bounded by
t/(pi ell^2), t/(pi d^2), t/(pi d^2); the other three basis factors by1/(pi t).
The Jacobian is1/|det B|<=1, and remaining26 variables integrate to probability1.
There is no conditioning on X=x, no event renormalization and no cancellation.
Summing gives an a.e. density upper bound; continuity of p_t and then limits in
R,delta give the pointwise bound, including equalities:

    p_t(x) <= A0/(R^2 delta^4), all t>0,
    A0=16 N 32^6 6^3/pi^6.                                    (2)

Use the rational outward constants

    Cbar=1/(160*3^6), Abar=16*14880*32^6*216/3^6,

since pi>3. For every T>0,

    G(x) <= Abar T/(R^2 delta^4) + Cbar/(5T^5).                 (3)

This proves a uniform tail at fixed positive S2-distance, rather than using
homogeneity on a compact angular set that approaches a singular direction.
For c>0 rational, the following explicit integer choices suffice:

    T=max(1,ceil(4 Cbar/(5c))),
    R=max(1,ceil_sqrt(4 Abar T/(c delta^4))), L=12R.             (4)

Here ceil_sqrt(q) means the least positive integer whose square is>=q.
T^5>=T and (3) imply G(x)<=c/2 whenever |x|>=R and dist(x,S2)>=delta.
The intentionally crude constants are finite and directly computable.

## 3. Explicit even rational-box bump and cylinder geometry

Suppose g0>0 is rational and g0<=essinf_Q G. Set c=g0/6,

    p=(6/5,1,4/5,3/5,2/5,1/5), hbox=1/100, delta=1/10,
    eta(u)=exp(-1/(1-u^2)) for |u|<1, eta(u)=0 otherwise,
    b(x)=product_i eta((x_i-p_i)/hbox)
         +product_i eta((x_i+p_i)/hbox).                       (5)

b is even, nonzero, nonnegative and C_c^infinity. Its two closed support boxes
are disjoint and lie outside Q because |x1|>=119/100>1. Every support point has
six distinct absolute coordinates, with consecutive gaps>=1/5-2/100=9/50.
Any vector in S3 has at most four distinct absolute coordinate values: the
eight sign patterns of three directions identify in ± pairs. Therefore any
y in S3 has two equal absolute coordinates, and for x in either support box,

    |x-y| >= (9/50)/sqrt(2) > 1/10=delta.                     (6)

Let C_b=union_(v in V)(supp(b)+R v), a finite union of closed cylinders. For
x=z+sv, z in supp(b), and any body pair plane W,

    dist(x,W)>=dist(z,W+R v)>=dist(z,S3)>=delta.               (7)

Thus the ENTIRE cylinders, not just supp(b), avoid true S2 uniformly. Moreover
Ab=0 off C_b; all nonlocal jump-in tails are retained within C_b.

The cube essential lower bound implies G(y)>=g0 at each off-S2 boundary point:
use actual off-S2 continuity and approach through interior points outside the
null exceptional set. Lower semicontinuity alone would not justify this step.
Homogeneity therefore gives G(x)>=g0/||x||_infinity^5 off S2 outside Q.
On supp(b), ||x||_infinity<=121/100 and (121/100)^5<3, so

    G(x)>g0/3=2c, K0(x)>c, and 0<=b<=1<=K0/c.                (8)

The compact support lies strictly inside Omega={G>c}, an open set by lower
semicontinuity. In particular b vanishes on a neighborhood of the closed set
E=Omega^c. With C0=1/c, b<=C0 K0 globally.

## 4. Genuine source density lower bound where it is needed

h=min(G,c) is bounded continuous: off S2 this uses G continuity; near each S2
point lower semicontinuity of G=+infinity makes h identically c. On Omega it is
locally constant c. Testing nu0=-Ah on compact subsets of Omega, or applying
Tonelli to d=(c-G)+=c-h, gives its full density there:

    m0(x)=(1/pi) sum_v integral_R (c-G(x+sv))+/s^2 ds,
    nu0|Omega=m0(x) dx.                                     (9)

The local small-s integrand is zero, and the tail is bounded by c/s^2. This
also proves local continuity and finiteness of m0 on Omega. Formula(9) includes
the original singular strata; it is not asserted positive everywhere.

By (3)-(7), Omega intersect C_b is contained in {|x|<R,dist(x,S2)>=delta}.
Fix any body direction v0. A pair plane W containing v0 keeps
dist(x+sv0,W)=dist(x,W)>=delta. For a pair plane not containing v0, the Gram
determinant formula and integral columns give

    dist(v0,W)^2=det Gram(v0,w1,w2)/det Gram(w1,w2)>=1/36.

For all x in the above bounded bad region and s in [L,2L], L=12R,

    dist(x+sv0,W)>=s/6-|x|>=R>=delta

for the other planes, and |x+sv0|>=s-|x|>=11R>=R. Hence G(x+sv0)<=c/2
uniformly on this SAME parameter segment. The single v0 term in (9) yields

    m0(x)>= (c/(2pi)) integral_L^(2L) s^-2 ds
           =c/(4pi L)>=c/(16L)=:mbar.                        (10)

This does not assume source positivity on axes/S2. Those sets cannot intersect
the cylinders by(7), and no new correction source is supported there.

## 5. Explicit Ab bound and whole-space measure domination

Writing t=1/(1-u^2), the derivative formulas within(-1,1) are

    eta'=-2u t^2 e^-t,
    eta''=(4u^2 t^4-8u^2 t^3-2t^2)e^-t.

For t>=0, t^k e^-t<=k! follows from the exponential series. Thus
||eta||infinity<=1, ||eta'||infinity<=4, ||eta''||infinity<=148.
The two supports are disjoint, so the six-dimensional Hessian row bound is

    ||D^2b||op,infinity <=(148+5*16)hbox^-2=2280000=:H.

In(1), split at |s|=1. On the inner part bound the second difference by
s^2 |v|^2 H; on the outer part by4||b||infinity. Therefore

    ||Ab||infinity <= (192H+128)/pi
                     <= (192H+128)/3 <145920043=:M.           (11)

The same estimate with L1 norms proves Ab in L1. Translation invariance/Fubini
gives integral Ab=0. Standard dominated convergence in(1) gives Ab continuous.
There are no new singular surface measures, and b is smooth across its support
boundary.

On Omega, (Ab)+ is zero off C_b and at most M on C_b. By(10),

    (Ab)+ dx|Omega <= (16LM/c) nu0|Omega.                     (12)

We separately justify the entire CLOSED complement E, including the level
interface, without assuming a boundary density. Define the positive probability
measures

    nu_t=[P_t K0-K0+U_t]/t dx, U_t=integral_0^t p_s ds.

Their positivity is the original semigroup positive-part inequality; total mass
is1, and test-generator convergence gives nu_t vaguely to nu0. On E, K0=0 a.e.
and b=0, so b<=C0 K0 gives

    nu_t(x)>=P_tb(x)/(C0 t), for a.e. x in E.

For any nonnegative compactly supported continuous phi, phi*1_E is upper
semicontinuous. The CLOSED-set Portmanteau direction is

    integral_E phi dnu0 >= limsup_(t down0) integral_E phi dnu_t
        >= (1/C0) integral_E phi Ab dx.                      (13)

The last limit uses uniform (P_tb-b)/t -> Ab and b=0 on E; b vanishes on a
neighborhood there, so Ab>=0 on E. The use of limsup with this direction is
essential. Vague convergence suffices for the compactly supported upper
semicontinuous test; alternatively approximate it from above by continuous
compactly supported functions. This proves a measure inequality on E, retaining
any old singular/interface source. Combining(12)-(13) gives the explicit global
bound

    (Ab)+ dx <= C1 nu0, C1=16LM/c >= C0.                     (14)

No positivity sampling or unidentified compact minimum appears in C1.

## 6. Actual admissible correction and transferred mass gain

Choose the rational positive amplitude

    epsilon=c/[32L(M+1)].                                    (15)

Then epsilon C0<1/2 and epsilon C1=M/[2(M+1)]<1/2. Consequently

    K_epsilon=K0-epsilon b >= K0/2 >=0,
    nu_epsilon=nu0-epsilon(Ab)dx >=nu0/2 >=0,
    A K_epsilon=nu_epsilon-delta0,
    nu_epsilon(R6)=1, integral K_epsilon=integral K0-epsilon integral b.

K_epsilon is even and L1. The original negative atom stays exactly-delta0;
the original positive singular part is unchanged, since only the smooth density
Ab dx is subtracted. No interface is discarded, and no new one is created.
On Q, b=0, hence cube domination is exactly preserved.

On the inner half of each bump box, eta>=e^-2>1/9 (e<3), so

    integral b >= 2 hbox^6/9^6 =2/(100^6*9^6)=:Bmin>0.        (16)

Normalize with N0=64(g0-c)=320c. Then K_epsilon/N0>=1_Q/64 and the general
finite-source comparison gives its valid all-L1/all-positive-common-scale
weak bound equal to its mass. Relative to the same baseline normalized R_c,
the exact mass decrement is epsilon integral b/N0, at least

    Delta_min=Bmin/[10240L(M+1)] >0.                          (17)

Using I=integral_[0,1]^5 G(1,y)^(6/5)dy, the baseline mass is
6^(6/5) I/g0^(6/5), so the correction yields that expression minus the positive
decrement. A strict improvement of this auxiliary bound is not a proof that it
is below64. For a certified I_U the cost-slack change is

    S_new >= S_old+(g0/6)^(6/5) Delta_min,
    S_old=64(g0/6)^(6/5)-I_U.

This tiny explicit construction does not establish S_new>0.

## 7. A fully rational noncompetitive instantiation, without a body evaluation

For existence without an unpaid g0 premise one may take

    gstar=1/[160 * (4*(4+122^2))^6 * 10^26].                  (18)

Proof: use the already accepted determinant160 body basis. Each inverse entry
has absolute value<=5!/160=3/4 by the cofactor formula. For x in Q, restrict the
remaining26 Cauchy parameters to[-1,1]; then basis coordinates have absolute
value<=6*(3/4)*27=243/2<122. For1<=t<=2, each basis density is at least
1/[4*(4+122^2)], using pi<4. Each remaining parameter has probability at least
2/[4*(4+1)]=1/10. The positive coarea density, with Jacobian1/160, is therefore
at least gstar. Integration over[1,2] gives G>=gstar everywhere on Q.

Equations(4)-(18) with g0=gstar specify every constant rationally/with integer
ceilings. The finite exact checker only checks these algebraic bounds and bump
derivatives; it does not evaluate a body kernel or source integral. This weak
gstar is an existence device, not a competitive global floor proposal. The
general construction can use any later certified same-body32 g0, with all
constants recomputed from that exact input under a new finite registration.

## 8. Scope of proof and remaining work

The uniform line tail, needed source density lower bound, exterior measure
domination and finite-source transfer have been given above, not left as
unverified domination assumptions. The author regards this as a complete proof
of a positive but very small admissible subtraction. Independent mathematical
review remains necessary. The proof does not pay a competitive I_U or the final
S>0 condition, optimize the bump/amplitude, combine different generators, prove
global positivity from samples, or certify new C6 values.
