# Two specified parameters: an integrated third-distance bound

Use the original body32 Cauchy frame, all parameters at the same semigroup
time s>0. Write v=(1,1,1,1,1,1), w=(1,1,1,-1,-1,-1), U0=span(v,w),
A={1,2,3}, B={4,5,6}, and x(u)=(1,1,1,u,u,u). Fix only

    D={x=(1,y), y in[0,1]^5: dist(x,U0)<rho,
                              mean_B(x)=u in[1/4,3/4]},
    rho=1/64, epsilon=1/1000.

Distances and quotient densities use Euclidean orthogonal measure;
the face norm uses dy. This proves a fully paid short-time L6/5 bound on
D. It is not a Green asymptotic or a complete I, g or C6 certificate.

## 1. Two large parameters and an actual disjoint complement

Partition the original parameter probability space into

    E={|a_v|>=1/2, |a_w|>=1/16},
    E_v={|a_v|<1/2},
    E_w={|a_v|>=1/2, |a_w|<1/16}.

Let P_s,B_v,s,B_w,s be their pushforward densities and B_s=B_v,s+B_w,s.
Then p_s=P_s+B_s, all terms are nonnegative, and both P_s and B_s are
pointwise at most p_s. To bound B_w,s we may enlarge E_w to the event
|a_w|<1/16, but the original B_s remains the disjoint complement. This
distinction permits one L1 bound for the entire B_s, rather than paying
two whole-density L1 bounds. No other parameters are assumed small.

On E the two densities are at most

    s^2/[pi^2(1/4+s^2)(1/256+s^2)].

Push the resulting positive product measure forward before evaluating
its density. Integration along the two removed parameters has coarea
Jacobian1/J2, where J2=sqrt(det Gram(v,w))=6. The remaining quotient is
exactly the original Q4 of pair type k=3. Hence

    P_s(x)<= R2/(16pi^6) * s^-2/
                         [(1/4+s^2)(1/256+s^2)]
            <=C2 s^-2,
    R2=2051/47038464, C2=64R2/pi^6.                       (1)

For pointwise rigor, the remaining30-column full frame has Gram
32I-vv^T-ww^T>=26I, so its density and the truncated convolution P_s
are continuous. The removed-parameter pushforward density is the
continuous quotient density divided by J2: choose four remaining
projected columns forming a basis of U0perp and use their linear change
of variables together with the two free parameters. Tonelli proves the
measure identity, and its density inequality extends everywhere by
continuity. No zero-probability fiber conditioning is used.

## 2. A missing-one quotient peak needs a new comparison

The original symbol psi is convex and invariant under coordinate signs
and permutations. For nonnegative coordinates with sum S, averaging all
coordinate permutations and using convexity gives

    psi(xi)>=psi((S/6)v)=10S,

because psi(v)=60. Thus for every body direction b,

    |b.xi|<=||xi||1<=psi(xi)/10,
    psi(xi)-|b.xi| >=(9/10)psi(xi).                        (2)

This holds after restriction to any orthogonal quotient W. Fourier
inversion and the triangle inequality therefore bound the density of
the quotient with parameter b omitted by

    sup_z q^(omit b)_(W,s)(z)
       <=(10/9)^dim(W) Q_W s^-dim(W).                    (3)

No new peak is assumed or numerically evaluated. In particular, at the
third deletion used below the loss is(10/9)^3. If the omitted direction
already projects to zero, (3) is still valid, though not sharp. This is
a comparison of densities at the actual common time s, not a change to
the parameter law in any event decomposition.

## 3. Geometry: two positive denominators, not three

For x in D the034 coordinates give |x-x(u)|<=eta=sqrt(3)rho. Indeed, the
A and B zero-sum residuals e,z satisfy dist(x,U0)^2=|e|^2+|z|^2, while
|x-x(u)|^2=dist(x,U0)^2+3e1^2<=3dist(x,U0)^2.

On the event |a_v|<1/2, the exact center residual x(u)-a_v v has constant
values p=1-a_v on A and q=u-a_v on B. Here p>=1/2. If q>=0 then
p-|q|=1-u>=1/4; if q<0 then p-|q|=1+u-2a_v>=u>=1/4.
The nearest body line has squared distance(3/2)(p-|q|)^2. Consequently
the distance to the origin and to every body line is bounded below by

    h_v0=sqrt(3)/2-eta,
    h_v1=sqrt(3/32)-eta.                                 (4)

On |a_w|<1/16 the center residual has p=1-a_w>=15/16 and
q=u+a_w>=3/16, with p-q=1-u-2a_w>=1/8. The corresponding bounds are

    h_w0=15sqrt(3)/16-eta,
    h_w1=sqrt(3/128)-eta.                                (5)

All four numbers are positive. They apply to the first two positive
deletions, with parameters b=v or b=w omitted, respectively.

A uniform positive third distance along the parameter fiber is FALSE.
For example u=a_v=1/4 gives x(u)-a_v v=(3/4)1_A. This lies in the plane
of (1,1,1,1,1,-1) and(1,1,1,-1,-1,1), neither equal to v. Such zeros are
retained and integrated in the next section.

At the original x, however, every body pair plane U other than U0 has
distance at least

    g=1/(4sqrt(2))-eta>0.                                (6)

For x(u), a different unsigned partition mixes a coordinate1 and u in
some block, giving distance at least(1-u)/sqrt(2). A signed plane with
the same A/B partition but a different span has an opposite-sign pair
inside a block, giving distance at least sqrt(2)min(1,u). Other signed
partitions dominate their unsigned distance by034. This proves the
uniform1/(4sqrt(2)) center bound; subtract eta for D. The only body lines
in U0 are v and w (a body vector constant on A and B has both constants
equal to +/-1). Thus no pair avoiding v, or avoiding w, can span U0.

## 4. The true Cauchy average pays isolated distance zeros

For z,b in a Euclidean quotient, lambda>0 and s>0, define

    J_s(z,b)=integral_R f_s(a)/(|z-a b|^2+lambda s^2) da,
    f_s(a)=s/[pi(a^2+s^2)].

For b!=0 put L=|b|, a*=z.b/L^2,
gamma=sqrt(|z|^2-L^2 a*^2+lambda s^2)/L>0. The one-dimensional Cauchy
convolution identity gives exactly

    J_s=(s+gamma)/
         [L^2 gamma(a*^2+(s+gamma)^2)]
       =(1+s/gamma)/
         [|z|^2+(lambda+L^2)s^2+2L^2 s gamma].             (7)

The convolution identity follows, for example, by multiplying the
characteristic functions exp(-s|xi|) and exp(-gamma|xi|); all densities
are integrable positive Cauchy densities. Since gamma>=sqrt(lambda)s/L,

    J_s(z,b)<=(1+L/sqrt(lambda))/|z|^2   when z!=0.        (8)

If b=0 the result follows directly with numerator1. No lower bound on
inf_a|z-a b| was used. Formula(7) includes exactly the a values where
that distance vanishes. There is no asymptotic approximation.

## 5. All complementary events are bounded

Fix the retained parameter b=v or w, with the relevant small threshold.
For U spanned by k=0,1,2 already deleted non-b directions, the remaining
31-k projected columns have Gram

    32I_(Uperp)-(P_(Uperp)b)(P_(Uperp)b)^T,

between26I and32I. Thus the positive deletion identity applies using the
actual counts31,30,29, and the time terms992s^2,960s^2,928s^2. Its product
proof is the inequality |z|^2<=32sum a_i^2; each deleted column has its
true Gram-Schmidt coarea factor. Every three distinct body directions
are independent. Quotient densities remain continuous by the lower Gram
bound, so the positive test-function proof extends pointwise.

After three deletions, retain the actual integral of f_s(a_b). The first
two denominators are bounded by (4) or(5), but the third is kept as

    |P_(span(i,j))perp(x-a_b b)|^2+928s^2.

The remaining dimension3 quotient omits b, so apply(3), not the original
Q3 without its loss. Extend the positive a_b integral to all R and use
(8), with z=P_(span(i,j))perp x and |P_(span(i,j))perp b|<=sqrt(6).
All these pairs avoid b and hence (6) applies to their original x.

Let C_rem=44743079/(2970pi^6) be the complete ordered-triple coefficient
over directions other than one fixed body line. The4495 unordered
counts are435,1740,435,580,1305, hence26970 ordered triples; by frame
symmetry this is the same coefficient for b=v and b=w. Every one is
retained in each respective positive bound. Set

    L=(10/9)^3 C_rem (1+sqrt(6/928))/g^2,
    M_v=L/(h_v0^2 h_v1^2),
    M_w=L/(h_w0^2 h_w1^2), M=M_v+M_w.                    (9)

Then B_v,s<=M_v and B_w,s<=M_w on D for all s>0, with all non-v/non-w
large parameters and cancellation events included. Before dropping time
terms, the sharper valid denominators are(h_b0^2+992s^2)(h_b1^2+960s^2).
No Q2 or fourth deletion is used. We have proved the new pointwise bound

    p_s(x)<=C2 s^-2+M on D.                              (10)

The first estimate in(1) may replace C2 s^-2 when useful. Unlike a
uniform third-distance estimate, (9) stays valid across its fiber zeros.

## 6. Finite short-time norm, with one L1 budget per actual branch

The original full density satisfies the positive-face slice bound

    integral_F p_s dy<=s/[pi(1+1024s^2)]<=s/pi.

Since P_s<=p_s and the actual disjoint complement B_s<=p_s, interpolation
at q=6/5 and Minkowski imply

    ||integral_0^epsilon p_s ds||_Lq(D)
      <=(2/3) C2^(1/6)pi^(-5/6)epsilon^(3/2)
         +(6/11) M^(1/6)pi^(-5/6)epsilon^(11/6).           (11)

In the second term we interpolate the combined actual B_s using its
L-infinity upper bound M_v+M_w. We do not assert that the enlarged
one-threshold events are disjoint, nor assign the slice symmetry to
either event. Positive time truncation and monotone convergence justify
the norm statement even at singular zero-measure sets.

Only rho=1/64, u in[1/4,3/4], epsilon=1/1000 and the two stated thresholds
are evaluated. Other plane types, the spatial endpoints, times above
epsilon and the rest of the face are not paid by this result.
