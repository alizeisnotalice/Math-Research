# Positive-band simplex Taylor lower certificates

Fix the true body32 Cauchy semigroup p_t, K_T=integral_0^T p_t, and a finite
time partition 0<a=t0<...<tn=T. Put

    B=B_(a,T)=integral_a^T p_t=H_a-H_T,
    J2_B(x)=integral_a^T t^-2 p_t(x)dt.

All quantities are finite and smooth on all space because a>0. In particular
there is no difference of infinite Green values. Positivity gives K_T>=B.

## 1. Whole-cell moment upper bound with no infinite tail

For a fixed center c let Delta_i=H_(ti)(c)-H_(ti+1)(c)>=0. The real Cauchy
parameter translation comparison is

    p_t(c+h)<=Phi(||h||_infinity/t)p_t(c),
    Phi(u)=exp(12 asinh(u/12)).

For completeness, one factor has exact ratio bound
f_t(z+b)/f_t(z)<=exp(2asinh(|b|/(2t))). The cube gauge represents h using
at most six signed body directions with coefficient absolute sum||h||inf:
h/||h||inf is on a five-dimensional cube facet, and its convex combination
uses at most six vertices. Parameter translation and concavity of asinh
therefore give the displayed comparison. It is a positive pushforward
argument; a zero-probability fiber is never used as a conditional law.

Let C=conv{v0,...,vm} be a convex polytope, c in C, and

    r=max_j ||vj-c||infinity,
    D=max_j ||vj-c||2^2.

The norm convexity implies ||x-c||infinity<=r for every x in C. On each
time band, t^-2 Phi(r/t)<=ti^-2 Phi(r/ti). Hence everywhere in C,

    J2_B(x)<=U2
       :=sum_(i=0)^(n-1) ti^-2 Phi(r/ti) upper(Delta_i).          (1)

There is no [T,infinity) term. For interval time values, obtain the mass
upper by outward subtraction followed by the upper endpoint; exact mass
positivity permits intersecting its lower endpoint with zero. If the
subtraction's upper endpoint is negative, reject the inconsistent input.
All weights in (1) are positive.

## 2. The correct lower Hessian and its Taylor coefficient

Let A have the32 body directions as columns. AA^T=32I and its orthogonal
right inverse R=A^T/32 satisfies ||Rh||2^2=||h||2^2/32. In the coarea product
F=product_v f_t(z_v), the logarithmic second derivative of a Cauchy factor
is at least-2/t^2. Therefore

    D_h^2 F/F=(D_h log F)^2+D_h^2 log F
       >=-2||Rh||2^2/t^2=-||h||2^2/(16t^2).

Positive coarea/time integration gives

    D_h^2 B(x)>=-||h||2^2 J2_B(x)/16.                           (2)

The derivative domination is legitimate locally: real Cauchy shifts have
finite multiplicative bounds and the first two factor derivatives are
bounded relative to the factor by constants times t^-1,t^-2; a<=t<=T.
This uses the negative Hessian lower bound, not a Hessian upper bound.

For x=c+h in C, convexity keeps the entire segment c+s h inside C.
Taylor's integral remainder and integral_0^1(1-s)ds=1/2 give

    B(c+h)>=B(c)+gradB(c).h-||h||2^2 U2/32.                    (3)

Both the value and gradient must belong to exactly the same a,T,center and
frame. In particular gradB(c)=gradH_a(c)-gradH_T(c).

The linear expression attains its minimum over a finite convex hull at a
vertex. The squared norm is convex, so ||x-c||2^2<=D. Consequently the
actual cell-wide certificate is

    K_T(x)>=B(x)
      >= L_B+min_j lower(gradB(c).(vj-c))-D*upper(U2)/32.         (4)

Here L_B is a lower endpoint of H_a(c)-H_T(c), every vertex dot product
uses outward interval arithmetic and a lower endpoint, and U2 uses upper
mass endpoints. Taking an exact rational minimum of the endpoint bounds is
safe. Formula (4) is valid for a simplex or another convex cell; it assumes
no convexity of B or K_T. A rational barycentric witness c=sum lambda_j vj,
lambda_j>=0, sum lambda_j=1 verifies the center/segment condition.

## 3. Augmented scalar-Gram extension, without transferring old values

For A_beta=[body32,beta e1,...,beta e6], beta>0,
AA^T=(32+beta^2)I. Thus ||Rh||2^2=||h||2^2/(32+beta^2), and the lower
Hessian in (2) has coefficient2/(32+beta^2) instead of1/16. Its Taylor
loss in (3)-(4) is

    D*upper(U2)/(32+beta^2).

Phi remains valid: represent the translation using the original six sparse
body columns and leave the extra axis variables unchanged. All H and
gradient data must nevertheless be for the augmented frame. The body32
numbers and the body32 A009 singular barrier do not automatically become
augmented-kernel certificates.

## 4. Exact global face tree and a same-plane tube rule

The ordered positive face is the closed5-simplex with six vertices

    wj=(1,1,...,1,0,...,0),  j=0,...,5,

where j free coordinates after the first are1. Its barycentric coordinates
are 1-y1,y1-y2,...,y4-y5,y5. This includes every tie and face boundary.

Split a longest edge (chosen with exact rational squared lengths and a
deterministic tie rule) at its exact midpoint. Replace one endpoint by the
midpoint in each child. The two closed children have exactly the parent
as union; their interiors are disjoint. In parent barycentrics, compare the
two endpoint weights: if lambda_i>=lambda_j, use weights lambda_i-lambda_j
on vi and2lambda_j on the midpoint, with all other weights unchanged;
the reverse inequality gives the other child. Thus every finite recorded
split tree preserves full coverage without a floating geometric tolerance.

For the base body32 model, the five A012 canonical planes in this chamber
are U_k=span(1_(1..k),1_(k+1..6)). A cell may be marked TUBE only if there
is one fixed k such that dist(vj,U_k)<=delta for every vertex. Projection
is linear and Euclidean norm is convex, so that same plane then contains
the entire cell's delta tube. It is invalid to choose a different closest
plane at each vertex. For example every root vertex is in S2 while the
root barycenter has squared distance1/9 from S2.

At T48=6/125, positivity gives K_T48>=K_(4/125). The proved A009 body32
barrier yields

    rho<=1/100000 => K_T48>=3/250000.

Since 3/250000> T48/4096=3/256000, this is a valid TUBE lower value. The
barrier is used at its actual distance, not as a claim that tube unions
are disjoint. For a BAND leaf, use (4) with a verified center-data source.
Every terminal leaf must have a rigorous lower bound at least the chosen
m_star>T/4096. An unresolved leaf, missing child, unproved region or tiny
uncovered volume invalidates a global certificate.

## 5. From a completely paid ordered face to the cube

For body32 the kernel is invariant under signed permutations and obeys
p_t(r y)=r^-6 p_(t/r)(y). Thus for0<r<=1,

    K_T(r y)=r^-5 K_(T/r)(y)>=K_T(y).

For each nonzero cube point choose r=||x||infinity and normalize/sort its
absolute coordinates to y in the ordered positive face. A completely paid
root tree therefore proves K_T(x)>=m_star on the whole cube. At zero,
K_T(0)=+infinity. This is the required fixed-T transfer, not a homogeneity
claim about H_a at fixed a. With the separately proved source interface,
m_T>=m_star>T/4096 would imply a strict cost below64.

This turn supplies the acceptance interface and a local calibration only.
It does not supply a paid root tree or any new global m_T. In particular
T=.032 remains refuted by the independent B010 counterexample elsewhere.
The later B015 counterexample also refutes the fixed body32 T=.048 target.
The T48 tube inequality above remains true on its stated neighborhood;
it is not evidence for a viable full-domain T48 certificate. The theorem
can be reused at a new T or frame only with matching newly verified data
and complete coverage; this turn proposes and evaluates no such new case.
