# R081 independent s=4 coordinate lemma — ACCEPT with explicit domain

This is a fresh s=4 derivation. No s=3/2 numerical constant is imported.
It is an independent temporary directory, not a Git clone. The original
repository and author files are read only; nothing was committed/pushed.

Fix t>0, h in [0,1/1000], epsilon>=0 and a FIXED target coordinate w.
When differentiating in q, t,h,epsilon,w are fixed. On any q interval
assume

    1/4<=q<=99/100, rho=tq>=12,
    R=4q-h in [1,4], a0=(1-epsilon)q,
    g=R-a0, 0<=g/R<=3/4,
    epsilon*q<=h/4, |w|<=r<=R+2h.

The assumptions imply epsilon<=h<=1/1000, a0>=R/4>0,
g<=3 and m=(|w|-a0)_+<=g+2h<4. These explicit sufficient
conditions replace the imprecise phrases “h small” and “rho sufficiently
large.” They apply to every connected matching interval on which the
assumptions hold. If w or t were allowed to vary with q, the derivative
claims would require their additional derivatives.

Set

    ell(z)=min(1,(R+rho+a0-|z|)/(2R)), |z|<=rho,
    capture={|z-w|<rho+a0}, phi=-log ell,
    J_j=(2rho)^-1 integral_capture phi(z)^j/ell(z) dz,
    J=J_0,
    m=(|w|-a0)_+, v=(g-m)_+,
    Q(R,u)=-log(1-u/(2R)), F(R,u)=2R Q(R,u)-u.

Then the formulas in the question are EXACT:

    J=1+[F(R,g)+F(R,v)-m]/(2rho),
    J1=R[Q(R,g)^2+Q(R,v)^2]/(2rho),
    J2=R[Q(R,g)^3+Q(R,v)^3]/(3rho).

They and log J are absolutely continuous in q, with a.e. bounds

    |J-1|<=4/rho<=6/rho,
    |d_q log J|<=711/(20rho)<100/rho,
    0<=J1<=1/rho,
    |d_q J1|<=5601/(1000rho)<40/rho,
    0<=J2<=1/(3rho).

An additional independently derived transport bound, useful for a new
second-profile argument, is |d_q J2|<3/rho a.e.

## 1. Exact edge integrals and central deletion

By evenness of ell we may reflect w and take w>=0. On each edge
layer of width g, parameterize u=g-y, with y the distance inward
from its source endpoint. Then ell=1-u/(2R), 0<=u<=g. In the
remaining central interval [-rho+g,rho-g], ell=1 and phi=0.
The layers do not overlap because rho>=12>g.

The complete inverse integral on either layer is g+F(R,g). More
generally the layer integrals of phi/ell and phi^2/ell are

    R Q(R,g)^2, (2R/3)Q(R,g)^3,

by the substitution Q=-log(1-u/(2R)), dQ=du/(2R-u).
The full cube therefore has inverse integral 2rho+2F(R,g) and two
equal nonzero edge-profile contributions.

If w<=a0, the whole source coordinate interval [-rho,rho] is captured,
m=0 and v=g; the formulas follow. If w>a0, the intersection with
capture deletes exactly the left interval [-rho,-rho+m), length m.
It cannot reach the right edge: m<4<2rho-g. When m<=g it deletes
part of the left edge, whose remaining length is v=g-m. The inverse
mass removed is m+F(R,g)-F(R,v). When m>g, it deletes the entire
left edge plus central length m-g; the latter has ell=1, phi=0.
The inverse mass removed is again m+F(R,g), whereas the nonzero
profile contribution removed is exactly the whole left edge.
With v=0 both cases are the displayed formulas. This explicitly
covers |w|>R, which is equivalent to m>g. The subtraction -m is
essential for J in this central-deletion regime. Boundary points have
zero coordinate Lebesgue mass.

## 2. Elementary uniform bounds

For 0<=u<=3R/4, let z=u/(2R)<=3/8. Then

    0<=Q<=log(8/5)<1/2,
    0<=F=2RQ-u<=2RQ<=R,
    0<=F_u=u/(2R-u)<=3/5,
    F_R=2Q-2u/(2R-u) in [-6/5,0],
    |Q_R|=u/[R(2R-u)]<=3/(5R),
    Q_u=1/(2R-u)<=4/(5R).

Here F_R<=0 follows from -log(1-z)<=z/(1-z), obtained by
integrating 1/(1-v)<=1/(1-z) from 0 to z. The strict log bound
needs no floating value: exp(1/2)>1+1/2+1/8=13/8>8/5.

Write A=F(R,g)+F(R,v)-m. Since m<4 and each F lies in [0,R],
-4<=A<=8 and |A|<=8. Therefore |J-1|<=4/rho and J>=2/3
for rho>=12. (The asymmetric lower bound gives an even stronger
J>=5/6, but is not needed.) The profile value bounds are

    J1<=R*(1/2)/(2rho)=R/(4rho)<=1/rho,
    J2<=R*(1/4)/(3rho)=R/(12rho)<=1/(3rho).

## 3. AC, positive-part kinks, and derivative bounds

R,a0,g and rho are affine in q; R and rho are bounded away from
zero. The functions m=(|w|-a0)_+ and v=(g-m)_+ are Lipschitz.
Q and F are smooth on a neighborhood of the compact relative domain
0<=u<=3R/4, R>=1. Their compositions and the displayed J profiles
are AC. Since J>=2/3, log J is AC as well. These statements permit
integration of all the a.e. derivative bounds across onset kinks.

Almost everywhere,

    R'=4, g'=3+epsilon in [3,4],
    m'=0 or -(1-epsilon), |m'|<=1,
    v'=0, g', or g'-m', so 0<=v'<=4.

At a kink no pointwise derivative claim is necessary. The ordinary
chain rule a.e. and the uniform partial bounds above give

    |A'|<=2*4*(6/5)+2*4*(3/5)+1=77/5.

Since rho'/rho=1/q and 1/q<=4,

    |J'|=|A'-A/q|/(2rho)
       <=(77/5+32)/(2rho)=237/(10rho),
    |(log J)'|<=237/(10rho)/(2/3)=711/(20rho).

For either u=g or u=v,

    |d_q Q(R,u)|<=4*3/(5R)+4*4/(5R)=28/(5R).

The prefactor has a useful EXACT derivative:

    (R/rho)'=(4-R/q)/rho=h/(q*rho).

Differentiating BOTH squared Q terms and the prefactor,

    |J1'|<=h/(4q*rho)+28/(5rho)
       <=[1/1000+28/5]/rho=5601/(1000rho).

Similarly differentiation of BOTH cubic Q terms and R/(3rho) yields

    |J2'|<=h/(12q*rho)+14/(5rho)
       <=[1/3000+14/5]/rho<3/rho.

These constants are uniformly valid in the central-deletion region
v=0 as well: its derivative vanishes a.e., while the -m term in A
was retained throughout. No s=3/2 derivative or moment constant was
used.

## 4. Actual computation and limits of evidence

The independent `check.py` ran with local Python/SciPy, PID 28899.
Exact Fraction arithmetic verifies the elementary rational constant
comparisons above. Direct quadrature checked all three integral
formulas at 24 fixed points spanning full capture, edge deletion and
central deletion, including negative w. The maximum absolute formula
error was 1.13e-16. Finite differences were diagnostic only; their
maximum rho-scaled magnitudes for log J,J1,J2 were about 2.039,
0.321,0.089. This is bounded sanity evidence, not proof of a uniform
bound. The analytic argument in Sections 1--3 supplies that proof.
