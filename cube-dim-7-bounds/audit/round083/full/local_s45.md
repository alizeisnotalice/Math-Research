# R083 independent s=9/2 coordinate lemma — ACCEPT

This is a fresh derivation on the s=9/2 domain. No s=4 value or derivative
constant is imported. The directory is independent staging, not a clone.
All author and original-repository sources are read only.

Fix t,h,epsilon,w while differentiating in q. Sufficient local conditions
on a connected q interval are

    2/9<=q<=9999/10000, 0<h<=1/1000, rho=tq>=16,
    1<=R=(9/2)q-h<=9/2, a0=(1-epsilon)q,
    g=R-a0, 0<=g<=7R/9, epsilon*q<=2h/9,
    |w|<=r<=min(9/2,R+2h).

These imply 0<=epsilon<=h, a0>=2R/9>0, g<=7/2, and
m=(|w|-a0)_+<=9/2. They hold throughout every actual R083 matching
interval for all sufficiently large L, after adding rho>=16 and h<=.001
to the finitely many eventual guards.

Set ell(z)=min(1,(R+rho+a0-|z|)/(2R)), phi=-log ell and
capture={|z-w|<rho+a0}, for -rho<=z<=rho. Then ell>=11/18.
For j=0,1,2 put J_j=(2rho)^-1 integral_capture phi^j/ell dz,
with J=J_0. Let Q(R,u)=-log(1-u/(2R)), F(R,u)=2R Q(R,u)-u,
v=(g-m)_+. The exact formulas are

    J=1+[F(R,g)+F(R,v)-m]/(2rho),
    J1=R[Q(R,g)^2+Q(R,v)^2]/(2rho),
    J2=R[Q(R,g)^3+Q(R,v)^3]/(3rho).

By reflection take w>=0. The two nonconstant edge layers each have
width g and are disjoint, since rho>=16>g. At distance y inward from
an endpoint, u=g-y and ell=1-u/(2R). Direct substitution gives the
edge integrals

    integral 1/ell = g+F(R,g),
    integral phi/ell = R Q(R,g)^2,
    integral phi^2/ell = (2R/3)Q(R,g)^3.

If w<=a0, no coordinate interval is deleted and m=0, v=g. If w>a0,
capture deletes exactly the left interval of length m. Since m<=9/2
and 2rho-g>=57/2, this deletion does not reach the opposite edge.
For m<=g its removed inverse integral is m+F(R,g)-F(R,g-m).
For m>g it deletes the whole near edge and central length m-g,
where ell=1 and phi=0. The removed inverse integral is then
m+F(R,g), while the removed J1,J2 contributions are the whole near
edge. Thus the formulas hold also in the central-deletion regime
|w|>R. The -m term in J is indispensable there. Endpoints have
zero coordinate Lebesgue mass.

For 0<=u<=7R/9 one has

    0<=Q<=log(18/11)<1/2,
    |Q_R|<=7/(11R), Q_u<=9/(11R),
    0<=F<=R, 0<=F_u<=7/11, -14/11<=F_R<=0.

The strict log bound follows from exp(1/2)>79/48>18/11.
For F_R=2Q-2u/(2R-u), the upper sign follows by integrating
1/(1-z)<=1/(1-z0) for z<=z0; the lower bound follows from Q>=0.
These bounds are independent of any old local lemma.

Writing A=F(R,g)+F(R,v)-m gives -9/2<=A<=9 and |A|<=9.
Hence |J-1|<=9/(2rho) and J>=23/32>2/3 when rho>=16.
The new moment bounds are

    0<=J1<=R/(4rho)<=9/(8rho),
    0<=J2<=R/(12rho)<=3/(8rho).

All affine functions have denominators bounded away from zero. The
positive parts m,v are Lipschitz. Q,F are smooth on a neighborhood
of the closed relative domain. Therefore J,J1,J2 and log J are AC
in q. Kinks need only a.e. chain rules, and all resulting bounds can
be integrated through them.

Almost everywhere, R'=9/2, g'=7/2+epsilon<=9/2,
m'=0 or -(1-epsilon), and v'=0,g',or g'-m', so 0<=v'<=9/2.
Thus

    |A'|<=2*(9/2)*(14/11)+2*(9/2)*(7/11)+1=200/11,
    |J'|<=[(200/11)+9*(9/2)]/(2rho)=1291/(44rho),
    |(log J)'|<=3873/(88rho)<90/rho.

For u=g or v, |d_q Q(R,u)|<=72/(11R). The exact prefactor
derivative is (R/rho)'=h/(q*rho). Differentiating both squared
Q terms therefore gives

    |J1'|<=h/(4q*rho)+72/(11rho)
          <=(9/8000+72/11)/rho<40/rho.

For completeness, differentiating both cubed terms gives
|J2'|<=(3/8000+36/11)/rho; no J2 transport is needed by R083.
All bounds hold when v=0 and in the central-deletion region.

The bounds used in the author proof, |J-1|<7/rho,
|d_q log J|<=90/rho, J1<=9/(8rho), |J1'|<40/rho and
J2<=3/(8rho), are therefore justified on the actual full matching
domain with explicit sufficient eventual guards.

Actual supplementary computation: one independent local Python/SciPy
run, PID 29940, checked all three integrals at twelve fixed
configurations spanning full capture, edge deletion, central deletion,
and both signs. Maximum absolute error was 2.22e-16. The maximum
rho-scaled finite-difference magnitudes for log J,J1,J2 were about
2.2552,.4161,.1106. These are diagnostics, not a proof of uniformity.
The analytic derivation above supplies the uniform lemma.
