# R085 independent uniform local lemma for 4<=s<=5 — ACCEPT

Fresh derivation; no s=4 or s=9/2 constant is assumed. Fix s,t,h,
epsilon,w during q differentiation. Sufficient local conditions are

    4<=s<=5, 1/s<=q<1, 0<h<=1/1000, rho=tq>=20,
    1<=R=sq-h<=5, a0=(1-epsilon)q,
    g=R-a0, 0<=g<=4R/5, epsilon*q<=h/s,
    |w|<=r<=min(5,R+2h).

The actual matching domain eventually meets these conditions; its stronger
relative bound g/R<=1-1/s can also be retained. They imply epsilon<=h,
a0>0, g<=4 and m=(|w|-a0)_+<=5. Put v=(g-m)_+ and

    ell(z)=min(1,(R+rho+a0-|z|)/(2R)), phi=-log ell,
    Q(R,u)=-log(1-u/(2R)), F(R,u)=2R Q(R,u)-u,
    J_j=(2rho)^-1 integral_{|z-w|<rho+a0, |z|<=rho} phi^j/ell dz.

Then ell>=3/5. The two boundary layers have width g, are disjoint,
and have the exact integrals

    integral_edge 1/ell = g+F(R,g),
    integral_edge phi/ell = R Q(R,g)^2,
    integral_edge phi^2/ell = (2R/3)Q(R,g)^3.

Reflection allows w>=0. When w<=a0 the full source interval is
captured. When w>a0 its left interval of length m is deleted.
This deletion cannot reach the opposite edge because m<=5 and
2rho-g>=36. If m<=g it removes only part of the near layer.
If m>g, it removes that entire layer and central length m-g,
where ell=1 and phi=0. Thus in both cases the exact formulas are

    J=1+[F(R,g)+F(R,v)-m]/(2rho),
    J1=R[Q(R,g)^2+Q(R,v)^2]/(2rho),
    J2=R[Q(R,g)^3+Q(R,v)^3]/(3rho).

In particular these formulas cover |w|>R, with v=0 and the necessary
-m contribution to J. Boundary points have zero coordinate mass.

For 0<=u<=4R/5, write z=u/(2R)<=2/5. The elementary bounds are

    0<=Q<=log(5/3)<3/5,
    |Q_R|<=2/(3R), Q_u<=5/(6R),
    0<=F<=6R/5<=6, 0<=F_u<=2/3,
    -4/3<=F_R=2Q-2u/(2R-u)<=0.

The strict logarithm bound uses exp(3/5)>1+3/5+9/50>5/3.
The sign F_R<=0 follows from -log(1-z)<=z/(1-z), proved by
integrating 1/(1-y)<=1/(1-z) on 0<=y<=z. The lower bound follows
from Q>=0. No old relative-gap estimate is used.

For A=F(R,g)+F(R,v)-m, -5<=A<=12, hence |A|<=12. Therefore

    |J-1|<=6/rho, J>=7/10>2/3,
    0<=J1<=9/(5rho), 0<=J2<=18/(25rho).

The affine functions R,a0,g,rho have denominators bounded away from
zero; the positive parts m,v are Lipschitz. Q,F are smooth near the
closed relative domain. Thus all three profiles and log J are AC.
The chain rule holds a.e. and its uniform estimates may be integrated
through all onset kinks.

Almost everywhere,

    R'=s<=5, g'=s-1+epsilon<=s,
    m'=0 or -(1-epsilon),
    v'=0,g',or g'-m', so 0<=v'<=s<=5.

Consequently

    |A'|<=2*5*(4/3)+2*5*(2/3)+1=21,
    |J'|<=[21+12*5]/(2rho)=81/(2rho),
    |(log J)'|<=243/(4rho)<120/rho.

For u=g or v, |d_q Q(R,u)|<=15/(2R). The prefactor has the
exact derivative (R/rho)'=h/(q*rho). Differentiating both squared
terms gives

    |J1'|<=9h/(25q*rho)+9/rho
          <=(9+9/5000)/rho<60/rho.

The value and derivative bounds used by the author are therefore
valid throughout the whole matching domain and uniformly in s in
[4,5]. In central deletion v=0 a.e.; the -m term remains in J.

One independent registered local Python/SciPy run (PID 31041) checked
the three exact integrals at 24 fixed configurations: s=4,4.5,5,
q=.3,.8, and full capture, edge deletion, positive/negative central
deletion. Maximum absolute discrepancy was 2.22e-16. Maximum
rho-scaled finite-difference diagnostics for log J,J1,J2 were
2.5064,.3531,.0616. These finite diagnostics support consistency;
the proof above supplies the uniform local theorem.
