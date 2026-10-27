# R066 · Typical-target first profile and bounded inverse-source payment

Status: full proposed operator mechanism, pending independent review. This
source adds a first-profile condition to the ambient target set. It does not
deduce an operator theorem from a numerical model, or from the scalar
certificate alone. All parameters are fixed, and the conclusion below is
only for a common sufficiently large dimension, with no explicit threshold.

## 1. Fixed tuple and the proposed complete conclusion

Write L=log n, T=n log n−n log log n, and set

    s=3/2, qmin=2/3, qmax=917/1000,
    delta=2^-30, rlo=1+delta, rhi=2751/2000-delta,
    sigma=83/1000, chi=9/200,
    x=(1-sigma)L, t=n/x, beta=1+1/t,
    zeta=n^-3, epsilon=zeta t, kappa=beta-zeta,
    h=2^-80/L, theta=exp(-chi x), k=nh theta/rhi.

Use e=2^-20 for eps_J, eps_N, eps_S, eps_T, eps_D and eps_R (the
relative radial-width error). Set eps_Ret=1/4096 and eps_A=1/65536.
Define the exact rational constants

    P0=1-e-2^-60,
    Ulo=P0(1-e)/(1+e), Uhi=(1+e)/[P0(1-e)],
    M=2(1+e)Uhi,
    b=2(1-eps_Ret)(1-e)Ulo-eps_S-eps_T,
    mL=637/2000, Gamma=mL b/M, gamma=317/1000.

The proposed theorem is, for every sufficiently large integer n,

    C_n <= F(A_n(t),B_n(t)-gamma n) < T-(9/40)n,          (1)
    F(A,B)=(sqrt A+sqrt B)^2.

Here A_n,B_n and C_n are precisely the original cover/load quantities
and centered common-radius cube weak-(1,1) constant. Arbitrary L1 input,
all levels and every positive witness radius are retained. Sections 8–9
include the complete original cover factor; a smaller model load by itself
is not being called a C_n improvement.

The new exact Fraction calculation in verify.py proves

    (1-eps_A)log(rhi/rlo)>mL,
    (1-eps_D)Gamma>gamma,
    gamma-sigma/(1-sigma)=207689/917000,
    gamma-sigma/(1-sigma)-9/40=341/229250>0.              (2)

The logarithm enclosures in that check are rational positive atanh series
with explicit geometric tails. There is no production search or sampling.

## 2. Full geometry, bounded kernel, and active reference radius

Use the original limiting corner measure

    D={||z||_infinity<t: some z_i<1-beta m(z)},
    m(z)=max(1,||z||_infinity),
    d lambda=2^-n m(z)^-n 1_D dz.

Before cap restriction let lambda0 denote the same density on the full
annuli. At a shell radius rho>1 each signed maximal face has measure
d rho/(2rho), with n−1 independent Uniform[-rho,rho] free coordinates.

Let A0={w:r=||w||_infinity in [rlo,rhi]}, let the FULL source band be
P0set={z:q=||z||_infinity/t in [qmin,qmax]}, and put P=D intersect P0set.
The notation P0set is distinct from the scalar product bound P0 above.
Matching is |r-sq|<h. For the entire source band, before deleting D's
complement, define

    rho=tq, R=sq-h, a0=(kappa-1)rho=(1-epsilon)q,
    g=R-a0=(s-1+epsilon)q-h,
    ell_R(z_i)=min(1,(R+kappa rho-|z_i|)/(2R)),
    G(z)=product_i ell_R(z_i),
    b_soft(z)=min(1,theta/G(z)).                         (3)

Eventually R>0 and a0>0 throughout P0set, so G>0 and 0<b_soft<=1.
Extend b_soft=0 off this band. The target radial factor is
a_rad(w)=r/rhi on A0 and zero off A0, so 0<=a_rad<=1.
After A is specified below, the final Borel kernel is exactly

    K(w,z)=1_A(w)1_P(z)1_{|r-sq|<h}
       1_{||z-w||_infinity<kappa rho} a_rad(w)b_soft(z). (4)

Write R_lambda(w)=integral K(w,z)d lambda(z) and
C_lambda(z)=integral K(w,z)d lambda(w).

For a target r the full matching interval has endpoints

    q_-=(r-h)/s, q_+=(r+h)/s, q_b=q_+.

The fixed trims and eventual h<delta/2 imply it lies in [qmin,qmax].
At the baseline q_b, R=r exactly; over matching R is in [r−2h,r].
An active source is one matching some target in A0. Its R exceeds
rlo−2h>1. Some unmatched first-band sources can have R<1; they have
zero final column and no local estimate requiring R>=1 is used on them.

The following are uniform eventual guards, now checked for qmax=.917,
not inferred from the older .909 endpoint:

    q>=2/3, 1<R<=r<3/2, rho>=t qmin>=4rhi>4,
    kappa>1, beta<2, a0>0, 0<=epsilon<1/2,
    h<min(delta/2,1/1000), epsilon qmax<=2h/3.

They imply, on every active match,

    0<g, g-R/3=epsilon q-2h/3<=0.                     (5)

Thus the relevant relative gap is g/R<=1/3; g itself may exceed 1/3.
The old absolute-gap condition is not reused. Targets and sources are
disjoint point sets since r<3/2 and rho>4, and both lie inside norm<t.

## 3. Coordinate normalizers and the new first profile

For |v|<=r define the full-source inverse-capture normalizer and profile

    J(v;q)=(2rho)^-1 integral_[-rho,rho]
                    1_{|z-v|<kappa rho}/ell_R(z) dz,
    Jphi(v;q)=(2rho)^-1 integral_[-rho,rho]
           1_{|z-v|<kappa rho}[-log ell_R(z)]/ell_R(z) dz.

Both integrals are taken BEFORE restriction to D. Set

    F(R,u)=-2R log(1-u/(2R))-u,
    m=(|v|-a0)_+, u0=(g-m)_+,
    Q(R,u)=-log(1-u/(2R)), H(R,u)=R Q(R,u)^2.

Direct tail integration gives

    J=1+[F(R,g)+F(R,u0)-m]/(2rho),
    Jphi=[H(R,g)+H(R,u0)]/(2rho).                       (6)

If m>g a missing central segment has profile zero, so the second formula
still holds with u0=0. Strict endpoints have no effect on these integrals.

For u<=g<=R/3, Q<=log(6/5)<1/5. One may use the safe bounds

    0<=H<=R/25, 0<=H_u<=6/25,
    |H_R|<=3/25.

For example Q_u=1/(2R-u) and Q_R=−u/[R(2R-u)], which give these
inequalities directly. On matching, R'=3/2, g'=1/2+epsilon<1;
where m>0, u0=(R−|v|)_+ has derivative 0 or 3/2, and where m=0
u0=g. Thus each H term has a.e. derivative at most 27/50 in absolute
value. Since (H(R,g)+H(R,u0))/q<=2R/(25q)<3/25, differentiation
of (6), including rho'=rho/q, gives the conservative bound

    0<=Jphi<=1/rho, |partial_q Jphi|<=8/rho.             (7)

Positive-part onsets give locally Lipschitz functions, so the a.e.
derivative integrates by absolute continuity. No everywhere-smoothness
claim is needed at capture onsets. In fact the constants above imply a
much smaller derivative, but no such sharpening is used.

The relative J-coordinate estimates are those proved in R064's local
lemma, and can also be obtained by differentiating the first formula (6):

    |partial_q log J|<=256/rho,
    |J-1|<1/(3rho),                                    (8)

uniformly on active matches. At q_b where R=r, if W is Uniform[-r,r],
Fubini gives E_W J=1, and the baseline bounds are

    -6g^3/rho^2<=E log J<=0,
    Var(log J)<=12g^3/rho^2.                            (9)

The pinned target coordinate contributes a deterministic log term in
[−g/rho,0]. These relative estimates require only R>=1, g/R<=1/3,
rho>=4r and the displayed matching guards. The larger qmax only changes
fixed compact endpoints and does not invalidate them.

### Exact target mean of Jphi

At the BASELINE R=r, ell_R(z) is exactly the target-coordinate capture
fraction E_W 1_{|z-W|<kappa rho}. Fubini cancels the inverse factor:

    E_W Jphi(W;q_b)=E_Z[-log ell_R(Z)]
                  =(2R/rho)[1-a+a log a],
    a=ell_pin=1-g/(2R)>=5/6.                           (10)

The bracket decreases with a on (0,1). Since R/q=s at the baseline
up to the harmless negative h/q term,

    n E_W Jphi<=mu0 x,
    mu0=3[1/6+(5/6)log(5/6)]
       =1/2-(5/2)log(6/5)<221/5000.                    (11)

This differs from dropping the capture indicator in each source
coordinate. Equation (10) is an ambient TARGET-law identity, not an
identity for an actual unknown increment measure or for a source law.

## 4. Five ambient target conditions, probability and profile transport

Let N_i(q)=1+1_{|w_i|<a0(q)}, Nsum=sum_i N_i, and
ell_pin=(1+a0/R)/2. At either endpoint q±, under a full target-face law,

    E Nsum(q±)=n+(n−1)a0(q±)/r, Var Nsum(q±)<=n/4.

Define A to be A0 intersected with the following five conditions:

1. |sum_i log J(w_i;q_b)|<=e.
2. |Nsum(q_-)-E Nsum(q_-)|<=e n/4.
3. |Nsum(q_+)-E Nsum(q_+)|<=e n/4.
4. w belongs to D.
5. sum_i Jphi(w_i;q_b)<=(221/5000+1/100000)x.

All failure probabilities below are uniform in the compact target
radius interval and in every signed maximal target face. The baseline
log mean and variance from (9) are O(L^2/n)+O(L/n), respectively
O(L^2/n); hence the first failure probability is o(1). Each count
failure is at most 4/(e² n). The cap complement on a positive pinned
target face is at most (3/4)^(n−1) eventually (rhi<3/2 and t tends
to infinity); on a negative face it is zero.

For condition 5 the n−1 free target coordinates are independent. By
(7), their total variance is at most n/rho_b²=O(L²/n). The pinned
term is at most 1/rho_b. Equation (11) bounds the mean of their sum
by mu0 x+1/rho_b. The threshold in condition 5 exceeds this by a
fixed positive multiple of x eventually. Chebyshev therefore gives
uniform failure o(1). No actual μ concentration is being assumed.

The union of five failures is at most eps_A eventually. Consequently

    lambda(A)>=(1-eps_A)n log(rhi/rlo)>mL n.             (12)

For the two count events, Nsum increases while

    partial_q ell_pin=−(1−epsilon)h/(2R²)<0.

The ratio Nsum/(2n ell_pin) is therefore nondecreasing. At q+ its
mean ratio is at most one; at q− its shortfall is at most 2h+1/n.
The endpoint fluctuation changes the ratio by at most e/4. Thus when
2h+1/n<=3e/4, on the ENTIRE matching interval,

    1-e<=Nsum/(2n ell_pin)<=1+e.                        (13)

The q-width is 2h/s. By (8) and n/rho<=x/qmin,

    |sum log J(q)-sum log J(q_b)|<=512Lh<2^-60.

The log baseline and e^v<=1/(1−v), e^-v>=1−v now give
P0<=product J_i<=1/P0. Also |J_i−1|<1/(3rho) implies each
J_i belongs to [1−e,1+e] eventually, pointwise and uniformly; no
union of n coordinate exceptional sets is used.

Similarly (7) gives

    |sum Jphi(q)-sum Jphi(q_b)|<=16 x h<=x/100000

eventually. Combining with condition 5 yields the uniform profile

    sum_i Jphi(w_i;q)<= (2211/50000)x.                  (14)

The matching bridge is deterministic on A. The target law was used
only to show that A has large ambient mass.

## 5. Every inverse source face: improved mean, unchanged safe variance

The complete inverse-capture coefficient per d rho/rho is

    U(w,rho)=[product_i J_i/(2ell_pin)] sum_i N_i/J_i.

Equations (13) and the product/individual bounds imply

    Ulo n<=U(w,rho)<=Uhi n.                             (15)

Fix ONE captured signed source face i. Its pin contributes Qpin
deterministically to V=−log G. Its free coordinates are independent
under densities 1_capture/(ell_R J_l). Their means are exactly
Jphi(w_l;q)/J_l. By nonnegativity, (14), and J_l>=1−e,

    E_inv(V | face i)
       =Qpin+sum_(l!=i) Jphi(w_l;q)/J_l
       <=Qpin+[2211/(50000(1−e))]x
       <Qpin+(223/5000)x.                              (16)

This is the new mean estimate: the target profile controls every
source face, without replacing capture by the whole source interval.

For the variance, retaining the old conservative fullsource bound is
valid. For a uniform FULLSOURCE free coordinate,

    E_full(phi²/ell_R)=2R Qpin³/(3rho),
    Qpin<=log(6/5)<Q5=34187/187500<1/5.

Dropping capture only in this nonnegative second moment and dividing
by J_l>=1−e gives

    Var_inv(V | face i)
       <=x [2s Q5³/(3(1−e))]<x/150.                   (17)

The exact Fraction gates verify (16) and (17). Eventually Qpin/x is
at most 1/5000. Since

    chi−223/5000−1/5000=1/5000,

Chebyshev on EACH captured face gives

    P_inv(V>chi x | face i)<=5000²/(150x)->0.            (18)

This is eventually at most eps_Ret. Only now mix these uniformly
bounded probabilities across source faces. No mixture variance is
identified with the sum of within-face variances; no source-cap-
conditioned independence or unrelated worst target is used.

On V<=chi x, G>=theta, so the ratio of b_soft to theta/G equals
one; everywhere that ratio lies in [0,1]. Thus the full-source SOFT
shell coefficient is between

    (1-eps_Ret)theta Ulo n and theta Uhi n.              (19)

## 6. Calibrated full rows and exact complete columns

For fixed r, the matching SOURCE log-width is exactly
log((r+h)/(r−h)). With the target factor r/rhi and u=h/r,

    1<=log((1+u)/(1−u))/(2u)<=1/(1−u²).

The last bound tends uniformly to one. Therefore

    (1−e)2h/rhi <=(r/rhi)log((r+h)/(r−h))
                                  <=(1+e)2h/rhi.

The factor s cancels in the logarithmic width; its ordinary q-width
is still 2h/s. Combining with (19) yields raw full-source rows

    2(1-eps_Ret)(1−e)Ulo k <=R_raw(w)<=Mk.              (20)

For the complete column fix z and an active target radius r>=R. Put
a_i=kappa rho−|z_i|>0. Its free target-coordinate fraction is
ell_i(r)=min(1,(r+a_i)/(2r)). If ell_i<1 exactly one target sign
is captured, so its sign count T_i=1; if ell_i=1 then T_i<=2.
Consequently T_i/(2ell_i)<=1, including tangencies in the upper
direction. Summing ALL pinned target faces gives shell coefficient

    (1/2)sum_i T_i product_(l!=i)ell_l(r)
        <=n product_i ell_i(r)<=n G(z).

The last inequality is monotonicity in r>=R. The radial weight now
cancels the Jacobian exactly: (dr/r)a_rad=dr/rhi. The matched
ordinary TARGET interval has width at most 2h. Hence

    C_raw(z)<=2nh G(z)/rhi,
    C_lambda(z)<=2nh theta/rhi=2k.                     (21)

An inactive source has zero column. Good-target and cap restrictions
can only decrease (21). Since M>2, Mk bounds both final degrees.

## 7. Source cap and one global selection threshold

The whole source band has lambda0 mass at most Csrc n, with
Csrc=(qmax−qmin)/qmin. Its corner-cap complement has negative pinned
faces excluded, and each free coordinate has fraction
1−(1−q)/(2rho). Thus its mass is at most

    Csrc n exp[−(1−1/n)eta x],
    eta=(1−qmax)/(2qmax)=83/1834.

Because a_rad,b_soft<=1, deleting this complement costs at most its
unweighted mass in every row. Relative to the actual normalization k,
the loss is bounded by

    (Csrc rhi/h)exp[−((1−1/n)eta−chi)x]->0,             (22)
    eta−chi=47/183400>0.

It is therefore at most eps_S k eventually. The cap step occurs only
after the fullsource per-face law, not inside its independence claim.

Choose the one original selection threshold

    tau=sqrt[A_n(t)(B_n(t)−gamma n)].

The original scalar bounds give B_n>=nL/2 eventually, A_n<=e^x=
n e^-sigma L, B_n<=nL. Thus the square root is positive and

    tau/k <=(rhi sqrt L/h)
       exp[−(sigma/2−chi(1−sigma))L]->0,
    sigma/2−chi(1−sigma)=47/200000>0.                  (23)

Eventually tau<=eps_T k. By (20)–(23), the final row lower bound
minus this SAME global tau is at least b k, both ambient degrees
are at most Mk, and the mass is at least mL n.

## 8. Real finite histories and genuine missing cap mass

This section retains the original largest-eligible-radius selection and
its real disjoint increments I_i=Q_i minus previous cubes. Use the
algorithmic dilation mathfrak a>=8n^4, which is independent of the
relative-profile parameters. At a positive load support point normalize
the last minimum positive-contribution cube radius to one and reflect
its center into the positive orthant. Let

    alpha=mathfrak a/(mathfrak a+1),
    ccap=(mathfrak a beta−1)/(mathfrak a+1),
    u=last center/(mathfrak a+1).

The actual uncaptured-increment measure, restricted to the limiting
cap D, is dominated by lambda_a=alpha^-n lambda there. The inclusion
of D in the finite cap follows from

    alpha−ccap m−(1−beta m)
       =[(beta+1)m−1]/(mathfrak a+1)>0.

The total finite cap deficit is at least this restricted missing mass.
For an actual target item j meeting A, its uncaptured-center inequality
gives ccap r_j<=||u−w||<alpha+rhi<3; eventually ccap>1, so r_j<3.
For an actual source point z=y/(mathfrak a+1), the common load support
0 in mathfrak a Q_i and y in Q_i give ||z||<r_i. Thus r_i>=rho>4,
so this source precedes the same target j in the original history.

For every mathfrak a>=8n^4 the displacement guard holds:

    alpha beta−kappa−2/(mathfrak a+1)
       =zeta−(beta+2)/(mathfrak a+1)>0.

If K(w,z)>0, the scaled actual centers differ by less than
kappa rho+(r_i+r_j)/(mathfrak a+1), which is less than
alpha beta r_i. Hence the source item really captures that target
center at its own eligibility test. Integrating included increment
pieces costs no more than the original item weight. Both new weights
a_rad,b_soft are at most one, so the actual weighted target row is
at most a_rad(w)tau<=tau on actual-mu_a-almost every target point.
No distribution for actual mu_a is assumed.

Let E be the limiting ambient energy of the final K, Dmax=Mk,
and ca=alpha^-n>=1. The difference identity for nu_a=lambda_a−mu_a
and disjoint A,P gives

    E_lambda_a,lambda_a−E_mu_a,mu_a
       =E_nu_a,lambda_a+E_mu_a,nu_a
       <=ca Mk [nu_a(A)+nu_a(P)].

Actual energy is at most tau ca lambda(A), while ambient energy is
ca² E. Therefore the total finite deficit satisfies

    Delta_a >=[ca² E−tau ca lambda(A)]/[ca Mk]
            >=[E−tau lambda(A)]/[Mk]
            >=lambda(A)b/M>=Gamma n>gamma n.          (24)

This is uniform in the algorithmic dilation for each sufficiently
large fixed n. It uses the same actual history and actual eligibility
throughout. Independent worst positions are not being added.

## 9. General L1 inputs, all positive radii, and full scalar cost

Insert (24) into the original finite cover/load/witness formula, with
its original captured budget tau and last self weight at most
mathfrak a^-n. It gives

    C_n <=[1+(beta+1/mathfrak a)^n/tau]
      [tau+mathfrak a^-n+alpha^-n J_mathfrak a(t)−gamma n]. (25)

For arbitrary L1 input use its absolute value. Every point of a strict
maximal level set has a true positive witness radius. Its open small
cube at radius R_x/mathfrak a covers the point; finitely many cover
each compact subset. All algorithmic dilates remain the original true
witnesses. The above normalization applies to every positive terminal
radius. Inner regularity then gives the entire strict level set.
Neither radius truncation nor an independent witness family is used.

For fixed n let mathfrak a tend to infinity only in the uniform NUMERIC
bounds (25). No selectors or histories need converge. At the chosen tau,

    C_n <=(1+A_n/tau)(tau+B_n−gamma n)
         =F(A_n,B_n−gamma n).                         (26)

The original scalar estimates also give

    B_n<=T−n log(1−sigma)+1
        <=T+[sigma/(1−sigma)]n+1,
    A_n+2sqrt(A_n B_n)
        <=n e^-sigma L+2n sqrt L e^-sigma L/2=o(n).

Equations (2) and (26) now bound the full root expression by

    T−(207689/917000)n
        +n[e^-L+e^-sigma L+2sqrt L e^-sigma L/2].

The bracket tends to zero and is eventually less than 341/229250.
Thus (1) follows, conditional on independent acceptance of this complete
mechanism. The endpoint coefficient 207689/917000 is not asserted.

## 10. Uniform eventual scope and actual evidence

All fixed tuples and target tolerances are independent of n. The guards
have one common threshold: h tends to zero, epsilon qmax/h=O(n^-2),
t tends to infinity, baseline log means/variances tend to zero, the new
target-profile variance is O(L²/n), count failures are O(1/n), perface
retention failure is O(1/x), and cap/tau divided by k are polynomials
in L times exponentials with strictly negative rates. The scalar root
remainder tends to zero. No explicit threshold from R056, R061 or R064
is transferred to this different tuple.

`verify.py` and `verification.json` record an actual new exact Fraction
run of 16 fixed gates, including the Fubini mean constant, target mass,
perface source mean and variance, positive rate gaps, energy budget and
complete scalar slack. These checks are neither a geometry proof nor a
search. They use no old passed grid, MATLAB, shared production run or
lower-bound computation. Only this owned round066 directory is edited.

The relative J input (8)–(9) is the accepted R064 coordinate lemma at
`audit/round064/relative_coordinate/LEMMA.md`, frozen SHA256
`c60d48d116863c9c9cac177e19a9fa4158f12fd8dcef0d5239a5ad0a5487d9f1`.
The new first-profile calculation (6)–(11), its target concentration
and sourceface mean interface are written explicitly here for a fresh
independent local audit. The full actual-history and original operator
transfer remain author responsibilities and require whole-chain review.
