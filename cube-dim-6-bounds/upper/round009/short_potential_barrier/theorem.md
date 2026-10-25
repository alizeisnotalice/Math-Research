# A positive coarea lower barrier for the body32 short potential

Use V={(1,s2,...,s6):sj=+/-1}, independent Cauchy variables Z_u of scale
t, f_t(z)=t/[pi(t^2+z^2)], and psi(xi)=sum_(u in V)|u.xi|. The joint
pushforward density of sum u Z_u is p_t. Define K_T=integral_0^T p_t and
S2=union_(v!=w) span(v,w). This concerns the alpha=1 model only.

## 1. Uniform quotient Fourier bounds, proved from moments

For independent Rademacher signs let X=sum epsilon_j xi_j. Then
EX^2=||xi||^2 and

    EX^4=3||xi||^4-2 sum xi_j^4 <=(8/3)||xi||^4.

Hölder gives EX^2<=(E|X|)^(2/3)(EX^4)^(1/3), whence

    8sqrt(6)||xi|| <= psi(xi)=32 E|X| <=32||xi||.

The upper inequality is Cauchy-Schwarz. These estimates hold on all of R6,
so they also hold on any subspace with its orthogonal Euclidean metric.

Fix two body directions v,w and put U=span(v,w), W=U^perp, Q=P_W.
Let Y=sum_(u!=v,w)u Z_u. On W, v.xi=w.xi=0, so the symbol of QY is
exactly psi restricted to W. Thus its four-dimensional density q_t obeys

    q_t(y) >= a t^-4-b||y||^2 t^-6,
    a=3/(4 pi^2 32^4),
    b=15/[8 pi^2 (8sqrt(6))^6],
    b/a=5/108.

This follows from cos(s)>=1-s^2/2: the positive Fourier integral uses
psi<=32||xi||, and the negative second moment uses psi>=8sqrt(6)||xi||.
Here |S3|=2pi^2, Fourier normalization is (2pi)^-4, and the spherical
second moment is ||y||^2/4. All fixed-time integrals converge absolutely.

Deleting one additional variable j gives, on W,

    psi_rest(xi)=psi(xi)-|u_j.xi|
       >=7sqrt(6)||xi||.

Therefore every such 29-variable quotient has a bounded continuous density
q_rest with

    sup q_rest <= C_rest t^-4,
    C_rest=3/[4 pi^2 (7sqrt(6))^4].

This also proves the required quotient spanning; no rank assumption is
silently added. The 30-variable Y itself has full-six-dimensional symbol
at least 6sqrt(6)||xi||, so its full density exists as well.

## 2. Truncating the plane component by positive density domination

For a fixed physical L>0 retain only ||P_UY||<=L. If this fails, then some
of the 30 variables satisfies

    |Z_j|>R=L/(30sqrt(6)),

because ||P_Uu||<=sqrt(6) and the triangle inequality applies before any
projection to W. The Cauchy tail satisfies

    P(|Z_j|>R)=(2/pi)atan(t/R)<=2t/(pi R).

The positive pushforward submeasure for this event has W-density

    integral_(|z|>R) f_t(z) q_rest(y-Q u_j z) dz
       <= C_rest t^-4 *2t/(pi R).

Summing the 30 positive submeasures gives an upper bound for the removed
W-density, not a conditional probability at a fiber:

    q_removed(y)<= (C_tail/L)t^-3,
    C_tail=(1800sqrt(6)/pi)C_rest.

All equalities and inequalities at this stage may be read as identities of
positive pushforward measures and their Lebesgue densities. Tonelli makes
the convolution formula valid. In particular no event of probability zero
has been used as a conditioning event.

Put kappa=sqrt(5/27)=2sqrt(b/a) and

    tau=min(T,aL/(4C_tail)).

When kappa||y||<=t<=tau, the Fourier second term is at most (a/4)t^-4,
and the removed density is also at most (a/4)t^-4. The retained W-density
therefore satisfies

    q_retained(y)>= (a/2)t^-4.                         (1)

## 3. Pair parameters and the exact coarea Jacobian

After a signed-coordinate symmetry take v=(1,...,1), and let w differ in
k of the last five coordinates, k=1,...,5. Then

    Gram(v,w)=[[6,6-2k],[6-2k,6]],
    D_k=det Gram=4k(6-k),
    (Gram^-1)_(11)=(Gram^-1)_(22)=6/D_k.

The pair map (z1,z2)->v z1+w z2 has two-dimensional Jacobian sqrt(D_k).
It is not one. In orthogonal U and W coordinates, positive coarea yields

    p_t(x)=(1/sqrt(D_k)) integral_U
              f_t(z1(x_U-u)) f_t(z2(x_U-u)) r_t(u,x_W) du,            (2)

initially as an almost-everywhere density identity. Here r_t is the joint
30-variable Y density. Formula (2) follows directly by changing the pair
variables in a nonnegative test-function integral and applying Tonelli.
It can equivalently be derived by solving four independent rest columns
in W: the full six-column determinant factors as sqrt(D_k) times the
orthogonal W determinant, leaving exactly the displayed Jacobian.

For ||x||<=sqrt(6) and the retained ||u||<=L, the pseudoinverse row norm gives

    |z_i(x_U-u)| <= B_k(L):=sqrt(6/D_k)(sqrt(6)+L),   i=1,2.

Thus for t<=tau the product in (2) is at least
t^2/[pi^2(tau^2+B_k(L)^2)^2]. Combining with (1) proves

    p_t(x)>=d_k(L)t^-2,
    d_k(L)=a/[2 pi^2 sqrt(D_k)(tau^2+B_k(L)^2)^2],                  (3)

whenever ||x||<=sqrt(6), kappa dist(x,U)<=t<=tau.

The positive-measure argument establishes (3) almost everywhere on this
closed intersection of a Euclidean ball and a tube. For every fixed t>0,
p_t is continuous because its Fourier multiplier is integrable. The
intersection is the closure of its interior (scale toward zero, then take
small interior perturbations). Consequently (3) holds at every point of
the intersection. This fixed-time continuity step includes the plane
itself; it is not an extension of G across its singular set.

## 4. Integration, all 496 planes, intersections and S1

Let d(L)=min_(k=1,...,5)d_k(L). For rho=dist(x,U)>0 with kappa rho<tau,

    K_T(x)>=d(L)[1/(kappa rho)-1/tau].                            (4)

For the present T=4/125, the minimum is exactly at D=20 (k=1 or5).
Indeed with A=sqrt(6)+L the coefficient as a function of D is proportional
to D^(3/2)/(tau^2 D+6A^2)^2. Its logarithmic derivative is positive when
tau^2 D<18A^2. Here tau^2 D<=36(4/125)^2<108<=18A^2, so this holds on
the entire interval 20<=D<=36. Thus only one scalar d(L) is needed per L.

For rho=0, integrate (3) from any positive lower cutoff to tau and send
the cutoff to zero. This proves K_T(x)=+infinity on the plane, including
all S1 lines and every plane intersection inside the cube.

Every pair is carried to one of the five canonical types by an orthogonal
signed permutation. For a fixed direction there are binom(5,k) partners,
so the unordered pair counts are 16 binom(5,k): 80,160,160,80,16, summing
to496. There is no omitted pair type. With

    delta(x)=dist(x,S2)=min_(v!=w)dist(x,span(v,w)),

select a closest plane in this finite union. The uniform bound becomes

    K_T(x)>= (d(L)/kappa)[1/delta(x)-kappa/tau]                 (5)

when delta>0 and the right side is positive, and is +infinity at delta=0.
Overlaps never require adding plane contributions or a factor496. This
is a lower bound; the closest-plane argument only selects one valid bound.

For theta=T/4096, define

    rho_star(L)=1/[kappa(1/tau+theta/d(L))].

Every rational 0<delta0<rho_star(L) gives the strict certificate

    K_T(x)>theta whenever ||x||_infinity<=1 and dist(x,S2)<=delta0.

The remaining set {x in [-1,1]^6:dist(x,S2)>=delta0} is compact and avoids
S2. It still requires a separate lower certificate. No global m_T or C6
improvement follows from the singular-neighborhood certificate alone.

## 5. A simple verified rational corollary at T=4/125

For L=3/4, the scalar certificate gives

    d/kappa >1/5000000000,
    kappa/tau <40000.

Consequently, for every cube point off S2,

    K_(4/125)(x) >= (1/5000000000)
                    *(1/dist(x,S2)-40000)_+.                    (6)

When the displayed positive part is nonzero, its distance is less than
1/40000<tau/kappa, so (5) applies with both rounded constants in the correct
direction. Outside that range, (6) is just positivity. On S2, K=+infinity.
In particular,

    dist(x,S2)<=1/80000  =>  K_(4/125)(x)>=1/125000>1/128000.

Both constant comparisons are evaluated with outward Arb384 enclosures and
independently reconstructed from the exact formulas in this proof. They
are well separated from equality; no serialized wide-ball status is used.
