# Self-contained radial Fourier lower bounds for alpha=3/2 body32

Let V be the 32 unoriented cube-body representatives in R6, and set

    psi(xi)=sum_(v in V)|v.xi|^(3/2),
    p_t(x)=(2pi)^-6 integral exp(-t psi(xi)) cos(xi.x) dxi,
    G(x)=integral_0^infinity p_t(x) dt.

This is the density of the sum of the 32 independent symmetric alpha-stable
variables along V; in particular p_t>=0. The stable exponent may equivalently
be constructed without assuming a density comparison: put

    K=[integral_R (1-cos u)|u|^(-5/2)du]^-1 >0,
    |s|^(3/2)=K integral_R (1-cos(us))|u|^(-5/2)du.

The defining integral converges at zero and infinity. Truncating the jump
variable away from zero and infinity gives symmetric compound-Poisson
probability measures for each of the 32 directions. Their characteristic
functions converge to exp(-t psi). Against any nonnegative Schwartz test
function, Fourier inversion and domination by its integrable Fourier
transform show that the limiting distribution is nonnegative. The continuous
inverse Fourier function below therefore satisfies p_t>=0 everywhere.
The coercivity below makes the characteristic function integrable for every
t>0, so the inverse Fourier integral is an absolutely convergent continuous
density. No pointwise lower comparison is inferred merely from ordering two
Fourier multipliers; every negative cosine term is paid separately.

## Spectral envelopes and exact angular mean

For X=sum_(j=1)^6 epsilon_j xi_j, Rademacher symmetry gives
psi=32 E|X|^(3/2). The second and fourth moments are

    EX^2=||xi||^2,
    EX^4=3||xi||^4-2 sum xi_j^4 <= (8/3)||xi||^4.

Hölder with exponents 5/4 and 5 implies
EX^2 <= (E|X|^(3/2))^(4/5)(EX^4)^(1/5), and thus

    c6 ||xi||^(3/2) <= psi(xi) <= 32 ||xi||^(3/2),
    c6=32(8/3)^(-1/4).

The upper inequality is the probability-space L^(3/2)-L2 inequality.
At xi=0 both inequalities are trivial.

There is a useful improvement for the positive zeroth-order integral. If
omega is uniform on S5, direct Gaussian/polar integration gives

    E_omega |v.omega|^(3/2)
      =6^(3/4)*2 Gamma(5/4)/(sqrt(pi) Gamma(15/4)),
    mu=E_omega psi(omega)
      =64*6^(3/4) Gamma(5/4)/(sqrt(pi) Gamma(15/4)).

For example, integrate exp(-||z||^2)|v.z|^(3/2) first in coordinates aligned
with v, and then in polar coordinates; the two radial integrals give exactly
this identity. Convexity of exp(-s) yields

    E_omega exp[-t r^(3/2) psi(omega)]
       >= exp[-mu t r^(3/2)].

This is an integral lower bound for p_t(0), not a claim psi<=mu pointwise.
Arb verifies c6<mu<32.

## The degree-2 bound and exact optimal tail

The elementary inequality cos(s)>=1-s^2/2, the two spectral envelopes, and
the angular identity E(omega.x)^2=||x||^2/6 give, for every x and t>0,

    p_t(x) >= a t^-4 - b ||x||^2 t^(-16/3),
    a=Gamma(4)/(96 pi^3 32^4),
    b=Gamma(16/3)/(1152 pi^3 c6^(16/3)).

The constant 1/(96 pi^3) includes |S5|=pi^3 and alpha=3/2. Replacing the
zeroth-order integral by Jensen is valid with
a_J=Gamma(4)/(96 pi^3 mu^4); b is unchanged.

On the whole cube [-1,1]^6, ||x||^2<=6. Since p_t is nonnegative,

    G(x) >= integral_T^infinity p_t(x)dt
         >= a/(3T^3) - (18b/13) T^(-13/3).

The right side has its unique positive maximum at

    T_star=(6b/a)^(3/4),
    g_2=(4a/39)*(a/(6b))^(9/4).

The same formulas with a_J give the Jensen degree-2 certificate. This is
an actual all-cube lower bound, not a point value, a Monte Carlo minimum,
or a unit-sphere substitute for the cube radius.

## Two finite higher-order variants

For K=1,3,5, define

    P_(2K)(s)=sum_(k=0)^K (-1)^k s^(2k)/(2k)!.

The degrees 2,6,10 are global lower bounds for cos. One proof starts with
cos>=1-s^2/2 and inducts through alternating upper/lower polynomials:
P_(2K)''=-P_(2K-2), and integrating the sign of the second derivative from
zero twice gives the next bound. Thus no small-frequency cutoff or
uncontrolled Taylor remainder is used.

For a unit vector n,

    E_omega (omega.n)^(2k)/(2k)!
      =1/[4^k k! (3)_k].

Let d_0=32 for the basic scheme and d_0=mu for the three Jensen schemes;
for positive even k>=2 put d_k=32, and for negative odd k put d_k=c6. Set

    A_k=Gamma(4+4k/3)/[96 pi^3 4^k k! (3)_k d_k^(4+4k/3)].

Separate positive and negative terms before applying the radial envelopes.
This proves the simultaneous all-space/all-positive-time lower inequality

    p_t(x) >= sum_(k=0)^K (-1)^k A_k ||x||^(2k) t^(-4-4k/3).

All Fourier integrals used here are absolutely convergent because c6>0.
All finite polynomial-tail integrals over [T,infinity) converge. For r=||x||>0,
choose T=r^(3/2) z^(-3/4), z>0. Positive integration gives

    G(x) >= r^(-9/2) F_K(z),
    F_K(z)=z^(9/4) sum_(k=0)^K (-1)^k A_k z^k/(3+4k/3).

Whenever F_K(z)>0, r<=sqrt(6) yields the all-cube lower bound

    g_K >= 6^(-9/4) F_K(z).

The lower bound is legitimate even where G=+infinity. At x=0, scaling gives
p_t(0)=t^-4 p_1(0), so G(0)=+infinity. In particular the statement holds
pointwise as an extended-valued inequality and hence for the cube essential
infimum. No singular layer or spatial boundary is omitted.

Tail optimization reduces exactly to a polynomial root because

    F_K'(z)=(3/4)z^(5/4) Q_K(z),
    Q_K(z)=sum_(k=0)^K (-1)^k A_k z^k.

For K=1 this recovers the closed degree-2 formula. The script verifies
A_2^2<3 A_1 A_3, implying Q_3'<0 everywhere. For K=5 it verifies
A_2^2<A_1 A_3 and 2 A_4^2<5 A_3 A_5; then

    Q_5'=(-A_1+2A_2z-A_3z^2)
            - z^2(2A_3-4A_4z+5A_5z^2)<0.

Thus each polynomial has exactly one positive root: it is positive at zero
and has negative leading coefficient. F_K increases to that root and then
decreases. Strict Arb brackets enclose the unique global optimum. A rational
point in the bracket separately supplies the constructive lower certificate;
the root enclosure also bounds the best value available to that scheme.

## What a noncompetitive budget means

For the same alpha=3/2 model the existing operator interface has cost
4^(4/3)I/g^(4/3)=(4B/g)^(4/3), where B=||G(1,.)||_(4/3).
With the fixed F008 norm-envelope value U, paying 64 requires the allocated
lower bound to exceed U/(4sqrt(2)). The script compares the optimal value
of each of these four lower-bound schemes with that threshold.

If the optimal scheme value is below the threshold, this excludes closure
by these specific lower certificates paired with that fixed norm budget.
It does not upper-bound the actual essential g or lower-bound the actual
operator/model cost, and it does not exclude sharper estimates for either
side. No alpha=1 g certificate is used anywhere in this construction.
