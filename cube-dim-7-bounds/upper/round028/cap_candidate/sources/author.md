# R027: cube-specific cap bound for the weighted selection load

Date: 2026-10-27. Native Sol, pure theory; zero mathematical numerical runs, MATLAB, parameter scans, commits or pushes. Only this assigned new directory is authored. The parent supplied the mechanism; this proof independently reconstructs every geometric inequality and is frozen for fresh audit. No approximate scalar bound is certified and no ledger is changed.

## 1. Statement

Let n>=1 be an integer and let M_n be the centered axis-parallel cube maximal operator on Lebesgue R^n, using all positive common radii and general L1 input. Fix t>1, tau>0 and a>2t. Put

    b=a(1+1/t),       alpha=a/(a+1),
    c=(b-1)/(a+1)>1, P=((a+1)/a)^n,
    A_a=(1+1/t+1/a)^n.

Define

    J_a(t)=1-((1+c-alpha)/2)^n+n log t
           -(n/2) integral_1^t
               ((1+c-alpha/rho)/2)^(n-1) d rho/rho,
    E_a(t)=a^(-n)+P J_a(t).

Then

    C_n <= (1+A_a/tau)(tau+E_a(t)).                         (1)

In particular define the positive quantity

    B_cube(t)=1-((1+1/t)/2)^n+n log t
              -(n/2) integral_1^t
                   [1+1/(2t)-1/(2rho)]^(n-1) d rho/rho.

Taking only an infimum of already proved scalar bounds as a tends to infinity yields

    C_n <= (1+(1+1/t)^n/tau)(tau+B_cube(t)).                 (2)

Optimizing the positive scalar tau gives

    C_n <= [(1+1/t)^(n/2)+sqrt(B_cube(t))]^2.                (3)

The square root is legitimate: B_cube(t)>0 for t>1, as proved below. All these are operator weak-(1,1) bounds, not costs of an unconnected model integral. Their independent review and ledger adoption remain the parent's responsibility.

## 2. Finite selection and unchanged covering factor

For any finite family of open cubes choose a subfamily recursively. Write B_j=B(y_j,r_j), and at selection time freeze

    I_j=B_j \ union_(i<j) B_i,
    v_j=|I_j|,
    w_j=v_j/(2a r_j)^n.

Select a maximal-radius remaining cube whose center y obeys

    sum_(j already selected) w_j 1_(bB_j)(y)<=tau.           (4)

The small-cube increments I_j are disjoint. Eligibility can only be lost, so selected radii are nonincreasing. For any unselected cube, its center must already have rejection sum >tau from selected cubes with at least its radius; otherwise the procedure could not select a smaller cube or terminate. Enlarging bB_j to (b+1)B_j covers every point of that rejected cube. Integration therefore gives

    |union(all supplied cubes)|
      <= (1+A_a/tau) sum_j v_j.                            (5)

This is the same finite cover proof with the actual new parameters, and requires no change to the source increments.

## 3. Paying captured previous terms and controlling ties

Fix x with positive W(x)=sum_j w_j 1_(aB_j)(x). Among positive contributing cubes of smallest radius, choose B_k with the **last selected index**. Translate x to 0 and normalize r_k=1. Weights remain invariant under this common transformation. Every other positive contributing term has index j<k: larger radii were selected earlier, and the last-index tie rule handles equal radii.

Call such a previous term captured if y_k is in bB_j. All captured terms together are a sub-sum of the actual selection condition (4) at y_k, so they cost at most tau, irrespective of their radii. The self-term is not in (4); it is separately bounded by

    w_k<=|B_k|/|aB_k|=a^(-n).                               (6)

For an uncaptured contributing term,

    ||y_k-y_j||_infinity>=b r_j,
    ||y_j||_infinity<a r_j,
    ||y_k||_infinity<a.

If r_j>=t, the last two inequalities instead imply

    ||y_k-y_j||_infinity<a+a r_j<=b r_j,

contradiction. Thus every unpaid term has 1<=r_j<t. The strict t>1 ensures the earlier-index payment argument never confuses the self-term with a large-radius term.

## 4. From unpaid increments to a geometric cap integral

For y in I_j of an uncaptured contributing term, set

    z=y/(a+1),    u=y_k/(a+1),
    m(z)=max(1,||z||_infinity).

Since y is in B_j and 0 is in aB_j,

    ||z||_infinity<r_j<t,
    r_j>=m(z),
    ||u||_infinity<alpha.

Uncapturedness and the reverse triangle inequality give

    ||u-z||_infinity>c r_j>=c m(z).                         (7)

The strict inequality may be enlarged to >= when defining an upper-bound domain. Also r_j^(-n)<=m(z)^(-n). Disjointness of the small increments and the Jacobian dy=(a+1)^n dz now imply

    sum_(uncaptured contributing j) w_j
       <=P H(u;c,t),

where

    H(u;c,t)=2^(-n) integral_
       {||z||infinity<t, ||u-z||infinity>=c m(z)}
          m(z)^(-n) dz.                                   (8)

Every unpaid increment is inside this domain. The proof does not assume that its original cube or star is disjoint from other cubes; it uses only the actual I_j disjointness. Together with (6), this proves

    W(0)<=tau+a^(-n)+P H(u;c,t).                           (9)

## 5. Exact coordinate monotonicity of H

By reflecting z and u in a coordinate axis, H depends only on the coordinate absolute values of u. It therefore suffices to take 0<=u_i<=alpha<1. We prove that H is nondecreasing in each u_i on this box, at fixed c>1 and t>1.

Inside ||z||infinity<1, m(z)=1. The captured region is the intersection of [-1,1]^n with u+[-c,c]^n. Its ith side length is

    L_i(1)=min(2,1+c-u_i).

Indeed its upper endpoint is 1 because u_i+c>1, and its lower endpoint is max(-1,u_i-c). The side lengths are positive and nonincreasing in u_i. Thus the uncaptured ball contribution, 1-2^(-n) product_i L_i(1), is nondecreasing in every u_i.

For the shell ||z||infinity=rho, 1<=rho<t, write

    L_i(rho)=min(2rho,(1+c)rho-u_i).

These are the captured side lengths in the other coordinates of a face. The positive ith face z_i=+rho is always captured in coordinate i: u_i<1<=rho and c>1 imply |u_i-rho|<c rho. The negative ith face z_i=-rho is captured exactly when

    u_i<(c-1)rho.

Consequently the total captured face measure at that radius is

    S(rho,u)=sum_i N_i(rho,u_i) product_(l!=i) L_l(rho),
    N_i(rho,u_i)=1+1_{u_i<(c-1)rho}.                        (10)

For each coordinate u_j, every summand in (10) is nonincreasing: either it contains the nonincreasing count N_j, or it contains the nonincreasing length L_j. All remaining factors are nonnegative. Therefore S is coordinatewise nonincreasing and the uncaptured shell contribution is coordinatewise nondecreasing.

This is a face-by-face integral identity, **not** differentiation of the volume of a moving captured intersection. Apart from sets where two coordinates tie for maximal absolute value, which have n-dimensional measure zero, the shell is parameterized by z_i=+/-rho and the other n-1 coordinates in [-rho,rho]. The volume Jacobian is exactly 1 in (rho,z_other). Thus

    H(u;c,t)=1-2^(-n) product_i L_i(1)
       +2^(-n) integral_1^t rho^(-n)
          [2n(2rho)^(n-1)-S(rho,u)] d rho.                 (11)

The same formula holds in dimension 1, interpreting an empty product as 1. Equality surfaces of a cap affect at most a null set of radii; the indicated strict captured condition matches the >= uncaptured domain. Formula (11) proves

    H(u;c,t)<=H(alpha*1;c,t).                              (12)

The endpoint u=alpha*1 need not be attainable by an actual chosen center. It is a valid maximum on the closed comparison box.

## 6. Evaluation at the comparison corner

For the stated parameters,

    c-1=(a/t-2)/(a+1)>0,
    alpha/(c-1)=a t/(a-2t)>t.

Hence alpha>(c-1)rho for every 1<=rho<=t. At the corner u=alpha*1 there is exactly one captured face per coordinate, and

    L(rho)=(1+c)rho-alpha,
    ell_a(rho)=L(rho)/(2rho)=(1+c-alpha/rho)/2.

One has 0<ell_a(rho)<1 on this whole interval. The ball term of (11) becomes 1-ell_a(1)^n. The total shell measure gives n d rho/rho after normalization, while the n captured faces give

    (n/2) ell_a(rho)^(n-1) d rho/rho.

Therefore H(alpha*1;c,t)=J_a(t) exactly as defined in Section1. In particular the factor n/2, the normalized volume 2^(-n), and the affine Jacobian P are all accounted for. Equations (9) and (12) prove

    W(x)<=tau+E_a(t)  for every x.                         (13)

Both J_a and its limit are positive: the ball term is positive, and the shell integrand n[1-ell_a^(n-1)/2]/rho is at least n/(2rho). This also covers n=1.

## 7. General L1, all scales, and the scalar limit

For f in L1 and lambda>0, E={M_n f>lambda} is open by L1 translation continuity of every fixed-radius cube average. For compact K contained in E, at each x choose an actual bad radius R_x>0 and supply the small cube B(x,R_x/a). Its a-star is the original centered common-radius bad cube. The small cubes form an open cover; take a finite subcover and apply the selection above. Every selected star still has average >lambda, so

    lambda sum_j v_j<=integral |f| W<=(tau+E_a(t))||f||_1.

Combine with (5) and then use inner regularity to obtain (1) for all of E. There is no radius restriction, finite-source approximation, product assumption, or regularity requirement beyond L1. Open and closed cube averages agree because their boundaries have zero Lebesgue measure.

Now hold n,t,tau fixed and let a tend to infinity through admissible a>2t. Then

    alpha->1, c->1+1/t, P->1, a^(-n)->0,
    A_a->(1+1/t)^n.

The integrands defining J_a converge uniformly on the fixed compact interval [1,t], since they are integer powers of continuous affine expressions in 1/rho with convergent coefficients. The ball terms also converge. Thus J_a tends to B_cube(t). The finite bound is already proved for the **same** operator and for every f,lambda at every finite admissible a. Taking its scalar constant infimum/limit proves (2). Finite subcovers, selections, and increments may change with a and need not converge.

In the limit ell(rho)=1+1/(2t)-1/(2rho) is in (0,1] on [1,t], with ell(t)=1. The ball term is positive for t>1, and again the shell part is at least (n/2)log t. Hence B_cube(t)>0. Optimizing tau>0 in (2) by tau=sqrt((1+1/t)^n B_cube(t)) gives (3).

The new coefficient is strictly smaller than the R025 coefficient at the same t:

    0<B_cube(t)<1+n log t,

because the subtracted ball-cap term and the captured-shell integral are positive. Therefore this is a strict symbolic improvement of that valid operator-bound family, without a numerical claim about its optimum.

## 8. Exact rational-plus-log form in n=7

Put

    q=1+1/(2t),       h=(1+1/t)/2.

For fixed t>1 the sixth-power integral in B_cube is exactly

    I_6(t)= integral_1^t [q-1/(2rho)]^6 d rho/rho
      =q^6 log t
       -3 q^5(1-t^(-1))
       +(15/8) q^4(1-t^(-2))
       -(5/6) q^3(1-t^(-3))
       +(15/64) q^2(1-t^(-4))
       -(3/80) q(1-t^(-5))
       +(1/384)(1-t^(-6)).                                 (14)

This follows term by term from the binomial theorem and integral_1^t rho^(-j-1)d rho=(1-t^(-j))/j, j>=1. Hence

    B_cube(t)=1-h^7+7 log t-(7/2) I_6(t)                    (15)

is exactly a rational function of t plus a rational function times log t. In compact form, its log coefficient is 7-(7/2)q^6 and its remaining part is

    1-h^7-(7/2) sum_(j=1)^6
      binom(6,j) q^(6-j)(-1/2)^j (1-t^(-j))/j.

The log coefficient need not be positive for every t>1; any future scalar certificate must retain its sign. For rational t all coefficients in (14)-(15) are rational. No cancellation was numerically evaluated here.

The finite-a integral has the same exact form with q_a=(1+c)/2 and p_a=alpha/2:

    I_(6,a)=q_a^6 log t+sum_(j=1)^6
       binom(6,j)q_a^(6-j)(-p_a)^j(1-t^(-j))/j,
    J_a=1-(q_a-p_a)^7+7 log t-(7/2)I_(6,a).

Remember the separate self-term a^(-7) in E_a=a^(-7)+P J_a. It is not part of this integral and cannot be omitted at finite a.

## Frozen scope

This branch considered only the proposed cube-cap mechanism. It did not evaluate t=5 or any other scalar, did not optimize t, and did not launch a numerical certificate. The original Stein–Strömberg method is the established baseline; the cube-specific cap deduction above is an internal mathematical result, not represented as a printed public theorem or a novel published record. The parent will coordinate fresh independent audit and any finite scalar task.
