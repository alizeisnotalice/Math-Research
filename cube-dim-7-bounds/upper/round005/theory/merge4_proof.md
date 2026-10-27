# R005: four-block merging with a two-moment Radau correction

2026-10-27. Sol proof only; this owner runs no numerical script, quadrature, MC or MATLAB. The fixed original body64/common-time Cauchy model and disjoint Hadamard partition are unchanged. Terra owns the separate numerical certificate. No numerical improvement, g, or C7 result is asserted here before that certificate.

## 1. Notation and unchanged scalar interface

For one7x8 Hadamard block put T=det(Sigma_block)^(1/7), and let T_i be iid copies. The R002 exact moment, valid for Re z>0, is

    m(z)=E T^(-z)
        =2^(-4z/7)Gamma(4+z)/Gamma(4+8z/7)
                       [Gamma(1/2+z/7)/Gamma(1/2)]^8.              (1)

Put U=T1+T2 and q(r)=E U^(-r), r>=0. Thus **q(0)=1 exactly** by probability normalisation. It is not evaluated by the Mellin integral at r=0, where Gamma(r) has a pole.

For r>0, R004 established

    q(r)=1/[pi Gamma(r)] integral_0^infinity f_r(t)dt,
    f_r(t)=Gamma(a_r+it)Gamma(a_r-it)m(a_r+it)m(a_r-it)>=0,
    a_r=r/2.                                                       (2)

Its entire derivation applies to every positive r, not just1/8. The already certified **complete** q(1/8) interval is reused from R004, with its Cauchy error, full tail and exact dyadic endpoints; it must not be re-evaluated in R005.

Let V=U1+U2=T1+...+T4, and q_four(s)=E V^(-s). Determinant Minkowski followed by AM–GM only between the two independent four-block sums gives

    M=E det(Sigma)^(-1/14)
       <=2^(-1/2) q_four(1/4)^2.                                  (3)

The original complete-G boundary interface is

    I<=F M^(7/6), F=1/[105·2^(25/6)pi^(13/12)].                     (4)

This is the complete time/full-boundary budget. An improved I expression alone supplies no global g or C7 result.

## 2. Tilted log-ratio: the exact four-block identity

Fix

    s=1/4, a=s/2=1/8, tau=1/16,
    dP_a(u)=u^(-a)dP_U(u)/q(a).

For independent U1,U2 under this tilted probability law, put

    Delta=log U1-log U2,
    W=cosh(tau Delta)-1>=0.

The elementary identity

    (U1+U2)^(-s)=2^(-s)(U1U2)^(-a)
                               cosh(Delta/2)^(-s)

implies

    q_four(s)=2^(-s)q(a)^2 E_a[g(W)],
    g(w)=[T8(1+w)]^(-1/4), w>=0.                                  (5)

Here T8 is the ordinary Chebyshev polynomial of the first kind. Because1/(2tau)=8,

    cosh(Delta/2)=cosh(8 tau Delta)=T8(cosh(tau Delta))=T8(1+W).

There is no missing2 in this substitution.

The two moments of W use only **real** pair moments:

    R1=E_a cosh(tau Delta)
       =q(a-tau)q(a+tau)/q(a)^2
       =q(1/16)q(3/16)/q(1/8)^2,
    R2=E_a cosh(2tau Delta)
       =q(a-2tau)q(a+2tau)/q(a)^2
       =q(0)q(1/4)/q(1/8)^2=q(1/4)/q(1/8)^2,
    mu=E_a W=R1-1,
    nu=E_a W²=R2/2-2R1+3/2.                                       (6)

The cosh² identity gives the coefficients in nu. All moments are finite: q(r) is finite for every r>0 by pair AM–GM and(1); the only boundary moment is q0=1.

U has a nondegenerate continuous law, so log U does too. Under the positive tilt this remains true; therefore mu>0 and nu>=mu²>0. Equivalently strict Holder gives R1>1. These strict analytical statements do not replace the need to separate their finite intervals in the numerical certificate.

## 3. Complete monotonicity and the global Hermite majorant

The roots of T8 are cos((2j-1)pi/16), j=1,...,8. Thus

    T8(1+w)=2^7 product_(j=1)^8(w+c_j),
    c_j=1-cos((2j-1)pi/16)>0.

Every factor(w+c_j)^(-1/4) is completely monotone on[0,infinity): its kth derivative has sign(-1)^k by the rising-factorial formula. Products preserve complete monotonicity by Leibniz' rule. Hence g>0, g'<0, g''>0, **g'''<0**, and g(0)=1.

Set d=nu/mu>0. Let P be the unique quadratic Hermite polynomial matching

    P(0)=g(0), P(d)=g(d), P'(d)=g'(d).

For any w>=0 other than the interpolation nodes, the standard repeated-node remainder gives

    g(w)-P(w)=g'''(xi) w(w-d)^2/6<=0,

for some xi in the interval spanned by0,d,w. At0,d equality holds. Thus P is a global upper majorant of g on the entire nonnegative half-line, not only between the two nodes.

The two-point measure with mass1-lambda at0 and lambda at d, where lambda=mu²/nu, matches the first two moments of W: lambda d=mu, lambda d²=nu. Since P is quadratic, its expectation depends only on those moments. Therefore

    E_a g(W)<=E_a P(W)
       =(1-lambda)g(0)+lambda g(d)
       =1-delta(mu,nu),
    delta(mu,nu)=mu²/nu ·[1-g(nu/mu)]>0.                            (7)

Also0<delta<1 because0<lambda<=1 and0<g(d)<1. This is the Radau two-moment correction. It needs no complex-parameter moment q(z), no multidimensional quadrature and no empirical variance assumption.

Combining(3),(5),(7) proves the new rigorous expression

    M<=M_Radau:=[q(1/8)^4/2]·[1-delta(mu,nu)]²,                    (8)
    I<=I_Radau:=F M_Radau^(7/6).

Relative to the exact R004 pair expression, I_Radau=I_pair[1-delta]^(7/3). The correction is analytically strictly positive. Its magnitude and a new rational budget endpoint have not been computed by this owner.

## 4. Monotonicity and safe interval assembly

Let k(d)=[1-g(d)]/d for d>0. Since g is decreasing convex and g(0)=1, 1-g is positive, concave, zero at0. Its secant slope k is positive nonincreasing. Thus

    delta(mu,nu)=mu·k(nu/mu)

is increasing in mu>0 at fixed nu, and nonincreasing in nu>0 at fixed mu. One may also differentiate: g convex gives1-g(d)>=d[-g'(d)]. Then

    partial_mu delta=(mu/nu)[2(1-g(d))+d g'(d)]>0,
    partial_nu delta=-(mu²/nu²)[(1-g(d))+d g'(d)]<=0.

These derivatives independently confirm the directions. Consequently using a lower mu and an upper nu is safe.

Suppose the complete certified pair-moment intervals have positive endpoints

    q(1/8) in[A-,A+] from R004,
    q(1/16) in[B-,B+], q(3/16) in[C-,C+], q(1/4) in[D-,D+].

All divisions below require A->0. Compute outward endpoints

    R1_L=B- C-/(A+)²,
    R2_U=D+/(A-)²,
    mu_L=R1_L-1, rounded downward,
    nu_U=R2_U/2-2R1_L+3/2, rounded upward.                          (9)

Then true mu>=mu_L and true nu<=nu_U. The numerical acceptance gate is mu_L>0, nu_U finite positive. Since true nu>=true mu², the correct enclosures necessarily have nu_U>=mu_L²; verify this with exact endpoint/rational arithmetic. A violated inequality is an inconsistency, not a reason to clip nu.

Compute a downward bound delta_L for

    delta(mu_L,nu_U)
      =mu_L²/nu_U ·[1-T8(1+nu_U/mu_L)^(-1/4)].                    (10)

Here d=nu_U/mu_L>0. The expression lies strictly between0 and1. Require0<delta_L<1; otherwise report unresolved separation rather than publish a midpoint. These conditions ensure the factor1-delta_L is positive and every subsequent power is legitimate.

Safe final upper expressions are

    M_upper=(A+)^4/2·(1-delta_L)²,
    I_upper=F M_upper^(7/6),                                      (11)

with all whole-expression endpoint operations outward. To certify quantitative improvement against R004's exact pair expression, compare(11)'s upper endpoint with its complete R004 **lower budget-expression endpoint**, e.g.(A-)^4/2 for M. A lower endpoint for a deterministic budget expression is not a lower bound on actual M or I.

For stable evaluation of T8, use T2(x)=2x²-1, T4=2T2²-1, T8=2T4²-1 for x>=1. To avoid cancellation in1-T8^(-1/4), with R=T8^(1/4) use

    1-T8^(-1/4)=(T8-1)/[R(R³+R²+R+1)].

The nonnegative T8-1 can itself be built from T2-1=2(x-1)(x+1), T4-1=2(T2-1)(T2+1), T8-1=2(T4-1)(T4+1). Positive interval arithmetic then makes the intended sign explicit.

## 5. Three new positive1D certificates: complete constants

For each new r>0 let a_r=r/2. Use exactly the R004 analytic integrand(2), with the following parameter-specific choices:

|pair exponent r|a_r|panel halfwidth h_r|Cauchy radius rho_r|panels on[0,32]|a_r-rho_r|a_r+rho_r|
|---|---|---|---|---:|---|---|
|1/16|1/32|1/64|3/128|1024|1/128|7/128|
|3/16|3/32|1/32|3/64|512|3/64|9/64|
|1/4|1/8|1/16|3/32|256|1/32|7/32|

Each has h_r/rho_r=2/3 and rho_r<a_r, so all Cauchy disks lie in |Im t|<a_r. The exact panel centre is(2j+1)h_r, j=0,...,n_r-1. The strip bound follows exactly as in R004 from Gamma(v)m(v) log-convexity:

    M_star(r)=Gamma(a_r-rho_r)Gamma(a_r+rho_r)
                         m(a_r-rho_r)m(a_r+rho_r).

Degree128, series order129,512-bit outward complex arithmetic give the **raw full-segment** error

    E_raw(r)=32 M_star(r)(2/3)^129/[1-2/3].                         (12)

It is32 times the envelope because all three intervals have length32. It is not n_r times the R004 error and is not normalised until multiplication by1/[pi Gamma(r)].

For Gamma recurrence tail order16,

    E_tail_q(r)=m(a_r)^2Gamma(a_r+16)^2/
                            [pi Gamma(r)·31·32^31].               (13)

This is the full normalised q tail t>=32. The upper Gamma argument a_r+16 changes with r and must not be copied from R004's257/16. Integrate the even Taylor coefficients on each panel with2h_r^(2k+1)/(2k+1). Reconstruct

    q_lower=max(0,[sum Re(P_j)-E_raw(r)]/[pi Gamma(r)]),
    q_upper=[sum Re(P_j)+E_raw(r)]/[pi Gamma(r)]+E_tail_q(r),         (14)

using the proper whole-expression outward endpoints. The polynomial may have negative coefficients; none is dropped. Actual coefficients are real and the enclosing imaginary intervals should include0.

There are1792 new panels total. R004 q(1/8) is reused, not recomputed. The interface is fully1D and uses the same already exercised Gamma-series kernel; this has substantially lower implementation complexity than nesting a complex q inside an outer Taylor series.

## 6. Frozen finite acceptance order

The separate `registration.md` gives the proposed computational contract. Briefly: probe only the already known parameterisation, lock sources/config, run four panels for each r as a12-panel pilot, and reuse them. If all pass,15 single-thread workers compute the remaining1780 panels. Record1792 unique(r,index) keys, all terminal intervals, CPU/wall/RSS and actual occupancy; do not infer full occupancy merely from requesting15 workers.

Assemble the three complete q intervals, load the exact complete R004 q interval, apply(9)-(14), certify0<delta_L<1 and strict improvement, and round a new I upper endpoint outward to a declared rational grid. Stop at this finite batch/main review. If interval separation fails, retain the analytic result and report the unresolved numerical comparison; no automatic finer grid, higher degree, extra exponent or MC is authorised.

This proof establishes the correction and safe error/endpoint directions. It does not execute the proposed batch or claim its numerical output. A whole-face g certificate remains independent work.
