# R092 one-sided counted-coordinate cap lemma

Status: independently proved geometric/probability lemma. This file changes
only the count used to pay source-cap deletion. The actual source inverse
law, exact coordinate profiles, weights and complete column are unchanged.
No n-dependent-s operator theorem is asserted here.

## 1. Entire lower strip is captured under a one-sided condition

Let 0<q<1, rho>0, a0>0, |w|<=r and 2rho>r+1. The source
coordinate capture interval is

    C_w={z: |z-w|<rho+a0}.

If w<a0, the ENTIRE strip

    S_q=[-rho,-rho+(1-q))

belongs to C_w. Indeed the lower capture endpoint is
w-rho-a0<-rho. The upper strip endpoint is <=-rho+1, while
w+rho+a0>=-r+rho+a0>-rho+1 by 2rho>r+1. Thus both endpoints
are paid with strict margin. In particular w can be far below -a0;
the symmetric restriction |w|<a0 is unnecessary. The assumptions
rho>4(s+1), r<s from the prior geometry more than suffice.

The original source corner-cap complement at rho=tq is exactly all
z_i>=-rho+(1-q). A captured negative source pin z_j=-rho is
outside this complement whenever q<1. A positive source pin can
omit at most one counted free coordinate.

## 2. Pre-cap inverse-law deletion bound

Keep the exact original profiles

    J_i=(2rho)^-1 integral_capture 1/ell_R(z) dz,
    ell_R(z)=min(1,(R+rho+a0-|z|)/(2R)),

and assume J_i<=1+e on the actual matching domain. For each counted
coordinate w_i<a0 the free normalized inverse-capture probability
of S_q is at least

    (1-q)/[2rho(1+e)],

since the whole strip is captured and 1/ell_R>=1. If there are
M_-(w,q)=#{i:w_i<a0(q)} counted coordinates, independence of
the FREE inverse coordinates before cap/soft conditioning gives

    P_inv(D^c | captured face)
       <=exp[-(M_--1)(1-q)/(2rho(1+e))].              (1)

For a negative pin the left side is zero. If M_->=nu_one*n and
q<=qmax<1, then with x=n/t and eta=(1-qmax)/(2qmax),

    P_inv(D^c | face)
       <=exp[-(nu_one-1/n)*eta*x/(1+e)].              (2)

This is a law only of the pre-cap inverse coordinate product. It is
not assigned to actual algorithmic increments and is not asserted
after conditioning.

Multiplication by bsoft<=theta/G, summation of all source face
coefficients and the original matched radial integration yield the
same weighted loss bound as before, with nu_one replacing nu:

    caploss(w)/k
       <=M_degree exp[-(nu_one-1/n)*eta*x/(1+e)].      (3)

The theta factors cancel exactly. No new row or column factor is
introduced. At fixed x, (3) need not tend to zero; this fact matters
in the separate critical-scaling diagnostic.

## 3. One additional full-target event supplies nu_one near 1/2

On a full signed target face of radius r, assume 0<a0(q_-)<r.
The n-1 free coordinates are Uniform[-r,r]. The one-sided excess
indicators 1_{w_i<a0(q_-)} are independent Bernoulli

    p_-=(r+a0(q_-))/(2r)>1/2.

The true target pin contributes one on a negative face and zero on a
positive face. Write this deterministic value as b_pin in {0,1} and
let E_-=b_pin+(n-1)p_-. Then E_->=(n-1)/2 and the free count
variance is <=(n-1)/4.

Add precisely one sixth target event

    |M_-(w,q_-)-E_-|<=e*n/4.                          (4)

Chebyshev bounds its failure by 4/(e^2 n), uniformly in r and
face. Because a0(q) increases, M_-(w,q) is nondecreasing in q.
For every matching q, (4) therefore implies

    M_-(w,q)>=n/2-1/2-e*n/4
             >=(1/2-e/2)n=:nu_one*n                 (5)

when n>=2/e. The number nu_one=(1-e)/2 stays near 1/2 independently
of s, unlike the old symmetric 1/s count. The pin omission in (1)
still costs at most one coordinate.

If each of the original five target-event failures and this new
failure is eventually <=epsA/8, their union costs <=6epsA/8<epsA.
Thus the same target-mass conclusion
lambda(A)>=(1-epsA)n log(rhi/rlo) is retained, after enlarging only
the threshold. Independence between target events is unnecessary.
For dimension-dependent error budgets one must separately verify
4/(e_n^2 n) and the integrated target failure; this lemma does not
silently assume that verification.

Adding (4) restricts A. It does not alter source capture, ell_R,
J,J1,J2 or their exact integral formulas. The complete column upper
bound only decreases under a restriction of A. Both kernel weights
remain <=1 and the original finite-history test is unchanged.

## 4. Exact strip mass useful for the critical diagnosis

For d_q=1-q, g=R-a0>=0 and rho>1 (as in the original source
band), the inverse integral of the full captured strip is independent
of w whenever w<a0:

    S(R,a0,d_q)=integral_0^d_q
                 1/min(1,(R+a0+y)/(2R)) dy
      =2R log[1+min(d_q,g)/(R+a0)]+(d_q-g)_+.

The exact inverse-law strip probability is S/(2rho J_i), not merely
the coarser lower bound in (1). When R->infinity with a0,d_q bounded,
S->2d_q. This factor two is needed for the critical cap-cost model;
discarding it can obscure cancellation against the scalar cap integral.

The independent bounded run in receipt.json checked both capture
endpoints at two exact rational examples, including w=-r outside
the symmetric count. The analytic inequalities above prove the full
lemma. No old source or paper was changed.
