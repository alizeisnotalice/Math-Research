# R092 one critical-scaling diagnostic — cap cancellation obstruction

Status: one conditional asymptotic model, not a new operator bound.
The proved one-sided lemma is in LEMMA.md. This file does not close
the n-dependent-s full chain or declare any C_n improvement. It
identifies a specific cost that the naive I/2 gain must pay.

Fix x>0, Q=qmax in (0,1), lambda>0 and let

    L=log n, t=n/x, s_n~lambda*L/x,
    a=1-log2, a*lambda<1.

One may formally choose theta=n^-kappa with
a*lambda<kappa<1, tau=n^alpha with 0<alpha<1-kappa,
and polynomial matching width h=L^-p. The full scalar H_tau
remainders then vanish: tau/n->0, A_n/n->0 and
A_n B_n/(tau*n)=O(L/n^alpha)->0. Also
tau/k=O(L^(p+1)n^(alpha-1+kappa))->0. These statements do not
by themselves prove the n-dependent geometric construction.

## 1. Exact scalar cap-integral limit

The original integral, with no deletion of its cap term, is

    I_n(t)=integral_1^t
       [1+1/(2t)-1/(2rho)]^(n-1) d rho/rho
      =integral_(x/n)^1
       [1-x(1/q-1)/(2n)]^(n-1) dq/q.

Extend the last integrand by zero below x/n. It converges pointwise
for q>0 to exp[-x(1/q-1)/2]/q. For n>=2 it is bounded by
exp[-x(1/q-1)/4]/q, an integrable function on (0,1). Dominated
convergence gives the exact limit

    I(x)=integral_0^1 exp[-x(1/q-1)/2] dq/q
        =exp(x/2) E1(x/2),
    E1(z)=integral_z^infinity exp(-v) dv/v.            (1)

Since the nonnegative core of B_n is O(1),

    B_n/n=L-log x-I(x)/2+o(1), A_n->exp x.            (2)

If the entire target log-mass log(s_n Q) could be extracted with
zero cap cost, the apparent coefficient would be
log lambda+log Q+I(x)/2. Source cap cannot simply be omitted.

## 2. True variance and retention scale

The local profile value bounds should use R/rho, not s/rho.
From Q(R,u)<=log2<1 and the exact formulas,

    J1<=R/rho<=s_n/t,
    J2<=2R/(3rho)<=2s_n/(3t).

Also |J-1|=O(R/rho)=O(L/n) uniformly on active matching,
using m_w<=r<=R+2h and R>=1. Thus the true centered source-face
variance, with the deterministic pin omitted, obeys

    Var V<=sum_(i!=j)J2_i/J_i
          <=2s_n*x/[3(1-o(1))]=O(s_n*x)=O(L).        (3)

The pin loss is O(1). The mean scale is

    mu(s_n)x=a*lambda*L-x*log2+o(1).

A soft threshold kappa L with kappa>a*lambda has gap of order L.
If the target-profile and transport estimates can be closed at this
scale, true-variance Chebyshev gives soft failure O(1/L). Its
integrated loss over target log-mass O(log L) is O(log L/L)=o(1).
Retention need not itself cost a positive constant coefficient in
this strict-slack regime. Conversely, a fixed Ret<1 would lose order
(1-Ret)log L and cannot preserve the full -n log log n term.

This is a sufficient conditional model for inexpensive retention,
not a completed proof of all n-dependent target and matching guards.

## 3. Exact pointwise cap limit under explicit typical-count conditions

Fix q in (0,Q]. On active matching R~s_n q->infinity,
a0~q and rho=nq/x. The exact strip mass from LEMMA.md satisfies
S(R,a0,1-q)->2(1-q). With J_i=1+O(L/n), every counted FREE
coordinate has strip probability

    p_i=(1-q)/(rho)*(1+o(1))
        =x(1-q)/(nq)*(1+o(1)).

Suppose the retained target arrays satisfy one-sided source-pin count
calibration N_-/n->1/2 and N_+/n->1/2 at this q, where
N_-=#{w_i<a0}, N_+=#{w_i>-a0}. These are the full-target
expectations in this scaling, since a0/r=O(1/s_n). The sixth event
in LEMMA.md, together with the old two-sided total count, can provide
such calibration if its n-dependent errors are independently paid.
For the SHARP cap limit additionally require that the transition count
#{i:a0<=w_i<a0+(1-q)} is o(n). This does not follow from the
one-sided count alone. It is the full-target typical behavior, since
that interval has probability <=(1-q)/(2r)=O(1/s_n) at fixed q;
paying its empirical exception in the retained target set is one
of the unclosed n-dependent interfaces of this conditional model.

For a captured positive source pin, the counted free coordinates give
factors (1-p_i). Coordinates in the transition interval can capture
part of the lower strip, contributing additional factors; there are
o(n) of them under the stated condition, and each has strip probability
O(1/n), so their aggregate logarithmic contribution is o(1). Other
uncounted coordinates capture none of that strip. The omitted pin
changes the count by at most one. As n p_i^2->0,

    P_inv(D^c | positive pin)
       ->exp[-x(1-q)/(2q)].                           (4)

Negative source pins have zero cap-complement probability for every
finite q<1. Since J_i->1, ellpin->1/2 and N_±/n->1/2, the
positive pins carry half of the complete inverse face coefficient.
Thus the relative UNCLIPPED weighted cap loss tends to

    f_x(q)=(1/2)exp[-x(1/q-1)/2].                    (5)

This factor 1/2 is a pin-face weight, not an extra independent
coordinate law. The exponential's factor two relative to the coarse
one-sided bound comes from the exact inverse strip mass S~2(1-q).
Mixing signs before computing the true face law would miss this.

Under the inexpensive-retention regime of Section 2, clipped and
unclipped cap losses differ by at most the total soft loss, which
is integrated o(1). Thus soft clipping does not remove the cost
(5) while preserving almost all the total log-mass.

## 4. Explicit conditional integrated model and its cancellation

To turn the pointwise calculation into a complete coefficient model,
assume all of the following are proved separately:

* raw radial rows are 2k up to errors whose target-log integrals are o(1),
  and both degrees are <=(2+o(1))*k with degree error times log L=o(1);
* target mass is n log(s_n Q)+o(n), with lower trim rlo->1;
* typical N_±, J_i and transition-interval calibration above holds,
  and q near zero has an integrable cap envelope, e.g. C exp(-c_x/q),
  before dq/q;
* total target-selection, transport, soft-retention and tau losses,
  measured in the normalized deficit, are o(n).

These conditions are not asserted to be a full proof here. Under them,
integrating (5) against the calibrated target dq/q gives the NECESSARY
model cap payment

    C_cap(x,Q)=(1/2) integral_0^Q exp[-x(1/q-1)/2] dq/q
              =(1/2)exp(x/2) E1(x/(2Q)).             (6)

The near-zero integral is finite and its envelope follows plausibly
from the one-sided bound with nu_one bounded away from zero; it must
still be incorporated in a full n-dependent proof. If rlo is kept
fixed above one, an additional -log rlo term lowers the coefficient.

The optimistic PAID coefficient is therefore

    C_model=log lambda+log Q+I(x)/2-C_cap(x,Q).        (7)

At Q increasing to one, C_cap increases to I(x)/2: the exact B cap
gain is canceled. This is a limit from Q<1, not a finite source band
including q=1. The endpoint q=1 has no effect on the integrals.

For every fixed 0<Q<1, let
g_x(Q)=log Q+[I(x)-exp(x/2)E1(x/(2Q))]/2. Then

    g_x'(Q)=[1-(1/2)exp(-x(1/Q-1)/2)]/Q>0,
    g_x(1)=0.

Consequently g_x(Q)<0 and C_model<log lambda for every Q<1.
Choosing Q<1 cannot improve the net balance: its log Q penalty
exceeds the saved cap cost. Since lambda<(1-log2)^-1,
log lambda<-log(1-log2). The conditional model thus gives no gain
beyond the R090 limiting coefficient. Extra nonnegative costs can
only lower it.

There is an even earlier obstruction if one uses only the UNIFORM
cap bound (3) of LEMMA.md: x and Q are fixed, so its exponent is
only a constant. Paying a fixed fraction of every target row would
lose order log L in the normalized deficit, destroying the target
T_n scale. One-sided counting repairs the old vanishing 1/s rate,
but does not make the cap probability vanish at fixed x. A radial
integrated treatment is needed even to reach the optimistic model.

The cancellation is a diagnosis of this calibrated inverse-kernel
model. It is not an impossibility theorem for different kernels,
cap-aware weights or a different covering argument.

## 5. One actual registered critical case and stopping scope

Exactly one case was registered and run locally:
x=2, lambda=3, Q=3/4, kappa=19/20 and alpha=1/40.
Independent log24 rational bounds verify a*lambda<kappa and
alpha+kappa<1. SciPy quadrature and its E1 implementation agree on
both (1) and (6). At the single n=1000000 diagnostic, I_n-I was
1.05479e-7. This is a finite consistency check, not proof of the limit.

The actual model values are

    I(2)=0.596347362323194,
    naive coefficient=1.109103897377926,
    positive-pin cap cost=0.174872920131172,
    paid conditional coefficient=0.934230977246754,
    log lambda=1.098612288668110.

The auxiliary run returned PASS with PID 33817; receipt.json records
exact rates, quadrature outputs and identities. No grid, old baseline
or author script was run. Sources were read only and their hashes
were checked before and after execution. No C_n improvement, explicit
dimension threshold, or completed n-dependent-s theorem is claimed.
This branch stops after this one lemma and this one scaling diagnosis.
