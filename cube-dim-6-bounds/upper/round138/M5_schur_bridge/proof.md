# C138: a body16 Schur bridge (conditional; not a new C6 bound)

Let V5={(1,epsilon2,...,epsilon5)} and V6={(u,epsilon6):u in V5}.
All stable variables below have Laplace transform exp(-sqrt(lambda)),
equivalently S=1/(4Y), Y~Gamma(1/2,rate 1). Write
Sigma_d=sum_{v in Vd} S_v vv^T and M_d=E(det Sigma_d)^(-1/(2d)).
Thus these are the concrete 16-line and 32-line laws, without a C5
operator-constant inference.

## Algebraic proof

Block Sigma6 as [[A,b],[b^T,c]], with
A=sum_u(Su+ + Su-)uu^T. Independence across pairs and the stable
Laplace transform show A has the law of 4 Sigma5. For each real z,

  z^T A z-2b^Tz+c
    =sum_u[Su+(u.z-1)^2+Su-(u.z+1)^2]
    >=sum_u 4Su+Su-/(Su++Su-)=H.

Minimizing in z gives c-b^T A^-1 b>=H. Every matrix is positive
definite almost surely. Since 4Su+Su-/(Su++Su-)=1/(Yu++Yu-),
H=sum_{j=1}^16 E_j^-1 with independent E_j~Exp(rate 1). Consequently

  M6 <= E[(det A)^(-1/12) H^(-1/12)]
     <= [E(det A)^(-1/10)]^(5/6) [E H^(-1/2)]^(1/6)
     = 2^(-5/6) M5^(5/6) B^(1/6),                 (1)

where B=E H^(-1/2). This uses Holder exponents 6/5 and 6 and makes
no claim that A and H are independent. Finiteness follows, for example,
from the pinned concrete M5 certificate below and
H^-1/2 <= (sum E_j)^1/2/16 (Cauchy--Schwarz), whose expectation is finite.

## Scalar identity, including its constant

For t>0, independence and nonnegative integration give

  E exp(-tH)=[int_0^infty exp(-x-t/x) dx]^16
            =[2 sqrt(t) K1(2 sqrt(t))]^16.

Here K1 can equivalently be defined by the displayed positive integral;
no numerical special-function identity is assumed in the derivation.
The Mellin identity x^-1/2=pi^-1/2 int_0^infty t^-1/2 exp(-tx)dt
and u=2sqrt(t), for which t^-1/2 dt=du, prove exactly

  B=pi^-1/2 int_0^infty [u K1(u)]^16 du.          (2)

In particular there is no missing factor 2 in (2).

## Pinned concrete M5 source (read-only)

The 5D astra_upper014/report.md defines the exact same V5 and Laplace
normalization in section 1: V={(1,+/-1,+/-1,+/-1,+/-1)}, independent S_v,
Laplace transform exp(-lambda^beta), and A=sum S_v vv^T. For beta=1/2,
this is exactly the law S=1/(4Y), Y~Gamma(1/2,rate1), since direct
integration gives E exp(-lambda/(4Y))=exp(-sqrt(lambda)). The source's
conditional Gaussian covariance 2t^2 A is not itself being identified
with Sigma5; its underlying A is. Its section 3 equation (4) states
J_(1/2)=16 Gamma(3/2) M5=8 sqrt(pi) M5. Section 4 proves

  J_(1/2) <= 4 Gamma(3/2)
    [Gamma(6/5)/Gamma(11/10)]^5 (2661/64)^(-1/10).

Its budget.json case alpha=1 beta=1/2 records the upward rational
Jcap=2044690556653767/10^15. Hence the concrete certified input is

  M5 <= Jcap/(8 sqrt(pi)).                        (3)

upper_alpha1_029/check_results.json independently re-encloses this
expression and says contained_in_locked_source_ball=true. The 5D
STATE/ledger was used only to locate sources; no current/live 5D
computation, state, process, or checkpoint was read or modified.
Source hashes and exact keys are in source_manifest.json. This input
is an explicit determinant-moment bound, not a consequence of a C5
operator upper bound. No stronger concrete M5 certificate was located
in the bounded search of the ledger and pinned determinant branches.

## One actual exact comparison already excludes this budget

Let J_U=2044690556653767/10^15 and M5_U=J_U/(8sqrt(pi)). The specific
upper-budget expression supplied by (1) has sixth power

  Q_U=2^-5 M5_U^5 B=J_U^5 B/(2^20 pi^(5/2)).

Set E_min=min_j E_j. Pointwise H<=16/E_min, while E_min~Exp(rate16).
It follows that B>=E sqrt(E_min)/4=sqrt(pi)/32. Thus

  Q_U>=J_U^5/(2^25 pi^2)
      >J_U^5/[2^25 (22/7)^2]=L.                 (4)

One new standard-library Fraction batch ran at nice10 in one process,
with no time cap. It proved L>1/10^7 and all three targets below 1/10^7:

| target | display only | L/target, display only |
|---|---:|---:|
| Theta(g0), g0=1/65000 | 5.863551139121912107973e-8 | 1.838984486755113363687 |
| Theta(gstar), gstar=1/63500 | 6.745257062289188138183e-8 | 1.598601726007602999221 |
| round128 complete saved M_U^6 | 6.711146303358518834960e-8 | 1.606726942719792467237 |

Here Theta(g)=(8/9)^6*24576^5*g^6, and the g0 value is checked equal
to the exact threshold saved by round128. The exact L, source hashes
and three true Fraction comparisons are in candidate_certificate.json.
Runtime was about 0.01455 seconds; there were zero Bessel evaluations,
zero optimization calls and zero old baseline reruns. Reproduce with

  nice -n 10 /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/venv/bin/python -B check_budget.py

Conclusion: this current concrete M5 certificate plus this Holder bridge
cannot beat the saved round128 moment upper budget or meet either
threshold. This is a strict lower bound for an upper-budget expression,
not a lower bound for actual M6, its optimal upper budget, or C6. A
stronger concrete M5 certificate could still make the bridge useful.
No global g claim or new C6 bound is supplied.

## Optional scalar tail lemma (not evaluated)

The following establishes a finite 1D strict packet if a future sharper
M5 input warrants one. The actual installed python-flint 0.9.0 API was
inspected: arb(u).bessel_k(1), with u the argument and 1 the order. The
prepared Bessel execution script was removed after (4) separated the
budget; no such scientific batch was launched.

For f(u)=[uK1(u)]^16, its positive-integral representation above gives
0<f(u)<=1, and f is decreasing: exp(-u^2/(4x)) decreases in u. Thus
the omitted [0,a] integral is in [0,a], and right/left endpoint sums
rigorously enclose the integral over [a,R]. Their exact difference is
h(f(a)-f(R))<=h, where h=(R-a)/N. No quadrature remainder estimate is used.

For u>=1, the standard positive K1 representation follows by substituting
x=(u/2)exp(t) in the defining integral and then y=cosh(t):

  K1(u)=u int_1^infty exp(-uy) sqrt(y^2-1)dy
       <=exp(-u)[u^-1+sqrt(pi/(2u))].

The last step uses sqrt(s^2+2s)<=s+sqrt(2s), followed by elementary
Gamma integrals. Since pi<4 and u>=1,
1+sqrt(pi*u/2)<=3u. The omitted infinity tail is therefore at most

  T=3^16 exp(-16R) sum_{k=0}^16
        [16! R^k/(k! 16^(17-k))].                (5)

This is the exact evaluation of 3^16 int_R^infty u^16 exp(-16u)du;
it can be enclosed with Arb. Both tail budgets are established
without launching a Bessel batch.

If B is enclosed in [BL,BU], equation (1) gives an upper certificate
(M6_bridge)^6 <= M5cap^5 BU/32. For a target m, sufficient M5 threshold
is (32*m^6/BU)^(1/5); the necessary threshold for this particular
scalar factor is at most (32*m^6/BL)^(1/5). No numerical B enclosure or
threshold is claimed in this round.
