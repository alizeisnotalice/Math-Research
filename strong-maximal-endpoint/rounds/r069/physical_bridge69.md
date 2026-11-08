# A physical interpretation of the qualified band capture

## Local cost, rather than global efficiency

Fix the deterministic shifted reference c_N=mu-sqrt(mu)/2 sqrt(log N/N), or c=mu. Let Q_ref=1/J_{*,2} be its capped scalar benchmark. Its logarithmic coefficient kappa is respectively alpha=(9-mu)/24 or beta=(3-mu)/6. All kernels below retain the full orthant labels.

In the global logarithmic centroid chart of r068, let 0<=F<=2 be the actual physical input. At a receiver (u,X,sigma), put

\[
 h(s)=[K_c^sF_{u+s}](X,\sigma),\quad
 j_u^c=\int_0^\infty e^{-s}[K_c^s\Phi_N(F_{u+s})](X,\sigma)ds.
\]

This j is a **local conditional cost**. It is not the physical global modular and need not have an integrable full-s envelope. The actual receiver average at this fixed c is integral f_(N,c)(s)h(s)ds.

Suppose the receiver is successful and, for some fixed 0<=epsilon<2K,

\[
 j_u^c\le J_{*,2}N^\epsilon,\qquad K>\kappa.
 \tag{1}
\]

Conditional Jensen gives J(h)<=j_u^c. Hence the scalar capture theorem applies to h, without requiring h to be monotone. For I=[q0-K log N,q0+K log N],

\[
 j_{u,I}^c:=\int_I e^{-s}[K_c^s\Phi_N(F_{u+s})](X,\sigma)ds
 \ge J_{*,2}N^{\kappa-K+o(1)}.
 \tag{2}
\]

The asymptotic error is uniform over inputs satisfying these fixed-cap, fixed-slack hypotheses. In particular K=kappa+delta gives a loss N^(delta+o(1)). This is a local efficiency-qualified physical source-cost statement, with no angular averaging baseline assumed active.

## The true KKT benchmark concentrates its cost in the band

For the unnotched optimizer itself, a*=q0+kappa log N+o(log N), the height-one plateau has length log N+o(log N), and soft-plus-cap cost is O(J*/log N). When K>kappa, its cutoff lies inside I. The plateau cost above the upper endpoint is at most exp(-(q0+K log N)). Therefore

\[\frac{J_I(h_*)}{J_*}\ge1-O(1/\log N)-N^{-(K-\kappa)+o(1)}.\]

This explains why a linear ramp is an inappropriate test of capture at true scalar optimum: the optimal corner plateau carries the dominant Orlicz cost. The weaker N^(kappa-K+o(1)) bound is uniform over all qualified terminals, not just the optimizer.

## The exact remaining reference transfer

Write A_u^I=e^(-u)j_(u,I)^c. Formula (2) gives

\[
 e^{-u}\le Q_{\rm ref}N^{K-\kappa+o(1)} A_u^I.
 \tag{3}
\]

For a measurable set E of successful receivers satisfying (1), fix an angular-label state and let E_X be its radial section. Regardless of monotonicity,

\[
 \int_{E_X}e^{-u}du\le\sup_{u\in E_X}e^{-u}
 \le Q_{\rm ref}N^{K-\kappa+o(1)}\sup_u A_u^I.
\]

The empty-section value is zero. No measurable choice of a boundary graph is needed. Integrating and using r068's fixed-action band transport gives

\[
 |E|\le Q_{\rm ref}N^{K-\kappa+o(1)} Z_I\int\Phi_N(f).
 \tag{4}
\]

The physical coarea factor 2^N/sqrt(N) appears on both sides and cancels. Since Z_I=O(log N), the shifted reference and K=alpha+delta give Q_N N^(delta+o(1)). This is a consistency check of a capture-to-transport route; the exact fixed-action upper bound r028 is already stronger than (4).

The same route tolerates a reference receiver mean at least 1-eta_N with eta_N=o(log N/N). Indeed scaling by (1-eta_N)^(-1) costs at most (1-eta_N)^(-N)=N^o(1), using Phi_N(a x)<=a^N Phi_N(x). Its height remains below a fixed cap. The cap-two/three deletion argument in the scalar proof works for any fixed initial cap by using a larger fixed repair cap. All fixed caps have the same logarithmic calibration. Consequently the band capture and (3)-(4) retain their exponents.

What remains unproved is that actual variable-action successful receivers can be transferred to this **one reference action** with eta_N=o(log N/N), while retaining (1) for a common actual source and controlling the discarded halo. Moving physical anchors additionally require a compatible chart allocation. Neither the local cost hypothesis nor this action comparison follows merely from the global modular. This note does not prove a general C_N bound.

## Why a pointwise capture hypothesis cannot be inferred from global efficiency alone

For fixed K, the critical band satisfies P(S in I)=O((log N)^2/N). Choose a fixed A>1 and sufficiently large N. There is a bounded compact increment set J contained in I^c with P(S in J)>1/A. At a fixed nonzero receiver x, put height A on precisely the source points reached by its fixed-ratio rectangle with total logarithmic increment in J, and zero on the other increment layers. The physical source is compact and bounded: all increments are nonnegative with bounded total, so its nonzero coordinates stay away from zero. The receiver average is A P(S in J)>1, but j_(u,I)^c=0 exactly. Strictness gives a positive neighboring receiver set by continuity of fixed-ratio averages. This does not satisfy the local near-optimal cost qualification (1).

The example can coexist with arbitrarily good global efficiency. Start from a compact near-Q_ref input with coordinate magnitudes bounded below by a positive number, using the compact approximation already justified in r067. Dilate the new packet about the common origin so that all its coordinate magnitudes are smaller than that lower bound and its modular is at most eta times the main modular. Their supports are disjoint; the height cap is preserved. The main halo remains successful by positivity, so the total score is at least the main score/(1+eta). At the small receiver the main source is unreachable because the fixed-ratio increments only decrease coordinate magnitudes, so j_(u,I)^c remains zero there. Thus global near-Q_ref efficiency does not imply capture at every successful point. This is a qualification counterexample, not a new endpoint lower bound.

## Why positive action envelopes still lose a power

The full labelled conditional c-flux is A_(N,c)(s)(K_c^s-Face_c^s), with both kernels positive and mass one. Taking total variation before integrating sources only gives a positive maximal bound of size 1+2 integral_C A_(N,c)(s)dc. In a critical action window of width rho=sqrt(log N/N), A_(N,c)(s)~N p_mu, so this is O(sqrt(N log N)). The known sharp positive c-envelope obstruction persists for capped Orlicz inputs: height-one small angular-label packets have Phi_N(F)=F, and concentrating them in a small source-w layer recovers the fixed-s envelope. Adding a radial band integral alone does not cure it. This does not refute an efficiency-qualified nonpositive response or source-packing argument.

## Prior truncation used only in the Python calculation

The Python program retains binomial counts |K-Np|<=sqrt(110N), without renormalization. Hoeffding gives an omitted probability d<=2exp(-220). This is a subprobability increment law with density f_tilde<=f. The omitted mean of any height-A terminal is at most A d. For its exact capped KKT minimizer at multiplier lambda_tilde, write J_tilde for its cost. It is feasible for the full-law problem, whereas the full-law dual at the same multiplier differs by at most A lambda_tilde d: the capped conjugate is A-Lipschitz and integral e^-s(R-R_tilde)ds=d. Consequently

\[0\le J_{\rm tilde}-J_{\rm full}\le A\lambda_{\rm tilde}d.\]

The code records A lambda_tilde d/J_tilde in logarithmically stable form for A=2. This controls the ideal optimization error due to prior truncation. It does not cover floating arithmetic, finite quadrature, finite root tolerances, or the numerical approximation to that ideal KKT minimizer. Those are checked separately against full finite sums in Mathematica and tightened tolerances, with no directed-rounding certificate claimed. The centered Stirling series in the Python Gamma log-density includes terms through N^-9; its next-term bound controls this separate analytic approximation.
