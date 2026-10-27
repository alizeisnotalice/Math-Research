# R032 — joint witness-mass and covering-flow interfaces

## Outcome and scope

This report gives two finite, computable interfaces that retain witness information. Every original instance maps explicitly into them, and their bound directions are proved. Neither interface alone proves a new uniform constant. No seven-dimensional optimization, operator lower bound, MATLAB run, commit, push, fetch, or rebase is used.

The operator remains the centered axis-parallel common-radius cube maximal operator on Lebesgue measure, with arbitrary L1 input and all positive radii. All constants below refer to that same operator. Source checkout: `f65a5bf92911ebed4cf803071bd7693978e5e582`; source hashes and actual small checks are in `verification.json`.

The main results are:

* An exact finite witness-mass LP and its packing dual. Its geometric supremum K_{n,a} obeys C_n <= K_{n,a} <= (1+1/a)^n C_n, so inf_a K_{n,a}=C_n.
* An adaptive-weight selection followed by a capacitated covering-flow LP. Its load includes selected and rejected volume in one coefficient, before taking a pointwise maximum. Its dual consists of a probability on witness incidence cells and independent fractional-knapsack problems on covering cells.
* A precise limitation: witness margins alone have no uniform positive gap. Every finite geometric family can have all margins simultaneously equal to 1+epsilon. Full joint incidence or flow information is therefore essential.

## 1. Exact witness-mass LP

Fix a>1 and a finite family Q_i=Q(y_i,r_i), i=1,...,N. Write V_i=aQ_i and U=union_i Q_i. Normalize lambda=1 for this section. Let E_sigma be the measurable incidence cells

    E_sigma = {x: {i:x in V_i}=sigma},    empty != sigma subset {1,...,N}.

Ignore cells of zero Lebesgue measure. All relevant cells have finite measure. The number of nonempty cells is finite and can be computed from the coordinate endpoints of the V_i. Define

    M(V)=min Σ_sigma m_sigma
         subject to m_sigma>=0,
         Σ_{sigma contains i} m_sigma >= |V_i|   for every i.        (P)

### Lemma 1: exactness and dual

Every original f with integrals over V_i greater than lambda|V_i| gives a feasible vector

    m_sigma=(1/lambda) integral_{E_sigma}|f|,
    Σ m_sigma <= ||f||_1/lambda.

Consequently

    lambda |U|/||f||_1 <= |U|/M(V).                              (1)

Conversely, every feasible finite vector m is realized by a nonnegative L1 function

    g(x)=Σ_sigma [m_sigma/|E_sigma|] 1_{E_sigma}(x).

It has integral Σm and integral over V_i at least |V_i|. Multiplication by 1+epsilon makes all witnesses strict. Thus (P) loses no information about minimal L1 mass for this finite witness family; replacing cells by a finite grid or dropping constraints would require a separate relaxation-direction argument.

The finite LP is feasible (g=1 on union V_i), has positive finite optimum, and its dual is

    M(V)=max Σ_i z_i |V_i|
         subject to z_i>=0,
         Σ_{i in sigma}z_i<=1    for every positive-volume cell.   (D)

Equivalently Σ_i z_i 1_{V_i}<=1 almost everywhere. Both optima are attained by elementary finite LP duality. In particular a geometric packing certificate with objective at least |U|/K proves (1) <=K for every original witness input. The inequality direction is important: a feasible primal supplies an upper bound on M and is insufficient for an upper bound on C_n; a feasible dual supplies a lower bound on M and has the correct direction.

### Lemma 2: exact variational sandwich

Let

    K_{n,a}=sup over all finite cube families |union Q_i|/M({aQ_i}).

Then

    C_n <= K_{n,a} <= (1+1/a)^n C_n.                              (2)

Proof of the left inequality: take a compact subset F of an original strict level set. For each x in F choose any true witness radius R_x. The open cubes Q(x,R_x/a) finitely cover F. The original f maps to a feasible vector in (P), so lambda|F| <=K_{n,a}||f||_1. Lebesgue inner regularity proves the assertion for the full level set. No scale truncation or compact-support assumption on f occurs.

Proof of the right inequality: choose an optimal g for (P). If x in Q_i, then V_i is contained in the centered cube Q(x,(a+1)r_i), so

    M_n g(x) >= (a/(a+1))^n.

For any 0<delta<1 the whole U lies in the strict level set at (1-delta)(a/(a+1))^n. Apply the definition of C_n and then let delta decrease to zero. This gives |U|/M(V)<=(1+1/a)^n C_n. Supremizing proves (2).

Since the old proof already supplies finite C_n, (2) also proves

    inf_{a>1} K_{n,a}=lim_{a to infinity} K_{n,a}=C_n.              (3)

No monotonicity in a is asserted. Equations (2)-(3) are an exact geometric LP reformulation, not an explicit improved estimate. A finite collection of numerical families bounds neither the supremum in K_{n,a} nor C_n from above.

## 2. Adaptive witness budgets and legal selection

Return to an original input f and lambda>0, with every V_i=aQ_i a true witness. Put

    M_i=integral_{V_i}|f|,
    |V_i| <= h_i <= M_i/lambda.                                 (4)

The budgets h_i are fixed for the input family before selection; they may depend on the actual witness masses. Two useful endpoints are h_i=|V_i| and h_i=M_i/lambda. Interpolations are also allowed. Fix any b>0 and tau>0. There is no need for b/a=1+1/t or a>2t in the following covering-flow lemma.

Select a largest-radius remaining cube satisfying

    Σ_{j already selected} p_j 1_{bQ_j}(y_i)<=tau.

After selecting Q_j define, once and for all,

    I_j=Q_j minus union_{ell<j}Q_ell,
    v_j=|I_j|,
    p_j=v_j/h_j <= v_j/|V_j|.                                   (5)

The budgets h_j are fixed; increments and weights do not change afterwards. Their nonnegativity makes eligibility monotone and selected radii nonincreasing. Zero increments simply give p_j=0.

Let S be the selected union, D=U minus S, and E_j=(b+1)Q_j. The paper's first-rejection-prefix argument, with p in place of w, proves

    P(x):=Σ_j p_j 1_{E_j}(x)>tau   for almost every x in D.        (6)

Indeed the rejected cube already fails eligibility before any smaller selected radius, or at termination. The rejecting prefix therefore has radii at least the rejected radius. Translating a point inside the rejected cube increases the relevant normalized distance by less than one, placing it in E_j. This proof does not use the old load estimate and remains valid for the full budget interval (4).

## 3. A joint covering-flow LP

Partition D by the E_j incidence pattern. Let D_sigma have positive volume d_sigma and let sigma={j:x in E_j}. Then P_sigma=Σ_{j in sigma}p_j>tau. Introduce variables z_{sigma,j} only for j in sigma, subject to

    Σ_{j in sigma} z_{sigma,j}=d_sigma,
    0<=z_{sigma,j}<=d_sigma p_j/tau.                             (7)

These constraints always have the explicit original-instance mapping

    z_{sigma,j}=d_sigma p_j/P_sigma.                             (8)

Put c_j=Σ_sigma z_{sigma,j}. Thus Σ_j c_j=|D|. The flow assigns every unit of rejected volume to an allowed selected witness, with capacity inherited from the actual first-rejection weights.

Partition the positive-volume part of union_j V_j by its witness incidence rho. The following is a finite LP in z and L:

    L_*(Q,h,b,tau)=min L
    subject to (7),
    Σ_{j in rho} (v_j+c_j)/h_j <= L for every witness cell rho.   (CF)

### Lemma 3: joint certificate and correct objective direction

Every original instance satisfies

    lambda |U|/||f||_1 <= L_*(Q,h,b,tau).                        (9)

For every feasible flow, write

    G_z(y)=Σ_j [(v_j+c_j)/h_j]1_{V_j}(y).

The same original witness masses, not independently chosen masses, give

    integral |f|G_z =Σ_j [(v_j+c_j)/h_j]M_j
                  >=lambda Σ_j(v_j+c_j)=lambda |U|.             (10)

The LP bounds the essential maximum of this **combined** load. Therefore (9) follows. If h_j=M_j/lambda, the inequality in (10) is equality. If h_j=|V_j|, (10) retains all true-witness inequalities without requiring f-dependent selection.

The targets in (P) and (CF) have opposite optimization directions: a feasible (D) dual gives a useful bound through M; a feasible (CF) primal gives a useful upper bound on L_*. A feasible (CF) dual is only a lower bound on the best certificate L_* and cannot certify a new upper bound.

### Relation to the old product

The capacities imply c_j<=p_j|E_j intersect D|/tau<=p_j|E_j|/tau. Since p_j=v_j/h_j and |E_j|/h_j<=|E_j|/|V_j|=((b+1)/a)^n=:A,

    G_z(y)<=Σ_j p_j[1+|E_j|/(tau h_j)]1_{V_j}(y)
            <=(1+A/tau)Σ_j p_j1_{V_j}(y).                      (11)

Thus any old pointwise load bound for the same adaptive selection recovers its product bound. The flow LP keeps the actual D, exact cover capacities and common witness incidence before maximizing, so it is no weaker than that recovery. It does not assume that the independently worst cover and load configurations are jointly realizable.

At h=M/lambda, define eta_j=M_j/(lambda|V_j|)>1. The loose recovery uses |E_j|/h_j=A/eta_j, making the single joint load Σ_j p_j(1+A/(tau eta_j))1_{V_j}. Replacing eta_j by one recovers the old loss. Keeping eta_j separately is legitimate, but needs a uniform joint estimate to improve asymptotics.

## 4. Exact dual and local fractional-knapsack cost

For a probability vector y_rho on the positive-volume witness cells, put

    kappa_j=(1/h_j)Σ_{rho contains j}y_rho.

For a cover cell sigma define the local cost

    chi_sigma(kappa)=min Σ_{j in sigma}kappa_j x_j
       subject to Σx_j=1, 0<=x_j<=p_j/tau.                     (12)

It is finite by (6), and is computed by sorting kappa_j increasingly and filling the available p_j/tau capacities until total one is reached. Finite LP minimax duality gives the exact identity

    L_*=max_{y>=0, Σy=1}
             [Σ_j kappa_j v_j+Σ_sigma d_sigma chi_sigma(kappa)]. (CFD)

One proof is to write max_rho load as max over the cell probability simplex, then interchange the minimum over the nonempty compact flow polytope with that maximum; the bilinear form separates by cover cells. This is also ordinary finite LP strong duality. The y variables are dual probes of pointwise load, not a replacement for the original f in (4) or (10).

Equations (CF)-(CFD) give a small Bellman-like local interface: a witness-cell probability sets the charge prices kappa, and each rejected cell pays its exact cheapest capacity-respecting local price. The missing global theorem is a dimension-dependent bound for this coupled game over every geometric input and every legal original witness instance. A diagnostic maximizing or minimizing only finitely many cells/families does not supply that theorem.

## 5. A simple margin-only obstruction

### Lemma 4: no universal strict witness-margin gap

For any finite geometric family and epsilon>0 choose

    f_epsilon=(1+epsilon)lambda 1_{union_i V_i}.

This is an L1 function and every original witness mass is exactly

    M_i=(1+epsilon)lambda|V_i|,   eta_i=1+epsilon.               (13)

Thus all margins can be arbitrarily close to one, simultaneously and with exactly consistent common f. Witness validity alone cannot yield eta_min>=1+delta_n for any positive delta_n. This rules out arguments whose entire improvement comes from such an assumed positive margin. It does not rule out incidence-based mass constraints, covering-flow savings, or a combined argument.

There is also an exact scaling identity at the full f-dependent endpoint. Under (13), compare with h_i=|V_i| and threshold tau_0. Use h_i=(1+epsilon)|V_i| and tau=tau_0/(1+epsilon). Then every selection decision is identical, the selected increments are identical, the capacities p_j/tau are identical, and

    L_*^{margin}=L_*^{geometric}/(1+epsilon).                   (14)

The fractional improvement disappears as epsilon decreases to zero. This identity applies even if the original selection has equality ties, because both sides of the eligibility inequality scale by exactly the same factor. It concerns the certificate family, not a lower bound for C_n.

## 6. Minimal computations, evidence and open work

`verify.py` performs just two finite examples, in dimensions one and two, with rational cube data, plus the elementary scaling check. It computes all incidence cells from coordinate cuts; no sampled points stand in for positive-volume cells. It verifies the witness-mass primal/dual direction, maps the uniform original f of (13), runs adaptive selection, checks every covering-cell capacity and the equality in (10), and solves the joint flow LP. The returned witness-cell dual is checked against the fractional-knapsack formula (CFD). Rationalized primal/dual certificates are checked exactly where recorded.

Actual results: the one-dimensional example has M(V)=48/5 and flow load 193/330; the two-dimensional example has M(V)=1224/25 and flow load 3666673/14434200. The mass-LP primal and dual, the flow primal and its knapsack dual, and the margin scaling identity were all checked as exact rational equalities after the small LP solves. These numbers only check the interfaces on the stated finite examples.

These are interface checks only. They do not establish a new C_n upper bound or diagnose high-dimensional optimal parameters. The report performs no scan of the excluded R030 soft-mixture family.

Remaining proof obligations for an improvement are explicit:

1. Construct a packing dual (D) uniformly with objective >=|U|/K_n, or a feasible covering flow (CF) uniformly with load <=K_n, for **all** finite geometric families/original inputs.
2. Exploit full witness incidence, actual rejected-volume allocation, or a scale-local potential to obtain a fixed proportional saving in the n log n regime. A bounded additive saving in a cap integral is insufficient.
3. Only then use Lemma 2 or Lemma 3 plus compact-level-set coverage to conclude a general-L1/all-positive-scales operator estimate. Every finite diagnostic remains separate until this uniform step is proved.

Dependencies: the operator/level-set regularity and first-rejection prefix in `paper/main.tex`; elementary finite LP strong duality; cube containment and Lebesgue incidence partitioning. No unverified external mathematical result is needed.
