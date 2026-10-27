# R039: arbitrary-pair payment optimization and the limits of strict merit descent

## Scope and outcome

This note proves an exact one-dimensional payment optimization for **every** witness pair, adjacent or nonadjacent. It has at most four affine load branches, at most eight rational endpoint/intersection candidates, and an exact plateau rule for lexicographic merit (L,sum q,E). It replaces the premature inference that stopping nonedge moves requires a multi-vertex mechanism.

A separate, smaller P3 example proves that even the revised all-pair merit can be stationary with a genuine payment-LP gap. Its optimal improvement is available either as one three-coordinate move or as two pair moves with a temporary energy increase on a load plateau. Thus one-step merit stationarity must not be confused with global payment optimality or with impossibility of all pair paths.

All statements concern finite witness-payment upper certificates for the same Lebesgue centered common-radius cube operator, general L1 input and every positive scale. No new C_n bound, operator lower bound, dimension-seven case, geometry/atom recomputation, old-LP re-solve, or 79-case production rerun is performed here. The separate authorized production diagnostic is limited to three cached final states times two pair domains, one step per domain, owned by Terra.

## 1. Fixed payment and true load

Let V_i=aQ_i be finitely many true witnesses of one fixed input f and lambda>0, and let 0<h_i<=integral_{V_i}|f|/lambda be fixed budgets. For q_i>=0 with

    sum_i h_i q_i=|U|,       U=union_i Q_i,

payment conservation gives

    lambda|U|<=integral |f| sum_i q_i1_{V_i}
                <=||f||_1 L(q).

For axis-parallel open boxes, coordinate Helly identifies the true load as

    L(q)=max_{C in complete maximal witness cliques}sum_{i in C}q_i.

The allowed graph is all original witnesses, including currently zero coefficients which a pair may activate. A complete maximal-clique list remains sufficient after every nonnegative reassignment. The mass-LP optimum |U|/M(V) uses h_i=|V_i|; general budgets in this section must not be silently substituted into that particular mass comparison.

## 2. Four clique classes for an arbitrary pair

Fix distinct i,j, not necessarily nonadjacent, and hold all outside coefficients fixed. Write

    P=h_iq_i+h_jq_j,
    x=q'_i in [0,P/h_i],
    y=q'_j=(P-h_i x)/h_j.

Partition the complete maximal-clique list by membership of i,j. For each nonempty class define its outside load maximum:

    A00=max_{C contains neither}sum_{k in C}q_k,
    A10=max_{C contains i only}sum_{k in C minus i}q_k,
    A01=max_{C contains j only}sum_{k in C minus j}q_k,
    A11=max_{C contains both}sum_{k in C minus {i,j}}q_k.

If a class is empty, omit its branch; do not interpret an unproved empty-class value as an actual clique. The exact new load is

    Phi(x)=max of the present branches
      A00,
      A10+x,
      A01+P/h_j-(h_i/h_j)x,
      A11+P/h_j+(1-h_i/h_j)x.                               (1)

This follows by partitioning every term in the true maximal-clique maximum. Nonadjacency is exactly the special case with no 11 branch. An adjacent pair may have no 10 or no 01 branch, and a complete graph may have only the 11 branch.

### Proposition 1: finite exact load minimum

Include domain endpoints 0,P/h_i. For every two present branches with different slopes, include their intersection abscissa if it belongs to the domain. There are at most six such intersections and at most eight candidates in total. The minimum m of Phi is the minimum of its exact rational values on those candidates.

Indeed a maximum of finitely many affine functions is convex and piecewise affine on the compact domain. Its minimum occurs at an endpoint or a breakpoint, unless it is constant on a segment; in that case a segment endpoint or a domain endpoint also attains the minimum. Parallel/coincident branches require no additional candidate. When P=0 the only state is unchanged.

No LP or floating root search is required. Recompute current branch prices after every change; stale prices do not optimize the new state.

## 3. Exact plateau interval and the revised merit

Write the present branches as a_r*x+b_r. Having obtained m, start with [ell,u]=[0,P/h_i] and impose every inequality a_r*x+b_r<=m:

    a_r>0:  u=min(u,(m-b_r)/a_r),
    a_r<0:  ell=max(ell,(m-b_r)/a_r),
    a_r=0:  b_r<=m is required.

The resulting nonempty closed interval is exactly the load-minimizing set.

Use global merit

    (L(q), S(q), E(q)),
    S(q)=sum_i q_i,   E(q)=sum_i h_iq_i^2,

in that lexicographic order. Outside contributions to S and E are fixed for the pair. The pair sum is

    x+y=P/h_j+(1-h_i/h_j)x.                                (2)

### Proposition 2: exact pair lexicographic minimizer

Within the load-minimizing interval:

* if h_i>h_j, choose x=u;
* if h_i<h_j, choose x=ell;
* if h_i=h_j, choose x=clip_[ell,u](P/(h_i+h_j)).

Then y=(P-h_i*x)/h_j. This exactly minimizes (L,S,E) for the pair. The first two cases follow from the strict sign of the slope in (2), leaving a unique S-minimizing endpoint. In the equal-budget case S is constant, and

    d/dx [h_i*x^2+h_j*((P-h_i*x)/h_j)^2]
       =(2h_i/h_j)*((h_i+h_j)x-P),

so clipping its unique unconstrained minimum gives the E tie rule. Every step conserves payment and stays nonnegative. The old pair is feasible, so the optimized merit never increases.

For a deterministic scan compare candidates by `(newL,newS,newE,i,j)` and accept only strict merit decrease. Differences in S/E may equivalently replace global values if the same current state is used for all candidates. Strict merit descent does not imply finite termination or global payment optimality; an iterated experiment would still need a finite move budget. The currently authorized diagnostic makes no iteration.

The added S level is substantive. In a load plateau it permits moving payment toward a larger budget, even if E increases. This makes R035 concentration-compatible plateaus available rather than blocking every such move with the old energy tie rule.

## 4. What active peak cliques can certify

Let A(q) be the current maximal cliques whose load equals L(q).

### Proposition 3: a necessary donor condition for immediate pair load decrease

Every nontrivial payment-conserving pair move increases one coefficient and decreases the other. Let j be the decreasing donor. If an active peak clique does not contain j, its load is unchanged or increases, depending on whether it contains the receiver. Thus a strict decrease of L requires

    j belongs to every C in A(q).                           (3)

In particular, an empty intersection of current active peak cliques rules out **every finite one-pair strict load decrease**, not merely a gradient direction and not merely nonedge pairs.

This does not rule out plateau pair moves. If an independently certified global payment-LP gap remains, convexity gives some strict-load direction toward a better feasible coefficient vector. Under the empty-intersection condition such an immediately strict direction changes at least three coordinates. Nevertheless a sequence of pair moves on a plateau can sometimes reach it; Proposition 3 does not exclude that path.

### Why the old E tie can produce a false barrier

For total payment P_total=sum h_iq_i,

    E(q)>=P_total^2/(sum_i h_i),

with equality exactly when all q_i are equal. This is weighted Cauchy--Schwarz or the weighted variance identity. Hence a constant-q state forbids every changed load-plateau state under lex(L,E), regardless of its future usefulness. The new S level may open such a state when budgets differ; it cannot help when all h are equal, because S is then fixed by payment.

## 5. A three-cube exact stationary example for the revised merit

Fix arbitrary n>=1 and a=4. Let original Q_0,Q_-,Q_+ all have radius 1 and centers 0,-5e_1,+5e_1. They are pairwise disjoint. The witnesses have radius 4, common budget h=8^n, and graph P3: the center witness intersects each leaf, but the two leaves are disjoint. Put

    v=2^n,  p=v/h=4^(-n),  |U|=3v.

Choose a full-selection threshold tau=3p. Every previously selected weight is p and any eligibility load is at most 2p, so all inputs are eligible for any fixed b; a legal nonincreasing-radius selection has increments v and baseline q_0=q_-=q_+=p. Its true load is

    L0=2p.                                                  (4)

The witness-mass primal assigns mass h to each of the two maximal edges, total 2h. Each leaf constraint is met and the center receives 2h>=h. The packing dual gives each leaf coefficient one and the center zero; it is feasible because no clique contains both leaves. Thus

    M(V)=2h,     Lopt=|U|/M(V)=3p/2.                        (5)

### Proposition 4: all-pair lex(L,S,E) stationarity with a positive LP gap

The constant vector q=p is stationary for every exact pair minimizer of the revised merit.

For a center-leaf pair, equal budgets make its pair sum 2p invariant, so the corresponding edge stays at load 2p. For the leaf-leaf pair, their sum stays 2p while center p is fixed; the larger leaf is at least p, so the true load is at least 2p. Therefore no pair can strictly lower L. All budgets equal h, so S is globally fixed at 3p by payment. Finally q=p is the unique global E minimum on that payment simplex. Any changed pair with unchanged L and S increases E. Thus no strict pair merit decrease exists.

The payment-LP improvement is nevertheless p/2, exactly 25 percent of the initial load. The three-coordinate reassignment

    q'_0=0,    q'_-=q'_+=3p/2

conserves payment and attains (5). This is a certificate-family obstruction, not a lower bound for C_n.

### A plateau pair path still exists

First transfer coefficient p/2 from the center to one leaf. The state becomes `(p/2,3p/2,p)` up to leaf order. Its load is still 2p, S remains 3p, but E rises from 3hp^2 to (7/2)hp^2. Then transfer the remaining center coefficient p/2 to the other leaf. The final load is 3p/2 and E=(9/2)hp^2.

Thus even this genuine one-pair merit stationary state is not evidence that all pair sequences fail. It distinguishes an immediately strict multi-coordinate direction from a pair path which temporarily worsens a secondary potential while keeping load fixed.

For true common-input witnesses take f=(1+epsilon)1_{union V_i}, lambda=1. No claim about an extremizing maximal-operator input is needed.

## 6. Limited read-only observations on the three R038 cached states

Before the new all-pair production diagnostic, this author performed only small exact cache reads to determine current peak structure. No pair scan or old LP was run here.

| cached final state | complete cliques | active peak cliques | their intersection | q constant |
|---|---:|---:|---|---|
| C/2/5 full_disjoint | 87 | 3 | empty | yes |
| C/1/6 full_disjoint | 71 | 24 | empty | no |
| A/3/15 partial | 37 | 24 | empty | no |

Thus Proposition 3 already rules out an immediate single-pair strict load decrease for all three current states. It does not determine whether a plateau with lower S, or lower E at fixed S, exists. That is the purpose of the authorized new-merit diagnostic.

The source cache identities read were:

    phasea_C_2_05_full_disjoint.json: aa259d09c5cee0e974f8d2a8f52b3bfd8b0fe1355e44007551749dc46f3726de
    phasea_C_1_06_full_disjoint.json: 2c755add961ccddc1e71d482c0d0491844a6358c272094bf53b9ab5c7f1a3649
    phasea_A_3_15_partial.json: 5582abceeb789961dd24714b10c336304b1a34f81f7817b8a612ce1b42c0b814

Production is restricted to these three final q_after vectors, with one complete nonedge-domain scan and one complete all-pair-domain scan for each: at most six scans, no iteration. It uses the frozen exact h/complete cliques and emits branch maxima, rational minimizer candidates, optimal interval, selected pair, new q, conserved payment and true new merit. It does not recompute the 79 baselines, union atoms, selection or old mass LP. New mass-LP certification belongs to the separately assigned numerical audit.

Possible conclusions must be kept distinct: immediate pair load decrease; load plateau with S decrease; load/S plateau with E decrease; or stationarity for this fixed one-step merit and tested pair domain. None alone proves a uniform C_n estimate or global payment optimality.
