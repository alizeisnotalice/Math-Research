# R037: a concentration obstruction and a legal independent-splitting interface

## Status, scope and outcome

This is a mathematical post-batch analysis and a new local payment interface, not a production rerun or a new C_n estimate. The actual 48-cube diagnostic motivated the obstruction, but the counterexample below is self-contained: four cubes in every dimension n>=2, an exact legal first-rejection history, exact covering-flow payments, and an exactly solved mass LP.

It proves that the R035/R037 rule requiring a positive saving from clique-to-largest-budget concentration can stall even when a positive fraction of the best payment-certificate improvement is available. A complementary independent-module splitting operation removes this obstruction. An arbitrary independent set has an explicit boundary bill. For a nonedge pair, an even stronger exact closed-form coordinate minimizer avoids the conservative bill entirely.

All statements concern the original fixed Lebesgue centered common-radius cube operator, general L1 inputs and all positive radii. They are finite-family certificate statements. They imply an operator bound only after a uniform estimate for every finite original witness instance, which is not proved here. No operator lower bound or dimension-seven specialization is asserted.

## 1. Conserved-payment framework

Let Q_i be finitely many original cubes, V_i=aQ_i true witnesses for a single f and lambda>0, and h_i=|V_i|. More general fixed positive budgets 0<h_i<=integral_{V_i}|f|/lambda are also permitted in Propositions 1-4, but the mass-LP comparison in the counterexample uses h_i=|V_i|.

For nonnegative load coefficients q assume

    sum_i h_i q_i=|U|,       U=union_i Q_i.

Then

    integral |f| sum_i q_i 1_{V_i} >=lambda |U|.                (1)

Every nonnegative reassignment preserving this scalar payment remains valid, even when it does not preserve the old covering flow's point-to-cube support. For axis-parallel open boxes, coordinate Helly gives

    L(q)=esssup sum_i q_i1_{V_i}
        =max over cliques C of the witness graph G sum_C q_i.  (2)

Every finite clique has positive-volume common intersection. Thus a complete maximal-clique list computes the true load.

## 2. Structural obstruction to all concentration candidates

### Proposition 1: maximal-budget support gives zero concentration saving

Suppose q_i>0 implies h_i=H=max_j h_j. For an arbitrary clique C choose k with h_k=max_{i in C}h_i and replace its charge by

    q'_k=sum_{i in C}h_i q_i/h_k,   q'_i=0 for i in C minus k.

Then the R035 saving

    s_C=sum_{i in C}q_i-q'_k

is zero for every C. If C has positive charge, its maximal budget is H and every positive term in it has budget H, proving q'_k=sum_Cq_i. If C has no positive charge both sides vanish.

Consequently any algorithm requiring s_C>0 makes no move. This conclusion does not depend on the maximal-clique search budget, the candidate pool, radius thresholds, or whether C contains initially zero coefficients. Expanding the same concentration pool cannot resolve this obstruction.

This does not prove that every s_C=0 reassignment is useless. It identifies the failure of the specified strict-saving concentration rule. In particular splitting payment between nonadjacent equal-budget witnesses is a different operation.

## 3. Four-cube, all-dimensional exact first-rejection example

Fix n>=2, a=4. Let all original Q_i have radius 1, with centers

    y_0=4e_1,  y_1=-4e_1,  y_2=4e_2,  y_3=-4e_2.

The small cubes are pairwise disjoint, each has volume v=2^n, and |U|=4v. Every witness has radius 4 and volume H=8^n. Its intersection graph is K_{2,2}: opposite pairs {0,1} and {2,3} are disjoint, while all four cross pairs intersect in positive volume. Put p=v/H=4^(-n).

Use the legal largest-radius tie order 0,2,3,1, the eligibility dilation b=7 and threshold tau=3p/2. The selected increments are v, since all Q_i are disjoint.

* Cube 0 is selected first.
* Cube 2 receives eligibility load p from 0, so is eligible.
* Cube 3 receives p from 0 and none from 2, because ||y_3-y_2||=8>7. It is eligible.
* Cube 1 receives p from each of 2 and 3, and none from 0. Its load 2p exceeds tau. It is rejected and selection stops.

Thus S=Q_0 union Q_2 union Q_3 and D=Q_1. The transported covering boxes are E_j=(b+1)Q_j=8Q_j. Both E_2 and E_3 contain all of D. The box E_0 contains exactly its half x_1>-4; the other half has x_1<-4. Boundaries are null sets. Hence the two rejected atom types have volumes v/2 and E-memberships {0,2,3} and {2,3}.

The exact proportional covering flow assigns

    c_0=(v/2)/3=v/6,
    c_2=c_3=(v/2)/3+(v/2)/2=5v/12.

The genuine legacy joint load is therefore

    q_0=7p/6,   q_1=0,   q_2=q_3=17p/12.                    (3)

It satisfies sum h_iq_i=H*(7/6+17/12+17/12)p=4v=|U| exactly. The true load is

    L0=q_0+q_2=31p/12.                                     (4)

All budgets equal H, so Proposition 1 makes s_C=0 for every clique. Algorithm D makes zero accepted moves, not merely zero moves in its truncated candidate pool.

### Exact witness-mass optimum

The mass LP has one variable for each of the four maximal edges. Assign mass H/2 to each edge. Every vertex receives total H and total mass is 2H, proving M(V)<=2H. Conversely the packing dual z_2=z_3=1, z_0=z_1=0 is feasible because V_2,V_3 are disjoint; its objective is 2H. Thus

    M(V)=2H,      Lopt=|U|/M(V)=2p.                          (5)

The available certificate gap is 7p/12, a fraction 7/31 of L0. It exceeds 10^(-4)L0, while the strict-saving concentration algorithm captures zero of it. The same half-gap implication fails on this standalone finite instance outside the registered batch; this example alone makes no claim about registered batch counts or case identities. This is not a claim of a C_n lower bound: the object being compared is a family of upper certificates.

Finally, split q_0 equally between vertices 0 and 1:

    q'_0=q'_1=7p/12,   q'_2=q'_3=17p/12.

Payment is conserved and every maximal edge has load 2p. The split attains the exact optimum (5).

These are true common-input witnesses: f=(1+epsilon)1_{union_i V_i}, lambda=1, makes every integral equal to (1+epsilon)H. All calculations above therefore apply to one consistent general-L1 input. No extremal input for the maximal operator is asserted.

## 4. Independent-module splitting

Let F be an independent set of witness vertices. Define

    P_F=sum_{i in F}h_iq_i,
    t_F=P_F/(sum_{i in F}h_i),
    q'_i=t_F for all i in F,   q'_i=q_i outside F.            (6)

This conserves payment, remains nonnegative, and satisfies t_F<=q_max:=max_{i in F}q_i.

### Proposition 2: exact module load formula

Assume every outside vertex with positive current q is adjacent either to all of F or to none of F. Let G_F be the common outside neighbors of F, L_good their maximum clique load, and L_rest the maximum clique load on vertices outside F. Then

    L(q)=max(L_rest,q_max+L_good),
    L(q')=max(L_rest,t_F+L_good)<=L(q).                       (7)

Each clique contains at most one F vertex. The active outside part of any clique meeting F consists only of common neighbors, and every common-neighbor clique can be joined to every member of F. This proves both equalities. Strict decrease holds if q_max>t_F, the affected branch is maximal, and the unaffected branch is strictly smaller.

Unlike clique concentration, independent splitting uses the **sum** of the eligible budgets to finance one common coefficient threshold. It remains effective when all positive coefficients have the same largest budget.

Support scope is essential. For this fixed-q postprocessing step, unchanged outside zero coefficients contribute nothing and may be omitted. If the rule is incorporated into an LP or a multi-step operation that allows such coefficients to become positive without recomputing the module condition, use the full allowed-positive graph, or constrain every omitted coefficient to remain zero. In a sequential algorithm the condition must be rechecked at each current state.

## 5. General independent set: an explicit boundary bill

For arbitrary independent F put

    B_i=N(i) minus [F union common outside neighbors of F],
    Delta_F=max_{i in F} max clique load in G[B_i],
    s_F=q_max-t_F>=0.                                      (8)

All boundary loads are computed from the current q; outside zero coefficients contribute zero even if their graph neighborhoods differ.

### Proposition 3: legal splitting with certified boundary control

The reassignment (6) satisfies

    L(q')<=max(L_rest,t_F+L_good+Delta_F)
           <=max(L_rest,L(q)-s_F+Delta_F).                   (9)

At any point at most one F witness is active. If it is i, partition its active outside neighbors into common neighbors and B_i; their respective loads are bounded by L_good and Delta_F. Outside union F the load is at most L_rest. Also q_max+L_good<=L(q), giving the second inequality.

Thus s_F>=Delta_F certifies nonincrease, with strict improvement when s_F>Delta_F and the unaffected branch is smaller. This is the splitting counterpart of the R035 concentration boundary bill. It is a local certificate only: a uniform high-dimensional bound for the aggregate bills is still missing.

For the four-cube example, F={0,1}, q_max=7p/6, t_F=7p/12, L_good=L_rest=17p/12 and Delta_F=0. Equation (7) gives the exact improvement from 31p/12 to 2p.

## 6. Exact nonedge-pair optimization without a conservative boundary estimate

For any nonadjacent i,j, keep their total payment P=h_iq_i+h_jq_j and every outside coefficient fixed. Define

    A_i=max_{C contains i}(sum_{k in C}q_k-q_i),
    A_j=max_{C contains j}(sum_{k in C}q_k-q_j),
    A_rest=max_{C contains neither i nor j}sum_{k in C}q_k,

using the complete maximal-clique list and the value zero for an empty last class. Because i,j are nonadjacent, no clique contains both. For x=q'_i in [0,P/h_i] and q'_j=(P-h_ix)/h_j the exact objective is

    Phi(x)=max(A_rest,A_i+x,A_j+(P-h_ix)/h_j).              (10)

### Proposition 4: closed-form optimal pair splitting

Let

    x_star=clip_[0,P/h_i]((P+h_j*(A_j-A_i))/(h_i+h_j)),
    q'_i=x_star,   q'_j=(P-h_i*x_star)/h_j.                 (11)

Then Phi(x_star) is the exact minimum over all nonnegative reassignments of this pair with payment P. The increasing and decreasing affine branches in (10) meet at the unclipped value in (11); if their meeting point lies outside the interval, the nearer endpoint minimizes their maximum. The constant A_rest does not change the validity of that minimizer. The old pair is feasible, so

    L(q')=Phi(x_star)<=Phi(q_i)=L(q).                        (12)

All quantities are rational for rational cubes and charges. No new LP, union calculation, or uncertain boundary estimate is needed after the complete clique data are cached. When the outside prices A_i=A_j, (11) reduces to the equal coefficient t=P/(h_i+h_j) of independent-module splitting.

Pair-coordinate optimality need not imply global payment-LP optimality. A future finite experiment may test how much gap repeated exact pair splits capture, but no such computation is run or result claimed in this note. The required reusable data are fixed geometry/allowed graph, complete maximal cliques, h, U, initial q, and verified mass-LP primal/dual bounds, each with source hashes. Recomputing the original geometry or baseline LP is unnecessary.

## 7. Production registration and evidence identities

R037 production's initially reported implementation differs from the written SPEC in tie-order seed domains/tie selection and in duplicate-sign handling for C. Those differences change finite case identities, not the legality of the exact original witness/flow certificates. They must be recorded; raw batch labels must not be called rigorous statistics for an identical pre-registered case set. Only independently audited actual rational instances support exact failure claims.

The all-n four-cube proof above does not depend on the stochastic generator, candidate counts, a floating solver, the 48-cube atomization, or a posteriori case labeling. The original SPEC is retained as the registered experimental target, and this file is explicitly a subsequent mathematical analysis.

Dependencies are R035 payment conservation and box Helly, elementary max-of-affine optimization, and finite mass-LP duality. No production computation or previous three-toy verification is rerun to create this proof. The general-L1/all-scales passage remains the original compact-level-set finite-family interface; a uniform all-family estimate would still be required before any new C_n conclusion.
