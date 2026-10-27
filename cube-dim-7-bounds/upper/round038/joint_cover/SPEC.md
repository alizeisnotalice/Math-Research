# R038 Phase A: frozen actual-cache nonedge-pair descent

Status: main coordinator authorized this phase. Production is owned by Terra. The SPEC freezes the mathematical rule and execution scope; actual production source hashes must be frozen before dispatch. Phase B/new geometry generation is not authorized by this document.

## 1. Authorized scope and identity

Input is exactly the manifest of **79 cached actual R037 failure patterns**, starting from each saved `q_before`. This is not a new realization of the originally written R037 seed registration. R037 implementation differences are already documented; the actual geometry, original-run generator hash and pattern IDs determine case identity.

Use 15 workers in a dynamic per-pattern queue. Every manifest item is attempted once and receives a terminal status. Do not regenerate 240 geometries or refill any missing items. Only missing geometry, h and complete maximal cliques may be restored, using the frozen **actual R037 implementation**, with its hash and deterministic environment recorded. Do not run original union atomization, selection, baseline payment construction, or baseline LP again. Saved U and q_before are the old production data, not newly measured quantities.

The baseline budgets are fixed to h_i=|aQ_i|>0. The allowed graph includes every original V_i=aQ_i, including vertices whose current coefficient is zero. All original witnesses share a legal common input as in R037. No dimension-seven case, MATLAB, lower-bound work, new stochastic geometry, or Git write operation is involved.

Before dispatch freeze: the ordered manifest/hash, actual-case cache hashes, all production code hashes, Python/library versions, and this SPEC hash. Changed code is a new version, not a silent continuation of the same run.

## 2. Admission checks

For each item restore or read geometry/h/adjacency and a **complete** maximal-clique list. If enumeration must be restored, retain the R037 2^20-node/65,536-clique limits; a cap is censored and must not produce a true-load claim. Verify all stored graph/clique identities against the frozen actual data and record any missing provenance.

Check rationally q_before>=0, h_i>0 and sum_i h_i q_before_i=U. Compute the initial true load from the complete clique list and match the cached L0 exactly. Missing/mismatched data is a cache failure, not a scientific descent failure. Do not repair it by rebuilding the original atoms/selection/LP within Phase A.

Original scalar optimum bounds are **historical diagnostics only**. They are not new certified oracle values for this phase; no global-optimality or half-gap certification is emitted from a scalar bound alone. At most three final candidates receive a new complete rational primal/dual LP certificate under section 7.

## 3. Fixed objective and pair feasible region

For the current nonnegative q define

    L(q)=max_{C in complete maximal cliques}sum_{k in C}q_k,
    E(q)=sum_i h_i q_i^2.

The mathematical objective is lexicographic (L,E): load first, energy second. The energy only resolves a load plateau; an energy increase is allowed when load strictly decreases. This is a payment-certificate objective, not a C_n value.

For a nonedge i<j, let P=h_iq_i+h_jq_j. Outside coefficients remain fixed. If P=0 the pair is unchanged and may be skipped. Let

    A_i=max_{C contains i}(sum_Cq-q_i),
    A_j=max_{C contains j}(sum_Cq-q_j),
    A_rest=max_{C contains neither i nor j}sum_Cq,

where the last empty maximum is zero. No clique contains both nonadjacent vertices. Therefore for x=q'_i, y=q'_j=(P-h_i x)/h_j, the exact pair objective is

    Phi(x)=max(A_rest,A_i+x,A_j+(P-h_i x)/h_j),
    0<=x<=P/h_i.                                             (1)

Recompute these current prices from the complete list at every accepted state. There is no support deletion: zero coefficients can become positive through a pair.

## 4. Deterministic root, plateau and tie rules

All following calculations are exact Fractions.

First balance the increasing/decreasing affine branches:

    x0=clip_[0,P/h_i]((P+h_j*(A_j-A_i))/(h_i+h_j)),
    m=Phi(x0).

The full minimum-load interval is

    lower=max(0,(P-h_j*(m-A_j))/h_i),
    upper=min(P/h_i,m-A_i).                                  (2)

It is nonempty. On this interval the pair energy h_i x^2+h_j y^2 has its unique minimum at the clipping of the equal-coefficient point:

    x=clip_[lower,upper](P/(h_i+h_j)),
    y=(P-h_i x)/h_j.                                         (3)

Thus (3) minimizes the pair's lexicographic (true load,global energy). Indeed Phi<=m gives exactly (2), and the pair-energy derivative is

    (2h_i/h_j)*((h_i+h_j)x-P).

Strict convexity proves the plateau tie rule. Outside energy is fixed. No numerical polynomial root or floating tie tolerance is used.

This is the formula currently implemented in Terra's `pair_balance.py`. Its candidate key `(newL,deltaE,i,j)` is equivalent to `(newL,newGlobalE,i,j)` since the current global E is common to all candidates. The admission test `(newL,deltaE)<(currentL,0)` is likewise the exact strict lexicographic improvement test.

## 5. Full-pair deterministic descent and stopping

At each accepted state enumerate **every nonedge pair** i<j (at most 4560 for N<=96), omitting only P=0 pairs which cannot change the state. No heuristic candidate pool, random subsampling, or radius cutoff is used.

For each pair compute (1)-(3). Among candidates whose `(newL,deltaE)<(currentL,0)`, select the lexicographically smallest `(newL,deltaE,i,j)`. Update only i,j. Verify rationally:

    x,y>=0,
    h_i x+h_j y=P,
    L(new q)=newL,
    E(new q)=oldE+deltaE,
    (newL,newE)<(oldL,oldE).

Recompute current clique prices for the next step. Strict lexicographic descent does not prove finite termination or global payment optimality; the execution has a fixed maximum of **128 accepted moves** per pattern.

Return `stop_pair_coordinate` only when a complete scan of all nonedges finds no strict lexicographic improvement. This certifies stationarity for these pair moves, not global LP optimality. If 128 moves are accepted, return `iteration_cap`, even if the next unperformed scan might have stopped. There is no wall-clock cap. Internal exact-invariant failures stop the affected worker and notify the coordinator; they are not ordinary scientific failures.

Clique merging is not part of the Phase A descent rule. A terminal old-R037-pool merge probe may be prepared only as an explicitly separate diagnostic and must not be inserted into the primary sequence or labeled pair descent.

## 6. Output and scientific interpretation

For every input item retain manifest/cache/code identity, initial and final rational q, U/h and complete clique provenance, initial/final L and E, accepted pair steps with old/new values and conserved pair payments, pair-evaluation count, terminal status, and worker/runtime receipts. One cache may serve multiple patterns; do not duplicate large geometry indiscriminately.

Report exact true-load reduction and the number of load-decreasing versus load-plateau/energy-decreasing moves. Historical half-gap indicators may be computed from old scalar bounds but must be labeled `historical_oracle_diagnostic`, not a newly certified PASS/FAIL. Resource censors, missing caches, invariant errors, iteration caps and coordinate stationarity remain distinct.

A remaining gap at a pair-coordinate stationary state would motivate a multi-vertex payment direction; it does not prove a lower bound for C_n. A large finite reduction does not prove a uniform-dimensional upper bound.

## 7. At most three new LP certifications

Only after the 79-pattern descent phase, select at most three distinct complete candidates deterministically, without rerunning geometry/atoms/selection:

1. The coordinate-stationary item with the largest positive residual relative to its historical optimum upper bound, ranked by `(Lfinal-Khi_historical)/Lfinal`; ties use the manifest order. Omit if none qualifies.
2. The remaining item with greatest exact relative true-load reduction `(L0-Lfinal)/L0`; ties use manifest order. Omit if no positive reduction remains.
3. The remaining item with largest positive residual relative to the historical optimum upper bound, using the same residual ratio regardless of coordinate/iteration terminal status. Omit if none qualifies.

Historical bounds only choose a candidate; they do not certify its quality. Each selected item gets the full mass LP on all original witnesses, h=|V_i| and complete cliques:

    min sum_C m_C,   m_C>=0, sum_{C contains i}m_C>=h_i,
    max sum_i h_i z_i, z_i>=0, sum_{i in C}z_i<=1.

Store actual rational feasible m and z, independently verify every constraint and both objectives, and provide Mlo<=M<=Mhi. Record any floating proposal/rational repair explicitly. No scalar-only claim is enough. With Klo=U/Mhi<=Lopt<=Khi=U/Mlo, a fresh half-gap statement is certified only if its exact inequalities have the appropriate worst-endpoint direction. Unresolved intervals stay unresolved.

This cap is three new LP certifications total for Phase A, not three per worker or template. No new geometry Phase B starts without a separate authorization.
