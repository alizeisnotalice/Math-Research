# R037-JP-v1: registered finite joint-payment experiment

Status: authorized by the main coordinator; production belongs to Terra. This document registers the model, finite cases, exact acceptance rule and stopping counts. The author of this SPEC does not execute the production batch. This is not a C_n bound or an operator lower-bound search. No dimension-seven specialization, MATLAB, or Git write operation is involved.

## 1. Single finite hypothesis and scientific target

For each eligible finite case, let L0 be the **true** essential maximum of a real legacy joint-payment load. Let Lopt=|U|/M(V), where M(V) is the exact witness-mass LP below, with budgets fixed throughout to h_i=|V_i|=|aQ_i|. Let Lalg be the final load of the deterministic R035 boundary-safe compression algorithm D.

The single hypothesis H is:

    if the witness graph is not one clique and
       L0-Lopt >= 10^(-4) L0,
    then Lalg <= (L0+Lopt)/2.                                 (H)

H concerns this explicitly bounded candidate algorithm on this finite batch. A failure is not a theorem that every possible clique move fails, and success is not a uniform upper bound. The purpose is to identify whether the certified boundary bill captures a useful fraction of the available *payment-certificate* improvement, or whether its failure is primarily a boundary, support, candidate-pool, or genuinely absent-LP-gap obstruction.

## 2. Fixed finite batch and process ownership

There are 15 worker strata: three geometry templates A/B/C times five parameter strata. Each worker handles exactly 16 geometry seeds, serially, with case index j=0,...,15. Total: 240 geometries. A and B each have two legacy payment patterns, full-selected and partial first-reject; C has full-selected only. Nominal total: 400 load patterns. Censored/unavailable cases remain in the ledger and are not silently replaced or replenished.

The parameter strata (n,a,N,L) are:

    S0=(8,4,48,6)
    S1=(16,8,64,8)
    S2=(32,16,80,8)
    S3=(64,32,96,12)
    S4=(32,8,96,12).

Set BLAS/OpenMP threads per process to one; use 15 processes rather than 15 independently multithreaded solvers. No wall-clock cutoff is prescribed. There are finite node, case, generator and algorithm-iteration budgets below; no automatic extension of the grid or memory-filling work.

Seeds are unique per geometry and per stochastic subdomain:

    SHA256('R037-JP-v1|template|stratum|case-index|domain').

Use the full digest as a reproducible integer seed and record it. Domain strings include `geometry`, `retry/k`, and `tie/o/radius/exact-radius`. No geometry seed is reused across parameter strata. PRNG implementation/version is recorded. All generated coordinates/radii are rational; random draws only choose integers from finite sets.

For case j define

    delta=[1/64,1/16,1/8,1/4][j mod 4],
    contact=[1/4,1/2,1,2][floor(j/4)],
    D=a*contact,       r_l=4^(-l),  l=0,...,L-1.

Original cubes are always Q_i=Q(center_i,radius_i), and witnesses V_i=aQ_i. Budgets never use margin-dependent h in this batch.

## 3. Exact geometry generators

### A. Four nested branches with annular contact probes

All centers have only their first two coordinates nonzero. Let u range over e1,e2,-e1,-e2; let v be u rotated counterclockwise by 90 degrees; put z_u=D*u.

Create 4L core cubes, one for each (u,l):

    center=z_u+2*r_l*u,     radius=r_l.

There are 8L distinct probe slots (u,l,s), s in {-1,+1}. Shuffle these slots by the geometry seed and take the first N-4L. For each chosen slot create

    center=z_u+(2-a+2*delta)*r_l*u+s*delta*r_l*v/2,
    radius=delta*r_l/a.

These probes are at the negative annular face of the corresponding witness branch, so they can distinguish nested contacts while carrying actual small-cube volume. No common-star assertion is assumed; the exact resulting graph decides the structure.

### B. Two stars with selective bridges

Again only two active coordinates. For s in {-1,+1} let z_s=(sD,0). Create 4L core cubes for (s,u,l), u in {e2,-e2}:

    center=z_s+2*r_l*u,     radius=r_l.

There are 6L bridge slots (l,t,s), t in {1/4,1/2,3/4}, s in {-1,+1}. Shuffle by the geometry seed and take the first N-4L. Each gives

    center=((1-2*t)*D, s*(2+delta/4)*r_l),
    radius=r_l.

Large-radius bridges can join star histories; smaller ones may remain separate. Whether a bridge is useful, isolated or creates a boundary obstruction is measured rather than imposed.

### C. Full-dimensional coded near-contact witnesses

This template has no active-coordinate reduction. It is admitted only if the original small Q_i are pairwise disjoint, permitting exact increments and union volume without a high-dimensional union approximation.

For vertex i=0,...,N-1 define

    l=i mod L,
    cluster=floor(i/L) mod 2,
    z_0=0, z_1=D*e1,
    nu=[1/n,1/4,1/2,3/4][floor(j/4)].

For each cluster draw a base sign vector in {-1,+1}^n. For each vertex independently flip each base bit with probability 1/8, implemented by an exact integer draw `randbelow(8)==0`. A sign vector must not duplicate an earlier vector at the same (l,cluster); redraw its bit choices deterministically from that geometry stream until it differs (at most 256 attempts for that vertex, otherwise reject this family candidate).

For each coordinate independently draw an outer mask with rational probability nu, using exact integer sampling. Put

    theta_k=sign_k*(1+delta) if outer, else sign_k*(1-delta),
    center_i=z_cluster+a*r_l*theta,
    radius_i=r_l.

Check every pair of original Q_i with strict open-box overlap. If any overlap, regenerate the whole family using the next `retry/k` seed, up to 32 family candidates. If none is pairwise disjoint, record `geometry_unavailable`; do not replace the case with another seed. Accepted C cases have v_i=|Q_i|, |U|=sum_i |Q_i| and full-selected q_i=a^(-n) exactly.

## 4. Exact union and increment algorithm for A/B

No Monte Carlo or sampled union volume is allowed. For any needed collection of Q/E boxes all inactive centers are zero. Let m=n-2. Form radius cuts {0} plus all inactive half-widths R of the boxes in that collection. On a slab

    rho_lo < ||x_inactive||_infinity < rho_hi,

the active boxes are precisely those with R>=rho_hi. Its inactive volume is exactly

    W=(2*rho_hi)^m-(2*rho_lo)^m.

Partition the first two coordinates by all active rectangle endpoints. Each open elementary 2D cell has constant Q/E membership. Multiply its exact rational area by W and aggregate by the membership bitmasks. Endpoints have zero Lebesgue measure. Streaming/sweep-line implementations and exact bitmask aggregation are permitted; the resulting rational atoms must be identical to this definition.

For Q-only atoms, |U| is their total positive-membership volume. Given an order, an atom's entire volume is assigned to its earliest selected member, giving the exact increments. During partial selection compute v_i=|Q_i minus union(previous selected Q)| from the same atoms. For D/E use the same construction on all input Q and final selected E_j. There is no need to atomize V: true witness loads use the complete maximal-clique list.

An elementary positive slab-cell budget of 8,000,000 is imposed per geometry/required atomization. Exceeding it records `geometry_censored`, with counts, rather than substituting approximate volume. No large array need be filled solely to reach this bound.

## 5. Two real legacy payment patterns

For A/B generate up to eight tie orders. Radii are always nonincreasing. Within each equal-radius class use an independent `tie/o/radius/exact-radius` domain, o=0,...,7. Deduplicate identical orders and identical resulting increment vectors. Using the complete true-load calculation below, select the order with largest full-selected L0; ties choose lexicographically smallest index order. This is an adversarial choice among legal same-radius histories, not an arbitrary payment vector.

### Full-selected

Use p_i=v_i/h_i with h_i=(2*a*r_i)^n and tau=N/a^n. Since every p_i<=a^(-n), all inputs are eligible regardless of center tests. Thus every chosen tie order is a legal largest-radius selection and D is empty. The joint load is q_i=p_i, and sum h_i q_i=|U| exactly. Template C uses this pattern only; its order does not affect disjoint increments.

### Partial first-reject (A/B only)

Use the chosen legal tie order within equal radii. Set

    K=[1/2,1,3/2,2][j mod 4],
    tau=K/a^n,       b=a+3,       E_j=(a+4)Q_j.

At each step a candidate i is eligible if

    sum over already selected j of p_j*1_{bQ_j}(center_i)<=tau.

Choose a largest-radius eligible cube, breaking equal radii by the fixed tie order. Define its increment and p_j=v_j/h_j once; do not later update it. Stop when no candidate is eligible.

Let S be the selected union and D=U minus S. Aggregate D by selected E incidence sigma, volume d_sigma. For every positive atom check exactly

    P_sigma=sum_{j in sigma}p_j>tau.

Failure is a geometry/selection-proof bug with an exact witness, not an ordinary negative scientific result. The explicit legal proportional covering flow is

    z_{sigma,j}=d_sigma*p_j/P_sigma,
    c_j=sum_sigma z_{sigma,j},
    q_j=(v_j+c_j)/h_j for selected j, and q_i=0 otherwise.

Check sum_i h_i q_i=|U| as an exact Fraction equality. This baseline retains true rejected volume and first-rejection capacities. It is a specific feasible joint flow, not a claim to have solved R032's flow LP optimally.

All original witnesses are simultaneously realizable by the common L1 input f=(11/10)*1_{union all V_i}, lambda=1. This gives M_i=(11/10)h_i. Selection, load, and the global LP are all on this same original finite family.

## 6. Complete witness graph and mass-LP oracle

The allowed graph is fixed to **all original V_i**, including initially unselected vertices with q_i=0. Reassignment may reactivate them; module tests therefore use this full allowed graph. A pair is adjacent iff all coordinates have strictly overlapping open intervals. Finite box Helly makes every clique a positive-volume common witness intersection.

Use deterministic-pivot Bron--Kerbosch with deterministic index tie breaking. Per geometry stop at either 2^20 recursive nodes or 65,536 output maximal cliques. Reaching either cap records `graph_censored`. An incomplete maximal-clique list must not be called a true maximum load or used as a complete LP certificate.

With the complete maximal-clique list C, calculate true load exactly by

    L(q)=max_{C} sum_{i in C}q_i.

For h_i=|V_i| define

    M(V)=min sum_C m_C,
         m_C>=0, sum_{C contains i}m_C>=h_i,

        =max sum_i h_i z_i,
         z_i>=0, sum_{i in C}z_i<=1 for every maximal clique C.

The optimum load at fixed total payment |U| is exactly Lopt=|U|/M(V). This identity is not used with a different h budget. It may activate any original witness.

Solve numerically after scaling h by max h, then rationalize and verify. A permissible repair is: clamp/rationalize nonnegative dual z and divide by max(1,max_C sum_C z), obtaining an exactly feasible dual; clamp/rationalize primal m and add each exact uncovered h_i residual to the lexicographically first maximal clique containing i. This obtains exact Mlo<=M(V)<=Mhi. Record both repairs and verify all inequalities with rational arithmetic.

If the oracle interval is too wide for H, at most one deterministic higher-precision/rationalization or alternate primal/dual solve is allowed; otherwise record `oracle_unresolved`. No LP failure or repaired wide interval is silently counted as H failure or success. Suggested objective interval target is (Mhi-Mlo)/Mhi<=2^(-20). Iteration/resource solver limits, if used, are recorded as resource censoring rather than wall-clock failure.

## 7. Deterministic compression algorithm D

Keep h and the full allowed graph fixed. Let q be the current exact rational coefficient vector. For any clique C of size at least two define

    anchor k = smallest index attaining max_{i in C}h_i,
    T_C=sum_C q_i,
    q'_k=sum_C h_i q_i/h_k, q'_i=0 for i in C minus {k},
    q'_i=q_i outside C,
    s_C=T_C-q'_k.

One iteration constructs this bounded pool:

1. Take the first 32 maximal cliques ordered by (-T_C,Ctuple), and first 32 ordered by (-s_C,-T_C,Ctuple); their union is the basis.
2. From every basis clique include itself and both intersections with {h_i<=t} and {h_i>=t} for every distinct original budget t. Remove duplicate subsets and size<2 subsets.
3. Retain the first 64 resulting ordinary candidates ordered by (-s_C,-T_C,Ctuple). Add every full-graph true-twin class of size at least two. There are at most 48 such disjoint classes for N<=96, so the pool has at most 112 candidates.

For each pool candidate with s_C>0 form

    B(C,k)=N(k) minus [C union common outside neighbors of C],
    Delta=max_{D in maximal cliques}sum_{i in D intersect B(C,k)}q_i,
    Lrest=max_{D in maximal cliques}sum_{i in D minus C}q_i.

This intersection formula is the exact induced-boundary clique load: every boundary clique lies in a full maximal clique. Candidate is **boundary-safe** iff s_C>=Delta. Construct q' and compute true L(q') from the complete list. Among all safe candidates select the lexicographically smallest tuple

    ( L(q'), sum_i q'_i, -s_C, Ctuple, anchor ).

Accept it, verifying rationally payment conservation, Delta, s, L(q')<=max(Lrest,L(q)-s+Delta)<=L(q). Zero coefficients may be reactivated only via anchors inside such a candidate clique; the full allowed graph never changes. Numerical pre-ranking/vectorization is permitted, but the stated pool and safe-candidate comparisons must ultimately have the same exact-rational result; do not discard an eligible candidate merely because of a floating tolerance.

Accept at most 24 moves per pattern. Stop earlier if no safe pool candidate exists. A plateau in L is allowed if s>0 because total sum q strictly decreases; all saved states and comparisons remain exact. Record `iteration_cap` or `stop_candidate_pool` distinctly.

At a stop with substantial residual LP gap, do one diagnostic pass over the **same final pool**: test whether any unsafe s<Delta candidate nevertheless gives a strictly smaller true load L(q'). If yes, record `boundary_bound_conservative` and its exact candidate; do not iterate a second algorithm. Also record whether any full-graph true-twin module exists, whether the stalled baseline involves initially zero anchors, and whether every accepted improvement used a maximal clique or a proper radius-threshold subclique.

## 8. Exact H classification and output certificate

Let

    Klo=|U|/Mhi,   Khi=|U|/Mlo,
    Klo<=Lopt<=Khi.

An eligible-gap case is certified by L0-Khi>=10^(-4)L0. A negligible-gap case is certified by L0-Klo<10^(-4)L0. Otherwise classify `gap_unresolved` or perform the one allowed oracle refinement.

For eligible cases:

    PASS if 2*Lalg<=L0+Klo,
    FAIL if 2*Lalg>L0+Khi,
    otherwise oracle/H unresolved.

All comparisons use exact rationals. Whole-star graphs and negligible-gap cases are controls outside H's antecedent. They remain useful diagnostics and are not removed from counts. Per-pattern compression stop budgets remain attached to any FAIL: it is a failure of this registered algorithm, not of all possible payment optimization.

Each output case includes: seed/domain and all rational centers/radii; complete graph adjacency and maximal cliques or a resource-censor reason; exact U/increments/partial proportional flows; full allowed h; initial and final rational q; every accepted (C,k,s,Delta,L,L') and exact conserved payment; Mlo/Mhi primal/dual certificates; classification and failure subtype. Normalized L0=1 numerical working units are allowed, but raw rational identities and exact ratios must be restored in the certificate.

Batch summary reports all 240 geometry statuses and all nominal 400 pattern statuses, separates complete from censored cases, and lists the worst certified improvement fraction and failure patterns. Retain only minimal reproducible rational candidates plus aggregate statistics; do not archive every redundant intermediate numerical array. No finite result is labeled a new C_n constant.
