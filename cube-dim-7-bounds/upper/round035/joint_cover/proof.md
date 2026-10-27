# R035 — conservation of witness payments and clique compression

## Outcome and fixed scope

This report adds an actual load-reducing operation to the R032 joint cover-flow interface. Cover charges may be reassigned after the flow is constructed, provided their total witness payment is preserved. For a clique module of the witness intersection graph, reassignment to the largest witness budget does not increase the true maximum load; strict decrease requires the conditions in Proposition 2. For an arbitrary clique there is an explicit boundary penalty; a four-interval example makes that penalty sharp and shows that blind compression can increase load.

The nested-cube example works in every dimension: a local own-witness payment scheme has load N/a^n while the compressed load is bounded independently of N. This is an obstruction to the local scheme, not an operator lower bound. The new operation strictly improves the R032 fixed-selected-increment certificate on finite examples. It still does not give a uniform improvement of the n log n coefficient: controlling boundary penalties over general geometric families remains open here.

The underlying operator is unchanged: Lebesgue measure, centered axis-parallel common-radius cubes, every positive scale, arbitrary L1 input. No n=7 special case, lower-bound work, MATLAB session, commit, push, fetch or rebase is used. The only new files are in this directory. The source checkout and hashes are recorded in `verification.json`.

## 1. Payment conservation permits more than spatial covering flows

Consider any finite family Q_i, witnesses V_i=aQ_i, a>1, and an original f,lambda with

    M_i=integral_{V_i}|f|>lambda |V_i|.

Choose budgets |V_i|<=h_i<=M_i/lambda. Suppose a previous valid covering argument, for example R032's flow, has assigned nonnegative charges d_i satisfying

    Σ_i d_i=|U|,        U=union_i Q_i.

Set q_i=d_i/h_i. Then

    integral |f| Σ_i q_i 1_{V_i} >=lambda Σ_i h_i q_i=lambda|U|.   (1)

### Lemma 1: charge-reassignment certificate

Every nonnegative q' satisfying

    Σ_i h_i q'_i=Σ_i h_i q_i=|U|                              (2)

is also legal in (1). Therefore

    lambda |U|/||f||_1 <= esssup_y Σ_i q'_i1_{V_i}(y).          (3)

The reassigned charge need not retain the old point-to-cube flow support or its first-rejection capacities. Those constraints were sufficient for constructing the first charge vector, not necessary for (1). Dropping them is justified by (1)-(2), rather than by assuming that a different covering flow exists. If h_i=M_i/lambda, (1) is an equality for every reassignment.

This is where fixed selected increments can be improved: they become transferable scalar payments, not permanently attached terms in the final load.

## 2. Exact graph description of true load

Let G be the graph with vertices i and an edge ij exactly when the open boxes V_i,V_j intersect. For an axis-parallel open box family, every clique has a positive-volume common intersection: in each coordinate pairwise overlap implies max lower endpoint < min upper endpoint. Thus

    L(q):=esssup Σ_i q_i1_{V_i}
         =max over cliques C of G Σ_{i in C}q_i,       q_i>=0.    (4)

It is enough to use maximal cliques. This is the true pointwise load, not an independent sum of pairwise intersection costs.

The R032 witness-mass LP can likewise use one mass variable for each maximal clique. Any mass supported on a smaller incidence set can be moved to a point in a containing maximal clique; the total mass is unchanged and no witness constraint decreases. Conversely every maximal clique has positive-volume intersection and is a realizable incidence cell: any additional box containing a point of that intersection would enlarge the clique. Hence this graph compression is exact, not a sampled-cell approximation.

In particular, if all witnesses share a common point,

    M(V)=max_i |V_i|.                                          (5)

For such a common-star family let r_* be the largest radius and z a common witness point. Then Q_i is contained in Q(z,(a+1)r_i), giving

    |U|/M(V)<=(1+1/a)^n.                                     (6)

Consequently a large load assembled entirely at one witness star does not by itself identify a difficult mass budget. The old load may reflect many fixed small payments which can all be financed by one large witness mass.

## 3. Exact compression on a clique module

A set C is a clique module if it is a clique and every outside vertex is adjacent either to all of C or to none of C. Equivalently its vertices have identical closed neighborhoods. The property only needs to hold in the support graph of the current positive coefficients; zero-charge outside vertices may be omitted.

Choose k in C with h_k=max_{i in C}h_i. Define

    q'_k=(Σ_{i in C}h_i q_i)/h_k,
    q'_i=0 for i in C minus {k}, and q'_i=q_i outside C.       (7)

Then (2) holds and q'_k<=Σ_{i in C}q_i. Let G_C be the common outside neighbors of C, and put

    L_rest=max clique weight in G outside C,
    L_good=max clique weight in G_C,
    T_C=Σ_{i in C}q_i,       T'_C=q'_k.

Empty graph maxima are zero.

### Proposition 2: exact true-load improvement

For a clique module,

    L(q)=max(L_rest,T_C+L_good),
    L(q')=max(L_rest,T'_C+L_good)<=L(q).                       (8)

Indeed a maximal clique containing any C vertex can be extended to contain all of C, while every common-neighbor clique extends with C. After compression the same statement holds with the single anchor k. The scalar decrease in the affected branch is

    s_C=T_C-T'_C=Σ_{i in C}q_i(1-h_i/h_k)>=0.                 (9)

If that branch determines the old maximum, h_i vary on positive charge, and the unaffected branch is strictly smaller, the improvement is strict. This is a rigorously monotone postprocessing operation on any existing joint cover-load certificate, including R032's coefficients (v_i+c_i)/h_i.

All disjoint true-twin classes may be compressed simultaneously. Their quotient graph has the same outside adjacencies, and the payment of each class is retained at its largest-budget representative. An LP minimizing the combined flow load can incorporate this operation by replacing every class's separate loads with its total charge divided by its representative budget. For an LP, the module property must hold in the full graph of vertices whose coefficients are allowed to become positive; a module computed only in a current support graph is valid only if every omitted coefficient is constrained to remain zero. With that support restriction, it is no weaker than the uncompressed flow LP and requires no change to the original selection or the general-L1 argument.

## 4. Arbitrary clique: savings versus a computable boundary bill

The module hypothesis cannot simply be removed. For any clique C and anchor k with largest h define

    B(C,k)={j outside C: j adjacent k but not adjacent to all of C},
    Delta(C,k)=max clique weight in the induced graph on B(C,k). (10)

Retain L_rest and L_good as above, where G_C now means outside vertices adjacent to every member of C. The same reassignment (7) is legal by payment conservation, but its load satisfies

    L(q')<=max(L_rest,T'_C+L_good+Delta(C,k))
           <=max(L_rest,L(q)-s_C+Delta(C,k)).                  (11)

Proof: at a point outside V_k only outside-C terms remain, bounded by L_rest. At a point inside V_k partition the active outside vertices into common neighbors and boundary vertices. Each part is a clique. The first costs at most L_good and the second at most Delta. Also C together with a common-neighbor clique is a clique, so T_C+L_good<=L(q).

In particular, s_C>=Delta is a sufficient condition for nonincrease. If s_C>Delta and the unaffected load is smaller, (11) gives a strict reduction. No independent coverage coefficient is multiplied into this saving: it is a direct improvement of the combined witness-payment load in (3).

This is a two-part bill: internal scalar payment savings pay for changed external incidences. A general high-dimensional theorem would need to find enough cliques with positive net savings, or sum their boundary bills in a scale-local potential. This report supplies the exact local inequality but no global bill bound.

## 5. An all-dimensional obstruction to local own-witness payments

Fix n>=1, a>2, N>=1, and for i=0,...,N-1 set

    r_i=4^(-i),       y_i=(2r_i,0,...,0),
    Q_i=Q(y_i,r_i),   V_i=aQ_i.

The Q_i are disjoint in the first coordinate: the next upper endpoint 3r_i/4 is less than the previous lower endpoint r_i. The V_i are nested, since each coordinate interval is a fixed interval containing zero multiplied by r_i. They have a common positive-volume intersection for every finite N.

If covering volume may only be paid to its own Q_i, each Q_i must pay its entire volume. Thus

    d_i=|Q_i|,  h_i=|V_i|,  q_i=a^(-n),
    L_own=N/a^n.                                             (12)

This is the unique own-witness flow, not merely one poor feasible flow. Yet the clique-module compression transfers all charges to V_0 and gives

    L_compressed=|U|/|V_0|
      =[1-4^(-nN)]/[a^n(1-4^(-n))].                         (13)

Equation (5) shows that (13) is exactly the geometric witness-mass LP certificate. The ratio L_own/L_compressed grows linearly in N. If a threshold is chosen large enough to select all these disjoint small cubes, R032's fixed-selected-increment flow has D empty and is exactly (12); it too is strictly improved by compression.

This example does not produce an operator lower bound: every witness is financed by the largest witness mass and the correct certificate stays bounded. It shows why a spatially constrained payment rule or permanently fixed selected increments may be unsuitable for proving a uniform upper bound.

## 6. A nontrivial module and a sharp boundary obstruction

The following intervals are witness boxes in dimension one, with a=3 and Q_i obtained by dividing each witness radius by three. Use a nonincreasing-radius increment order; the resulting charges are as stated.

    V_0=(-1,5),      h_0=6,
    V_1=(-1/4,5/4),  h_1=3/2,
    V_2=(1,6),       h_2=5,
    V_3=(11/2,8),    h_3=5/2.

The small-cube union has volume 14/3. Its increment charges give

    q=(1/3,1/3,4/15,1/3),    Σh_iq_i=14/3.

C={0,1} is a clique module: vertex 2 meets both and vertex 3 meets neither. The old true maximum is 14/15. Compression gives q'_0=5/12, q'_1=0 and the improved true maximum

    L(q')=41/60 <14/15=L(q).                               (14)

The whole graph is not a clique; another maximal clique {2,3} remains. Thus this is more than the common-star special case.

Now move only V_3 to (9/2,7), preserving its volume, its small-cube volume, all four charge coefficients, and total coverage volume. It meets V_0 and V_2 but misses V_1. The module property fails. For the same C and anchor,

    s_C=1/4,  Delta=1/3,
    L(q)=14/15,     L(q')=61/60=L(q)+Delta-s_C.             (15)

Thus the boundary penalty in (11) can be attained exactly. A clique alone does not permit monotone compression. Ignoring a small number of changed incidences can reverse the improvement.

To realize both examples by actual witnesses take f=(1+epsilon)1_{union V_i}, lambda=1. It is L1 and all witness margins equal 1+epsilon. Every certificate and charge reassignment above then applies to a single consistent original input. Alternatively use h_i=(1+epsilon)|V_i| and divide all q by 1+epsilon; the joint integral in (1) is then exactly |U|.

## 7. Evidence, reproducibility, and what is still missing

`verify.py` checks only three rational families: one two-dimensional nested family with N=6, the nontrivial one-dimensional module, and its sharp-boundary modification. It computes cube-union volumes from coordinate cells, enumerates maximal cliques, solves the compressed graph mass LP, and checks all displayed charges, payments, clique loads, module decisions and boundary identities as rational equalities. LP solutions are rationalized and verified against both primal and dual constraints. This is a small interface diagnostic, with no parameter scan.

The new actionable interface is post-flow payment conservation plus clique compression, with exact load recomputation or the certified boundary estimate (11). It can reduce a valid old product or R032 joint certificate before integrating against the original f. No finite diagnostic is presented as a uniform C_n bound.

The remaining general-dimension obligations are concrete: identify geometrically forced compressible or approximately compressible witness cliques; control the aggregate Delta bills without another independent worst-case factor; and obtain a uniform saving of order n log n rather than merely an additive cap improvement. Families with no useful clique modules are not excluded by cube geometry. No claim is made that the local operation alone supplies those obligations.

Dependencies: R032 payment identity and original-instance mapping, source `paper/main.tex` level-set regularity/finite witnesses, elementary coordinate Helly for boxes (proved above), and finite LP duality. Source hashes distinguish this report from the earlier soft-mixture route; that route is not rescanned.
