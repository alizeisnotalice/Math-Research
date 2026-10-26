# Sparse4 joint scalar preparation and strict release boundaries

This is the new frame body32+(1/80)sparse4, not the original body32 model. The152 unoriented directions comprise32 full-support cube lines and120 four-support lines (15 supports,8 sign representatives per support). Their common weights are invariant under all signed coordinate permutations. In particular the first-coordinate marginal rate is32+80/80=33. The proposed exact time32/825 satisfies33T=32/25; it is a marginal normalization, not an optimized time. **No positive time is substituted or physically evaluated in this preparation.**

## Fan adapter

Read only the frozen050 global_fan.json.gz after verifying both compressed and decompressed hashes. The full fan has478 indexed primitive integer rays and13229 simplex entries, not21 entries. The adapter checks all records, positive rates, six distinct ray IDs per simplex and every six-dimensional absolute determinant. This is structural adapter verification; it does not independently prove fan completeness, which remains gated on F050.

For ray r with q=p/d>0 in lowest terms, scale its whole column by d:

    R=d r, Q=p,  R/Q=r/q.

For a simplex with primitive determinant D, the scaled determinant isD'=D∏d_j. Consequently

    D'/∏Q_j=D/∏q_j.

Every actual entry was checked with exact arithmetic. The true rate was also rebuilt for each ray as the sum of152 weighted absolute integer forms, rather than trusting the stored denominator. Actual maxima are d≤40, |R_i|≤520, Q≤20773, D≤32 and D'≤2048000000. For the proposed point numerator bound50 and denominator50, a safe node numerator bound is156000 and denominator bound1038650. All fit int64. The group weight24 times D' exceeds32-bit range in the worst case; use int64 or Arb, never int32. Products of six Q values need not fit int64 and must be arbitrary precision or sequential Arb divisions.

The canonical complete converted-data digest, including every scaled ray, simplex and orbit representative, is
`bbe9274e5b235d5c7f94b8c226ca124bd2c1f778ba012a124440cdf711cd47b3`.
The complete ordered weight-stream digest is
`77f1584527cf64989a5cfdd62aeea08e51dc0cf7883462bda31cf789d44390f2`.
The derived large input is not duplicated in the archive: prepare_adapter.py reproduces it from the pinned external source. A future numerical driver must read all13229 entries with a count/stream gate and cannot invoke the old21-cone parser.

## Independent orbit multiplicities and global negation

For x0=(50,50,50,45,27,0)/50, a full exact loop through all720 permutations and64 signs was performed without any kernel. It yielded3840 distinct point images, each with multiplicity12. This also has an independent counting derivation:6 choices of the zero position,10 choices of the three repeated50 positions,2 arrangements of45/27, and32 nonzero-coordinate sign choices yield3840. The stabilizer is3! from equal50 coordinates times2 from the zero-coordinate sign, hence12.

No image is its own global negative. Pairing the3840 images yields1920 representatives, with exact combined weight12+12=24 per representative. No additional factor12 or2 may be applied later. With13229 simplices,

    actual point-orbit/cone terms=13229·1920=25399680,
    logical group/cone terms=13229·46080=609592320.

The complete group action was checked through a Counter, including every image multiplicity and its negative partner; these are not inferred by blindly dividing a term count.

For each cone, DD5 of the odd function atan(b/T) is even under simultaneous negation of all nodes: the function contributes−1 and order5 contributes(−1)^5. Thus H is even term by term. The full epsilon-dependent DD5 of atan(b/epsilon) has the same property. Every Laurent coefficient, including its constant, is therefore also even. This proves the negative pairing for both prepared scalar functionals, including zero/repeated nodes, before any physical finite-part interpretation.

## Scalar mathematics and compiled-core boundary

H uses ordinary Hermite DD5 of atan(b/T), with true normalized jets through order5. At zero, its odd normalized jets are(-1)^((k−1)/2)/(k T^k), not zero. The source reuses C010's already validated analytic-jet algorithm, with a persistent DD class available to a future driver.

The G candidate is **only the formal epsilon^0 coefficient** of DD5 atan(b/epsilon). At a nonzero node its value jet has constant(π/2)sign(b), all higher normalized jets have constant0; at zero every jet has constant0. The omitted negative powers are retained in the formal tests, not declared absent. For example an all-zero tuple has formal constant0 but a nonzero epsilon^-5 term: it is not a finite physical Green value.

Eight exact rational node multisets each provide two FORMAL checks (16 total): H's symbolic common-T derivative agrees with minus an independently expanded real resolvent product; the G Laurent vector (constant plus epsilon^-1,-3,-5) agrees between triangular DD and independent exact Hermite polynomial weights. Global-negative parity of the entire vector is checked. T and epsilon remain indeterminates; no value32/825 or other positive time is sampled.

The new C++ source compiles the Arb512 joint scalar algorithms and exports a function explicitly named FORMAL_constant_and_H. The executable only permits ABI mode; there is deliberately no physical scalar or full-frame CLI in this phase. One build and one ABI invocation passed. No compiled numeric H/G anchor was executed; the16 checks are formal algorithm checks, not a numerical new-frame or new-time certificate. A future full driver must keep its node cache across terms rather than use the scalar wrapper as a per-term allocator. An unbounded tuple cache of25million96-byte keys would be inappropriate; the simple ray/image node upper count is478·1920=917760 before rational coincidences, so resource accounting must be explicit.

## New true-S2 eligibility: separate gate

The sparse4 directions create a new S2. The body32 point's old eligibility is irrelevant. Before any physical G interpretation, separately register and execute exact_eligibility in s2_preflight_plan.py. It constructs all152 unique unoriented directions, checks every line and all11476 unordered pairs using exact rational Gram arithmetic, and rejects any zero distance. For pair(v,w), with a=v·v,b=v·w,c=w·w, xv=x·v,xw=x·w,

    dist(x,span(v,w))²
       =|x|²−(c xv²−2b xv xw+a xw²)/(ac−b²).

Distinct unoriented directions give ac−b²>0. Positive beta scaling does not change their spans. Reject a nonpositive Gram determinant, negative computed distance, missing pair, or any exact zero distance; no binary64 rank tolerance is permitted. The check has NOT been executed here. Its CLI explicitly refuses a run without a new registered stage. If the point lies on any new S2 plane, do not return the formal constant as finite G; stop and report ineligibility.

## Required future release sequence

Before a new full-frame numerical invocation, all of the following must be pinned and accepted: F050's complete fan audit; B052's finite-part theorem for this152-direction frame; the new exact off-true-S2 point preflight; scalar/native integration acceptance; and a separate main-agent release/finite registration. New full outer-driver code and any necessary build belong to that later release, not to this preparation's one compilation. Only then may one joint G/H full value at the single registered x0,T be attempted. K=G−H is merely a necessary point probe, never a global mT or Cn bound. No automatic release occurs when a dependency file appears.
