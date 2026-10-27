# R015 independent true-G point audit

2026-10-27. **Accepted, within the stated dependencies and evidence boundary.** This audit independently checks the mathematics, actual called source path and exact saved-interval accounting. It does not repeat the scientific phase/DD run. The point is the original exact R003 index2710 dyadic point, not a newly searched or rounded point. Parent retains responsibility for the overall proof and integration. No earlier sealed artifact was edited.

## 1. Dependencies and what was actually checked

The accepted R009/R010 structural certificates supply46rays,135ordered cells and212top simplices covering the ordered frequency cone, followed by the full signed permutation group B7. The R011 Fourier/DD identity supplies the positive-cutoff expression. These structural proofs are prior dependencies, not recertified by a new enumeration here. Current sources are hash locked in sources.json. In particular, the numerical manifest and all its members, the frozen R012 manifest and all members, the R015 theory manifest and members, the numerical pre-run hash list, the original R011 extensions, R001 proof and R003 exact point source were checked by the independent checker. Source reading is not claimed as independent phase execution.

The independent clone HEAD is14f3507f406185667f84146328f1db9da3f9949a. Only audit/round015/green_point was written. Registration preceded the local arithmetic run. The sole checker run completed exit0 on its first attempt, with empty stderr; no failure or repair occurred in this audit. The numerical producer's saved receipt says completed/error null,5040new Gchunks,256bits. It is separate producer evidence.

## 2. Why epsilon zero computes the original G

Write u_r=xi_r/q_r and b_(Q,r)(x)=x dot Q u_r. For epsilon>0 the accepted complete formula is

    H_epsilon(x)=(2pi)^(-7) sum_(Q,S) alpha_S
         integral_epsilon^infinity product_(r in S)(t-i b_(Q,r)(x))^(-1) dt,
    alpha_S=|det W_S|/product_(r in S)q_r.

It equals the integral of the actual nonnegative Cauchy density p_t. At the fixed point, every group/ray knot is nonzero. Thus each cone integrand is bounded near t0 and at most t^-7 in absolute value for t>=1. Each improper time integral is absolutely convergent, and the sum is finite. Dominated convergence for these cone integrals and monotone convergence H_epsilon↑integral_0^infinity p_t=G identify the same epsilon0 limit. One may not justify this solely by observing a small cutoff or an apparently generic point.

Taking real parts of the complex log limit yields DD6 log|b|. At any repeated nonzero b, its kth confluent coefficient is

    (log|b|)^(k)/k! = (-1)^(k-1)/(k b^k), k>=1.

Mixed knot signs cause no issue: distinct-node recurrence denominators remain nonzero, while local derivatives at repeated nonzero knots are legitimate. No analytic branch across zero is assumed. A fully repeated integer knot1 has DD6=-1/6; individual cone contributions may be negative. The real recurrence, positive alpha weights and full signed sum must preserve them.

For integer knots a_r/B, B=dx*lcm(q_r), order6 scaling gives B^6 DD_a log|a| because the additive constant -logB has zero positive-order DD. This B^6 appears exactly once. Opposite phases have conjugate contributions, so the322560half-phases get factor2/(2pi)^7 once. At this point all7absolute coordinates are nonzero and distinct, so its signed-permutation stabilizer is trivial. The5040permutation chunks contain64half-sign choices each; no further fibre factor is present.

## 3. Exact nonzero evidence, source changes and normalization

The new kernel is compared to the frozen R012 kernel by whole-module AST and textual diff: exactly three function bodies and three diff hunks change, namely dd_real, phase and generic_chunk. Manual source review confirms these are precisely logabs/nonzero Hermite coefficients, epsilon0 plus the46exact-integer nonzero guard, and the epsilon0 generic call. The real recurrence, alpha weights, knot compiler, B^6 and phase normalization are unchanged. The actual run.py imports this module and calls generic_chunk(j,256). Legacy positive-epsilon scalar/origin/main helpers remain in the module but are not called; they are not valid general epsilon0 entry points or evidence that such tests ran in R015.

Each new phase calculates its46exact integer projections and rejects any zero before DD. The wrapper checks the complete ordered64phase IDs and aggregates nonnegative zero counts. All5040terminal chunks are completed, range0..5039, without duplicates, and have zero zero-ray count. This is producer runtime evidence plus independent accounting/source verification, not an auditor recomputation of14,837,760ray projections. Opposite phases preserve nonzero status. The corresponding R012 all-phase zero receipt is checked independently as additional consistency, but it does not substitute for the new guards.

All5040Gchunks are new, including pilot0..3. The source refuses an existing terminal ledger, computes the first4as G, then dispatches4..5039. No R012 H interval enters the G sum. Perchunk CPU/wall fields are positive; the first4G interval sum is enclosed by the saved pilot and passes the registered planning width. The final total is a direct sum of already normalized signed intervals, without a second factor. This follows from the source and ledger; hash checks alone would not establish scientific completeness.

Counts reconstruct exactly:68,382,720top requests;659,175,346unique subDD computations;14,837,760log evaluations;0zero knots, repeated top nodes or Hermite evaluations in the actual generic batch. The latter does not replace the mathematical repeated-knot formula. There are2502strictly negative chunk intervals and0intervals crossing zero; all signed endpoints are retained. Same counts as R012 reflect the unchanged node/cache structure, not interval reuse.

The exact point is

    (1073741824,963735967,1028576497,992591807,36795719,796227913,992507474)/1073741824.

Its hex strings match R003's two index2710minimum records, the R012 saved result and the frozen testcase. First coordinate1, all others strictly between0 and1, and7distinct absolute values are checked exactly. A span of k<=3body directions has at most2^(k-1)<=4distinct absolute coordinate values, so it is offS3. The stronger nonzero-ray certificate above directly suffices for the present continuity argument.

## 4. Independent exact endpoints and H comparison

The auditor adds every saved chunk lower endpoint using Fraction, and separately every upper endpoint. Only these two fixed sums are converted outward with local Arb768. They are represented exactly at this precision, so this step loses nothing. The resulting independent interval is

    [229614710670979215991566231096334957101598057539195956428066673092728624708171,
     229614710670979215991566231096334957101598057539195956428066675099327124903095]
       /60708402882054033466233184588234965832575213720379360039119137804340758912662765568.

It lies within the reported256bit outward interval, is positive and has width<10^-12. This slightly tighter endpoint is **only arithmetic sharpening of a fixed saved certificate**, not a new point or kernel evaluation. It remains approximately3.7822558290172949438809082638268589e-6. Exactresults.json records the full rational width, positive H-to-G gap and every target gap. The frozen same-point H_(1/1024) upper is strictly less than the independent G lower. This is a consistency test; the epsilon0 limit argument in section2 is what identifies trueG.

## 5. Continuity and essential-floor quantifiers

Set delta=min_(Q,r)|b_(Q,r)(x0)|>0 and R=max_r||u_r||. Finiteness and nonzero knots give delta>0 without needing a numerical lower estimate. For ||x-x0||<delta/(2R), all |b_(Q,r)(x)|>=delta/2. Every cone integrand is then uniformly bounded by(2/delta)^7 for0<t<=1 and t^-7 for t>=1. This common integrable bound and continuity of the integrands prove continuity of each cone integral and the finite full sum throughout an ambient neighborhood; the same cone/MCT identity identifies that sum with finite trueG there.

For any legal cube essential floor g, G>=g a.e. on[-1,1]^7 implies g<=G(x0). Otherwise a neighborhood of x0 has G<g by continuity, and its inward intersection with the cube has positive7Dvolume, a contradiction. Although x0_1=1, the other coordinates lie in the face interior, so one can move the first coordinate slightly inward and fit a ball. For a face essential floor the relative neighborhood has positive6Dmeasure directly. Lower semicontinuity alone would not justify this obstruction; here actual continuity has been proved. No quantitative spatial cover is asserted.

The exact reported G upper and the independent sharper upper are strictly below1/219000,1/194000,1/179000. Therefore all three corresponding cube/face essential or global floors for this originalG are impossible. This is a genuine target exclusion, not just failure to certify a lower bound.

## 6. The fixed-budget exclusion is limited

Let I0=7723627498344/10^18 be the prescribed complete-I upper budget, as a fixed numeral. Any legal positive floor g satisfies g<=G(x0)<=G_U, hence

    [7^(7/6) I0/g^(7/6)]^6 >=7^7 I0^6/G_U^7 >128^6.

The independent exact Fraction recomputation uses the producer's wider G_U, reproduces assessment.json's exact rational sixth power and verifies a strictly positive residual. Thus even the best legal floor cannot make this fixed-numerator certificate expression<=128. The displayed sixth root approximately158.38947 is optional and was not recomputed by this auditor; the exact power test is sufficient.

Since I0 is an upper budget, this is not a lower bound on actualI or M, not a lower bound on C7, not a new C7upper near158.39, and not failure of the entire original Green model. A sharper numerator may change the certificate cost. No globalg, newI, operator bound, whole-face grid or new phase is obtained here.

## 7. Reproduction and stopping boundary

minchecker.py is independent, read-only and stdout-only. It imports no author generator/checker and computes no logarithm, phase or DD. Supply --numeric for the frozen green_point directory, --r012 for frozen kernel_full, --point-theory for the frozen R015 theory directory, --repository for the main repository containing R001/R003, and --extensions for R011 extensions.md. Run using a Python with python-flint0.9.0; redirect stdout externally if saving exactresults.json. All paths are explicit CLI inputs; no working directory or author's absolute path is required by the checker.

The local run used the supplied venv Python and python-flint, exited0, and wrote exactresults.json by shell redirection. It did not run the author checker, any scientific generator, MATLAB, MC or a new point. Local rational/interval skill guidance was used only for proof direction and endpoint workflow; no unverified generic interval API claim is borrowed from that skill. Parent can rerun this portable bookkeeping interface. Source hashes, actual receipt and audit manifest accompany the review. This audit is complete and stops here.
