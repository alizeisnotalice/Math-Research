# C147: strict dimensional reduction is available; one real prototype completed

For every current C130 pair the joint law of its two nonnegative linear forms can be represented as a positive two-dimensional pushforward. A constructive convolution proof produces an exact exponential-polynomial formula on a finite rational angular fan; it does not assume independent outputs or generic distinct columns. The original six-dimensional positive expectation becomes a genuine two-dimensional integral against this density. The full general proof, singular rank1/rank0 alternatives and normalization arguments are in `PROOF.md`.

A useful refinement keeps the geometry small: comparison lines for two ray-integration endpoints all annihilate the newly added column. Thus they give its own direction, not new arbitrary break directions. The final fan needs only the original column directions. With m active columns there are at most m-1 nonempty sectors (at most5 here), although the number of signed exponential-polynomial terms still needs measurement. Duplicate columns, equal directions, zero columns and exact resonance are handled directly; no problematic distinct-root partial fraction assumption is used.

Positivity belongs to the complete density because it is a convolution of probability measures. Signed polynomial/exponential terms must be recombined before interval bounding; they are not independent positive contributions. Total mass1 follows by Tonelli. A separate exact bivariate rational Laplace identity is designed to catch symbolic generation errors at each stage.

## One actual source-bound prototype, no k integral

The registered SymPy1.14.0 execution used only cone0/body pair(16,17), zero-based. Its four zero columns disappear with exact mass1, and

    A=E5/60+E6/30, B=E5/60,
    f(A,B)=1800 exp(-30A-30B), 0<B<A.

The symbolic run verifies determinant-1/1800, inverse(E5,E6)=(60B,30(A-B)), mass1 and exact joint Laplace

    1800/[(s+30)(s+t+60)].

It also verifies E A=1/20, E B=1/60 and covariance1/3600>0, so output independence is explicitly false. The source-to-density identity is exact and the original cone weight is retained.

By symmetry of the complete integrand, the pair expectation additionally equals

    900 [integral_0^infinity exp(-30u)k_h(u)du]^2.

This special simplification does not claim A and B are independent. The proof derives the closed expression for this one-dimensional transform and its lambda0/2 removable cases; no erf/erfc/exp value or k integral has been evaluated. The original C130 h and k_h definition remain unchanged. This prototype alone supplies no numerical improvement of mu2, M or C6.

`finite_registration.json`, `prototype_result.json` and `exit_receipt.json` preserve the exact source pins, one real case, executed symbolic jobs and success. Actual symbolic wall0.764189875s, CPU0.661367s;0 new k/special-function values,0 old6D grid,0 native or optimizer calls. Generic integration-by-parts and resonant monomial identities were checked through degree4.

## Next two-case generator: prepared, not launched

`generate_next.py` is a separate dual-gated symbolic packet. Fixed source cases are:

- cone0 pair(9,16):3 active columns,3 directions, including axis boundaries;
- cone0 pair(0,1):6 active columns,2 directions, including repeated columns/rates and parallel columns of different scales.

The compiler constructs the exact fan and density using the proved recurrence. It treats eta0 separately and checks the whole two-dimensional Laplace rational expression against the source product after every stage:7 checks across both cases. It proves mass1 only after whole-expression rational cancellation. Fixed caps are10000 nonzero terms per stage and100000 raw endpoint contributions total, with maximum polynomial degree4 and maximum5 sectors. It stops on structural/algebraic failure, preserves a failure receipt, and does not restart. It uses one process at nice>=10, one-thread environment, no time cap and no k/kernel/integration/native evaluation.

Exact case/proof/runner hashes and finite limits are in `next_registration.json` and `next_release_requirements.json`. Current runner SHA256 is `696e6bdaf629fdfc6404e76c525149fc101b798a83878b9dae43fcc6ba4c0daf`. Only its saved-only source/case preflight has run. There is no `next_parent_release.json`, `next_started.json` or `next_result.json`. Main and an independent reviewer must approve this new generator before execution; this author has not independently audited its own implementation.

Terra's source-pattern enumeration is pinned by hash and used only as external preparation evidence, not as a density proof or an excuse to stop. This packet does not generate densities for all10416 pairs or evaluate any of the9961 canonical keys. Whole-set expression size, cancellation stability, rigorous quadrature/tails and cost versus geometric width remain open empirical/engineering questions. No new M upper, globalg or C6 claim is made. Nothing was staged or pushed, and C146 was neither modified nor audited here.
