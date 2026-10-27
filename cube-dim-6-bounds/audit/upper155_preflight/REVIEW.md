# C155 v2 final independent static gate — PASS

The v2 source and Terra's concrete prepared packet resolve the initial review's accounting and geometry-schema issues. This is a source/mathematical gate, not a runtime PASS. Sol executed no prepare, case selection, polynomial verifier, FLINT/SymPy API, numerical preflight or fixture. Only read-only source text and SHA256/JSON metadata handling were used.

## Exact released candidate

Author source: `/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-pushforward147-audit/cube-dim-6-bounds/upper/round155/exact_laplace_verifier_v2`.

Terra target: `/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-brascamp108/cube-dim-6-bounds/upper/round155/exact_laplace_verifier_v2`.

- Author manifest: `58dbfb4dfe35d96ad540870f76c7d99d46e7fdb68ee80ffbeec3bfc5c017a675`.
- Execution whitelist: `f950fa5acb7cbc2328dd77725b184072d2c01d9cc6d7d21079270f2f7ac8f933`.
- Terra registration: `97c032d7c55612477407351ad5141be09e78326ec1d5cb854ad8afd7227eee9d`.
- Terra provenance: `5583403fa013a0d0c8e7b71e1dbad3e77bc756b584593b7b10cb64e99fb784f2`.
- Terra prepare receipt: `b0256c55c0c95e39e9819648a4a23e92bc3104950a335b6125146008f0bdb469`.
- Unchanged first packet: `6587f978228bc7478f25919425d0ed915ff24270786f4a9961027da17b3d495b`.
- Unchanged source manifest: `3012a7f9a22c838334158228ae2e44c8e2e4f84e10acdaf2fa6039b626e34d9b`.

All six source files match the whitelist in both clones and the new registration. The first packet and source manifest are byte-identical to the original Terra v1 preparation. Their source leaves match their original hashes and lengths. Cases remain `case_024, case_015, case_008, case_003, case_000`, one each for active-column counts 2 through 6, with exactly 15 saved stages. V2 preparation does not rescan149 or choose easier cases. Main reports its one actual prepare invocation exited zero; this prepares inputs without qualifying the FLINT API or the law identities at runtime. No v1 or v2 science/start/result evidence was present at this static snapshot.

## Mathematical certificate inherited and retained

The initial mathematical review remains applicable. Exact ray substitution and the binomial theorem produce the correct Jacobian/Gamma factors. Primitive signed-normalized affine triples retain their rational scalar; repeated equal factors add powers. Taking maximum factor powers gives a valid common denominator, and a positive scalar LCM clears all remaining rational coefficients. The constructed integer polynomial is the complete transformed-density-minus-source-law numerator. Its exact zero predicate is mandatory; no random/modular substitution, multivariate gcd, omitted term or unreduced-origin evaluation is used.

The independent rank-two source pivot gives an absolutely continuous source law; support fan boundaries carry no original mass. Nonnegative rays allow a strict positive tilt at which every signed polynomial-exponential term is absolutely integrable. Equality of finite signed-measure Laplace transforms identifies the tilted measures. Removing the positive weight identifies the complete density with the original probability law, so nonnegativity almost everywhere and mass one follow without assigning positivity to individual signed terms.

The finite case/stage, degree, factor, scalar/input bit and visible polynomial size caps remain unchanged. They can cause an honest finite failure or scaling obstruction. Success would verify saved laws only, with no new density, k integral, beta2/Cn bound, live149/151 change, baseline replay or same-workload speedup claim.

## Closed geometry-schema gap

The core now requires an integer active count, exactly six raw two-dimensional source columns, exactly k two-dimensional saved active vectors, two-dimensional support rays and exactly two lambda coordinates. Dot and determinant helpers independently reject any dimension other than two; `zip` cannot silently drop a covector coordinate. The 8192-bit guard applies directly to raw columns, saved active vectors/rays and lambda coordinates before dot products/target factors, in addition to the existing derived coefficient/factor guards. These checks preserve the law while rejecting malformed or oversized representations.

## Closed API-accounting gap and precise budget

Every explicit context setup, construction, polynomial addition/multiplication, length/zero/degree/coefficient query and fmpz-to-int conversion passes through `Ledger.invoke`. Context acquisition and the final zero test are now included. The counter categories are `setup`, `construction`, `arithmetic`, and `query`. Library import/bootstrap internals and private FLINT implementation instructions are explicitly outside observable-call accounting.

The exact registered budget is 200,000 **setup + construction + arithmetic invocation attempts** per stage and 3,000,000 total. It is not a cap on all API calls. Queries are completely counted and reported separately. Their number is finite under the existing caps: each polynomial inspection is charged to an already budgeted construction/arithmetic result, performs at most four fixed queries plus the visible-term number of coefficient conversions, and each stage adds one final zero query. Thus a conservative bound is `queries <= (native_polynomial_terms+4)*budgeted_attempts + attempted_stages`. No additional query budget or hidden time limit is asserted.

The ledger distinguishes attempted, completed, raised-exception and inflight-unknown state, by stage/category and globally. It records uncertainty before invocation bookkeeping. A returned call increments completed; a raised call increments raised-exception and retains its descriptor. Rejection at a spent budget is separately recorded and does not invoke the API. Quiescent exact totals require no inflight call, attempted=completed+raised-exception, category identities and category totals. An interrupted bookkeeping update therefore fails exactness rather than silently asserting zero.

The driver preserves an exclusive, fsynced per-stage intent with the previous exact snapshot before entering a new stage, and records completed stage snapshots independently. A caught failure preserves the observed stage/global ledger and identifies the failed/inflight stage. An unresolved call or inconsistent bookkeeping makes exact totals `None` and the status UNKNOWN. A process death may prevent a failure receipt entirely; its intent is not proof, and the operator must report uncheckpointed totals as unknown, retaining completed snapshots as lower bounds. No restart/reissue is authorized by a missing or failed receipt.

## Remaining execution gate

The prepared registration binds all six source hashes, byte-identical packet/source-manifest, provenance, caps and exact counter/state semantics. `run_prepared` requires a new exact-hash `155_V2` release with main and independent approvals, one process, five one-thread environment variables, nice at least 15, no time/CPU cap and no automatic restart. Original source/result pins are checked before the verifier and again at successful end. Original v1 preparation is preserved, and unexpected v1 science evidence stops v2 preparation.

There is no remaining static blocker for these exact bytes and this expressly limited budget meaning. Runtime/API compatibility and useful performance are still hypotheses: the first actual Terra stage qualifies the stub-described interface. Any API/size-cap failure stops without patching, replaying a primary/SymPy baseline, substituting cases or expanding scope. Main retains the final concrete release and proof responsibility; only Terra may execute the unique finite batch. This review changed no author/live source, index or remote.
