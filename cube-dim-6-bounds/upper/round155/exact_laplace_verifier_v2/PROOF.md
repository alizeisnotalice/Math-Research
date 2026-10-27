# C155 v2: same saved inputs, explicit API accounting and guarded dimensions

This is a separate unrun source package. The v1 source and Terra's first v1 preparation evidence remain unchanged. V2 must reuse the exact first case_packet bytes, SHA256 `6587f978228bc7478f25919425d0ed915ff24270786f4a9961027da17b3d495b`, and source_manifest bytes, SHA256 `3012a7f9a22c838334158228ae2e44c8e2e4f84e10acdaf2fa6039b626e34d9b`. The original v1 registration is SHA256 `c9324fcfac4a4ca1162446c7c1ea2faceff127a2d7961758ace6d90ff6a072df`. The selected source cases remain case_024, case_015, case_008, case_003, case_000, one each for m2..6, with exactly15 complete saved stages. Prepare.py never rescans149 or reselects cases based on subsequent progress.

## Mathematical certificate, unchanged

For a saved density term P(z)exp(-lambda.z) on the cone z=a u+b v, u,v>0, substitute the exact ray coordinates and expand P by the binomial theorem. With positive rates its transform is

    det(a,b) sum_ij p_ij i!j! / [L_a(S,T)^(i+1) L_b(S,T)^(j+1)],
    L_a=lambda.a+a_1 S+a_2 T.

Subtract the product of original active input factors (1+b_j S+c_j T)^-1. Normalize each rational affine factor to a rational scalar times a signed-normalized primitive integer coefficient triple. Equal factor powers add. Absorb each scalar power into the rational coefficient. If E_F is the largest denominator power of each primitive factor and B clears all rational coefficient denominators, the complete rational difference multiplied by B*product_F F^E_F is an integer polynomial N(S,T).

Only exact fmpz_mpoly setup/construction/addition/multiplication and exact zero testing are used. N identically zero is a strict rational-function identity, not a modular/random sample, and needs no CRT height assumption or multivariate gcd. Common-denominator growth remains a practical risk, controlled by explicit finite caps and preserved failure rather than dropping terms or expanding the batch.

All support rays are nonnegative/nonzero. An exact integer tilt q makes every affine rate positive at S=T=q, hence on S,T>=q. The tilted signed density is a finite signed measure, termwise absolutely integrable. The rational identity gives the same transform as the tilted original positive exponential pushforward law. Finite-measure Laplace uniqueness identifies these tilted measures. Removing the strictly positive tilt identifies the complete density with the original probability law, yielding nonnegativity almost everywhere and mass1. Individual signed terms are never presumed positive, and removable origin poles are not evaluated separately.

Raw six-column vectors, saved active vectors, support rays and covectors are now explicitly exactly two-dimensional. Every raw-column and lambda rational numerator/denominator receives the8192-bit input guard before dot products or target factors. Rays/active vectors also receive it. The dot/determinant helpers themselves reject dimensions other than2, so zip cannot truncate an unchecked covector. Derived affine values and primitive integer coefficients retain their bit guards. These strengthen schema/resource validation without changing the law or the selected data.

## Complete explicit API accounting and budget meaning

The old packet's `native_arithmetic_per_stage=200000` and `native_arithmetic_total=3000000` numbers remain unchanged. In v2 their binding semantics are ALL explicit verifier setup, construction and arithmetic invocation attempts, including fmpz_mpoly_ctx.get, constant, from_dict and polynomial add/multiply. They are NOT a cap on all FLINT API calls. All explicit queries are separately counted in the same ledger: len, is_zero (including the final zero test), total_degree, coeffs, and every fmpz-to-int conversion. Context setup is no longer omitted. Library-import/bootstrap internals and FLINT's private implementation instructions are outside explicit-call counts; no claim is made that these are instrumented.

Every explicit invocation passes through api_ledger.Ledger.invoke, classified setup/construction/arithmetic/query. Stage and global ledgers contain attempted, completed, raised_exception and an INFLIGHT_UNKNOWN descriptor. The uncertainty marker is set before budget/counter bookkeeping or the actual call. Successful return increments completed; an invocation that raises increments raised_exception and retains its exception descriptor. No call that raises is counted as completed. A budget rejection occurs before invocation and is recorded separately. Quiescent states assert attempted=completed+raised_exception; the driver additionally checks per-category totals.

An ordinary caught stage validation/API failure preserves the observed ledger. If bookkeeping was interrupted or an invocation remains unresolved, API_totals_exact is false and the exact total is None. A process crash may leave no failure receipt at all: the durable stage intent marks that stage/API proof status UNKNOWN. Prior completed stage snapshots remain observed lower bounds, never invented final zero counts. Terra's terminal receipt must report final API counts None and list the in-flight/unknown stage if no valid exact terminal receipt exists. This accounting is deliberately observable-call accounting; it does not assert mathematical completion for an invocation whose process died.

Queries have no new independent resource ceiling; they are counted completely. Their number is finite under the existing construction/arithmetic and visible polynomial term caps, because each polynomial inspection performs finitely many queries and at most100000 coefficient conversions. No hidden query budget is folded into the3M claim. No time/CPU ceiling is introduced. These visible-boundary guards do not instrument private FLINT temporary allocations.

## Runtime scope and evidence

All original case-data caps, five cases and15 stage checks remain fixed. V2 is not permission to retry failed v1 science; no v1 verifier science has been claimed here, and any unexpected existing v1 run evidence must be returned to main before proceeding. V2 prepare-only freezes provenance and byte-identical inputs; science still needs a concrete main/independent exact-hash release and one Terra process, single-thread environment, nice15, no restart/retry. A runtime/API/cap failure stops and preserves evidence.

Only local API stubs/source and these new source files have been read/written/static-parsed by Sol. No preparation, FLINT call, numerical/symbolic checker, density regeneration or old baseline has run. V2 changes no live149/151 source and produces no k integral, beta2 or Cn bound. Saved producer whole-case CPU is contextual only; no same-workload speedup or live replacement claim follows from this prepared verifier.
