# C160 independent code-only preflight

Verdict: PASS for Terra preparation only. No mathematical/source blocker was found in the proposed normalized genuine-cap model or the two-phase implementation. This is not a science release, a validation of C156's producer BU, or a new M/global-g/C6 certificate. PhaseA and PhaseB remain separately unreleased. Full independent saved audit of the exact C156 result plus main acceptance is a hard prerequisite to either actual science phase.

The source package manifest SHA256 is `7719ad8534d9f3bb4048bb7bf2c0f4960ef76e7f6be6f8677816a39412e378ab` at `/Users/zhengzhihao/Desktop/T/Math-Research-cube6-sol-gaussian106/cube-dim-6-bounds/upper/round160/normalized_joint_with_true_beta2/PACKAGE_MANIFEST.json`.

## Exact source/cap identity

All seven packaged file hashes/byte counts and all17 external role-addressed source pins match. Syntax was parsed without importing/executing the scripts. The exact C156 result pin is f9a0689da14fb34d48c3ffb65b6abaef06193906ff71cc3cab77420074bee21f. The selected field is budgets.finite_best.valid_true_beta2_upper and is the same stored exact upper endpoint as budget_expression_interval[1], at h1/128 in both result and inputs. No display decimal, single-pair result or lower endpoint is consumed.

The code performs N=ceil(2^96 B) in integer arithmetic, checks (N-1)den(B)<2^96num(B)<=Nden(B), and uses betaD=N/2^96. This proves B<=betaD<B+2^-96, conditional only on the input B being the valid complete cap. It stores B's original exact rational and source selector/hash. The giant original denominator is not embedded into the feature cap; rounding is upward. This reviewer did not execute that quantization or generate its certificate.

Only beta_by_h[5][1] changes. The preflight verifies h_values[5]=1/128 and Phi_table[n][2]=n(n-1)/2 over the complete stored n0..32 domain; flattened index32 is correct. Thus it consumes the same true binomial second-count moment bounded by C130/C156, with no missing factor2. All remaining48 cap entries including beta0, the full Phi table, all8 h values, complete CDF intervals, cutoff/state/tail parameters and original128 threshold remain pinned/unchanged. The complete original128 dyadic49-weight candidate is the sole warm start. No C5 cap, optimistic bL or mixed old/new cap table is used. Acceptance of the complete base128 cap-table hash is explicit in each phase release.

## Exact normalized cost and conditional interface

For the exact dyadic nonnegative simplex x, x0>0 and sum x=1. Define a_i=x_i/cap_i using exact Fraction arithmetic. The original moment interfaces give E[d F_i]<=cap_i, including betaD only after full C156 saved acceptance. Hence E[dW]<=sum a_i cap_i=1, W=sum a_i F_i>0. Holder with exponents6 and6/5 gives M^6<=E[dW]*(E W^-1/5)^5<=H^5. No extra146 cost correction is needed for this newly normalized true-upper-cap model. The source records the cap-cost sum1 rather than asserting the actual unknown cost equals1.

Only x is quantized to2^-44. The resulting a remains generally rational and is not quantized. PhaseB's outward2^-192 W tables enclose that exact same a and alter neither coefficients nor moment cost. Positivity, root/cap and table agreement checks come from the pinned already-audited128 Joint interface; the123/128/146 base interfaces are reused here rather than re-proved or numerically rerun.

A strict M-slot comparison remains conditional on the unchanged determinant/whole-I interface. The complete-I expression (3/16)^6*H_U^5/384^5 has no missing cost factor in this model. No global-g/general-L1/C6 conclusion follows from an M-slot comparison alone. The source explicitly keeps global_g_proved=false/new_C6_bound=false and requires a later independent saved certificate audit.

## PhaseA selector and float-to-exact candidate

Exactly one SLSQP invocation uses the new midpoint whole-CDF model and original128 simplex warm start. The feature matrix has all735471 prefix states,16 positive tail-relaxation rows and49 normalized features. Midpoint category probabilities come from differences within one complete CDF vector. The float mass is checked, never normalized; ordinary probability log/exp calls are registered separately from forbidden kernel/erfc/native/density calls.

The objective returns its analytic gradient in the same counted call. Its guard allows at most149 optimizer calls, then one final call at the dyadic exact new candidate: total<=150. There are no retries or second starting points. Linear feasibility/gradient callbacks do not evaluate an unregistered functional. Sanitizing the retained numerical candidate and allocating the dyadic remainder to its intercept produces an exactly nonnegative sum-one integer simplex, checked before use. Cleaning/normalizing x is not probability renormalization. A local failure or budget exhaustion can produce a candidate only; it excludes no parameter family and proves no optimum.

All objective, Euler residual and old/threshold comparisons in A are expressly float diagnostics. The old128 exact bound is converted only for selection, not used to announce a strict M bound. No old H/root functional is replayed: even the warm-start evaluation uses the new cap model. PhaseA writes no PhaseB preparation/run automatically.

## PhaseB new exact certificate

B preparation requires the completed actual A result/exit hash, the exact same candidate/hash, numerical improvement over saved128 and a nonempty explicit main selection. B has a new registration, separate main/independent dual release, selected candidate/coefficient/job hashes and renewed source156 saved acceptance. A release is not inherited as permission for B.

Fifteen disjoint contiguous jobs exactly partition the735471 prefix states. Each spawned worker uses the new candidate with the pinned Joint code, one-thread environment inherited from the gated parent, nice>=10, version checks and its own65536-entry new root cache. The collector requests16 new tail roots; the total logical new requests are735487. Iterator skipping/re-enumeration for indexing/hash checks does not evaluate additional root functionals. No old functional root values or caches from another candidate are consumed.

The exact probability convention is correct: pU-weighted prefix gives the lower H enclosure, pL-weighted prefix plus all16 pL-binomial tail terms gives the upper enclosure. Coupling concerns each whole endpoint law before the tail relaxation. Tail W retains the positive intercept and final-h terms, so Wtail<=W and reciprocal fifth-root tail values bound upward. No monotonicity is assumed across the prefix/tail seam. Both endpoint prefix+tail probability sums are checked exactly equal to S^32. Weighted sums use S^33, with S=2^192, matching the inherited reciprocal-root certificate inequalities. No discarded positive mass or optimistic renormalization appears.

The table hash binds the same exact new a across all workers and the collector. Exact H_U^5 is compared with the pinned old128 M_power6_upper and exact Theta(1/65000), Theta(1/63500). Theta65000 is independently matched to the stored128 exact threshold. H_L is not treated as an M lower bound. Files contain the full new root records for a subsequent independent saved audit, and failures retain available output without automatic restart or precision/cutoff/cap expansion.

## Source156 acceptance and operational boundary

common.check_source156_acceptance requires the exact same156 result hash, independent_saved_review_passed=true, main_accepted=true, and nonempty hash-verified evidence rows. Both actual phase gates call it before science. The evidence must be the genuine full independent saved review/seal for C156's complete upper expression, not a producer result or code-only preflight renamed as an audit. Main/operator must semantically verify that fact: the code checks trusted acceptance assertions and evidence hashes, and intentionally does not infer proof validity from arbitrary prose. This division is explicit in the handoff. Until that full saved acceptance exists, producer B and derived betaD remain unaccepted candidates; code-only preparation may freeze them but cannot launch science.

The target root must be an independent existing clone with the exact project/round/root structure. P-relative output keeps each phase in its new160 root; old sources are absolute/pinned/read-only. Each phase has exclusive started/results/receipts and one allowed run. Pool exceptions do not authorize reissue, fallback or optimistic completion; all already submitted work remains within the one finite registered batch. There is no time/CPU ceiling. Main owns machine allocation and original process/terminal monitoring. This review does not create a release or run a preflight/fixture.

## Exact package pins and actions

| File | SHA256 |
|---|---|
| PROOF.md | 5b0e3a06bc86c284f0f7b78f906a46f882da1d379f705ea58a4a6ef4eba7fad7 |
| common.py | 98e46ba5675521f340625ab77de1edbbc00d60d014766886a53c7d329277ab0b |
| prepare_for_terra.py | b8ec18a3fcb662ccf4cb4ad1eeac9170a9889d2a96b203f9ba54036792010615 |
| run_phaseA.py | f8e18977cd58f25361c92ff3bb56d012d4a9fcc7a16656d795ff32c48ec5e6c4 |
| run_phaseB.py | e70f4bd5bf812e5e83d2d878363c5ce3a5bc806dfe67bd8f03e6af097c977584 |
| HANDOFF_TO_TERRA.md | 0ab43c62c2fbf9536ee2bc9fc8341c16a819bd5adf181473db3032da07a33c95 |
| STATIC_SOURCE_PINS.json | 23bd838863e24018b566e0539969098fa1d3c65e4141690b014241f3514201b0 |

Only source/proof reading, static JSON schema/field identity inspection, AST syntax parsing and hash/byte checks were performed. No rational-ceiling/Phi fixture, prepare, mathematical preflight, optimizer, Joint import, H/root evaluation, C156 replay, old functional, FLINT/native kernel or density run occurred. Only this clone's project audit/upper160_preflight was written. No recursive agent, stage, commit or push occurred.
