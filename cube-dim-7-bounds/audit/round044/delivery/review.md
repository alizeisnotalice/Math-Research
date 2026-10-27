# R044 delivery review

Read-only evidence review against the compute clone and the separate num-audit clone. No calculations were rerun. This review does not modify producer or main-clone files.

## Evidence and provenance

- Producer report: `/Users/zhengzhihao/Desktop/T/cube-r030-compute/cube-dim-7-bounds/upper/round044/compute/R044_REPORT.md` (SHA-256 `4fd3b655d12dae680fd9c528624922c0ed2af39cb8122d2497efe1d14dee0c83`).
- Main run registration: `registration.json` SHA-256 `6e610b50931701de3b56ba799e99e0d1dbe703de5d1d3836030f376cf2a35d59`; its summary records 76 cases and 2.908 s wall. It pins the input Phase A bundle/manifest and producer code hashes.
- Separate triple-retry, B, and C stages have their own registrations and summaries. The retry registration explicitly says no full-LP rerun; B fixes 15 cases, at most 512 new triples each; C fixes 15 workers and 50,000 branch nodes/case. Their summary registration hashes match their respective registration files.
- Independent audit: `/Users/zhengzhihao/Desktop/T/cube-r030-num-audit/cube-dim-7-bounds/audit/round044/README.md`, `lp_batch.json`, `audit_c.json`, and `candidates.json`. The stdlib code is a separate num-audit implementation; this is the independent mathematical check. Producer `verify_exact.py` is producer-written checking code and should not itself be described as an independent author.

Checked evidence hashes: `summary.json` `8e65894a82bc8ca15d82bc8ca9538895a6dbe87da47820f64a292cda27a34d49`; `triple_retry_summary.json` `a41d67b1dbb4e2e3c16ff823077da1a16166b94cf5e3065f6d8e8edacaef0b68`; `B_summary.json` `20aa3e0609811b6a73b934bf64b9e85240e62add7ba1ea34fbdab8091b466c78`; `C_summary.json` `c9e92983e202c33e4931fd1472b0e1034240b3e17ae708f85bb80ee9accde736`. Monitor hashes: initial `278239863f9f947740afcf53fc3a6621741b394db600cd67e102f822c57d00ca`, retry `4f3b0e2b13f824d8b38e61e68af197bb0d699158cedee6bbaac6182a2a780413`, B `9db01dfd5747f60fb82b1f2c6e3d0ee33f7d200068503c489abd8d303380085b`, C empty-file SHA-256 `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`. Independent audit summary hashes: `README.md` `a8989f66a3602a12b49abfbe4ae08a44e3c1a1814406cabb4a2ea81c3f7f0562`, `lp_batch.json` `939eb2e0f752d8183b92fac11b9ef945ff10e38d7377f793ef8002a419e42653`, `audit_c.json` `6735caab4dedf89683ec3cb23f1e15fbef74ff46a6d6a64cc833127eeaaa2d8d`, `candidates.json` `20873a60b9911fc8c3f5211cd51c52af8fd581672d90d892681095801cc8d400`.

## Exact certificate and candidate findings

- Independent exact checks verify all 76 recorded primal/dual vectors are nonnegative, satisfy every primal vertex and dual clique inequality, and reproduce their objectives. Hence each provides a valid bracket `M_lower ≤ M ≤ M_upper` and the corresponding payment bounds. **Every one of the 76 exact primal/dual objective gaps is strictly positive.** Floating-point display rounding to zero is not equality; these are feasible brackets, not exact LP optima.
- For the full global proposals, exact payment and full-clique checks show 67 of 76 strictly reduce the initial maximum load. The other nine do not give a strict improvement in these saved candidates.
- The independent audit confirms two saved R044-B candidates (`C/1/13`, `C/4/07`) are feasible strict improvements, each about 1/6 of initial load. These start at R038 coordinate stops; they are not certified all-pair/triple optima. B's 13 no-result states reached the finite 512-new-triple cap; this is not proof those triples cannot improve.
- The independent C audit verifies all 38 saved rational vectors, their payment, nonnegativity, full-clique loads, and changed-coordinate supports. They are feasible candidates, not minimum-support certificates: the underlying MILP status and reported numerical gap do not prove integer optimality. Of the 76 C records, 67 have numerical status 0; 38 contain exact verified vectors and 29 status-0 results have no exact vector. The other nine had no certified initial strict gap.
- The triple retry logged 2,418 rationalization-feasible candidates and 14 rejected rationalizations. The 14 failed vectors were not retained. The audit therefore cannot reconstruct their coordinates; these are failures to produce an exact candidate, not mathematical infeasibility or excluded solutions.

## Resource and runtime wording

- R044-B's sole monitor observation is at 0.547 s: 15 workers, 1096.6% aggregate CPU (about 10.97 cores), and 1,156,096 KiB aggregate RSS. This is one transient sample, not sustained saturation or memory exhaustion.
- The initial LP and triple-retry monitor each sampled before worker launch (zero workers); C finished before its monitor interval and has an empty monitor file. Their summed worker CPU and wall times are task receipts, not evidence of continuously full CPU use.
- Registrations bind code/input hashes and declared budgets. They do not establish a frozen Python/HiGHS environment or an independently verified runtime build. The report explicitly says the first run did not establish HiGHS's actual internal thread count; the retry used fresh processes with `threads=1` and did not repeat the full LP.

## Scope limits

R044 works on 76 cached R038 payment vectors, excludes the three cases with prior R038 mass-LP certificates, and does not rerun the original geometry selection, atoms, or baseline mass LP. It yields finite feasible payment redistributions only: **no new `C_n` bound, no exact LP-optimality assertion, no exhausted triple-search result, and no seven-dimensional lower-bound work.** The producer report's exact claims and the independent audit agree on these limits; no report correction is indicated from this review.
