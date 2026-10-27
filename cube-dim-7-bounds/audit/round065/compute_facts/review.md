# R065 compute facts review

**Result: FACTS ACCEPT, with scope/provenance limits below.** This is an artifact and execution-record audit only. No producer grid, quadrature, or integral was rerun; no MATLAB was used.

## Run identity and finite coverage

The registered producer is `upper/round065/compute/derivative.py`, SHA-256 `00c182a664e6dfd379e752e34b2ff7af256891139458611ba5a7a703dc846109`. The registration (`derivative_registration.json`, SHA-256 `a2febc3807ce058b036b569a3112e5f4eada91e6a4c0ce1431cf83c2c2ae964b`) records 128 tasks × 1,800 states = 230,400 states and 15 workers. Its `code_sha256` matches the current producer source. `derivative_summary.json` (SHA-256 `84f0ca241bd5bc4fbdddcefff8ce771a7466d36bf0de2c9e7433a259ecc14048`) binds that registration hash and reports 128 chunk receipts.

I checked each of the 128 receipt files against the SHA-256 in the summary: all matched; every chunk is `completed`, has 1,800 states, and has no recorded failure. Totals are 230,400 states and 1,280 independent checks. Parsing each check record gives 1,024 direct quadrature/one-sided derivative crosschecks plus 256 baseline-mean checks. The archived `original_chunks.json.gz` contains all 128 original UTF-8 chunk payloads byte-for-byte; their decoded content and hashes match the current chunk files. The archive has a post-run manifest and a post-run clean-copy verification receipt; those are preservation checks, not registration-time or original-run receipts.

Registration specifies the tested family as `q0=2/3+(909/1000−2/3)(2i+1)/256`, `s=3/2`, `r=s*q0` fixed during differentiation, `R=s*q−h`, `a0=(1−epsilon)q`, and `g=R−a0`, with `epsilon=factor*h/(s*q0)` fixed during differentiation. This parameterization and registration differ from R057, so R065 is not the same passed grid repeated. It remains a finite tested family. The point list includes the exact onset `±a0` and nearby points, but finite one-sided checks at/on a kink do not prove a global or almost-everywhere derivative bound. The reported numerical candidate/max is not a global supremum or a `C_n` bound.

## Precision and resources

Registration distinguishes the stable formula evaluated at 160 decimal digits from direct quadrature/difference crosschecks at 256 decimal digits (more than 512 bits). This is **not** 512 decimal digits. The summary reports wall time 19.7459 s. Summing per-chunk CPU fields gives 247.9498 s over 15 distinct worker PIDs. Monitor snapshots report a peak aggregate CPU sample of 1448.8% (about 14.488 cores) and peak aggregate RSS of 1,910,224 KiB (about 1.82 GiB). The CPU percentage is a sampled peak, not sustained utilization; cumulative worker CPU exceeding wall time is expected under parallel execution. These records do not support claims of full-run all-core saturation or memory exhaustion.

## Scalar candidate correction and provenance

The preserved initial scalar source `candidate_initial_sign_error.py` has SHA-256 `2f2b576dd8756c67b6833912b74d5d604638f9e3a02b58b14d23a8b48a492616`; its paired JSON records that same source hash and `tau_identity: false`. The source contains `tau=sigma/2−chi*(1+sigma)`. The corrected `candidate.py` has SHA-256 `653dec8dfc6d4d16734cde1728a6f7264660debcfcf2eb0a55cbc2f5c6fc2cc1`; its paired certificate records that source hash and `tau_identity: true`, with `tau=sigma/2−chi*(1−sigma)`. Both exact source/JSON pairs are preserved, so the false gate is evidenced as a transcription/sign error rather than a mathematical counterexample. These JSON records are scalar-candidate receipts only; their stated scope expressly excludes geometric guards and the operator theorem.

The producer certificate uses `m=3097/10000`. The separate R064 author-side value `m_L=1239/4000` is not present in the inspected R065 producer registration/certificate, so it must remain distinct provenance and must not be presented as the same receipt or parameter. I did not locate an R064 author source in the scoped producer/audit paths searched; this review therefore makes no claim about how that separate value enters the full proof. The 8 scalar gates certify only the producer's conditional scalar candidate, not the complete `C_n` theorem.

## Source hashes

- `derivative.py`: `00c182a664e6dfd379e752e34b2ff7af256891139458611ba5a7a703dc846109`
- `derivative_registration.json`: `a2febc3807ce058b036b569a3112e5f4eada91e6a4c0ce1431cf83c2c2ae964b`
- `derivative_summary.json`: `84f0ca241bd5bc4fbdddcefff8ce771a7466d36bf0de2c9e7433a259ecc14048`
- `archive/archive_manifest.json` (post-run): `b364b485f69ef315ceeab72586e14433d76e91d0f7cfb36915cd7a9d541a44c0`
- `archive/original_chunks.json.gz`: `6068f9ca479b1c3895f64a508c0899c827210855ba6d7949fdba015cbc15e711`
- `candidate_initial_sign_error.py` / `.json`: `2f2b576dd8756c67b6833912b74d5d604638f9e3a02b58b14d23a8b48a492616` / `31e90fbd9d232825ed2346cc1594b91be0421c266268b215be3e52cd01b3af1c`
- `candidate.py` / `candidate_certificate.json`: `653dec8dfc6d4d16734cde1728a6f7264660debcfcf2eb0a55cbc2f5c6fc2cc1` / `e81bdad4218cf039b9f8861691b49233a137d7b64afdb865945724816b900d8c`
