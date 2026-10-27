# R067 compute facts review

**Result: FACTS ACCEPT, limited to execution and artifact identity.** I inspected registration, producer source, run log, summary, all 256 chunk files, monitor samples, candidate certificate, and the separate ratio probe. I did not rerun producer mathematics or use MATLAB. This review does not certify the analytic inequalities or a `C_n` theorem.

## Registered production batch

`registration.json` SHA-256 is `9cdcd902b9aa10eac076289fde13a7f137ac357e91e146f0f91baadf39524e28`; it registers `phi_diagnostic.py` SHA-256 `b221843d9e8c0cf0a55d6f892684092eda640f4a2b01cf2eeb9e7b4b4f5f198a`. The actual source hash matches. The run log start record carries this registration hash, and the summary binds the same registration. The producer checks that the registered source hash matches its own file before starting; each chunk also records that source hash.

The registered finite grid has 256 q midpoints, five `h/r` values, six `rho/r` values, and three epsilon factors. The code evaluates 25 points per combination (including kink and nearby points), giving 5×6×3×25 = 2,250 states per chunk and 576,000 total. All 256 chunk files match their summary SHA-256; every chunk reports `completed`, 2,250 states, and no failure. All per-chunk source hashes match registration. The run log records 256 completions and zero failures.

The 2,560 independent-check records comprise 2,048 direct quadrature plus finite-difference derivative checks (four points at each of two h values in every chunk) and 512 baseline records (two h values per chunk). Each baseline record carries the mean comparison, second moment, and variance. The direct checks are registered and executed at 256 decimal digits (>512 bits); the closed-form diagnostic uses 160 decimal digits. At exact kink points the code computes both left and right derivative formulas and compares both to one-sided finite differences; the receipt stores the maximum derivative error for the pair, not separate left/right error fields. Baseline checks use the registered uniform target and source-mean formula and also record second moment/variance. These checks establish what was evaluated for this finite run, not a global pointwise/a.e. theorem.

## Timing and sampled resources

The summary and run log give wall time 70.4013 s. Summing chunk CPU fields gives 953.8347 s across 15 distinct worker PIDs, an average of about 13.55 worker-CPU seconds per wall second. Fourteen monitor snapshots show a sampled peak aggregate CPU of 1462.2% (about 14.622 cores) and peak RSS 3,765,840 KiB (about 3.59 GiB). These are sampled peaks; they do not imply continuous full-core use or memory exhaustion.

## Separate post-run ratio probe

`ratio_probe.json` SHA-256 is `a166556fecf122824ccf73d3aea205a43d0ca4413b8f13c4c7e3adb608198315`. Its embedded scope explicitly calls it a **post-run** finite H/J diagnostic, separate from the production batch. It records 1,536 rows across 256 q midpoints, `h/r=2^-48`, epsilon factor 1/2, `rho/r` in `{4,2^256}`, and `w` in `{0,a,r}`. The recorded maximum `rho_scaled_Jphi_over_J` is approximately `0.0456988229162`, at q index 255, `rho/r=2^256`, `w=0`. The probe source hash in the JSON matches `ratio_probe.py` (`4a97a66403f71a301986f1c6614465c50b218601d807e1062c38aeb4911bff3b`), and its `phi_source_sha256` matches the producer source. This statistic is a finite pointwise H/J ratio; it is not a typical-target sum mean and does not establish an all-target/source-face statistic.

R067's registered q grid begins from `917/1000`, and its `phi` observable/source and checks differ from R065's `909/1000` derivative batch. Thus it is a distinct finite diagnostic, not a rerun of the R065 passed grid. Neither the grid nor the ratio probe proves a global error/rate, sharpness, or `C_n` optimality.

## Conditional scalar candidate

`candidate_certificate.json` SHA-256 is `95d153705058bd41da86c697e2f638a9691e179278138ecea42ddf2933278464`; it binds `candidate.py` SHA-256 `680b71774a86e04016609e53514bba5d99023f908584428d49c9e063aa3f84b9`. Its nine listed gates are true for the stated scalar parameters, including `m=637/2000`. The certificate labels its scope “conditional scalar parameters; no target profile or operator theorem.” It cannot by itself activate a `C_n` conclusion. The separately accepted R066 full theory is a distinct proof artifact and is outside this compute-facts review; no theorem status is inferred from these numerical receipts.

## Hashes inspected

- `registration.json`: `9cdcd902b9aa10eac076289fde13a7f137ac357e91e146f0f91baadf39524e28`
- `phi_diagnostic.py`: `b221843d9e8c0cf0a55d6f892684092eda640f4a2b01cf2eeb9e7b4b4f5f198a`
- `summary.json`: `f63cc094eeb3df59845d980a78a4b6c4d955efe9882a5daa1b6a260c2d406517`
- `ratio_probe.py`: `4a97a66403f71a301986f1c6614465c50b218601d807e1062c38aeb4911bff3b`
- `ratio_probe.json`: `a166556fecf122824ccf73d3aea205a43d0ca4413b8f13c4c7e3adb608198315`
- `candidate.py`: `680b71774a86e04016609e53514bba5d99023f908584428d49c9e063aa3f84b9`
- `candidate_certificate.json`: `95d153705058bd41da86c697e2f638a9691e179278138ecea42ddf2933278464`
