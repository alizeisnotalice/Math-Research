# C061 cutoff-gradient pilot: frozen preparation, zero production calls

This handoff prepares the finite cutoff pilot only. The original unexecuted
`registration.json`, `input_manifest.json`, `jet_adapter.py` and
`check_preparation.py` were found in the old audit clone and retained unchanged.
There was no preparation result or production receipt and no observed process
running from this C061 path when continuation began. C055 was not run or changed.

The existing C016 Arb512 all-real-node band/physical-gradient binary is reused
without compilation or source changes. Its exact fan, source, binary and runtime
hashes match the old records. The six physical derivatives are computed by the
same call as the band value; the adapter selects physical coordinates 2 through
6 as the five face derivatives. Its metadata binds model, normalization, exact
center and exact time window. No full-G gradient is substituted.

The sole preparation invocation completed all 12 already registered finite
configurations. Configurations 1–6 read the saved C016 scalar/Hermite evidence;
they did not re-run scalar kernels. Configurations 7–10 checked the exact S3
witness, prescribed low-point source/integer bounds, symbolic normalized jets
and saved old-control prerequisite. Configurations 11–12 exercised synthetic
metadata adaptation and wrong-window rejection; synthetic jets are not kernel
evidence. Exact arithmetic used local SymPy 1.14.0 and python-flint 0.9.0 /
FLINT 3.6.0. The preparation body took 0.1354282090032939 seconds; its recorded
wrapper wall time was 0.3280997920082882 seconds, exit 0, empty stderr.
The receipt and all checkpoints are retained. There were 0 native invocations,
0 production kernels, 0 compiles and 0 old-seed recalculations. A later source
syntax-only inspection parsed the four Python files; it added no mathematical
configuration, native invocation or compilation.

`policy_amendment_180s.json` preserves the old registration while superseding its
inactive 600/1800-second fields with the user's 180-second per-program limit.
`prepare_once.py` cannot launch kernels and refuses a prior preparation receipt.
`run_case_once.py` is prepared but has never been invoked. It requires explicit
main release and the final F057 PASS evidence, rechecks all pins, records one
case attempt before launch, refuses existing case directories, allows only the
registered branch, and places its whole execution including the owned native
child under a 180-second alarm. On timeout/failure it retains output and receipt,
kills/reaps only its own child process group, and has no retry or compile path.

The exact low point is read from A038 `results.json`, `points[2]`, with denominator
100000000000000 and numerators
`[100000000000000,99999998999999,99999997999988,99999996999867,57255226550170,43659674865]`.
The second prescribed center is `(1,1,1,1,3/5,0)`, witnessed by three independent
body columns with nonzero coefficients. The branch is fixed:

1. Run low-primary at `a=1/250,T=1/2` once.
2. If its lower endpoint is strictly above `1/65000`, use S3-primary as the
   second and final native case. If its upper endpoint is strictly below the
   target, use low-fallback at `a=1/1000,T=1/2` as the second case, and allow
   S3-fallback as the third only when that fallback lower endpoint is strictly
   above the target.
3. A target-overlapping interval, invalid native result, timeout or other
   failure stops the batch. No precision increase, old failed-window rerun,
   extra point, extra cutoff or kernel replay is permitted.

For a center strictly above target, at most eight predeclared halfwidths
`h=a/2^k`, `k=1,...,8`, are tried in that order. Each box clips physical face
coordinates to `[0,1]` and uses the fixed A057 `Q=I5/32` theorem, all 32 corners,
and the same stored center jet. No box calculation calls the kernel. The first
certified halfwidth is retained. If all eight widths remain unresolved, the
local-radius outcome is reported separately. Main has explicitly clarified
that this does not prevent the already-budgeted same-window S3 case after a
strict center pass; its release must set
`allow_S3_after_center_pass_with_unresolved_local_radius=true`.

Main has accepted the A057 seal
`235b2d544a25d63bb88ab0dd48a15830799e808ba86a526e06dbd7b28837b664`
and F057 final audit seal
`73f8be649904c41215c1594699046e27e2fd04bd080676a405d4758be62abc55`.
For the runner's PASS field, the F057 artifact is `results.json`, status
`PASS_EIGHT_INDEPENDENT_AUDIT_CONFIGURATIONS`, SHA256
`03cadd552f6d6ddc3340440c089a5658413f6443bc309051799debb634f5dfa2`.
The preparation result's historical `A057_independently_accepted=false` predates
that main acceptance and is intentionally not rewritten. Acceptance of the
theory is not a production release. The supplied release template remains
`production_authorized=false`; main will write the actual release separately.

The inherited C016 limitation remains: no independent full-group gradient
cross-certificate is claimed. Scalar DD5/6 evidence, the physical chain formula
and Euler consistency support the implementation; Euler alone does not prove
all six derivative components. Each future native result will additionally
reconstruct Euler with exact rational interval endpoints, but this does not
close that independent-gradient gap.

There is no global g-floor certificate, no new C6 bound and no closure of the
main partial-I argument. No MATLAB, recursion/delegation, stage, commit or push
was used. Main retains release, proof and integration responsibility.
