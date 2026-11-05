# R004: a compact nonproportional kernel improvement

The certified kernel is

\[
K_4(u,v)=K_1(u,v)-\frac1{200}B(32u)B(32v),
\qquad C_2<3.590865312933.
\]

Here K1 is the project's R1 scaled kernel and B is the centered cubic B-spline. Only its `a["0,0"]` coefficient changes. The positive singular coefficient, alpha=47/40, R=33/20, h=1/32, and the other 1180 classes are unchanged. The half integral decreases by exactly 1/409600. This changes the shape, rather than proportionally scaling the kernel.

The result uses the Lin–Tomar comparison principle (Theorem 2.2 and Lemma 3.1 of the supplied paper): an even integrable kernel majorizing the unit diamond, with nonnegative coordinate-jump generator as a distribution outside the origin, bounds C2 by half its integral. It applies to general L1 inputs and all positive scales for centered axis-parallel squares and Lebesgue measure. The singular exponent is below two; the same finite-spline regularity and the complete C4 check give the required distributional condition away from zero. The full spline support beyond r=R is included; outside the entire kernel support, C1 gives the nonnegative jump integral. No integrability of AK at the origin or finite-signed-measure representation is assumed.

The mathematical proof and independent mathematical review are staged separately by the project without rewriting their bytes. Their original paths and hashes, together with the source-to-package mappings, are in `source-bindings.json`. This directory preserves the numerical inputs and full new running record needed by that proof. It claims neither a sharp constant nor an optimal spline direction, and does not perform a Lean build.

## Default verification

This directory belongs at `round004/kernel` in the published project. It reuses, without copying:

- `kernel-rescaling/vendor/hl_check.py`, `hl_core.py`, and their upstream identity/license;
- `kernel-rescaling/certificates/scaled.json`, used to verify the exact single-coefficient change.

Install the small Python dependency with `python -m pip install -r round004/kernel/requirements.txt`, or use an environment already containing python-flint 0.9.0. From the project root:

```sh
python -I -B round004/kernel/verify_result.py --out /tmp/round004-kernel-verification.json
```

The verifier locates all inputs relative to its own file. It rejects Python `-O`. It checks the static manifest, exact a00-only change and all other mathematical fields, the complete new C4 report and stdout, source hashes and historical worker provenance, then actually recomputes C1/C2 and the full-orbit direct mass at 192 bits. It also integrates the spline polynomial pieces exactly and rechecks the large-step negative control. Exact rational Arb endpoints, rather than printed decimal strings or floats, prove

\[
3.590865312932196537<\tfrac12\int K_4
<3.590865312932196538<3.590865312933.
\]

The default verifier **does not recompute C4 leaves**. The unchanged full report records a new complete 128-bit Arb run for K4: 1368 knot cells, 443757 accepted leaves, no failures, with 10 pool workers. It is a per-knot-cell running summary, not an independently verifiable per-leaf polynomial certificate. The verifier checks the expected cell set and its size, excluding missing cells or duplicate substitutions, every worker's provenance against its own historical main process, and each stdout cell's leaf/Taylor/depth fields against the report. Historical absolute paths and binary hashes remain metadata; they are not opened or required to match the current host. Current environment provenance is recorded separately.

`portable-receipt.json` and `portable-execution.json` record an actual default verification from a relocated project subset with working directory `/tmp`. They do not assert a second full C4 run. Static manifest inputs exclude these generated receipts, avoiding a circular hash dependency.

## Optional complete numerical replay

To replace the historical complete C4 run with a fresh one, from the project root run:

```sh
python -I -B kernel-rescaling/vendor/hl_check.py \
  round004/kernel/candidate.json --claim 3590865312933/1000000000000 \
  --only c1,c2,c3,c4 --procs 10 --maxdepth 8 --kdepth 10 \
  --out /tmp/round004-kernel-full-replay.json
```

Reduce `--procs` to suit the machine; this changes the scheduling, not the cell set or inequalities. The original checker uses 128-bit Arb. The separate default verifier's fresh cheap calculations use 192 bits. No equality of historical and current FLINT binaries or absolute paths is a mathematical premise of these fresh checks.

For the symbolic direction integral, `symbolic/direction_mass.wl` and its unchanged actual Mathematica 15.0 log are included. The default Python verifier does not require Mathematica. The negative control proves that epsilon>=3/200 is infeasible on the ray K1-epsilon B(32u)B(32v); convexity gives feasibility for 0<=epsilon<=1/200. The intermediate interval and any optimum are unresolved.

Only a subset of the research freeze is copied here. The 1294-point screening, process samples, monitor wrappers, intermediate checks and old C4 are omitted; their absence does not weaken the retained mathematical checks. The normalized execution summary and its source clarification are historical metadata, not independent proof of CPU usage. `source-bindings.json` distinguishes unchanged source bytes, the adapted verifier, reused project inputs, and separately staged proof/review; it does not claim that omitted source files were reverified by this package.
