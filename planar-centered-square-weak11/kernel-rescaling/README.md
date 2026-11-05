# Portable kernel-rescaling verification

This small package reproduces the exact positive rescaling of Tomar's square-maximal kernel. With the comparison theorem described in `PROOF.md`, the certified bound is

\[
C_2 < 3.590867754339.
\]

It applies to general Lebesgue L1 inputs and all positive scales for the centred axis-parallel square maximal operator. It does not identify the sharp constant, an optimal kernel, or an extremizer. No Lean build is claimed here.

## Quick verification

Tested with Python 3.12.13 and python-flint 0.9.0. From this directory:

```sh
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python -I -B verify.py --out runs/current-verification.json
```

The verifier also works from a different current directory: it resolves inputs relative to its own file. Do not use Python `-O`, which disables assertions; the script explicitly rejects that mode. Runtime environments and generated files are not distributed in this package.

The default run takes approximately seven seconds on the machine used for the independent review. It:

- verifies the package hashes and the explicitly copied subset of the old freeze;
- checks all 1181 coefficient classes satisfy `Knew = K/kappa`, with `kappa = 10000993/10000000`;
- checks the old C4 report covers the expected 1368 knot cells with no duplicates, 443754 accepted leaves, no reported failures, and internally consistent historical worker provenance;
- actually reruns old/new C1 and C2 at 192-bit Arb precision, comparing exact rational lower endpoints to `kappa` and `1`;
- independently expands coefficient orbits and integrates both kernels, comparing exact rational mass endpoints to the claim. No displayed decimal or float is a proof input.

The resulting JSON records the **actual current environment** separately from the old report's historical provenance. New machine paths or FLINT binary hashes need not match the old host. The vendored checker/core source contents must match those used by the historical report. The dependency pin is for repeatability, not an assertion that all compatible machines have identical binaries.

## What the default run does not recompute

The archived `records/original-full-report.json` is the complete output report of the original C1-C4 run. It contains per-knot-cell summaries, not each leaf's Bernstein coefficients and remainder proof. Therefore default verification of its coverage and hashes **does not numerically recheck every old C4 leaf**. It reuses the recorded full C4 run and transfers C4 to the new kernel by exact positive linearity, as explained in `PROOF.md`.

For a fresh, self-contained numerical C4 replay, run the unmodified vendored checker. The first command reruns all original conditions, including full C4. The second independently reruns all conditions for the scaled kernel, instead of using homogeneity:

```sh
.venv/bin/python -I -B vendor/hl_check.py certificates/original.json --claim 3.59123 --procs 15 --out runs/fresh-original-full.json
.venv/bin/python -I -B vendor/hl_check.py certificates/scaled.json --claim 3590867754339/1000000000000 --procs 15 --out runs/fresh-scaled-full.json
```

Use a smaller `--procs` if needed. These are optional larger runs and were not repeated during package validation. The original full run used 15 workers and completed with zero failures. A different host must rely on its new checker verdict and recorded provenance, not pretend to reproduce old timing or binary hashes.

## Contents and hash scope

`vendor/hl_check.py` and `vendor/hl_core.py` are byte-for-byte upstream copies, with Apache-2.0 license and source commit recorded in `provenance/UPSTREAM.json`. The wrapper `verify.py` is adapted from the independent review script; its original SHA and exact adaptation scope are recorded there.

`provenance/original-rescaling-freeze.json` is the unchanged original 13-file freeze manifest. Only five frozen artifacts are copied, mapped in `provenance/subset-bindings.json`; the original certificate and vendor sources have separate upstream bindings. The large process samples, absolute-path runner, old virtual environment, and unrelated outputs are omitted. The verifier checks this **subset**, not all 13 original artifacts. Historical absolute paths inside the unchanged reports are metadata only and are never opened by the portable verifier.

`package-manifest.json` binds static package inputs. `runs/` contains generated validation output and is excluded from that input manifest. `release-receipt.json` binds the delivered static manifest and the recorded portable validation run; it is an archival receipt, not an additional mathematical assumption. New runs should use another output filename.

## Sources

- T. S. Tomar, [paper](https://www.mefyi.com/papers/weak-type-squares/weak-type-squares.pdf), Theorem 2.2, Lemma 3.1, Sections 3-5 and Appendix A. Its Theorem 2.2 is Lin's comparison principle.
- [Upstream certificate/checker repository](https://github.com/tejasvi/square-maximal-weak-type-v4), commit `3af09fc359ce262f4a7f63d34b0d97e6d4d95d1a`.

The package supplies numerical verification and the short rescaling proof. It does not bundle or recompile the 2299-module Lean development, and does not claim to reprove the comparison principle or the correctness of Arb itself.
