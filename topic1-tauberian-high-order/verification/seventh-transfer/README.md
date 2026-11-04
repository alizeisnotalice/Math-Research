# R039 conditional seventh-coefficient transfer: computation archive

`coefficients.wl/log/json` records the exact conditional arithmetic for fixed dimensions `n=2,...,16`. For each dimension the stored coefficient is computed from the scalar-H condition via `2^(n+7) * (s_n - (n+26/9) K_n)`,
where `s_n` is the displayed coefficient from the expansion of `L_n`. The script checks the `L_n` series coefficient against its finite binomial formula, records the conditional coefficient, checks the sign expression for `n>=8`, and exports the 15 rows as JSON. The review-valid R039 proof retains condition `H`; this calculation does not prove `H` or turn the transfer into an unconditional geometric result.

The Wolfram log prints `zero dimensions: Missing[KeyAbsent,n]` because `Lookup` is applied to an empty selection of rows. This is an output-format artifact: it is neither a numerical zero nor an unresolved coefficient. The saved JSON has 15 rows with `zero: false`, and `summary.py/log` independently reports `zero dimensions: []`. The script explicitly checks `series match all: True` and `n>=8 sign check: True`; the `n>=8 sign proof components` line is a counterexample-set reduction whose result is `False` (empty set). The retained n=2 conditional coefficient is `-16015360/243`.

`summary.py/log` is a portable Python summary using the adjacent `coefficients.json`. It checks the row count, stored series flags, zero list, distinct kernel/process counts, and the n=2 regression. `verify_fractions.py/log` is an additional independent Fraction check of all 15 stored rows. The original proof-side verifier assumed the original directory layout; this archive contains Terra’s explicitly modified portable copy, whose only functional change is to use an adjacent default JSON path and support `--input`. It was actually run in `r039-portable-check` against the same JSON (exit status 0); its log is preserved. No old `j<=6` coefficient table or redundant precheck log is included.

## Scope

These files verify a finite set of exact coefficient calculations and the stored Wolfram output. They do not prove the scalar condition `H`, geometric supply, a global all-input seventh-order bound, or an unconditional Tauberian coefficient. The independent R039 review validates the frozen transfer within its stated conditional scope; condition `H` remains explicit.

## Reproduction

From this directory:

```sh
wolframscript -file coefficients.wl
python3 summary.py
python3 verify_fractions.py
```

The Wolfram script writes `coefficients.json` beside itself; run it in a disposable copy if preserving the archived JSON. Both Python scripts read the adjacent JSON. The `verify_fractions.log` records Terra’s actual run of the modified portable copy; the original layout-coupled verifier was not copied or rerun.

## Whitelist and SHA-256

- `coefficients.wl` — 1470 bytes — SHA-256 `e2d592678c7577888cb4b0d6750caae6d3b230d16fe1bcc551cbe7915bc23f3b`
- `coefficients.log` — 514 bytes — SHA-256 `c19013ed4488e2b34919f2680b528a38db9855c8bf240fdd686bf96ce7b4a3ea`
- `coefficients.json` — 2755 bytes — SHA-256 `1c17308eb384c7d0f4f19ae4fc7f35249e1f9ff8433fbd1995fb6dc38a8b4dd9`
- `summary.py` — 440 bytes — SHA-256 `3c773c8d2cb3f7f4a6a5feb6927cd99f250d74197d86420d27a4935b9137fd60`
- `summary.log` — 110 bytes — SHA-256 `ecf65cdfcab7bf44b7f78de1c3ab11314549e6805230c87ae072350a7864465f`
- `verify_fractions.py` — 1378 bytes — SHA-256 `345853683eb88561358aae8dc8ad43700f624cb6ccb92ee34bb63a02aec6601a`
- `verify_fractions.log` — 617 bytes — SHA-256 `e8c0cd45f94b314d3cff47c629ce60da4b48606079b54bfb7974aa91fe852b7e`
