# R25 computation archive: n-dimensional candidate lower family

This archive keeps the general symbolic calculation and the three exact-Fraction regression cases. `general_next.log` includes the general identities and independent subkernel regressions for dimensions 2, 3, 4, and 5, including the requested 4D check. The Fraction cases are `n=2,3,4` and use the supplied three-face geometry estimate as an input; they are finite regressions, not independent geometric proofs.

The computation and reviewed argument support only a one-sided lower construction. They do not determine the actual next asymptotic term of `C_n` or establish a matching upper bound. The additional triple_constant.wl/log pair was actually run locally: the normalized three-face leading constant and the residual slack identity are exact for integer n>=3, with n=2 evaluated separately as zero. This checks algebra only; the geometric three-face covering remains part of the independently reviewed proof.

Excluded: the known failed initial log (`Series::ivar`) and duplicate `general_next_results.txt` output. The retained `general_next.log` is the corrected successful run.

## Whitelist and SHA-256

- `general_next.wl` — 2682 bytes — SHA-256 `e2fcaba58ce99eb79aa22f5ade958139ac0366a4fae5a6c5366599bb9257e88e`
- `general_next.log` — 1324 bytes — SHA-256 `4b8dbe641ae3890b1f0040ab0014096fc6f873a1254d6dfae37f36748bc0e024`
- `finite_examples.py` — 2343 bytes — SHA-256 `e2313c997f45470cd40241b6b0067484105bd04bf1a0a18d4766faaab65e0416`
- `finite_examples.json` — 3639 bytes — SHA-256 `4a3c3bbc8de4e55fe873b75805b572d668826e635a6d31ad84071c84c4941ae1`
- `finite_examples.log` — 101 bytes — SHA-256 `e095696f87cb7f97ca9d2d80520941fff6d857723797f242083ac6231c4cd05b`

## Reproduction

Run WolframScript 15 ARM with wolframscript -file general_next.wl and wolframscript -file triple_constant.wl; run python3 finite_examples.py from this directory. The symbolic batch used four actual subkernels for dimensions 2 through 5. The three Fraction cases are periodic budget tests, not certified finite truncations. Their safety coefficient is 10^6, which need not equal the particular c_n in the general theorem.

- triple_constant.wl: 541 bytes, SHA-256 f0473816228f8743d8fc1d43d71ab9ec2876cea0998d0618097fe7ce9549cbdd

- triple_constant.log: 229 bytes, SHA-256 e9dc66c889a393dbb8cb1f493d108f831024412c530e00abf53910df22724fe4
