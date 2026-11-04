# R041 eighth-order feedback: computation archive

This package retains the exact two-arm moments, conditional denominator/budget algebra, all outer-box support branches, and the seven rational boundary inputs cited for the frozen lower-construction calculation. It keeps only the scripts and recorded outputs needed to reproduce those items.

## Exact moments and refund

`refund.wl/log` assumes `0<d<1/2` and `theta=d/(1-d)`. It records the two clipped-arm moment branches, the exact baseline `I0=(5-2d)d/3`, the feedback difference `I0-I1`, its expansion, and the triple integral `e^4/(4f)`. `boundary_fraction.py/json` independently integrates the same piecewise polynomial with exact Fractions at `d=1/4,1/10,1/100` and checks the exact closed form. Its other four rows check the inside/outside support identities and positivity for `n=3,4`, `u=1/16,1/64`. These are finite algebra checks; no actual geometric deleted-set measure is computed by them.

## Budget and exact original denominator

`budget.wl/log` performs the formal fixed-dimension coefficient optimization under the explicitly stated asymptotic-remainder and triple-bound assumptions. Its `k1=-544 n/81` is the coefficient of `u` when expanding `kappa` around the bare constant `8n/3`. `exact_original.wl` instead writes `kappa=(8n/3)g^2+(320n/81)g^3 u`; expanding the moving baseline `kappa_*=(8n/3)g^2` gives the same total coefficient `-544n/81`. `exact_original.wl/log` also checks the exact budget remainder through `u^9`, numerator series, equality of the closed forms for `J_n`, and the special `n=2` inverse-denominator correction and coefficient `15694/6561`. These are conditional algebraic calculations, not a proof that the geometric hypotheses hold.

The construction target is `J_n=(n-1)n^3(2609+1971n+81n^2)/26244 + (8/27) 1_{n=2}` (the polynomial formula alone applies for `n>=3`), with `J_2=15694/6561`. The specified finite-family halo ratio itself has the stated expansion through `J_n a v^8` with remainder `O_n(a v^9)`; passing to the supremum `C_n` yields only the corresponding one-sided lower limit after subtracting through `a v^7` and normalizing by `a v^8`. This is not an all-input eighth-order coefficient.

## Outer support branches

`support_branches.wl/log/json` covers fixed dimensions `n=3,...,17`, balanced even supports `j=2,4,...,n` when applicable: 15 parallel subkernels and 71 support rows. The log records kernel IDs `{1,...,15}` and actual PIDs `{92323,92324,92325,92326,92327,92328,92329,92330,92331,92332,92333,92334,92335,92336,92337}`; the two `j=2` identities are `True`. The JSON stores each branch's leading coefficient and expected order; for `j>=4`, the minimum inside/outside exponents are 13 and 15. This is an outer-box algebraic support calculation, not a geometric proof of the support decomposition or a measure bound for the actual deleted set.

## Seven exact rational inputs

The seven rows in `boundary_fraction.json` comprise three `d`-moment cases and four support cases (`n=3,4`, `u=1/16,1/64`). All recorded exact continuity, closed-form, identity, and positivity flags are true. The Fraction script prints the JSON to standard output; the JSON is the recorded output. It does not compute the asymptotic limit or replace the symbolic argument.

## Geometric and asymptotic scope

All computed coefficients depend on the scalar construction assumptions stated in the frozen proof, including its budget remainder and triple-intersection estimate. Geometry, strict halo inclusion, finite-box localization, and the proof that the chosen input realizes the stated coefficient remain mathematical proof obligations; these computations do not supply them. No claim about a matching upper bound or global eighth-order asymptotic follows from this archive.

## Reproduction

From this directory:

```sh
wolframscript -file refund.wl
wolframscript -file budget.wl
wolframscript -file exact_original.wl
wolframscript -file support_branches.wl
python3 boundary_fraction.py
```

The proof finite-box step uses `N=ceil(u^-10)` and records a truncation error `O_n(a/N)=O_n(a u^10)`; this normalized error comes from the geometric proof, not from the arithmetic scripts. The parallel support run used 15 kernels. Its JSON export writes beside the WL script, so use a disposable copy if preserving the archived output. The Fraction script prints its JSON rather than writing a file. No computations were rerun for this archive.

## Whitelist and SHA-256

- `refund.wl` — 699 bytes — SHA-256 `1e6f958e6ec03d4bed9ce619e910d7b918e3d65beb1ed46128079e0956e447da`
- `refund.log` — 387 bytes — SHA-256 `ed24a98984ad22298457b9de0f21654ad5bc39caad5f937edf20a2b03a59d2ac`
- `budget.wl` — 1063 bytes — SHA-256 `78c9e6f0a92c5b936574717b9868961a89952c38cd226850ab65d49ec94d0de8`
- `budget.log` — 650 bytes — SHA-256 `b494e377dde697d69d1e05627d0842253e07f8cad6d84b38c2b92638531f0a3b`
- `exact_original.wl` — 1333 bytes — SHA-256 `ffd27d042eec143d5ab86b67bd7ee3904120cd34acaddecb9e06d5e5ab5bfe61`
- `exact_original.log` — 520 bytes — SHA-256 `165548ff743de46835109c6b3a3a3ee643c593bf4a27fd2d6c2c2401ecd0ccb4`
- `support_branches.wl` — 1764 bytes — SHA-256 `2dfeb441f0e59281ff92cbb07626d2dfe7718a1941808e8a62006b03d2dd0587`
- `support_branches.log` — 477 bytes — SHA-256 `f4439b0edcfe27f6b7a7e24d7dc56336a878d1ff894c1797e23853891503cb7b`
- `support_branches.json` — 20288 bytes — SHA-256 `899c7b01df4c7133e3bd8ca4b672a85b892ebbab6fff9af939f3b09df77d5185`
- `boundary_fraction.py` — 1513 bytes — SHA-256 `5737ddb22879bb399c984a6e76d594367348a8af614f70e4e2ccb2c8195f48a2`
- `boundary_fraction.json` — 1102 bytes — SHA-256 `faa814d8d3d626f6938490b7033ffece8c56707fe105d3fcb704509d3c8b5907`
