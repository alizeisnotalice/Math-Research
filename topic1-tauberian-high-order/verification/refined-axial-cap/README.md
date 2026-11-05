# R051 refined axial-cap scalar kernel

This archive preserves the four successful computation files from `r051-computation` byte for byte. The frozen source is `scratch/r051-refined-axial-cap/refined_axial_cap.md` (SHA-256 `37bbb0ec63fac1dce33fb9883e14499488baef490341df32a6fdfb0da427c979`). The completed independent Sol R003 audit is archived in `reviews/051-refined-axial-cap.md`; this README records computation scope only.

## Wolfram check

Run `wolframscript -file axis_cap.wl`; `axis_cap.log` is the actual output. The general scalar cap inequality is proved in the source for arbitrary integer `n>=2` from monotonicity of the power derivative. The WL `Reduce` checks for normalized caps are only the cases `n=2,3` with `R=a=L=1`; they do not prove the all-`n` statement. Other scalar outputs check the positive-part identity, an affine majorant, low-branch power increments for `n=2,3`, and the normalization `a=v^n`.

## Fraction checks

Run `python3 axis_fraction.py`; `axis_fraction.json` is its actual output. It contains eight exact rational cases, all for `n=2,3`, checking the positive-part identity, cap inequality margin, and cases where `z` lies outside `L`. These examples do not establish the arbitrary-`n` derivative argument.

## Geometric scope

The scalar calculations do not prove existence of the actual geometric cap, indexing hypotheses, finite-overlap bounds, or non-overcounting of the multiplicities. In particular, the displayed estimate converting `sum_h (m_h/A_h) R_i^(n-1)` to `R_N/c` is conditional on the stated lower bounds for `A_h`, the identity `sum_h m_h=(a/v)R_N`, and the geometric indexing assumptions. An external `O(v^2)` weight then changes `R_N=O(v^8 S)` to `O(v^{10}S)`. These inputs must come from the proof, not from the finite scalar checks.

## Reproduction commands

```sh
wolframscript -file axis_cap.wl
python3 axis_fraction.py
```

No computation was rerun for this archive.

## Whitelist and SHA-256

- `axis_cap.wl` — 1060 bytes — SHA-256 `31a7672e72744b59aacfc3782b182fd46225a18e1b8a3fcb4a5b65ae077c7da1`
- `axis_cap.log` — 569 bytes — SHA-256 `05e7c9a440d2500bf25de8e36a0f9bb7a795d35d09d79054c994dceeefecc424`
- `axis_fraction.py` — 501 bytes — SHA-256 `58f08c631c7cb136a61df6dc90433e3a2fc9efb1914e997cd06141aa0a90e1cb`
- `axis_fraction.json` — 1439 bytes — SHA-256 `62abaa03953f0168cbd33f58b50dc06770fa181ea2718398d619bb3dc8fc8666`

