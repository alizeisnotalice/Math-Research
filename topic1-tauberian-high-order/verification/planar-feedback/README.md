# R28 planar feedback: minimal computation archive

This final package retains the clipped thin-arm degeneration check and the independent price-check Wolfram run. The earlier general `feedback.wl/log` and `exact_examples.py/json/log` remain byte-preserved in the sibling `planar-feedback-superseded/` directory outside this package; the local proof and these two focused checks made them redundant for this archive.

## Scope and evidence

`clipped.wl/log` computes the true clipped two-arm integral in the family `0<eta<1/2` and compares it with the formal untruncated kernel. The true normalized integral is
\[
\eta/2+\eta^3/3-\eta^4/6\longrightarrow 0,
\]
whereas the formal expression `1/3+3 eta^2/2-eta^6/6` tends to `1/3`. This is an analytic degeneration counterexample to uniform use of the formal thin-arm kernel.

`price_check.wl/log` checks the Tonelli kernel identity, the exact weighted two-arm bound, the scalar optimization and its upper-bound inequality, and the coefficient conversion. These are integral and scalar-algebra checks; the package alone does not prove the geometry or any global `C_2` result.

The two independent audits, `r028-review-rigidity/audit.md` and `r028-review-price/audit.md`, both report PROOF-VALID within their stated scopes. Together they validate the full planar local candidate. This is a local planar price result, not a global upper bound or a determination of the seventh-order coefficient of `C_2`.

## Additional recorded Wolfram check

The following rigidity check was run by the reviewer through `wolframscript -code`; no standalone script or log was saved. The command and output below are transcribed from that run and were not rerun for this archive:

```sh
wolframscript -code 'Clear[u,e,h,y]; Print[$Version]; Print[Factor[4/27-u^2(2-u)/8]]; cut=e h/(e+h); inter=Integrate[y,{y,0,cut}]+Integrate[h(1-y/e),{y,cut,e}]; Print[FullSimplify[inter,Assumptions->0<h<e]]; Print[FullSimplify[e h-2inter,Assumptions->0<h<e]]; Print[FullSimplify[e h-2inter<=h^2,Assumptions->0<h<e]]'
```

Recorded output: Wolfram 15.0.0 ARM; `((-4+3*u)^2*(2+3*u))/216`; `(e^2*h)/(2*(e+h))`; `(e*h^2)/(e+h)`; `True`; final `Null`, exit status 0.

## Reproduction

From this directory:

```sh
wolframscript -file clipped.wl
wolframscript -file price_check.wl
```

Compare each run with its `.log` file. The recorded `-code` check above has no archived script and is included only as an existing reviewer output record.

## Whitelist and SHA-256

- `clipped.wl` — 1644 bytes — SHA-256 `bcdb7a652443339745cad35081e06a814507b13e34522f6d2febbe3c16a51dbe` — source `r028-scan/`
- `clipped.log` — 622 bytes — SHA-256 `ceb7b95510db657752a838d15f6c9f24da248d43576b4f166405da7ac1f83404` — source `r028-scan/`
- `price_check.wl` — 872 bytes — SHA-256 `5e8a200dd15f812ca180a409342412c783627c369477da1ff69e4cea3c443cf7` — source `r028-review-price/`
- `price_check.log` — 325 bytes — SHA-256 `f3fc93cfb3c41d6ace2f1bfb1037f65ab9b4b88017c724f80425bc053d21b6c3` — source `r028-review-price/`
