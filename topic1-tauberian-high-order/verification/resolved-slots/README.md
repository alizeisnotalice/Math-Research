# R035 resolved-slot helper probe: arithmetic archive

`helper_probe.wl/log` checks only the fixed probe domain `0<eps<1/100`, `0<=del<=eps/100`. It expands both candidate core-area expressions and asks Wolfram to reduce the sets where either is strictly below, or no greater than, `2 eps^2`. Both violation sets return `False`, so both expressions are strictly above the threshold throughout the stated domain. This is the local algebraic input for the resolved-slot argument.

The probe supports only the fixed-`tau` resolved conclusion (the resolved edges have `min(L_1,L_2) >= tau v^3 R_e`). It does not establish a rawSlot estimate, cover edges whose minimum intersection length degenerates below that scale, or prove a full n-dimensional seventh-order bound. The independent review reports PROOF-VALID within the fixed-`tau` resolved scope; these computations are not a geometric proof.

The `False` outputs are empty violation sets, not failed checks. No unresolved Wolfram result is presented as a successful proof. The geometric helper existence, uniqueness, and allocation are supplied by the proof and independent review, not by this probe.

## Reproduction

From this directory:

```sh
wolframscript -file helper_probe.wl
```

Compare the output with `helper_probe.log`. No computation was rerun for this archive.

## Whitelist and SHA-256

- `helper_probe.wl` — 424 bytes — SHA-256 `b8b2da058dd95364e6b8af2a96107b257f25fd791ee85fcabdf22fb32d3f5f55`
- `helper_probe.log` — 229 bytes — SHA-256 `bf983e602693fa8004db4ac5cdf7a3379ae36c54e6e9371d4448020846d0f7fc`
