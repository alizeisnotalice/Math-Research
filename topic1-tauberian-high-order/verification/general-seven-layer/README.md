# R038 general seventh-layer computation archive

`general_seven.wl/log` checks finite symbolic algebra for the proposed seventh-layer weighted sum under the stated fixed-parameter assumptions (`n>=2` fixed, `0<v<1`, `c>0`; formal expansions hold for fixed constants). It verifies the telescoping residual, endpoint term, truncated weights, weighted-series expansion, exact ratio factorization, and coefficient identity after writing `a=v^2 w`, `0<w<=1`.

The last identity gives conditional coefficient pair `K` and `-(n+26/9)K`; substituting `a=v^n` means `w=v^(n-2) in (0,1]` for fixed `n>=2`. This computation is conditional algebra. The R038 independent review validates the separate supplying argument that removes the condition for its stated result; the computation itself does not prove the geometric supply or replace that argument.

A display-encoding artifact in the log prints `vÂ²` in one label; the source expression is `v^2 w`. This does not change the Wolfram output `k - (k*(26+9*n)*v)/9` or the subsequent `True` coefficient check. The `True` is the result of that finite symbolic identity, not a geometric theorem.

## Reproduction

From this directory:

```sh
wolframscript -file general_seven.wl
```

Compare output with `general_seven.log`. No computation was rerun for this archive.

## Whitelist and SHA-256

- `general_seven.wl` — 1356 bytes — SHA-256 `2fe2e3c77207f7ed5bf459fd55e93961ac52ff463935db7d528780c48228bcf1`
- `general_seven.log` — 561 bytes — SHA-256 `98860cd23430ad0e50e7ee8a86754c1c7c7c9aefe40f90638ab98e3fe093cb2b`
