# R24 computation archive: 2D rational lower family

This archive contains only the calculation inputs and outputs retained for the final rational family. The older fitting, stationary-point, IFT, branch-table, and obsolete certificate materials were intentionally excluded; the rational family does not depend on them.

## Scope

`family_asymptotics.wl/log` records the symbolic series for the chosen rational family. `finite_family.py/json/log` records one exact-Fraction sample at `u=1/1024`; `finite_family_check.wl/log` is its high-precision Wolfram cross-check. These finite checks establish the listed sample guards and sample margins only. `extra_guards.py/log` adds exact cached-root checks for the stronger 18-fold double-corner guard, the slice right-hand side bounds, the prescribed `N`, and the rational width identity; it reuses the cached 400-step `h` interval and performs no new root search. It does not calculate a numerical value of `c`: it certifies that the mass equation has a unique solution in the stated slice, hence existence and uniqueness of the mass quantile. None of these finite checks is the all-small-parameter proof.

The reviewed proof separately establishes the analytic family for every sufficiently small real `u`, chooses the unique `z` from the mass equation `H(z)=D-2A`, and sets `c=w-z`. Thus `c` is a mass quantile selecting the actual hole of prescribed mass; it is not a second computed root of the budget equation. The proof's actual center-average estimate and convex-box argument control the full strict halo over all radii, which is distinct from the finite sample guard checks and from the periodic halo formula.

The result is a one-sided lower bound for `C_2` from a concrete family. It does not identify the true next coefficient or assert an upper bound/global optimum.

## Whitelist and SHA-256

- `family_asymptotics.wl` — 431 bytes — SHA-256 `1384c5f79f08578cfaeb59d14e271cce7e34fdc903cb6e3b8c5d6d3b9b38f9e3`
- `family_asymptotics.log` — 123 bytes — SHA-256 `1fa90248a0dbc8ed462da8a7b43c038ad8076ca45507da077cd9c9c5c6df90a4`
- `finite_family.py` — 2822 bytes — SHA-256 `918ff04fa64e5fb3c7a69258ce522a4a66bc66639fec8954ab5df46a14ff1396`
- `finite_family.json` — 2912 bytes — SHA-256 `ee848252d0cd3999a03676b75abfcc428535519ca01d166f89046f3aa59702fc`
- `finite_family.log` — 609 bytes — SHA-256 `021842730f8aedf24d5ed487c62f50d1d0baf31eb8eb0a44212d42f051070e07`
- `finite_family_check.wl` — 1673 bytes — SHA-256 `f2b0d93f6f1891f410021048190b67906b0e39aefc9b7ddd9e39e92278d63ede`
- `finite_family_check.log` — 435 bytes — SHA-256 `142a256652cb0f7e0d8dcbde310685fe94e5cbe34efeb009ad37d0b406fb5c3d`
- `extra_guards.py` — 962 bytes — SHA-256 `5b8acd650d2e462aebf30bb789d37eed91a67e468fa531b6f2dbeb8362a602a4`
- `extra_guards.log` — 226 bytes — SHA-256 `ed8157c36fbc3aec274c1dfa1eeef4e956ad20531e14c06849528b14f6b7ad51`
