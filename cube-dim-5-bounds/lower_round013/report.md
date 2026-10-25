# Lower round 013 execution report

Six new body-centered configurations were registered and actually evaluated. All finite-family upper endpoints lie strictly below the frozen complete algebraic C_cert. No new bound was found, and no ledger update or commit was made.

| Candidate | Exact inner ratio, decimal display | Nodes | Seconds |
|---|---:|---:|---:|
| anisotropic_9_8 | 1.46756612206647431165267631898537011772209059944583884799726 | 70225 | 2.744589833 |
| anisotropic_19_16 | 1.31661822276524337688090040793379268555290032678629359245254 | 14620 | 1.039026958 |
| anisotropic_5_4 | 1.15418314071507262321366922546364148928631377362787738511803 | 4010 | 0.693445458 |
| isotropic_11_8 | 1.17109487367772260542530483407095394005444154124935673215153 | 35887 | 1.399827375 |
| isotropic_23_16 | 1.10235223675498538710769957276902400309918871180452038303277 | 3543 | 0.845927708 |
| isotropic_3_2 | 1.06328855825815807185638488690411087566288013245824541925367 | 3370 | 0.792210667 |

For every row the exact equality-radius finite-family ratio is at most the displayed inner value plus `1280*102048/2^80`. This rational error was added and the resulting exact upper endpoint was actually compared to C_cert's lower isolating endpoint. The comparison is not against a rounded display or a safety threshold. Exact fractions and all hashes are in `completion_receipt.json`.

All six runs processed 102048 finite configurations. Direct distinct-atom enumeration and all ten geometry endpoint inequalities passed for every retained box. The recursive slicing and coordinate-mask union engines agreed exactly. All jobs stayed inside their own 120-second/200000-node caps.

`mechanism_and_scope.md` gives the point set, exact subset mass, single common radius, reflection/period normalization, finite-family error bound, and explicit width-two collar plus smooth L1 transfer. The omitted witness blocks and spacings remain unresolved. This is a finite six-point exclusion, not a global BCC obstruction.

Reproduce from the repository root:

```sh
python3 cube-dim-5-bounds/lower_round013/bcc_finite.py anisotropic_9_8
python3 cube-dim-5-bounds/lower_round013/bcc_finite.py anisotropic_19_16
python3 cube-dim-5-bounds/lower_round013/bcc_finite.py anisotropic_5_4
python3 cube-dim-5-bounds/lower_round013/bcc_finite.py isotropic_11_8
python3 cube-dim-5-bounds/lower_round013/bcc_finite.py isotropic_23_16
python3 cube-dim-5-bounds/lower_round013/bcc_finite.py isotropic_3_2
```

Source SHA256: `f4ae984271e4c404a77b348f51c795897069552966f09a77bff25486dfae562a`. Full-reference SHA256: `f385811e5460af63abc117804965461eed88d8ecb5b9d2cd30268387646f45d2`.

The next mechanism must gain substantial geometry or coverage; the best value here is about 1.46757, well separated from C_cert. A finer spacing mesh is not supported by this result and was not launched.
