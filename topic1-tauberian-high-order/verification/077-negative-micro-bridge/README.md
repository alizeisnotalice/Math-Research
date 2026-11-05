# R077 computation package

Frozen proof: `../../proofs/077-negative-micro-bridge.md` (SHA256 `42abc4975857d15e9bff0899e01b262243194cecd4e97d287caee7e78b9544dc`). Review: `../../reviews/077-negative-micro-bridge.md` (SHA256 `90472e67dd358b5cc9fc85f077c289f2b0e9b7f71d69b185119f2fddd3607fae`).

The two Wolfram files check scalar identities used by the negative-shear bridge. `joint.wl/log` checks the general-n budget derivative through u^9, micro-mass normalization 2K_n, and a determinant-at-M_new/mass-at-d_plus identity. `joint_oldmass.wl/log` checks the separate determinant-at-M_old/mass-at-d_minus ordering; at n=3 the normalized scalar contributions are 4 and 2. `joint_fraction.py/json` provides three positive rational denominator/ratio samples. These computations do not establish the signed-shear geometry, intrinsic-host stability, or construction.

For fixed n>=3, the proof compares actual finite inputs and obtains micro-minus-macro coefficient binom(n,2) kappa_0^3/(256n) at scale a v^(n+6); with kappa_0=8n/3 this equals 2K_n. In dimension three, comparison with the same-parameter positive pure-macro family gives (6+o(1)) a v^9. This is a lower comparison, not an absolute C_n ninth-order coefficient or an upper bound, and it does not show micro-hole optimality.

Reproduce the Wolfram runs with `wolframscript -file joint.wl` and `wolframscript -file joint_oldmass.wl`; reproduce the rational samples with `python3 joint_fraction.py`. The captured outputs are the adjacent `.log` and `.json` files. The initial failed Wolfram attempt is intentionally omitted; staged logs are from corrected successful runs.

The manifest is `manifest.json` in this directory. Its paths are relative to the project root (`../../` from this package).
