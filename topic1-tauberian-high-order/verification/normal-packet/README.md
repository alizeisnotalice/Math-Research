# R042 normal-packet and tail-obstruction: minimal computation archive

This package contains the two successful Wolfram checks cited by the frozen source: `coercivity.wl/log` and `three_cycle.wl/log`. The failed first `three_cycle` run using the protected Mathematica symbol `I` (`three_cycle-first-invalid.log`) is intentionally excluded. The manuscript’s displayed two-arm expression is a two-dimensional cross-sectional moment after proposing the target transverse width `T_j^perp`; it is not itself a full-dimensional measured slot, and the retained computations do not verify that geometric realization.

## Regular-packet scalar coercivity

`coercivity.wl/log` checks the stated compact domain `1/2<=b<=1`, `1/6<=x,y<=1`, `|z|<=1/12`, with `c=1/1000`. The displayed `False` values are empty violation sets: the positive-part domain constraints, the `z`-coefficient violation set, and the endpoint (`b=1/2,1`) coercivity violation sets. It also records the exact decomposition `G-G0=3(x+y)z^2` and confirms the remainder is affine in `b`. The manuscript gives a separate rational hand proof for this local packet inequality; these symbolic checks do not supply its geometry or prove the actual packet selection.

## Three-cycle calibration

`three_cycle.wl/log` explicitly assumes the scalar model `I/u^(15/2)=kappa^3`, `J/u^(15/2)=kappa^5/64`, `kappa>0`. Under that assumption, at `kappa=6` it returns `(I0,J0,I0/J0,I0-J0)=(216,243/2,16/9,189/2)`. The Young-square completion has residual zero and its violation set is empty (`False`). These `I0,J0` are safe determinant/budget scalar quantities; they are not the greedy or global quantities `I,J`, and the four numbers are not a global counterexample or an identification of the actual greedy packet sums. The source separately constructs a safe periodic hole family, but this computation does not verify that geometric construction.

## Scope boundary

The source discusses a regular-window packet lemma and leaves a common global budget interface unresolved. Its `N8` statement is conditional. These computations neither close that missing interface nor prove an unconditional global eighth-order upper bound. They are finite-domain scalar checks only, not geometric evidence that can replace the packet proof, actual halo inclusion, or tail accounting.

## Reproduction

From this directory:

```sh
wolframscript -file coercivity.wl
wolframscript -file three_cycle.wl
```

Compare output with the adjacent `.log` files. No computation was rerun for this archive.

## Whitelist and SHA-256

- `coercivity.wl` — 723 bytes — SHA-256 `c2184351ed2a6e2d72d28244338706a57e8c22dec8b565b3b8880669c444282b`
- `coercivity.log` — 283 bytes — SHA-256 `ca21c6dcc208a3dd19997cf616a2d749125dce7afe7ac806ce69d7e142a30b64`
- `three_cycle.wl` — 452 bytes — SHA-256 `46416d3c28cf89fd78427b09fc50c15e81543418ea033de6cbd50dd6851faec6`
- `three_cycle.log` — 241 bytes — SHA-256 `b8df0ed425f17102a3dd957f4909a7ab2407483f69c9d6e99dc33c5d03eb4558`

Frozen source `normal_packet_and_tail_obstruction.md` SHA-256: `dd291e7bc4294483571eaa8dcfb85e8e409266fbfb8eb8fb59ebcea6b9e2dcff`.
