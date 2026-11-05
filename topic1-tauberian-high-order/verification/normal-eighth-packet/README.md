# R053 normal eighth-packet: scalar and finite quadrature evidence

Frozen mathematical source: `scratch/r053-normal-eighth-packet/actual_normal_eighth_price.md`, SHA-256 `d83dc772fc2c1ca5698a39eb2bd15702395d5e8238b11209fd2b6e7765b1e642`. The formal reviews are `scratch/r053-review/audit.md` (SHA-256 `a6ebc52f1274e1d5cb00eadc7aad1ce8dd791d8713006a880f878e2ebee8cbf0`) and `scratch/r053-review-adversarial/audit.md` (SHA-256 `7ab3b0eec28e61df3f5324bbef8c0425277343d80379c59da13e84d0f846a0f0`); this package records computations only. The public mathematical source is `proofs/053-normal-eighth-packet.md`; reviews are `reviews/053-normal-eighth-packet.md` and `reviews/053-normal-adversarial.md`.

## Exact symbolic checks

`asymmetric_cap.wl/log` treats exactly the displayed two-branch model on `0<eps<1/10`, `0<=z<=1`, with `p=1-eps z`, `q=1+eps z` on the plus branch and reversed on the minus branch, `f=1/eps`. The recorded Wolfram 15 ARM run gives the two branch series through order eps^2, zero eps^2 coefficient in z for their sum, and `False` for the minus-branch ordering violation set. This is the stated narrow parameter range, not a uniform geometric estimate.

`local_optimum.wl/log` is a formal scalar model with bounded fields and fixed `n>=2`: it checks a square completion, the quadratic optimizer, the closed formal `J_n` identity for `n>=3`, and the separate `n=2` denominator correction. The `8/27` appearing as that R41 `J_2` correction is distinct from the additional R54 micro-hole `8/27`; the latter is not computed by this package.

`check_crosscap_identity.wl/log` is the independent exact normalization kernel. Its recorded residuals are zero for the sorted-fibre identity (8a), nominal opposite-root `F`, and all four source/refund scalar conversions; the double cross-hat difference is even and divisible by `c^2`. These are finite algebraic identities, not uniform geometry. They do not prove the general unequal-`P,Q` chart geometry or any uniform remainder.

## Numerical counterprobe and reproduction

`u1_counterprobe.py` defines isotropic square-hat inputs and computes both orientations; it does not cover general unequal `P,Q` end widths. `u1_counterprobe.json` records 18 float64 SciPy quadrature / polygon-clipping cases (`tau=0,1,2,5,10,20`; `eps=1e-3,2.5e-4,6.25e-5`), run with `ProcessPoolExecutor(max_workers=6)`. The logged worker PIDs are 16539–16544. Its `estimated_quad_error_scaled` fields are numerical estimates, not certified bounds or total strict remainders.

`u1_small_eps.json` records two additional cases, `tau=10,20`, `eps=1.5625e-5`, from a separately recorded call to the same `run` function (worker PID 16585); the original command record is not copied because it contains a private absolute path. From this directory, the portable equivalent is:

```sh
python3 -c 'import json; from u1_counterprobe import run; print(json.dumps([run((t,1.5625e-5)) for t in (10,20)], indent=2))'
```

The preserved JSON is the actual prior output; no computation was rerun for this archive. The 18+2 floating-point checks are exploratory only, not proof or certification.

## Scope boundary

Nothing here certifies a uniform geometric construction, cross-cap noninteraction, or a Tauberian conclusion. Exact symbolic identities and numerical quadrature have distinct roles; the latter's error estimates are not proof bounds. The package does not establish a global geometric upper bound or identify an exact coefficient beyond the source's stated scope.

## Reproduction commands

```sh
 wolframscript -file asymmetric_cap.wl
 wolframscript -file local_optimum.wl
 wolframscript -file check_crosscap_identity.wl
 python3 u1_counterprobe.py
 python3 -c 'import json; from u1_counterprobe import run; print(json.dumps([run((t,1.5625e-5)) for t in (10,20)], indent=2))'
```

## Files and SHA-256

- `asymmetric_cap.wl` — 1145 bytes — SHA-256 `5b8fe3c043f7ca3c1c8ed56d3645edea18193f675e6c0e012becfed6969d91eb`
- `asymmetric_cap.log` — 335 bytes — SHA-256 `dd33d78dc2faf3a728bc981f92376687854b6e445ecc50d10eb64a97152bc372`
- `local_optimum.wl` — 1066 bytes — SHA-256 `5d7e4c706600710e3caa925e400f036e864ea9ea84564beb9616d56de6322280`
- `local_optimum.log` — 507 bytes — SHA-256 `e412508dfda5d615604305d38541e4fcf12c0d4093769615ab90dac462ad6219`
- `u1_counterprobe.py` — 2072 bytes — SHA-256 `d662ab4a34a615d0facf247a5a499a655de3ff98e2e8b086e3679bce12e0e11d`
- `u1_counterprobe.json` — 6925 bytes — SHA-256 `94c4d527e191c5100467764275d12c9f7ca0ce3248ec6e34fb3d14a1d7b0b932`
- `u1_small_eps.json` — 648 bytes — SHA-256 `689a3f5abe48400990e07584f7fc0f54ab78c2f5cc7a4d4847fa72eb6d758d1d`
- `check_crosscap_identity.wl` — 1779 bytes — SHA-256 `ff4bf287f8e5d43881da291020106a43a29794add9d5e2444f9c14330321adf9`
- `check_crosscap_identity.log` — 449 bytes — SHA-256 `436fb558aaa449f74dd5c97ef4cf009b37422e48c2a5638b97580e62d54538ba`
