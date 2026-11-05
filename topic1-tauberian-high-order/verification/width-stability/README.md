# R059 width stability: local scalar computation evidence

Frozen source: `scratch/r059-width-stability/width_rms_from_actual_budget.md`, SHA-256 `e16380a6aa2192654a3f5ff51a7773275d4ee132dbe18d1a7c28699fcc03beed`.

Independent mathematical review: `reviews/059-width-stability.md`, SHA-256 `4825a4898789d31689a67d220136acc6048719bae08b3ca59c70bfb8008f5c81` (PROOF-VALID, same N8 scope).

## Recorded Wolfram checks

`width_amgm.wl/log` record Mathematica 15.0 ARM, PID 32488. The output gives zero for the stated Hessian residual, `False` for the violation set of the cross-multiplied inequality on `0<=lo<=xx<=up`, zero for its positive-term decomposition residual, and zero for the two-dimensional `cosh` gap residual.

Reproduce with `wolframscript -file width_amgm.wl`. The script and successful output are copied unchanged; the initial invalid run is excluded.

## Exact rational spot checks

`width_fraction.py/json` use Python `Fraction` arithmetic on three specified inputs. All assumptions are true, and the reported margins are `0`, `275/2304`, and `359/1600`; all pass. These finite examples check rational transcription only.

Reproduce with `python3 width_fraction.py`.

## Scope boundary

The computation checks local scalar identities and a few rational cases. It does not replace the arbitrary-dimensional geometric proof, show a global ninth-order Tauberian remainder, or promote this local width estimate to a global coefficient. All inputs, dimensions, and geometric interfaces must be handled by the mathematical source.

## Files and SHA-256
- `width_amgm.wl` — 688 bytes — SHA-256 `3457b2a4d32942e68522debb4aaf9284d4c44aacef451e2637ed5c1584ba60a4`
- `width_amgm.log` — 283 bytes — SHA-256 `4e2e73f3192d16996817b4ae076ea38b183a20ffdbe815e7c2dbf2ef1705628b`
- `width_fraction.py` — 679 bytes — SHA-256 `85f7c5f4f69771e1a1634b13416b6470c8faefa53cfca806a84e4c7a9d7fe25d`
- `width_fraction.json` — 715 bytes — SHA-256 `46f197b6df16a4742a1b3701ac87b6303b1f23b983c7ada053eafb7038d1651d`
