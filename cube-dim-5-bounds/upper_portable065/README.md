# upper_portable065 portable extraction

This directory contains a compact, reproducible extraction of the first 2048 accepted rows in the frozen `upper_cover061` state. The bundle carries each node, its four integer root coordinates, the author-saved box bounds, the matching cache point and five exact rational interval pairs (the H interval followed by four physical-gradient intervals), the stored precision bits, the original failed-certificate SHA-256 and old minimum, and the complete new certificate object.

The builder checks the three frozen input hashes, row-by-row equality with the final 2048 `passed` entries, and the cache point/key match. It copies rational strings as stored and does not recalculate H, gradients, logarithms, or margins. The underlying H ball remains a shared production input; this bundle is not an independent value certification and does not establish a global `g` or `C5` bound.

## Reproduce

Use Python 3 and the standard library. Supply the frozen state and certificate paths explicitly, and choose a new output directory; existing output directories are refused. The owner manifest is pinned to its recorded absolute path and hash in the builder.

```sh
python3 build_portable_bundle.py \
  --state /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_cover061/batch0001/state.json \
  --certificates /Users/zhengzhihao/Desktop/T/cube-finite-dim-scratch/upper_cover061/batch0001/certificates.json \
  --output-dir /path/to/new/output-directory
```

The gzip payload uses a zero modification-time header. A successful run writes `registration.json`, `records.json.gz`, and `manifest.json` to the requested output directory. A failed run preserves the registration and writes `failure.json` with the stage and reason.

## Frozen run and verification

- Registered and executed once as PID `60429`; elapsed time was `1.993213` seconds, below the registered 120-second budget.
- Extracted records: `2048`.
- Compressed bundle: `1,703,698` bytes; SHA-256 `e25c6cbcad0fabcdc8df5ade205de098c8c952e18f73a9d57a1d888cc40291f4`.
- gzip mtime header: `00000000`.
- Source hashes were unchanged after extraction: state `d0d62d99d57b3f15d3acb1f4a635b188e196b0de6e155b461460ccbd64a72d31`; certificates `a3331af02c0a3df5548977971b85e62f2dac3b39f591658289569ab4f402b151`; owner manifest `2e1c6a69d9fd55ff7c09d65cd356c88bf726151f07744c6b7311c4d43cda3504`.
- Read-back checks confirmed 2048 records, four root integers and box coordinates per row, five stored interval pairs, four directed gradient intervals per certificate, and the registered PID.

## File whitelist

Keep only these files for this extraction:

- `build_portable_bundle.py` — standard-library builder.
- `README.md` — scope, command, run evidence, limitations, and whitelist.
- `replay_001/registration.json` — pre-extraction registration, completed status, and PID.
- `replay_001/records.json.gz` — the compact 2048-row portable data bundle.
- `replay_001/manifest.json` — input hashes, builder hash, output hash/size, and scope.

The full state, cache, journal, and unused records are not copied into this directory.
