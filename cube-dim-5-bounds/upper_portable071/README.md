# upper_portable071 portable extraction

This package carries the 15,143 passed parent records from all 17 frozen `upper_cover066` blocks. Each record contains the node, four integer root coordinates, stored box bounds, its matching cache point and five exact rational interval pairs (H followed by four physical-gradient intervals), cache precision bits, the original failed-certificate SHA-256 and old minimum, and the complete new certificate object. The certificate object is copied intact, including `all16_corners`.

The builder checks the frozen final-state and chain-manifest hashes, all 17 block certificate hashes and cumulative counts, and exact node/certificate correspondence between every accepted block row and `final/state.json`'s `passed` list. It verifies each cache point matches its anchor key, then copies stored values verbatim. It does not recompute H, gradients, logarithms, or margins. The 18,895 unpassed rows are excluded; a failed margin certificate is not presented as an H counterexample.

The final cache (178,722 entries) is read-only and is not copied into the packet. The optional logical-frontier/tree extraction was omitted; this packet contains accepted parent evidence only. The underlying H ball remains shared production input, so this is not independent value certification and proves no global `g` or `C5` bound.

## Reproduce

Use Python 3 and its standard library. All four input/output locations are supplied explicitly. The output directory must be new; existing directories are refused. The chain manifest is pinned by the expected frozen owner hash in the builder. Block data paths are resolved from `--chain-root` using the manifest's 17-block ordering.

```sh
python3 -B build_portable_bundle.py \
  --state /path/to/upper_cover066/final/state.json \
  --chain-manifest /path/to/upper_cover066/manifest.json \
  --chain-root /path/to/upper_cover066 \
  --output-dir /path/to/new/output-directory
```

The total program wall limit is 180 seconds, including input hashes, loading, compression and post-extraction hash checks. A run registers its PID before loading data. On failure it preserves that registration and writes `failure.json` with the failing stage and reason.

## Frozen run and verification

- Registered and run once as PID `64567`; elapsed time `3.80305` seconds.
- Source counts: `34,038` tested, `15,143` passed and included, `18,895` not passed and excluded.
- Bundle: `12,971,805` bytes; SHA-256 `54fcdd5af1745b4b4e8e8c97791f5bd27f102cfe9876ece22c6fb455f871ece7`; gzip mtime header `00000000`.
- Frozen final-state SHA-256 remained `42c67156078cb0eb134910a5248821eb3c826c96981afb3d1d5b2fe3242a871b`.
- Frozen owner/chain-manifest SHA-256 remained `75d5754fde6813285fe59d36b68680fd31675c375ae46cfcbb261f2c6649e8a6`.
- The builder checked and recorded every block certificate SHA-256. Read-back confirmed 15,143 packet records and the expected field dimensions.

## File whitelist

- `build_portable_bundle.py` — standard-library extractor; SHA-256 `9a68cbb1759fb17ddd1086bbbe9ed58d1b93bd1eb5cb55af9ed00e55781041fc`.
- `README.md` — scope, reproduction, run evidence and whitelist.
- `replay_001/registration.json` — pre-extraction registration and completed PID record.
- `replay_001/records.json.gz` — compact 15,143-record packet.
- `replay_001/manifest.json` — source/block hashes, builder hash, packet hash/size and scope.

The full state, cache, journal, block rows that did not pass, and unused records are not copied into this directory.
