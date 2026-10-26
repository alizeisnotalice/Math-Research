# upper_portable096

## Completed portable extraction

The registered 098 extraction and independent checker completed successfully. The packet is in [`run098_v3/`](run098_v3/); its receipt records full input, code, manifest, chunk, and checker hashes. It contains 52 gzip chunks totaling 13,764,568 bytes.

| Check | Result |
|---|---:|
| Roots | 52,360 |
| Accepted leaves | 100,998 |
| Splits | 48,638 |
| Epsilon `1/64` / `1/128` | 100,804 / 194 |
| Unique `(epsilon, point)` keys | 100,998 |
| Selected outward quantization | 160 bits for all leaves |
| Independent checker | PASS |

Source locks: final state `73bc873f23f7b68f6f70132e590d28ebbe4cd619012ac9c81980f9fa500314af`; production manifest `b318e998a41d3b8d78ab663e8a457244e1b87b03c380b05f7c27664fdd0f557d`. Both hashes were unchanged after extraction.

Use `build_compact_packet_v3.py`, `compact_math_v3.py`, and `check_compact_packet_v3.py` for this source. The frozen v1/v2 drafts remain unchanged. The v3 adapter follows the actual 098 schema: `passed[].leaf_epsilon`, certificate `anchor` or `point_key`, and `cache_by_epsilon[epsilon][point]`. It retains the actual source mechanism, quadratic label, and penalty coefficient for per-leaf checking. This v3 adapter was needed because v2 expected a global cache and earlier certificate field names.

The checker is self-contained, reads only the packet, rebuilds the complete canonical root tree, checks unique leaf nodes and epsilon-point keys, verifies the source mode and coefficient for each epsilon, and recomputes directed log, gradient/H, and local margins with exact rational arithmetic and the 896-bit/224-term fixed-point log routine. It made zero H/gradient/kernel/Green/Arb calls. The source H balls remain shared production inputs.

## Reproduce the packet check

From the repository root, run the checker into a new output path:

```sh
python3 cube-dim-5-bounds/upper_portable096/check_compact_packet_v3.py \
  --bundle-dir cube-dim-5-bounds/upper_portable096/run098_v3 \
  --output /path/to/new/checker_summary.json
```

The extraction receipt is [`run098_v3/receipt.json`](run098_v3/receipt.json). Build registration and checker registration are also retained there.

This is portable local closed-root evidence only. It does not independently recompute H, establish a global `g` or `C5` bound, or replace the pending mathematical qualification and parent review.
