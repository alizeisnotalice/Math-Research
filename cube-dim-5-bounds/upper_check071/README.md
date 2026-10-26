# upper_check071 independent portable checker

This checker verifies the 15,143 local parent certificates in the `upper_portable071` packet using only Python's standard library and exact rational arithmetic. It reconstructs each box from its four integer root coordinates and interleaved path, checks saved box and anchor/cache correspondence, checks the five stored H/physical-gradient interval pairs and directed gradient/H bounds, and independently checks the saved `Delta_lower` against a fixed-point logarithm lower endpoint.

It then recomputes the old 16-corner signed-sum quadratic minimum and the new 81-point envelope
`psi(v) = sum(v_i^2)/4 - (sum(abs(v_i)))^2/36` exactly. The old minimum must equal the saved value and be nonpositive; the new minimum and worst offset must exactly match the certificate, with a positive minimum. All 15,143 nodes must be unique.

The fixed-point `atanh`/log routine is a standard-library port of the reviewed round010 implementation at SHA-256 `2ad20352bf6a013d9e861f39c18e72903893dd4ee986cb8c740fab38536f9341`. This checker uses 896-bit fixed point and 224 series terms. It makes no H, gradient, kernel, Arb, or special-function calls. It checks only the stored intervals and local certificate collection. The underlying H ball remains shared production input; the result is not a full-tree, global `g`, or `C5` certification.

## Reproduce

Use Python 3 and the standard library. `--bundle` must point to the frozen packet with SHA-256 `54fcdd5af1745b4b4e8e8c97791f5bd27f102cfe9876ece22c6fb455f871ece7`. `--output` must be a new JSON path; registration and prefix sidecars use the same name stem, and existing paths are refused.

```sh
python3 -B check_portable_certificates.py \
  --bundle /path/to/upper_portable071/replay_001/records.json.gz \
  --output /path/to/new/summary.json
```

The run registers its PID, checker hash, bundle hash and limits before loading the packet. It writes a durable prefix checkpoint every 64 rows and stops at 540 seconds. A 590-second SIGALRM is a hard fallback; total runtime is capped at 600 seconds. The complete run writes only the summary and small registration/prefix files; it does not write a duplicate per-row result file.

## Recorded run

- Checker SHA-256: `daceb23caffd6797f06a966e8838f7a1c11ef4d719f89f5865eb32c6e94d8075`.
- Registered PID: `66047`; complete run time: `32.249589` seconds.
- Records and unique nodes: `15,143`; saved Delta/log lower checks: `15,143`.
- Old 16-corner minimum matched and was `<= 0` for all `15,143` records.
- New 81-point minimum matched and was `> 0` for all `15,143` records.
- Bundle hash was checked against the frozen packet hash above.

See `summary.json` for exact worst minima and `summary.json.registration.json` and `summary.json.prefix.json` for the run registration and completed checkpoint.

## File whitelist

- `check_portable_certificates.py` — portable independent checker.
- `README.md` — method, reproduction command, scope, run evidence and whitelist.
- `summary.json` — exact aggregate result and limits.
- `summary.json.registration.json` — pre-load registration and PID.
- `summary.json.prefix.json` — final durable prefix checkpoint.
