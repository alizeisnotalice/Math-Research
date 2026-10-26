# upper_portable084 portable closed-root evidence

This package makes the `upper_cover070` closed-root certificate collection portable without copying its 680,238,201-byte final state or the rest of its cache. It contains all 52,360 integer roots, all 45,600 necessary internal split records (`node`, `axis`, `children`), and all 97,960 accepted leaf records. Each leaf keeps its complete original certificate and only the referenced cache `point`, five exact interval pairs (H and four physical gradients), and precision bits. Every certificate mode is preserved as recorded.

The package is split by contiguous root-index ranges into 26 gzip JSON chunks. Each root's tree is wholly contained in one chunk. Chunks use mtime zero and are SHA-256 indexed by the package manifest. No journal, unused cache values, attempts/geometry diagnostics, or failed split-node certificates are copied.

## Source locks and extraction

- Frozen final state: `8727af764ea088b1c30522f0204d6630ad6e0a4d804c7698c5ed2d646cedf343`.
- Frozen 070 manifest: `5f8caa2fe34775f1dfec518f484870828054966d10fc13575b28a211a8c5fdd2`.
- Owner report: `073809314163bd23af026b21e32ae866f1496d3047eb5ca8e79e603ff96346d2`.
- Builder PID `72041`; elapsed `12.587485` seconds. It rechecked all three source hashes after extraction.
- Compressed chunks total `70,601,199` bytes; largest chunk is `4,591,929` bytes. All 26 chunks are below 25 MB; total packet is below 2 GB.
- Certificate counts: 77,453 `Q=I4/4-J4/36`, 16,341 `signed_basis_81`, 4,166 `full_frame_16`.

Rebuild into a new output directory with Python 3 and the standard library:

```sh
python3 -B build_portable_tree.py \
  --state /path/to/upper_cover070/final/state.json \
  --manifest /path/to/upper_cover070/manifest.json \
  --source-report /path/to/upper_cover070/report.md \
  --output-dir /path/to/new/package-directory
```

The builder refuses an existing output directory and writes a PID/source-lock registration before loading the state. It performs only extraction and structural field checks; it makes no H, gradient, kernel, Green-function, Arb, or special-function calls.

## Independent check

Run the strengthened standard-library checker with the packaged directory and a new summary path:

```sh
python3 -B check_portable_tree_v2.py \
  --bundle-dir /path/to/package-directory \
  --output /path/to/new/check-summary.json
```

The checker requires no source-state path or private file. It verifies chunk hashes and zero gzip mtimes, requires ranges to partition exactly `[0,52360)` without overlap or gaps, and compares every `(rid, root)` pair against an independently generated nested catalog of all `31>=i>=j>=k>=l>=0` tuples. It builds actual global sets for leaf nodes, all tree nodes and leaf anchors, then checks exact unique-set cardinalities. It reconstructs each physical box from the root and four-bit interleaved path, verifies split parents/children and exact root coverage, checks every anchor/cache key and stored H/physical-gradient interval, and validates the directed gradient/H bounds. It independently checks every stored `Delta_lower` with the 896-bit, 224-term directed fixed-point `atanh` series ported from round010 checker SHA-256 `2ad20352bf6a013d9e861f39c18e72903893dd4ee986cb8c740fab38536f9341`.

Each leaf is checked under its actual saved curvature mode: the 77,453 single-basis leaves use their 16-corner `Q=I4/4-J4/36` certificate; the 16,341 signed-basis leaves use their saved 81-point psi certificate; and the 4,166 full-frame leaves use their 16-corner `Q_full=I4/16` certificate. It compares exact rational minima and worst corners/offsets with the stored certificate. It does not substitute one mode for another.

### Strengthened v2 run

- Checker SHA-256: `0e7e7db53c0e12db1dbeb9f5d8da866f50fc0e34ad40ea3733cb445b6ae21f58`.
- Registered PID `72868`; elapsed `97.528233` seconds.
- Before the full run, `--self-test-guards` rejected repeated/skipped/out-of-range root ranges, a valid tuple assigned to the wrong ID, and malformed/out-of-order root tuples.
- Result: `PASS_portable_closed_root_certificate_package`.
- Checked 52,360 closed roots, 97,960 leaves and 45,600 internal splits; all 97,960 leaf log lower bounds and local certificates passed exact checks.
- Explicit set checks found 97,960 distinct leaf nodes, 97,960 distinct anchors, and 143,560 distinct total tree nodes. Root ranges were contiguous and exact; all 52,360 roots matched the canonical tuple catalog by ID.
- The certificate-mode counts match the source inventory. The minimum saved margin by mode and minimum Delta are recorded exactly in `package_001/independent_check_v2_replay004.json`.

The prior `check_portable_tree.py` / `independent_check_replay003.json` record is retained unchanged as local history. That pass checked local chunk trees and certificates, but did not check catalog-by-ID correspondence, contiguous global root ranges, or actual global anchor uniqueness. The complete portable verification command and archive whitelist use `check_portable_tree_v2.py`; the v1 checker, v1 report and failed attempt sidecars are not needed in a concise archive.

Two initial checker attempts stopped on checker bookkeeping errors: one missing size constant before reading a chunk, and one zero-count certificate-mode comparison after checking the first chunk. Their registration/prefix sidecars are retained as diagnostics. The corrected run above checked all 26 chunks.

The H balls in the original state are shared historical production inputs. This portable package and its checker establish only a closed-root local certificate collection; they do not assert a global `g` or `Cn` bound.

## Essential file whitelist

- `build_portable_tree.py` — portable extractor.
- `check_portable_tree_v2.py` — strengthened portable independent checker.
- `README.md` — sources, method, commands, scope, run evidence and whitelist.
- `archive_whitelist.json` — exact hashes and sizes for every required archive file.
- `guard_selftest.json` — saved result of the finite root/range negative guard checks.
- `package_001/registration.json` and `package_001/manifest.json` — source locks and chunk index.
- `package_001/chunk01.json.gz` through `chunk26.json.gz` — complete portable root/tree/leaf evidence.
- `package_001/independent_check_v2_replay004.json` and its `.registration.json` / `.prefix.json` sidecars — strengthened successful full check evidence.

The v1 checker, `independent_check_replay003.json` and its sidecars, and the two failed v1 attempt sidecar pairs may remain in the local directory as historical diagnostics; they are excluded from the archive whitelist.
