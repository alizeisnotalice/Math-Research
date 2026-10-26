#!/usr/bin/env python3
"""Build a compact closed-root certificate packet from frozen upper_cover070."""
from __future__ import annotations

import argparse
import datetime as dt
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import sys
import time


LOCKS = {
    "state": "8727af764ea088b1c30522f0204d6630ad6e0a4d804c7698c5ed2d646cedf343",
    "manifest": "5f8caa2fe34775f1dfec518f484870828054966d10fc13575b28a211a8c5fdd2",
    "owner_report": "073809314163bd23af026b21e32ae866f1496d3047eb5ca8e79e603ff96346d2",
}
ROOTS = 52360
LEAVES = 97960
SPLITS = 45600
ANCHORS = 97960
ROOTS_PER_CHUNK = 2048
MAX_CHUNK_BYTES = 25_000_000
MAX_TOTAL_BYTES = 2_000_000_000


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def write_json(path: Path, obj: object) -> None:
    payload = json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True) + "\n"
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", encoding="utf-8") as f:
        f.write(payload)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--state", required=True, type=Path)
    ap.add_argument("--manifest", required=True, type=Path)
    ap.add_argument("--source-report", required=True, type=Path)
    ap.add_argument("--output-dir", required=True, type=Path)
    args = ap.parse_args()
    started = time.monotonic()
    state_path = args.state.expanduser().resolve(strict=True)
    manifest_path = args.manifest.expanduser().resolve(strict=True)
    report_path = args.source_report.expanduser().resolve(strict=True)
    out = args.output_dir.expanduser().resolve()
    if out.exists():
        raise SystemExit(f"refusing to overwrite existing output directory: {out}")
    out.mkdir(parents=True)
    registration = {
        "task": "upper_portable084",
        "status": "registered_before_state_load",
        "registered_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "pid": os.getpid(),
        "state_path": str(state_path),
        "manifest_path": str(manifest_path),
        "source_report_path": str(report_path),
        "output_dir": str(out),
        "source_locks": LOCKS,
        "expected_counts": {"roots": ROOTS, "leaves": LEAVES, "splits": SPLITS, "unique_anchors": ANCHORS},
        "partition": {"roots_per_chunk": ROOTS_PER_CHUNK, "maximum_chunk_bytes": MAX_CHUNK_BYTES},
        "scope": "portable evidence extraction only; no H/gradient/kernel/special-function recomputation",
    }
    write_json(out / "registration.json", registration)
    stage = "source_hashes"
    try:
        actual = {
            "state": sha256_file(state_path),
            "manifest": sha256_file(manifest_path),
            "owner_report": sha256_file(report_path),
        }
        if actual != LOCKS:
            raise ValueError(f"source lock mismatch: {actual!r}")
        with manifest_path.open(encoding="utf-8") as f:
            owner_manifest = json.load(f)
        if owner_manifest.get("task") != "upper_cover070" or owner_manifest.get("global_g_or_C5_claim") is not False:
            raise ValueError("unexpected task or scope in frozen owner manifest")
        if owner_manifest.get("external_files", {}).get(str(state_path)) != LOCKS["state"]:
            raise ValueError("owner manifest does not lock the requested final state")

        stage = "load_and_lock_state_schema"
        with state_path.open(encoding="utf-8") as f:
            state = json.load(f)
        if state.get("status070") != "closed_root_cover_candidate_waiting_main_review":
            raise ValueError("final state status does not match the frozen closed-root candidate")
        if state.get("epsilon") != "1/64" or state.get("target") != "1/10000":
            raise ValueError("epsilon/target mismatch")
        if state.get("queue") or state.get("inflight") or state.get("resume_jobs"):
            raise ValueError("frozen source is not closed")
        if [len(state[k]) for k in ("roots", "passed", "splits", "cache")] != [ROOTS, LEAVES, SPLITS, 180040]:
            raise ValueError("frozen root/leaf/split/cache counts do not match registered counts")

        stage = "partition_source_rows"
        chunk_count = (ROOTS + ROOTS_PER_CHUNK - 1) // ROOTS_PER_CHUNK
        leaves_by_chunk: list[list[dict]] = [[] for _ in range(chunk_count)]
        splits_by_chunk: list[list[dict]] = [[] for _ in range(chunk_count)]
        seen_nodes = set()
        seen_anchors = set()
        mode_counts: dict[str, int] = {}
        for row in state["passed"]:
            node = row["node"]
            rid = node[0]
            if type(rid) is not int or not 0 <= rid < ROOTS:
                raise ValueError(f"invalid leaf root index: {node!r}")
            key = (rid, node[1])
            if key in seen_nodes:
                raise ValueError(f"duplicate leaf node: {node!r}")
            seen_nodes.add(key)
            cert = row["certificate"]
            if cert.get("passed") is not True or not isinstance(cert.get("anchor"), str):
                raise ValueError(f"invalid passed leaf certificate at {node!r}")
            anchor = cert["anchor"]
            if anchor in seen_anchors:
                raise ValueError(f"duplicate accepted cache anchor: {anchor}")
            seen_anchors.add(anchor)
            cache = state["cache"][anchor]
            if ",".join(cache["point"]) != anchor or len(cache["values"]) != 5:
                raise ValueError(f"cache point/anchor or five-value schema mismatch at {node!r}")
            if any(len(pair) != 2 for pair in cache["values"]):
                raise ValueError(f"incomplete cache intervals at {node!r}")
            mode = cert.get("mechanism", cert.get("quadratic", "MISSING"))
            mode_counts[mode] = mode_counts.get(mode, 0) + 1
            leaves_by_chunk[rid // ROOTS_PER_CHUNK].append({
                "node": node,
                "certificate": cert,
                "cache_record": {"point": cache["point"], "values": cache["values"], "bits": cache["bits"]},
            })
        for row in state["splits"]:
            node = row["node"]
            rid = node[0]
            if type(rid) is not int or not 0 <= rid < ROOTS:
                raise ValueError(f"invalid split root index: {node!r}")
            splits_by_chunk[rid // ROOTS_PER_CHUNK].append({k: row[k] for k in ("node", "axis", "children")})
        if (len(seen_nodes), len(seen_anchors)) != (LEAVES, ANCHORS):
            raise ValueError("leaf or unique anchor count mismatch")
        if mode_counts != {"Q=I4/4-J4/36": 77453, "signed_basis_81": 16341, "full_frame_16": 4166}:
            raise ValueError(f"unexpected source certificate-mode counts: {mode_counts!r}")

        stage = "write_chunks"
        chunk_rows = []
        total_bytes = 0
        for chunk_index in range(chunk_count):
            root_start = chunk_index * ROOTS_PER_CHUNK
            root_stop = min(root_start + ROOTS_PER_CHUNK, ROOTS)
            roots = state["roots"][root_start:root_stop]
            leaves = leaves_by_chunk[chunk_index]
            splits = splits_by_chunk[chunk_index]
            chunk_modes: dict[str, int] = {}
            for leaf in leaves:
                cert = leaf["certificate"]
                mode = cert.get("mechanism", cert.get("quadratic", "MISSING"))
                chunk_modes[mode] = chunk_modes.get(mode, 0) + 1
            packet = {
                "format": "upper_portable084_closed_root_chunk_v1",
                "task": "upper_cover070",
                "epsilon": state["epsilon"],
                "target": state["target"],
                "root_range": [root_start, root_stop],
                "roots": roots,
                "splits": splits,
                "leaves": leaves,
            }
            name = f"chunk{chunk_index + 1:02d}.json.gz"
            chunk_path = out / name
            with chunk_path.open("wb") as raw:
                with gzip.GzipFile(filename="", mode="wb", compresslevel=9, fileobj=raw, mtime=0) as zipped:
                    with io.TextIOWrapper(zipped, encoding="utf-8", write_through=True) as text:
                        json.dump(packet, text, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
            size = chunk_path.stat().st_size
            digest = sha256_file(chunk_path)
            if size > MAX_CHUNK_BYTES:
                raise ValueError(f"{name} exceeds 25,000,000-byte limit: {size}")
            total_bytes += size
            chunk_rows.append({
                "file": name,
                "sha256": digest,
                "bytes": size,
                "root_range": [root_start, root_stop],
                "roots": len(roots),
                "leaves": len(leaves),
                "splits": len(splits),
                "certificate_modes": chunk_modes,
            })
            del packet, roots, leaves, splits
        if total_bytes > MAX_TOTAL_BYTES:
            raise ValueError(f"total gzip payload exceeds 2GB bound: {total_bytes}")

        stage = "verify_sources_unchanged"
        if sha256_file(state_path) != LOCKS["state"] or sha256_file(manifest_path) != LOCKS["manifest"] or sha256_file(report_path) != LOCKS["owner_report"]:
            raise ValueError("a frozen source hash changed during extraction")

        stage = "write_package_manifest"
        manifest = {
            "task": "upper_portable084",
            "status": "extracted_closed_root_certificate_package",
            "source_locks": LOCKS,
            "source_paths": {"state": str(state_path), "manifest": str(manifest_path), "owner_report": str(report_path)},
            "source_state_bytes": state_path.stat().st_size,
            "registered_pid": os.getpid(),
            "elapsed_seconds": round(time.monotonic() - started, 6),
            "source_counts": {"roots": ROOTS, "leaves": LEAVES, "splits": SPLITS, "source_cache_entries": 180040, "unique_leaf_anchors": ANCHORS},
            "certificate_modes": mode_counts,
            "partition": {"roots_per_chunk": ROOTS_PER_CHUNK, "chunk_count": chunk_count, "max_chunk_bytes": MAX_CHUNK_BYTES, "total_chunk_bytes": total_bytes},
            "chunks": chunk_rows,
            "fields": "root coordinates; split node/axis/children; each complete original leaf certificate; referenced cache point, five exact H/physical-gradient interval pairs, and bits",
            "omitted": "full state, unused cache entries, journal, attempts/geometry diagnostics, split-node failed certificates",
            "scope": "portable local closed-root certificate evidence; shared historical H balls; no global g or Cn claim",
            "builder_sha256": sha256_file(Path(__file__).resolve()),
        }
        write_json(out / "manifest.json", manifest)
        registration["status"] = "completed"
        registration["elapsed_seconds"] = round(time.monotonic() - started, 6)
        registration["package_manifest_sha256"] = sha256_file(out / "manifest.json")
        write_json(out / "registration.json", registration)
        print(json.dumps({"status": manifest["status"], "pid": os.getpid(), "elapsed_seconds": manifest["elapsed_seconds"], "chunks": chunk_count, "total_gzip_bytes": total_bytes, "certificate_modes": mode_counts}, indent=2))
        return 0
    except BaseException as exc:
        registration["status"] = "failed"
        registration["failure_stage"] = stage
        registration["failure_reason"] = f"{type(exc).__name__}: {exc}"
        registration["elapsed_seconds"] = round(time.monotonic() - started, 6)
        write_json(out / "registration.json", registration)
        write_json(out / "failure.json", {"stage": stage, "reason": f"{type(exc).__name__}: {exc}", "pid": os.getpid()})
        raise


if __name__ == "__main__":
    sys.exit(main())
