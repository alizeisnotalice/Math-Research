#!/usr/bin/env python3
"""Prepare a compact outward-quantized closed-root certificate packet.

Build the locked 098 mixed-epsilon packet from the exact frozen source.
"""
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
from fractions import Fraction as F

from compact_math_v3 import G, SUPPORTED_EPSILONS, root_catalog, quantize_values, assess_compact, mode_from_certificate


PRECISIONS = (160, 192, 256)
ROOTS_PER_CHUNK = 1024
MAX_CHUNK_BYTES = 25_000_000


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def write_json(path: Path, obj: object) -> None:
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", encoding="utf-8") as f:
        f.write(json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True) + "\n")
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def epsilon_for(row: dict, cert: dict, state: dict) -> F:
    """Resolve source leaf epsilon and require certificate agreement when present."""
    raw = row.get("leaf_epsilon")
    if raw is None:
        raise ValueError("source accepted leaf is missing leaf_epsilon")
    cert_eps = cert.get("leaf_epsilon", cert.get("epsilon"))
    if cert_eps is not None and str(F(str(cert_eps))) != str(F(str(raw))):
        raise ValueError("leaf and certificate epsilon fields disagree")
    epsilon = F(str(raw))
    if epsilon not in SUPPORTED_EPSILONS:
        raise ValueError(f"unsupported source epsilon {raw!r}")
    return epsilon


def cache_for(state: dict, epsilon: F, anchor: str) -> dict:
    """Resolve only the exact epsilon-point cache namespace; never cross-reuse."""
    cache = state["cache_by_epsilon"]
    eps = str(epsilon)
    subcache = cache.get(eps)
    if not isinstance(subcache, dict) or anchor not in subcache:
        raise KeyError(f"missing source cache record for ({eps}, {anchor})")
    record = subcache[anchor]
    tagged = record.get("epsilon", eps)
    if str(F(str(tagged))) != eps:
        raise ValueError("cache record epsilon differs from bucket epsilon")
    if record.get("point") is None or ",".join(record["point"]) != anchor:
        raise ValueError("cache point does not match canonical epsilon-point key")
    return record


def point_key_for(cert: dict) -> str:
    key = cert.get("anchor", cert.get("point_key"))
    if not isinstance(key, str):
        raise ValueError("source certificate has no canonical point key")
    return key


def source_mode_fields(cert: dict, epsilon: F) -> dict:
    mode_from_certificate(cert)
    mechanism, quadratic = cert["mechanism"], cert["quadratic"]
    penalty = cert.get("penalty_coefficient")
    expected_penalty = str(1 / (16 * epsilon * epsilon))
    if penalty is not None and str(F(str(penalty))) != expected_penalty:
        raise ValueError("saved penalty coefficient does not match epsilon")
    if quadratic == "Q_full=I4/16" and epsilon != F(1, 64):
        raise ValueError("Q_full=I4/16 source mode is only expected at epsilon=1/64")
    if quadratic == "I4/(16epsilon^2)" and epsilon != F(1, 128):
        raise ValueError("epsilon-scaled full-frame source mode is only expected at epsilon=1/128")
    return {"curvature_mode": "full_frame_16", "source_mechanism": mechanism,
            "source_quadratic": quadratic, "source_penalty_coefficient": expected_penalty}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--state", required=True, type=Path)
    ap.add_argument("--source-manifest", required=True, type=Path)
    ap.add_argument("--expected-state-sha256", required=True)
    ap.add_argument("--expected-manifest-sha256", required=True)
    ap.add_argument("--output-dir", required=True, type=Path)
    args = ap.parse_args()

    began = time.monotonic()
    state_path = args.state.expanduser().resolve(strict=True)
    manifest_path = args.source_manifest.expanduser().resolve(strict=True)
    output = args.output_dir.expanduser().resolve()
    if output.exists():
        raise SystemExit(f"refusing to overwrite output directory: {output}")
    output.mkdir(parents=True)
    source_locks = {
        "state": args.expected_state_sha256,
        "source_manifest": args.expected_manifest_sha256,
    }
    registration = {
        "task": "upper_portable096",
        "status": "registered_before_source_load",
        "registered_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "pid": os.getpid(),
        "source_paths": {"state": str(state_path), "manifest": str(manifest_path)},
        "source_locks": source_locks,
        "finite_scope": "all 100998 accepted leaves and closed tree in frozen upper_mixed098 final_state",
        "root_count_expected": 52360,
        "g": str(G),
        "epsilon_counts_expected": {"1/64": 100804, "1/128": 194},
        "epsilon_policy": {"supported": [str(e) for e in SUPPORTED_EPSILONS], "source_field": "passed[].leaf_epsilon"},
        "cache_identity": "state.cache_by_epsilon[epsilon][canonical_point_key]; cross-epsilon cache reuse forbidden",
        "builder_sha256": sha256_file(Path(__file__).resolve()),
        "checker_sha256": sha256_file(Path(__file__).with_name("check_compact_packet_v3.py")),
        "math_helper_sha256": sha256_file(Path(__file__).with_name("compact_math_v3.py")),
        "quantization_precision_order": list(PRECISIONS),
        "chunk_roots": ROOTS_PER_CHUNK,
        "chunk_max_bytes": MAX_CHUNK_BYTES,
        "algorithm": "exact Fraction outward quantization of stored rational interval endpoints; no H/gradient/kernel production calls",
    }
    write_json(output / "registration.json", registration)
    stage = "source_hashes"
    try:
        actual = {
            "state": sha256_file(state_path),
            "source_manifest": sha256_file(manifest_path),
        }
        if actual != source_locks:
            raise ValueError(f"frozen input hash mismatch: {actual!r}")
        source_manifest = json.load(manifest_path.open(encoding="utf-8"))
        locked_state = source_manifest.get("external_files", {}).get(str(state_path))
        if locked_state != source_locks["state"]:
            raise ValueError("source manifest does not lock the supplied final state hash")
        stage = "load_state"
        state = json.load(state_path.open(encoding="utf-8"))
        catalog = root_catalog()
        if len(catalog) != 52360 or state.get("roots") != [list(r) for r in catalog]:
            raise ValueError("source roots differ from the canonical 52,360-tuple catalog")
        for key in ("queue", "inflight", "resume_jobs"):
            if state.get(key):
                raise ValueError(f"source is not closed: {key} is nonempty")
        if not isinstance(state.get("passed"), list) or not isinstance(state.get("splits"), list) or not isinstance(state.get("cache_by_epsilon"), dict):
            raise ValueError("source state does not have mixed098 passed/splits/cache_by_epsilon schema")
        if state.get("target") != str(G) or state.get("cache_key_format") != "(canonical epsilon,canonical point)":
            raise ValueError("source target or epsilon-point cache key schema mismatch")

        stage = "project_and_quantize"
        chunk_count = (len(catalog) + ROOTS_PER_CHUNK - 1) // ROOTS_PER_CHUNK
        leaves_by_chunk = [[] for _ in range(chunk_count)]
        splits_by_chunk = [[] for _ in range(chunk_count)]
        seen_nodes = set()
        seen_cache_keys = set()
        epsilon_counts: dict[str, int] = {}
        mode_counts: dict[str, int] = {}
        precision_counts = {b: 0 for b in PRECISIONS}
        upgrade_attempts = 0
        unresolved = []
        for source_row in state["passed"]:
            node = source_row["node"]
            rid, path = node
            if type(rid) is not int or not 0 <= rid < len(catalog):
                raise ValueError(f"invalid source leaf node: {node!r}")
            key = (rid, path)
            if key in seen_nodes:
                raise ValueError(f"duplicate source leaf node: {node!r}")
            seen_nodes.add(key)
            cert = source_row["certificate"]
            if cert.get("passed") is not True:
                raise ValueError(f"nonaccepted row found in source passed list: {node!r}")
            epsilon = epsilon_for(source_row, cert, state)
            epsilon_text = str(epsilon)
            mode_fields = source_mode_fields(cert, epsilon)
            mode = mode_fields["curvature_mode"]
            mode_counts[mode] = mode_counts.get(mode, 0) + 1
            epsilon_counts[epsilon_text] = epsilon_counts.get(epsilon_text, 0) + 1
            anchor = point_key_for(cert)
            cache_key = (epsilon_text, anchor)
            if cache_key in seen_cache_keys:
                raise ValueError(f"duplicate source epsilon-point pair: {cache_key!r}")
            seen_cache_keys.add(cache_key)
            cache = cache_for(state, epsilon, anchor)
            point = cache["point"]
            if len(point) != 5 or point[0] != "1" or ",".join(point) != anchor:
                raise ValueError(f"source anchor/cache-key mismatch at {node!r}")
            cache_epsilon = cache.get("epsilon", cache.get("eps", epsilon_text))
            if str(F(str(cache_epsilon))) != epsilon_text:
                raise ValueError(f"cache epsilon mismatch at {node!r}")

            attempts = []
            selected = None
            for bits in PRECISIONS:
                pairs = quantize_values(cache, bits)
                candidate = assess_compact(node, list(catalog[rid]), point[1:], epsilon, bits, pairs, mode)
                if candidate["passed"]:
                    selected = {"node": node, "point_tail": point[1:], "epsilon": epsilon_text,
                                "cache_epsilon": epsilon_text, "quant_bits": bits,
                                "value_pairs": pairs, **mode_fields}
                    if attempts:
                        selected["precision_attempts"] = attempts
                    precision_counts[bits] += 1
                    upgrade_attempts += bool(attempts)
                    break
                attempts.append({"quant_bits": bits, "value_pairs": pairs, "reason": candidate["reason"] if "reason" in candidate else str(candidate["minimum"])})
            if selected is None:
                unresolved.append({"node": node, "attempted_bits": list(PRECISIONS), "anchor": anchor})
                continue
            leaves_by_chunk[rid // ROOTS_PER_CHUNK].append(selected)

        if unresolved:
            write_json(output / "unresolved_precision.json", {
                "status": "not_packaged_unseparated_leaf_exists",
                "count": len(unresolved),
                "leaves": unresolved,
                "precision_counts": precision_counts,
                "scope": "no leaf dropped; inspect source proof margin or permitted precision escalation",
            })
            raise ValueError(f"{len(unresolved)} leaves remain uncertified after 160/192/256-bit outward quantization")

        for split in state["splits"]:
            node = split["node"]
            rid, path = node
            if type(rid) is not int or not 0 <= rid < len(catalog):
                raise ValueError(f"invalid split node {node!r}")
            splits_by_chunk[rid // ROOTS_PER_CHUNK].append({"node": node, "axis": split["axis"], "children": split["children"]})

        if len(seen_nodes) != len(state["passed"]) or len(seen_cache_keys) != len(state["passed"]):
            raise ValueError("source leaf/epsilon-point unique-set cardinality mismatch")
        if len(state["passed"]) != 100998 or epsilon_counts != {"1/64": 100804, "1/128": 194}:
            raise ValueError(f"frozen source accepted-leaf/epsilon counts mismatch: {len(state['passed'])}, {epsilon_counts}")
        registration["source_counts"] = {
            "roots": len(catalog), "leaves": len(state["passed"]), "splits": len(state["splits"]),
            "source_cache_entries": sum(len(v) for v in state["cache_by_epsilon"].values()), "unique_epsilon_point_pairs": len(seen_cache_keys),
            "epsilon_counts": epsilon_counts,
        }
        registration["selected_precision_counts"] = {str(k): v for k, v in precision_counts.items()}

        stage = "write_compact_chunks"
        chunk_manifest = []
        total_bytes = 0
        for chunk_index in range(chunk_count):
            start = chunk_index * ROOTS_PER_CHUNK
            stop = min(start + ROOTS_PER_CHUNK, len(catalog))
            leaves = leaves_by_chunk[chunk_index]
            splits = splits_by_chunk[chunk_index]
            curvatures: dict[str, int] = {}
            quant_bits_hist: dict[str, int] = {}
            epsilon_hist: dict[str, int] = {}
            for leaf in leaves:
                curvatures[leaf["curvature_mode"]] = curvatures.get(leaf["curvature_mode"], 0) + 1
                q = str(leaf["quant_bits"])
                quant_bits_hist[q] = quant_bits_hist.get(q, 0) + 1
                e = leaf["epsilon"]
                epsilon_hist[e] = epsilon_hist.get(e, 0) + 1
            chunk = {
                "format": "upper_portable096_mixed_epsilon_root_chunk_v3",
                "root_range": [start, stop],
                "leaves": leaves,
                "splits": splits,
            }
            name = f"chunk{chunk_index + 1:03d}.json.gz"
            path = output / name
            with path.open("wb") as raw:
                with gzip.GzipFile(filename="", mode="wb", compresslevel=9, fileobj=raw, mtime=0) as zipped:
                    with io.TextIOWrapper(zipped, encoding="utf-8", write_through=True) as text:
                        json.dump(chunk, text, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
            size = path.stat().st_size
            if size > MAX_CHUNK_BYTES:
                raise ValueError(f"chunk {name} exceeds 25MB: {size} bytes")
            total_bytes += size
            chunk_manifest.append({
                "file": name, "sha256": sha256_file(path), "bytes": size,
                "root_range": [start, stop], "leaves": len(leaves), "splits": len(splits),
                "curvature_modes": curvatures, "quant_bits": quant_bits_hist, "epsilon_counts": epsilon_hist,
            })
            del chunk, leaves, splits

        stage = "source_hash_recheck"
        if sha256_file(state_path) != source_locks["state"] or sha256_file(manifest_path) != source_locks["source_manifest"]:
            raise ValueError("source hash changed during package extraction")
        package_manifest = {
            "task": "upper_portable096",
            "status": "portable_compact_certificate_packet",
            "source_locks": source_locks,
            "source_counts": registration["source_counts"],
            "curvature_modes": mode_counts,
            "quantization_bits": {str(k): v for k, v in precision_counts.items()},
            "adaptive_upgrade_leaf_count": upgrade_attempts,
            "g": str(G), "epsilon_counts": epsilon_counts,
            "supported_epsilons": [str(e) for e in SUPPORTED_EPSILONS],
            "interval_encoding": "each exact source endpoint pair [a,b] becomes [floor(a*2^q), ceil(b*2^q)-floor(a*2^q)] as canonical decimal integers; width is nonnegative",
            "precision_attempt_witness": "for each q>160 leaf, the packet retains every lower-precision outward pair attempt so checker can confirm it did not separate",
            "log_and_gradient": "not stored; recomputed from outward H/physical-gradient intervals by checker",
            "omitted": "source cache, Delta, directed gradient ratios, old margins, journals, and production attempts",
            "partition": {"roots_per_chunk": ROOTS_PER_CHUNK, "chunks": len(chunk_manifest), "max_chunk_bytes": MAX_CHUNK_BYTES, "actual_total_chunk_bytes": total_bytes, "under_25MB_total": total_bytes < 25_000_000},
            "chunks": chunk_manifest,
            "scope": "local accepted-leaf closed-root evidence; shared source H balls; no global g or C5 claim",
            "builder_sha256": sha256_file(Path(__file__).resolve()),
        }
        write_json(output / "manifest.json", package_manifest)
        registration["status"] = "completed"
        registration["elapsed_seconds"] = round(time.monotonic() - began, 6)
        registration["package_manifest_sha256"] = sha256_file(output / "manifest.json")
        write_json(output / "registration.json", registration)
        print(json.dumps({"status": package_manifest["status"], "pid": os.getpid(), "elapsed_seconds": registration["elapsed_seconds"], "source_counts": registration["source_counts"], "quantization_bits": package_manifest["quantization_bits"], "upgrade_leaf_count": upgrade_attempts, "total_chunk_bytes": total_bytes, "under_25MB_total": total_bytes < 25_000_000}, indent=2))
        return 0
    except BaseException as exc:
        registration["status"] = "failed"
        registration["failure_stage"] = stage
        registration["failure_reason"] = f"{type(exc).__name__}: {exc}"
        registration["elapsed_seconds"] = round(time.monotonic() - began, 6)
        write_json(output / "registration.json", registration)
        write_json(output / "failure.json", {"stage": stage, "reason": f"{type(exc).__name__}: {exc}", "pid": os.getpid()})
        raise


if __name__ == "__main__":
    sys.exit(main())
