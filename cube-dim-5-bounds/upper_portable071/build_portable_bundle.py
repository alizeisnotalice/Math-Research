#!/usr/bin/env python3
"""Package passed upper_cover066 records using only the Python standard library.

The builder verifies frozen byte hashes and exact stored-object correspondence.
It never recomputes H, gradients, logarithms, or certificate margins.
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


EXPECTED_STATE_SHA256 = "42c67156078cb0eb134910a5248821eb3c826c96981afb3d1d5b2fe3242a871b"
EXPECTED_CHAIN_MANIFEST_SHA256 = "75d5754fde6813285fe59d36b68680fd31675c375ae46cfcbb261f2c6649e8a6"
EXPECTED_TESTED = 34038
EXPECTED_ACCEPTED = 15143
EXPECTED_UNPASSED = 18895
WALL_LIMIT_SECONDS = 180


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, obj: object) -> None:
    path.write_text(json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True) + "\n", encoding="utf-8")


def assert_time(started: float, where: str) -> None:
    elapsed = time.monotonic() - started
    if elapsed > WALL_LIMIT_SECONDS:
        raise TimeoutError(f"180-second total wall limit exceeded at {where}: {elapsed:.3f}s")


def write_failure(out: Path, stage: str, exc: BaseException, started: float) -> None:
    write_json(out / "failure.json", {
        "status": "failed", "stage": stage, "exception_type": type(exc).__name__,
        "reason": str(exc), "pid": os.getpid(),
        "elapsed_seconds": round(time.monotonic() - started, 6),
    })


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--state", required=True, type=Path)
    parser.add_argument("--chain-manifest", required=True, type=Path)
    parser.add_argument("--chain-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    state_path = args.state.expanduser().resolve(strict=True)
    manifest_path = args.chain_manifest.expanduser().resolve(strict=True)
    chain_root = args.chain_root.expanduser().resolve(strict=True)
    output_dir = args.output_dir.expanduser().resolve()
    if output_dir.exists():
        raise SystemExit(f"refusing to overwrite existing output directory: {output_dir}")
    if state_path != (chain_root / "final" / "state.json").resolve():
        raise SystemExit("--state must identify final/state.json within --chain-root")

    output_dir.mkdir(parents=True)
    started = time.monotonic()
    registration = {
        "task": "upper_portable071", "status": "registered_before_extraction",
        "registered_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "pid": os.getpid(), "accepted_records": EXPECTED_ACCEPTED,
        "tested_records": EXPECTED_TESTED, "unpassed_records_not_packed": EXPECTED_UNPASSED,
        "total_wall_limit_seconds": WALL_LIMIT_SECONDS,
        "state_path": str(state_path), "chain_manifest_path": str(manifest_path),
        "chain_root": str(chain_root), "output_dir": str(output_dir),
        "expected_state_sha256": EXPECTED_STATE_SHA256,
        "expected_chain_manifest_sha256": EXPECTED_CHAIN_MANIFEST_SHA256,
        "scope": "literal portable extraction and correspondence checks; no numerical recomputation",
    }
    registration_path = output_dir / "registration.json"
    write_json(registration_path, registration)
    stage = "input_hashes_and_manifest"
    try:
        if sha256_file(state_path) != EXPECTED_STATE_SHA256:
            raise ValueError("final state SHA-256 does not match frozen lock")
        if sha256_file(manifest_path) != EXPECTED_CHAIN_MANIFEST_SHA256:
            raise ValueError("chain manifest SHA-256 does not match frozen owner lock")
        with manifest_path.open(encoding="utf-8") as stream:
            chain_manifest = json.load(stream)
        blocks = chain_manifest.get("blocks")
        if len(blocks) != 17 or chain_manifest.get("task") != "upper_cover066":
            raise ValueError("expected the frozen 17-block upper_cover066 chain")
        if chain_manifest.get("global_g_or_C5_claim") is not False:
            raise ValueError("unexpected global claim flag in chain manifest")
        assert_time(started, "manifest verification")

        stage = "verify_blocks_and_collect_passed_rows"
        accepted_rows = []
        source_file_hashes = []
        cumulative_tested = 0
        cumulative_accepted = 0
        cumulative_unpassed = 0
        for block_number, block in enumerate(blocks, start=1):
            if block.get("block") != block_number:
                raise ValueError(f"unexpected block sequence at {block_number}")
            cert_path = chain_root / f"batch{block_number:04d}" / "certificates.json"
            cert_hash = sha256_file(cert_path)
            if cert_hash != block.get("certificates_sha256"):
                raise ValueError(f"block {block_number} certificates hash mismatch")
            source_file_hashes.append((cert_path, cert_hash))
            with cert_path.open(encoding="utf-8") as stream:
                rows = json.load(stream)
            cumulative_tested += len(rows)
            for row in rows:
                passed = row["certificate"].get("passed")
                if not isinstance(passed, bool):
                    raise ValueError(f"block {block_number} has a non-boolean certificate passed flag")
                if passed:
                    accepted_rows.append(row)
                    cumulative_accepted += 1
                else:
                    cumulative_unpassed += 1
            if block.get("cumulative_tested") != cumulative_tested:
                raise ValueError(f"block {block_number} cumulative tested count mismatch")
            if block.get("cumulative_accepted") != cumulative_accepted:
                raise ValueError(f"block {block_number} cumulative accepted count mismatch")
            assert_time(started, f"block {block_number} processing")

        if (cumulative_tested, cumulative_accepted, cumulative_unpassed) != (
            EXPECTED_TESTED, EXPECTED_ACCEPTED, EXPECTED_UNPASSED
        ):
            raise ValueError(
                "unexpected chain totals: "
                f"tested={cumulative_tested}, passed={cumulative_accepted}, unpassed={cumulative_unpassed}"
            )

        stage = "load_final_state_and_match"
        with state_path.open(encoding="utf-8") as stream:
            state = json.load(stream)
        if len(state.get("cache", {})) != 178722:
            raise ValueError("final cache entry count differs from the frozen expected count 178722")
        if len(accepted_rows) != EXPECTED_ACCEPTED:
            raise ValueError("accepted row count differs from 15143")

        wanted = {}
        for row in accepted_rows:
            node = row["node"]
            key = (node[0], node[1])
            if key in wanted:
                raise ValueError(f"duplicate accepted node in block certificates: {node!r}")
            wanted[key] = row["certificate"]
        matched = set()
        for passed_row in state["passed"]:
            node = passed_row["node"]
            key = (node[0], node[1])
            if key in wanted:
                if key in matched:
                    raise ValueError(f"accepted node occurs more than once in final passed list: {node!r}")
                if passed_row["certificate"] != wanted[key]:
                    raise ValueError(f"certificate differs between block row and final passed list: {node!r}")
                matched.add(key)
        if len(matched) != EXPECTED_ACCEPTED:
            missing = EXPECTED_ACCEPTED - len(matched)
            raise ValueError(f"{missing} accepted block rows have no exact final passed match")

        stage = "extract_portable_records"
        output_rows = []
        required_cert = {"Delta_lower", "loggradient_directed", "minimum_log_margin", "worst_offset", "passed"}
        for index, row in enumerate(accepted_rows):
            node = row["node"]
            certificate = row["certificate"]
            if not required_cert.issubset(certificate) or not certificate["passed"]:
                raise ValueError(f"incomplete passed certificate at accepted row {index}")
            anchor = certificate["anchor"]
            cache_record = state["cache"][anchor]
            point = cache_record["point"]
            if not isinstance(point, list) or ",".join(point) != anchor:
                raise ValueError(f"cache point/key mismatch at accepted row {index}")
            values = cache_record["values"]
            if len(values) != 5 or any(len(pair) != 2 for pair in values):
                raise ValueError(f"cache record lacks five complete rational interval pairs at row {index}")
            bits = cache_record["bits"]
            if not isinstance(bits, int) or isinstance(bits, bool):
                raise ValueError(f"cache precision bits is not an integer at row {index}")
            root = state["roots"][node[0]]
            if len(root) != 4 or any(not isinstance(x, int) or isinstance(x, bool) for x in root):
                raise ValueError(f"root coordinates are not four integers at row {index}")
            if len(row["box_lower"]) != 4 or len(row["box_upper"]) != 4:
                raise ValueError(f"box bounds are not four-dimensional at row {index}")
            output_rows.append({
                "node": node,
                "root_integers": root,
                "box_lower": row["box_lower"],
                "box_upper": row["box_upper"],
                "cache_record": {"point": point, "values": values, "bits": bits},
                "source_failed_certificate_sha256": row["source_failed_certificate_sha256"],
                "old_minimum_log_margin": row["old_minimum_log_margin"],
                "certificate": certificate,
            })
        assert_time(started, "record extraction")

        stage = "write_compressed_packet"
        packet_path = output_dir / "records.json.gz"
        with packet_path.open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", compresslevel=9, fileobj=raw, mtime=0) as compressed:
                with io.TextIOWrapper(compressed, encoding="utf-8", write_through=True) as text:
                    json.dump({"format": "upper_portable071_records_v1", "task": "upper_cover066", "records": output_rows},
                              text, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
        assert_time(started, "packet compression")

        stage = "post_extraction_source_hashes"
        if sha256_file(state_path) != EXPECTED_STATE_SHA256:
            raise ValueError("final state hash changed after extraction")
        if sha256_file(manifest_path) != EXPECTED_CHAIN_MANIFEST_SHA256:
            raise ValueError("chain manifest hash changed after extraction")
        for cert_path, locked_hash in source_file_hashes:
            if sha256_file(cert_path) != locked_hash:
                raise ValueError(f"block certificate source changed after extraction: {cert_path}")
        assert_time(started, "post-extraction source hash check")

        stage = "write_packet_manifest"
        manifest = {
            "task": "upper_portable071", "status": "extracted",
            "records": len(output_rows), "tested_source_rows": cumulative_tested,
            "unpassed_rows_excluded": cumulative_unpassed,
            "final_cache_entries_read_only": len(state["cache"]),
            "pid": os.getpid(), "elapsed_seconds": round(time.monotonic() - started, 6),
            "inputs": {
                "state_sha256": EXPECTED_STATE_SHA256,
                "chain_manifest_sha256": EXPECTED_CHAIN_MANIFEST_SHA256,
                "block_certificate_sha256": [b["certificates_sha256"] for b in blocks],
            },
            "builder_sha256": sha256_file(Path(__file__).resolve()),
            "packet_file": packet_path.name, "packet_sha256": sha256_file(packet_path),
            "packet_bytes": packet_path.stat().st_size,
            "gzip_mtime": 0,
            "numeric_content": "stored rational strings copied verbatim; no H/gradient/log/margin recomputation",
            "scope": "15143 accepted parent records only; no frontier/tree materialization",
            "claim_limit": "underlying H ball is shared production input; not independent value certification; no global g or C5 bound",
        }
        assert_time(started, "manifest finalization")
        write_json(output_dir / "manifest.json", manifest)
        registration["status"] = "completed"
        registration["elapsed_seconds"] = round(time.monotonic() - started, 6)
        write_json(registration_path, registration)
        print(json.dumps(manifest, sort_keys=True, indent=2))
        return 0
    except BaseException as exc:
        write_failure(output_dir, stage, exc, started)
        raise


if __name__ == "__main__":
    sys.exit(main())
