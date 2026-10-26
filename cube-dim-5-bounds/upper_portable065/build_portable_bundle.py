#!/usr/bin/env python3
"""Extract the frozen upper_cover061 2048-row evidence into a portable bundle.

Standard library only. This program checks byte hashes and row correspondence;
it does not evaluate H, gradients, logarithms, or any mathematical margin.
"""
from __future__ import annotations

import argparse
import datetime as dt
import gzip
import hashlib
import json
import os
from pathlib import Path
import sys
import time


EXPECTED = {
    "state": "d0d62d99d57b3f15d3acb1f4a635b188e196b0de6e155b461460ccbd64a72d31",
    "certificates": "a3331af02c0a3df5548977971b85e62f2dac3b39f591658289569ab4f402b151",
    "owner_manifest": "2e1c6a69d9fd55ff7c09d65cd356c88bf726151f07744c6b7311c4d43cda3504",
}
OWNER_MANIFEST = Path(
    "/Users/zhengzhihao/Desktop/T/Math-Research-cube5-upper/"
    "cube-dim-5-bounds/upper_cover061/manifest.json"
)
RECORDS = 2048
TIME_BUDGET_SECONDS = 120


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, obj: object) -> None:
    payload = json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    path.write_text(payload + "\n", encoding="utf-8")


def fail(path: Path, stage: str, exc: BaseException, started: float) -> None:
    details = {
        "status": "failed",
        "stage": stage,
        "exception_type": type(exc).__name__,
        "reason": str(exc),
        "pid": os.getpid(),
        "elapsed_seconds": round(time.monotonic() - started, 6),
    }
    write_json(path / "failure.json", details)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--state", required=True, type=Path)
    parser.add_argument("--certificates", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    state_path = args.state.expanduser().resolve(strict=True)
    certificates_path = args.certificates.expanduser().resolve(strict=True)
    output_dir = args.output_dir.expanduser().resolve()
    if output_dir.exists():
        raise SystemExit(f"refusing to overwrite existing output directory: {output_dir}")
    output_dir.mkdir(parents=True)
    started = time.monotonic()

    registration = {
        "task": "upper_portable065",
        "status": "registered_before_extraction",
        "registered_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "pid": os.getpid(),
        "records": RECORDS,
        "time_budget_seconds": TIME_BUDGET_SECONDS,
        "state_path": str(state_path),
        "certificates_path": str(certificates_path),
        "owner_manifest_path": str(OWNER_MANIFEST),
        "output_dir": str(output_dir),
        "expected_input_sha256": EXPECTED,
        "scope": "literal extraction and correspondence checks; no H/gradient/log/margin recomputation",
    }
    registration_path = output_dir / "registration.json"
    write_json(registration_path, registration)
    stage = "input_hashes"
    try:
        actual_hashes = {
            "state": sha256_file(state_path),
            "certificates": sha256_file(certificates_path),
            "owner_manifest": sha256_file(OWNER_MANIFEST),
        }
        if actual_hashes != EXPECTED:
            raise ValueError(f"input hash mismatch: {actual_hashes!r}")

        stage = "load_frozen_json"
        with state_path.open(encoding="utf-8") as stream:
            state = json.load(stream)
        with certificates_path.open(encoding="utf-8") as stream:
            certificate_rows = json.load(stream)
        if time.monotonic() - started > TIME_BUDGET_SECONDS:
            raise TimeoutError("120-second extraction budget exceeded while loading frozen inputs")

        stage = "extract_2048_records"
        if len(certificate_rows) != RECORDS or len(state.get("passed", [])) < RECORDS:
            raise ValueError("expected exactly 2048 certificate rows and at least 2048 passed entries")
        tail = state["passed"][-RECORDS:]
        extracted = []
        required_certificate_fields = {
            "Delta_lower", "loggradient_directed", "minimum_log_margin",
            "worst_offset", "checked_points", "passed",
        }
        for index, (source_row, passed_row) in enumerate(zip(certificate_rows, tail)):
            node = source_row["node"]
            certificate = source_row["certificate"]
            if node != passed_row["node"] or certificate != passed_row["certificate"]:
                raise ValueError(f"certificate row {index} differs from corresponding state passed tail object")
            if not certificate.get("passed") or not required_certificate_fields.issubset(certificate):
                raise ValueError(f"certificate row {index} is not a complete passed certificate")

            anchor = certificate["anchor"]
            cache_record = state["cache"][anchor]
            point = cache_record["point"]
            if not isinstance(point, list) or ",".join(point) != anchor:
                raise ValueError(f"cache point/key mismatch for row {index}")
            values = cache_record["values"]
            if len(values) != 5 or any(len(pair) != 2 for pair in values):
                raise ValueError(f"cache record does not contain five complete interval pairs for row {index}")
            if not isinstance(cache_record["bits"], int) or isinstance(cache_record["bits"], bool):
                raise ValueError(f"cache bits is not an integer for row {index}")

            root_index = node[0]
            root = state["roots"][root_index]
            if len(root) != 4 or any(not isinstance(x, int) or isinstance(x, bool) for x in root):
                raise ValueError(f"root is not four integers for row {index}")
            if len(source_row["box_lower"]) != 4 or len(source_row["box_upper"]) != 4:
                raise ValueError(f"box bounds are not four-dimensional for row {index}")
            if source_row["certificate"] != passed_row["certificate"]:
                raise ValueError(f"certificate object mismatch at row {index}")

            extracted.append({
                "node": node,
                "root_integers": root,
                "box_lower": source_row["box_lower"],
                "box_upper": source_row["box_upper"],
                "cache_record": {"point": point, "values": values, "bits": cache_record["bits"]},
                "source_failed_certificate_sha256": source_row["source_failed_certificate_sha256"],
                "old_minimum_log_margin": source_row["old_minimum_log_margin"],
                "certificate": certificate,
            })
        if time.monotonic() - started > TIME_BUDGET_SECONDS:
            raise TimeoutError("120-second extraction budget exceeded while checking rows")

        stage = "write_bundle"
        bundle = {
            "format": "upper_portable065_records_v1",
            "task": "upper_cover061",
            "records": extracted,
        }
        bundle_path = output_dir / "records.json.gz"
        encoded = json.dumps(bundle, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
        with bundle_path.open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", compresslevel=9, fileobj=raw, mtime=0) as compressed:
                compressed.write(encoded)

        if sha256_file(state_path) != EXPECTED["state"] or sha256_file(certificates_path) != EXPECTED["certificates"]:
            raise ValueError("frozen input hash changed after extraction")
        stage = "write_manifest"
        manifest = {
            "task": "upper_portable065",
            "status": "extracted",
            "records": len(extracted),
            "pid": os.getpid(),
            "elapsed_seconds": round(time.monotonic() - started, 6),
            "inputs": {
                "state_sha256": EXPECTED["state"],
                "certificates_sha256": EXPECTED["certificates"],
                "owner_manifest_sha256": EXPECTED["owner_manifest"],
            },
            "builder_sha256": sha256_file(Path(__file__).resolve()),
            "bundle_file": bundle_path.name,
            "bundle_sha256": sha256_file(bundle_path),
            "bundle_bytes": bundle_path.stat().st_size,
            "numeric_content": "exact stored rational strings only; no values recomputed",
            "scope": "portable row evidence only; underlying H ball remains shared production input",
            "claim_limit": "not an independent value certification and not a global g or C5 bound",
        }
        write_json(output_dir / "manifest.json", manifest)
        registration["status"] = "completed"
        registration["elapsed_seconds"] = round(time.monotonic() - started, 6)
        write_json(registration_path, registration)
        print(json.dumps(manifest, sort_keys=True, indent=2))
        return 0
    except BaseException as exc:
        fail(output_dir, stage, exc, started)
        raise


if __name__ == "__main__":
    sys.exit(main())
