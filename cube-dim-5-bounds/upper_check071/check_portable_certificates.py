#!/usr/bin/env python3
"""Independent exact checker for the portable upper_cover066 certificate packet.

This checks stored H/gradient interval arithmetic and the local box certificates.
It makes no H, gradient, kernel, or special-function calls.
"""
from __future__ import annotations

import argparse
import datetime as dt
from fractions import Fraction as F
from itertools import product
import gzip
import hashlib
import json
import os
from pathlib import Path
import signal
import sys
import time


PACKET_SHA256 = "54fcdd5af1745b4b4e8e8c97791f5bd27f102cfe9876ece22c6fb455f871ece7"
SOURCE_ATANH_SHA256 = "2ad20352bf6a013d9e861f39c18e72903893dd4ee986cb8c740fab38536f9341"
RECORD_COUNT = 15143
SCALE_BITS = 896
SERIES_TERMS = 224
SCALE = 1 << SCALE_BITS
PREFIX_LIMIT_SECONDS = 540
TOTAL_LIMIT_SECONDS = 600
ALARM_LIMIT_SECONDS = 590
CHECKPOINT_EVERY = 64
TARGET = F(1, 10000)
PROGRAM_STARTED = time.monotonic()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def canonical_json(path: Path, obj: object) -> None:
    data = json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True) + "\n"
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", encoding="utf-8") as stream:
        stream.write(data)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(tmp, path)


def ceildiv(a: int, b: int) -> int:
    return -(-a // b)


def log_atanh_fixed(z: F) -> tuple[int, int]:
    """Fixed-point interval for log((1+z)/(1-z)), ported from pinned source."""
    if not F(0) <= z <= F(1, 3):
        raise ValueError("atanh argument outside [0, 1/3]")
    zl = z.numerator * SCALE // z.denominator
    zu = ceildiv(z.numerator * SCALE, z.denominator)
    z2l = zl * zl // SCALE
    z2u = ceildiv(zu * zu, SCALE)
    pl, pu = zl, zu
    lo = hi = 0
    for j in range(SERIES_TERMS):
        d = 2 * j + 1
        lo += 2 * pl // d
        hi += ceildiv(2 * pu, d)
        pl = pl * z2l // SCALE
        pu = ceildiv(pu * z2u, SCALE)
    hi += ceildiv(2 * pu * SCALE, (2 * SERIES_TERMS + 1) * (SCALE - z2u))
    return lo, hi


LOG2_FIXED = log_atanh_fixed(F(1, 3))


def log_fixed(x: F) -> tuple[int, int]:
    """Fixed-point interval [lo,hi]/SCALE enclosing log(x)."""
    if x <= 0:
        raise ValueError("log argument must be positive")
    exponent = x.numerator.bit_length() - x.denominator.bit_length()
    y = x / (F(2) ** exponent)
    if y < 1:
        exponent -= 1
        y *= 2
    if y >= 2:
        exponent += 1
        y /= 2
    if not F(1) <= y < F(2):
        raise ArithmeticError("power-of-two normalization failed")
    lo, hi = log_atanh_fixed((y - 1) / (y + 1))
    if exponent >= 0:
        return lo + exponent * LOG2_FIXED[0], hi + exponent * LOG2_FIXED[1]
    return lo + exponent * LOG2_FIXED[1], hi + exponent * LOG2_FIXED[0]


def fraction(text: str) -> F:
    if not isinstance(text, str):
        raise TypeError("expected an exact rational string")
    value = F(text)
    if str(value) != text:
        raise ValueError(f"non-canonical rational string: {text!r}")
    return value


def reconstructed_box(root: list[int], path: str) -> tuple[list[F], list[F]]:
    if len(root) != 4 or any(type(x) is not int for x in root):
        raise ValueError("root_integers must contain exactly four integers")
    if any(x < 0 or x > 31 for x in root) or any(root[i] < root[i + 1] for i in range(3)):
        raise ValueError("root integers must satisfy 31 >= x0 >= x1 >= x2 >= x3 >= 0")
    if not isinstance(path, str) or any(bit not in "01" for bit in path):
        raise ValueError("node path must contain only interleaved 0/1 bits")
    lo = [F(x, 32) for x in root]
    hi = [x + F(1, 32) for x in lo]
    for depth, bit in enumerate(path):
        axis = depth % 4
        mid = (lo[axis] + hi[axis]) / 2
        if bit == "0":
            hi[axis] = mid
        else:
            lo[axis] = mid
    return lo, hi


def linear_lower(offset: tuple[F, ...], directed: list[tuple[F, F]]) -> F:
    return sum((min(glo * v, ghi * v) for v, (glo, ghi) in zip(offset, directed)), F(0))


def old_q(offset: tuple[F, ...]) -> F:
    return sum((v * v for v in offset), F(0)) / 4 - sum(offset, F(0)) ** 2 / 36


def psi(offset: tuple[F, ...]) -> F:
    return sum((v * v for v in offset), F(0)) / 4 - sum((abs(v) for v in offset), F(0)) ** 2 / 36


def check_entry(row: dict) -> tuple[F, F, F, tuple[F, ...]]:
    node = row["node"]
    if not isinstance(node, list) or len(node) != 2 or type(node[0]) is not int or node[0] < 0:
        raise ValueError("node must be [integer root index, bit path]")
    root = row["root_integers"]
    lo, hi = reconstructed_box(root, node[1])
    saved_lo = [fraction(x) for x in row["box_lower"]]
    saved_hi = [fraction(x) for x in row["box_upper"]]
    if saved_lo != lo or saved_hi != hi:
        raise ValueError(f"root/path box reconstruction mismatch at {node!r}")

    cert = row["certificate"]
    anchor_text = cert["anchor"].split(",")
    point_text = row["cache_record"]["point"]
    if point_text != anchor_text:
        raise ValueError(f"cache point and certificate anchor differ at {node!r}")
    if len(point_text) != 5 or fraction(point_text[0]) != 1:
        raise ValueError(f"cache anchor is not (1, four coordinates) at {node!r}")
    point = [fraction(x) for x in point_text]
    if ",".join(point_text) != cert["anchor"]:
        raise ValueError(f"cache key does not match cache point at {node!r}")
    if any(not l <= x <= h for l, x, h in zip(lo, point[1:], hi)):
        raise ValueError(f"cache point falls outside reconstructed box at {node!r}")

    cache_record = row["cache_record"]
    values_text = cache_record["values"]
    if len(values_text) != 5 or any(len(pair) != 2 for pair in values_text):
        raise ValueError(f"cache record must contain five interval pairs at {node!r}")
    values = [tuple(fraction(x) for x in pair) for pair in values_text]
    if any(a > b for a, b in values):
        raise ValueError(f"reversed H/gradient interval at {node!r}")
    Hlo, Hhi = values[0]
    if not 0 < Hlo <= Hhi:
        raise ValueError(f"invalid positive H interval at {node!r}")
    if type(cache_record["bits"]) is not int or cache_record["bits"] < 1:
        raise ValueError(f"invalid cached precision bits at {node!r}")

    directed_text = cert["loggradient_directed"]
    if len(directed_text) != 4 or any(len(pair) != 2 for pair in directed_text):
        raise ValueError(f"four directed gradient intervals required at {node!r}")
    directed = [tuple(fraction(x) for x in pair) for pair in directed_text]
    for (glo, ghi), (dlo, dhi) in zip(directed, values[1:]):
        if glo > ghi:
            raise ValueError(f"reversed directed log-gradient interval at {node!r}")
        ratios = (dlo / Hlo, dlo / Hhi, dhi / Hlo, dhi / Hhi)
        if glo > min(ratios) or max(ratios) > ghi:
            raise ValueError(f"stored log-gradient interval misses derivative/H range at {node!r}")

    delta = fraction(cert["Delta_lower"])
    if delta <= 0:
        raise ValueError(f"nonpositive Delta lower endpoint at {node!r}")
    loglo, _ = log_fixed(Hlo / TARGET)
    if delta.numerator * SCALE > loglo * delta.denominator:
        raise ValueError(f"saved Delta is not below independent fixed-point log lower endpoint at {node!r}")

    lower = tuple(l - x for l, x in zip(lo, point[1:]))
    upper = tuple(h - x for h, x in zip(hi, point[1:]))
    old_values = []
    for bits in product((0, 1), repeat=4):
        v = tuple(upper[j] if bits[j] else lower[j] for j in range(4))
        old_values.append(delta + linear_lower(v, directed) - 4096 * old_q(v))
    old_min = min(old_values)
    saved_old_min = fraction(row["old_minimum_log_margin"])
    if old_min != saved_old_min or old_min > 0:
        raise ValueError(f"old 16-corner minimum mismatch or not <=0 at {node!r}")

    axes = [sorted({lower[j], F(0), upper[j]}) for j in range(4)]
    if any(len(axis) != 3 for axis in axes):
        raise ValueError(f"81-point certificate has a degenerate coordinate at {node!r}")
    new_values = []
    offsets = []
    for v in product(*axes):
        new_values.append(delta + linear_lower(v, directed) - 4096 * psi(v))
        offsets.append(v)
    new_min = min(new_values)
    new_index = new_values.index(new_min)
    new_worst_offset = tuple(fraction(x) for x in cert["worst_offset"])
    if new_min != fraction(cert["minimum_log_margin"]) or new_min <= 0:
        raise ValueError(f"new 81-point minimum mismatch or not >0 at {node!r}")
    if offsets[new_index] != new_worst_offset:
        raise ValueError(f"new 81-point worst offset mismatch at {node!r}")
    if cert["checked_points"] != 81 or cert["passed"] is not True:
        raise ValueError(f"certificate status/point count mismatch at {node!r}")
    if cert["mechanism"] != "signed_basis_81":
        raise ValueError(f"unexpected certificate mechanism at {node!r}")
    return delta, old_min, new_min, new_worst_offset


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bundle", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path, help="new summary JSON path; sidecars use this stem")
    ap.add_argument("--preflight-only", action="store_true", help="hash and schema check only; do not parse rationals or recompute margins")
    args = ap.parse_args()

    started = PROGRAM_STARTED
    bundle = args.bundle.expanduser().resolve(strict=True)
    output = args.output.expanduser().resolve()
    registration_path = output.with_name(output.name + ".registration.json")
    prefix_path = output.with_name(output.name + ".prefix.json")
    preflight_path = output.with_name(output.name + ".preflight.json")
    for path in (output, registration_path, prefix_path, preflight_path):
        if path.exists():
            raise SystemExit(f"refusing to overwrite existing output or sidecar: {path}")
    output.parent.mkdir(parents=True, exist_ok=True)

    source_code_hash = sha256_file(Path(__file__).resolve())
    registration = {
        "task": "upper_check071",
        "status": "registered_before_bundle_load",
        "registered_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "pid": os.getpid(),
        "bundle_path": str(bundle),
        "bundle_sha256": PACKET_SHA256,
        "records": RECORD_COUNT,
        "mode": "preflight_only" if args.preflight_only else "full_exact_check",
        "prefix_checkpoint_stop_seconds": PREFIX_LIMIT_SECONDS,
        "total_wall_limit_seconds": TOTAL_LIMIT_SECONDS,
        "checkpoint_every_records": CHECKPOINT_EVERY,
        "fixed_log_bits": SCALE_BITS,
        "fixed_log_terms": SERIES_TERMS,
        "fixed_log_source_sha256": SOURCE_ATANH_SHA256,
        "checker_sha256": source_code_hash,
        "scope": "portable local certificate collection only; shared H ball; no global g or C5 claim",
    }
    canonical_json(registration_path, registration)

    stage = "bundle_hash_and_load"
    checked = 0
    old_worst: F | None = None
    new_worst: F | None = None
    minimum_delta: F | None = None
    log_lower_checks = 0
    previous_handler = None
    try:
        def alarm_timeout(signum, frame):
            raise TimeoutError("590-second SIGALRM safety limit reached")

        previous_handler = signal.signal(signal.SIGALRM, alarm_timeout)
        signal.alarm(ALARM_LIMIT_SECONDS)
        if sha256_file(bundle) != PACKET_SHA256:
            raise ValueError("bundle SHA-256 does not match the reviewed frozen packet")
        with gzip.open(bundle, "rt", encoding="utf-8") as stream:
            packet = json.load(stream)
        if time.monotonic() - started > TOTAL_LIMIT_SECONDS:
            raise TimeoutError("180-second wall limit exceeded loading bundle")
        if packet.get("format") != "upper_portable071_records_v1" or packet.get("task") != "upper_cover066":
            raise ValueError("unexpected packet format or source task")
        rows = packet.get("records")
        if not isinstance(rows, list) or len(rows) != RECORD_COUNT:
            raise ValueError("packet must contain exactly 15,143 records")

        stage = "metadata_schema_preflight"
        seen = set()
        for i, row in enumerate(rows):
            node = row.get("node")
            if not isinstance(node, list) or len(node) != 2 or type(node[0]) is not int or not isinstance(node[1], str):
                raise ValueError(f"bad node schema at row {i}")
            key = (node[0], node[1])
            if key in seen:
                raise ValueError(f"duplicate node in packet: {node!r}")
            seen.add(key)
            root = row.get("root_integers", [])
            if len(root) != 4 or any(type(x) is not int for x in root):
                raise ValueError(f"bad root integer schema at row {i}")
            if any(x < 0 or x > 31 for x in root) or any(root[j] < root[j + 1] for j in range(3)):
                raise ValueError(f"root integers are not ordered in [0,31] at row {i}")
            if len(row.get("box_lower", [])) != 4 or len(row.get("box_upper", [])) != 4:
                raise ValueError(f"bad root/box dimensions at row {i}")
            cache = row.get("cache_record", {})
            cert = row.get("certificate", {})
            if len(cache.get("point", [])) != 5 or len(cache.get("values", [])) != 5:
                raise ValueError(f"bad cache point/value dimensions at row {i}")
            if len(cert.get("loggradient_directed", [])) != 4:
                raise ValueError(f"bad directed gradient schema at row {i}")
            if not all(isinstance(x, str) for x in cache.get("point", [])):
                raise ValueError(f"cache point must use string coordinates at row {i}")
            if any(len(pair) != 2 or not all(isinstance(x, str) for x in pair) for pair in cache["values"]):
                raise ValueError(f"cache values must be five string interval pairs at row {i}")
            if any(len(pair) != 2 or not all(isinstance(x, str) for x in pair) for pair in cert["loggradient_directed"]):
                raise ValueError(f"directed gradients must be four string interval pairs at row {i}")
        if args.preflight_only:
            result = {
                "status": "PREFLIGHT_SCHEMA_PASS",
                "bundle_sha256": PACKET_SHA256,
                "records": len(rows),
                "unique_nodes": len(seen),
                "numerical_recomputation": False,
                "elapsed_seconds": round(time.monotonic() - started, 6),
                "checker_sha256": source_code_hash,
            }
            canonical_json(preflight_path, result)
            registration["status"] = "preflight_complete"
            registration["elapsed_seconds"] = round(time.monotonic() - started, 6)
            canonical_json(registration_path, registration)
            print(json.dumps(result, sort_keys=True, indent=2))
            return 0

        stage = "full_exact_prefix_check"
        def prefix_status(status: str) -> dict:
            return {
                "status": status,
                "checked_records": checked,
                "target_records": RECORD_COUNT,
                "pid": os.getpid(),
                "elapsed_seconds": round(time.monotonic() - started, 6),
                "minimum_delta_seen": None if minimum_delta is None else str(minimum_delta),
                "minimum_old16_seen": None if old_worst is None else str(old_worst),
                "minimum_new81_seen": None if new_worst is None else str(new_worst),
            }

        if time.monotonic() >= started + PREFIX_LIMIT_SECONDS:
            canonical_json(prefix_path, prefix_status("stopped_at_prefix_limit_before_math"))
            registration["status"] = "stopped_at_prefix_limit"
            registration["elapsed_seconds"] = round(time.monotonic() - started, 6)
            canonical_json(registration_path, registration)
            return 2

        for row in rows:
            delta, old_min, new_min, _worst_offset = check_entry(row)
            if minimum_delta is None or delta < minimum_delta:
                minimum_delta = delta
            if old_worst is None or old_min < old_worst:
                old_worst = old_min
            if new_worst is None or new_min < new_worst:
                new_worst = new_min
            checked += 1
            log_lower_checks += 1
            elapsed = time.monotonic() - started
            if checked % CHECKPOINT_EVERY == 0:
                canonical_json(prefix_path, prefix_status("in_progress"))
            if elapsed >= PREFIX_LIMIT_SECONDS:
                canonical_json(prefix_path, prefix_status("stopped_at_prefix_limit"))
                registration["status"] = "stopped_at_prefix_limit"
                registration["elapsed_seconds"] = round(elapsed, 6)
                canonical_json(registration_path, registration)
                return 2
            if elapsed >= TOTAL_LIMIT_SECONDS:
                raise TimeoutError("600-second total wall limit exceeded during exact checking")
        if checked != RECORD_COUNT or log_lower_checks != RECORD_COUNT:
            raise ValueError("full prefix did not check all 15,143 records")
        if time.monotonic() - started >= TOTAL_LIMIT_SECONDS:
            raise TimeoutError("600-second total wall limit exceeded before final summary")

        stage = "write_full_summary"
        summary = {
            "status": "PASS_portable_local_certificate_collection",
            "bundle_sha256": PACKET_SHA256,
            "records": checked,
            "unique_nodes": len(seen),
            "old16_minimum_nonpositive_count": checked,
            "new81_minimum_positive_count": checked,
            "saved_Delta_independent_log_lower_checks": log_lower_checks,
            "old16_minimum_of_minima": str(old_worst),
            "new81_minimum_of_minima": str(new_worst),
            "minimum_saved_Delta": str(minimum_delta),
            "fixed_log_bits": SCALE_BITS,
            "fixed_log_terms": SERIES_TERMS,
            "fixed_log_algorithm_source_sha256": SOURCE_ATANH_SHA256,
            "checker_sha256": source_code_hash,
            "elapsed_seconds": round(time.monotonic() - started, 6),
            "scope": "local certificate collection only; shared H ball; not full-tree, global g, or C5 certification",
            "recomputations": "stored H/gradient interval arithmetic only; no H, gradient, kernel, Arb, or special-function calls",
        }
        canonical_json(output, summary)
        canonical_json(prefix_path, prefix_status("complete"))
        registration["status"] = "completed"
        registration["elapsed_seconds"] = round(time.monotonic() - started, 6)
        canonical_json(registration_path, registration)
        print(json.dumps(summary, sort_keys=True, indent=2))
        return 0
    except BaseException as exc:
        registration["status"] = "failed"
        registration["failure_stage"] = stage
        registration["failure_reason"] = f"{type(exc).__name__}: {exc}"
        registration["elapsed_seconds"] = round(time.monotonic() - started, 6)
        canonical_json(registration_path, registration)
        if time.monotonic() - started <= TOTAL_LIMIT_SECONDS:
            canonical_json(prefix_path, {"status": "failed", "stage": stage,
                                          "reason": f"{type(exc).__name__}: {exc}",
                                          "checked_records": checked,
                                          "elapsed_seconds": round(time.monotonic() - started, 6)})
        raise
    finally:
        signal.alarm(0)
        if previous_handler is not None:
            signal.signal(signal.SIGALRM, previous_handler)


if __name__ == "__main__":
    sys.exit(main())
