#!/usr/bin/env python3
"""Portable exact verifier for the upper_portable096 098 mixed-epsilon packet.

This checker is self-contained and does not read the original source state.
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
import sys
import time


ROOT_COUNT = 52360
PRECISIONS = (160, 192, 256)
G = F(21, 200000)
SUPPORTED_EPSILONS = {F(1, 64), F(1, 128)}
LOG_BITS = 896
LOG_TERMS = 224
LOG_SCALE = 1 << LOG_BITS
PROGRAM_STARTED = time.monotonic()


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


def ceildiv(a: int, b: int) -> int:
    return -(-a // b)


def rat(s: str) -> F:
    if not isinstance(s, str):
        raise TypeError("canonical rational strings required")
    x = F(s)
    if str(x) != s:
        raise ValueError(f"noncanonical rational {s!r}")
    return x


def log_atanh(z: F) -> tuple[int, int]:
    if not F(0) <= z <= F(1, 3):
        raise ValueError("atanh input outside [0,1/3]")
    zl = z.numerator * LOG_SCALE // z.denominator
    zu = ceildiv(z.numerator * LOG_SCALE, z.denominator)
    z2l = zl * zl // LOG_SCALE
    z2u = ceildiv(zu * zu, LOG_SCALE)
    pl, pu = zl, zu
    lo = hi = 0
    for j in range(LOG_TERMS):
        d = 2 * j + 1
        lo += 2 * pl // d
        hi += ceildiv(2 * pu, d)
        pl = pl * z2l // LOG_SCALE
        pu = ceildiv(pu * z2u, LOG_SCALE)
    hi += ceildiv(2 * pu * LOG_SCALE, (2 * LOG_TERMS + 1) * (LOG_SCALE - z2u))
    return lo, hi


LOG2_INTERVAL = log_atanh(F(1, 3))


def log_fixed(x: F) -> tuple[int, int]:
    if x <= 0:
        raise ValueError("log argument must be positive")
    e = x.numerator.bit_length() - x.denominator.bit_length()
    y = x / (F(2) ** e)
    if y < 1:
        e -= 1
        y *= 2
    if y >= 2:
        e += 1
        y /= 2
    if not F(1) <= y < F(2):
        raise ArithmeticError("log normalization failed")
    lo, hi = log_atanh((y - 1) / (y + 1))
    if e >= 0:
        return lo + e * LOG2_INTERVAL[0], hi + e * LOG2_INTERVAL[1]
    return lo + e * LOG2_INTERVAL[1], hi + e * LOG2_INTERVAL[0]


def root_catalog() -> list[tuple[int, int, int, int]]:
    return [(i, j, k, l)
            for i in range(32)
            for j in range(i + 1)
            for k in range(j + 1)
            for l in range(k + 1)]


def get_box(root: tuple[int, ...], path: str) -> tuple[list[F], list[F]]:
    if len(root) != 4 or any(type(x) is not int for x in root):
        raise ValueError("root catalog entry must have four integer coordinates")
    if not (31 >= root[0] >= root[1] >= root[2] >= root[3] >= 0):
        raise ValueError("root catalog entry is not nonincreasing in [0,31]")
    if not isinstance(path, str) or any(bit not in "01" for bit in path):
        raise ValueError("path must consist of binary split bits")
    lo = [F(x, 32) for x in root]
    hi = [x + F(1, 32) for x in lo]
    for depth, bit in enumerate(path):
        j = depth % 4
        mid = (lo[j] + hi[j]) / 2
        if bit == "0":
            hi[j] = mid
        else:
            lo[j] = mid
    return lo, hi


def decode_pairs(pairs: list, bits: int) -> list[tuple[F, F]]:
    if bits not in PRECISIONS or len(pairs) != 5:
        raise ValueError("precision must be 160/192/256 and exactly five interval pairs are required")
    scale = 1 << bits
    result = []
    for pair in pairs:
        if len(pair) != 2 or any(not isinstance(x, str) for x in pair):
            raise ValueError("compact pair must contain decimal integer strings [lower,width]")
        if any(str(int(x)) != x for x in pair):
            raise ValueError("lower/width integer strings must be canonical")
        lower_i, width_i = map(int, pair)
        if width_i < 0:
            raise ValueError("quantized interval width must be nonnegative")
        result.append((F(lower_i, scale), F(lower_i + width_i, scale)))
    return result


def evaluate(node: list, root: tuple[int, ...], point_tail: list[str], epsilon_text: str, bits: int, pairs: list, mode: str) -> dict:
    epsilon = rat(epsilon_text)
    if epsilon not in SUPPORTED_EPSILONS:
        raise ValueError(f"unsupported per-leaf epsilon {epsilon_text!r}")
    _, path = node
    lo, hi = get_box(root, path)
    if len(point_tail) != 4:
        raise ValueError("point_tail must contain four exact rational coordinates")
    point = tuple(rat(x) for x in point_tail)
    if any(not a <= p <= b for a, p, b in zip(lo, point, hi)):
        raise ValueError("point lies outside its reconstructed cell")
    values = decode_pairs(pairs, bits)
    Hlo, Hhi = values[0]
    if not F(0) < Hlo <= Hhi:
        return {"passed": False, "reason": "H lower endpoint is not strictly positive"}
    if Hlo <= G:
        return {"passed": False, "reason": "H lower endpoint does not exceed g"}
    delta_i, _ = log_fixed(Hlo / G)
    delta = F(delta_i, LOG_SCALE)
    gradients = []
    for dlo, dhi in values[1:]:
        ratios = (dlo / Hlo, dlo / Hhi, dhi / Hlo, dhi / Hhi)
        gradients.append((min(ratios), max(ratios)))
    lower = tuple(a - p for a, p in zip(lo, point))
    upper = tuple(b - p for b, p in zip(hi, point))

    def lin(v):
        return sum((min(gl * x, gu * x) for x, (gl, gu) in zip(v, gradients)), F(0))

    if mode == "Q=I4/4-J4/36":
        margins = []
        for bits4 in product((0, 1), repeat=4):
            v = tuple(upper[j] if bits4[j] else lower[j] for j in range(4))
            q = sum((x*x for x in v), F(0))/4 - sum(v, F(0))**2/36
            margins.append(delta + lin(v) - q/(epsilon*epsilon))
    elif mode == "signed_basis_81":
        axes = [sorted({lower[j], F(0), upper[j]}) for j in range(4)]
        if any(len(a) != 3 for a in axes):
            return {"passed": False, "reason": "signed 81-point grid degenerates"}
        margins = []
        for v in product(*axes):
            q = sum((x*x for x in v), F(0))/4 - sum((abs(x) for x in v), F(0))**2/36
            margins.append(delta + lin(v) - q/(epsilon*epsilon))
    elif mode == "full_frame_16":
        margins = []
        for bits4 in product((0, 1), repeat=4):
            v = tuple(upper[j] if bits4[j] else lower[j] for j in range(4))
            q = sum((x*x for x in v), F(0))/16
            margins.append(delta + lin(v) - q/(epsilon*epsilon))
    else:
        raise ValueError(f"unknown curvature mode {mode!r}")
    minimum = min(margins)
    return {"passed": minimum > 0, "minimum": minimum, "delta_lower": delta}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bundle-dir", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    args = ap.parse_args()
    start = PROGRAM_STARTED
    bundle = args.bundle_dir.expanduser().resolve(strict=True)
    output = args.output.expanduser().resolve()
    regpath = output.with_name(output.name + ".registration.json")
    prefixpath = output.with_name(output.name + ".prefix.json")
    for p in (output, regpath, prefixpath):
        if p.exists():
            raise SystemExit(f"refusing to overwrite output path {p}")
    output.parent.mkdir(parents=True, exist_ok=True)
    registration = {
        "task": "upper_portable096_checker",
        "status": "registered_before_packet_read",
        "registered_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "pid": os.getpid(),
        "bundle_dir": str(bundle),
        "fixed_log_bits": LOG_BITS,
        "fixed_log_terms": LOG_TERMS,
        "scope": "compact interval certificate verification only; no source-state dependency or H/gradient/kernel calls",
    }
    write_json(regpath, registration)
    totals = {"roots": 0, "leaves": 0, "splits": 0, "log_checks": 0}
    modes: dict[str, int] = {}
    epsilon_counts: dict[str, int] = {}
    quant_bits_counts: dict[str, int] = {str(p): 0 for p in PRECISIONS}
    global_leaf_nodes = set()
    global_leaf_cache_keys = set()
    global_tree_nodes = set()
    minimum_by_mode: dict[str, F] = {}
    minimum_delta = None
    stage = "read_manifest"
    try:
        manifest = json.load((bundle / "manifest.json").open(encoding="utf-8"))
        if manifest.get("task") != "upper_portable096" or manifest.get("status") != "portable_compact_certificate_packet":
            raise ValueError("unexpected packet manifest status")
        if manifest.get("g") != "21/200000" or manifest.get("supported_epsilons") != ["1/64", "1/128"]:
            raise ValueError("packet g or supported epsilon set mismatch")
        source_counts = manifest["source_counts"]
        if source_counts.get("roots") != ROOT_COUNT or source_counts.get("unique_epsilon_point_pairs") != source_counts.get("leaves"):
            raise ValueError("source root/epsilon-point count lock mismatch")
        roots_expected = root_catalog()
        if len(roots_expected) != ROOT_COUNT:
            raise ValueError("canonical root catalog has wrong size")
        chunks = manifest.get("chunks")
        if not chunks:
            raise ValueError("packet manifest has no chunks")
        next_root = 0
        modes_expected = manifest.get("curvature_modes", {})
        quant_expected = manifest.get("quantization_bits", {})
        lower_precision_counts = {str(p): 0 for p in PRECISIONS}
        upgraded_leaf_count = 0

        for chunk_index, meta in enumerate(chunks, start=1):
            stage = f"chunk_{chunk_index:03d}"
            path = bundle / meta["file"]
            if path.stat().st_size != meta["bytes"] or sha256_file(path) != meta["sha256"]:
                raise ValueError(f"chunk lock mismatch: {path.name}")
            with path.open("rb") as f:
                header = f.read(8)
            if header[:2] != b"\x1f\x8b" or header[4:8] != b"\0\0\0\0":
                raise ValueError(f"chunk gzip timestamp/header invalid: {path.name}")
            with gzip.open(path, "rt", encoding="utf-8") as f:
                packet = json.load(f)
            if packet.get("format") != "upper_portable096_mixed_epsilon_root_chunk_v3":
                raise ValueError(f"bad chunk format: {path.name}")
            root_start, root_stop = packet["root_range"]
            if root_start != next_root or not root_start < root_stop <= ROOT_COUNT:
                raise ValueError(f"root ranges are not contiguous: {packet['root_range']}")
            if packet["root_range"] != meta["root_range"]:
                raise ValueError(f"manifest/chunk root range mismatch: {path.name}")
            for rid in range(root_start, root_stop):
                # Canonical root tuple generation is the only root catalog source.
                if roots_expected[rid] != tuple(roots_expected[rid]):
                    raise ValueError("canonical root catalog inconsistency")
            root_nodes = {rid: {} for rid in range(root_start, root_stop)}
            split_map = {rid: {} for rid in range(root_start, root_stop)}
            chunk_modes: dict[str, int] = {}
            chunk_bits: dict[str, int] = {}
            chunk_epsilons: dict[str, int] = {}
            for split in packet["splits"]:
                rid, path_bits = split["node"]
                if type(rid) is not int or not root_start <= rid < root_stop:
                    raise ValueError(f"split outside chunk root range: {split['node']!r}")
                key = (rid, path_bits)
                if key in global_tree_nodes:
                    raise ValueError(f"duplicate global tree node {key!r}")
                global_tree_nodes.add(key)
                if split["axis"] != len(path_bits) % 4 or split["children"] != [[rid, path_bits + "0"], [rid, path_bits + "1"]]:
                    raise ValueError(f"split geometry mismatch at {key!r}")
                if path_bits in root_nodes[rid]:
                    raise ValueError(f"duplicate split path {key!r}")
                root_nodes[rid][path_bits] = "split"
                split_map[rid][path_bits] = split
            for leaf in packet["leaves"]:
                rid, path_bits = leaf["node"]
                if type(rid) is not int or not root_start <= rid < root_stop:
                    raise ValueError(f"leaf outside chunk root range: {leaf['node']!r}")
                key = (rid, path_bits)
                if key in global_tree_nodes:
                    raise ValueError(f"duplicate global tree node {key!r}")
                global_tree_nodes.add(key)
                if key in global_leaf_nodes:
                    raise ValueError(f"duplicate global leaf node {key!r}")
                global_leaf_nodes.add(key)
                root_nodes[rid][path_bits] = "leaf"
                point_tail = leaf["point_tail"]
                if len(point_tail) != 4:
                    raise ValueError(f"bad point tail at {key!r}")
                anchor = "1," + ",".join(point_tail)
                epsilon_text = leaf["epsilon"]
                if leaf.get("cache_epsilon") != epsilon_text:
                    raise ValueError(f"leaf/cache epsilon mismatch at {key!r}")
                if str(rat(epsilon_text)) != epsilon_text:
                    raise ValueError(f"noncanonical per-leaf epsilon at {key!r}")
                cache_key = (epsilon_text, anchor)
                if cache_key in global_leaf_cache_keys:
                    raise ValueError(f"duplicate global epsilon-point cache key {cache_key!r}")
                global_leaf_cache_keys.add(cache_key)
                epsilon_counts[epsilon_text] = epsilon_counts.get(epsilon_text, 0) + 1
                chunk_epsilons[epsilon_text] = chunk_epsilons.get(epsilon_text, 0) + 1
                bits = leaf["quant_bits"]
                if bits not in PRECISIONS:
                    raise ValueError(f"unsupported quantization precision at {key!r}")
                mode = leaf["curvature_mode"]
                if mode != "full_frame_16":
                    raise ValueError(f"unknown curvature mode at {key!r}")
                epsilon = rat(epsilon_text)
                mechanism = leaf.get("source_mechanism")
                quadratic = leaf.get("source_quadratic")
                penalty = leaf.get("source_penalty_coefficient")
                expected_penalty = str(1 / (16 * epsilon * epsilon))
                if penalty != expected_penalty:
                    raise ValueError(f"source penalty coefficient does not match epsilon at {key!r}")
                if epsilon == F(1, 64):
                    if mechanism not in {"full_frame_16", "full_frame_16_new_target_g105"} or quadratic != "Q_full=I4/16":
                        raise ValueError(f"invalid saved 1/64 full-frame source mode at {key!r}")
                elif epsilon == F(1, 128):
                    if mechanism != "075_full_frame_16" or quadratic != "I4/(16epsilon^2)":
                        raise ValueError(f"invalid saved 1/128 full-frame source mode at {key!r}")
                attempts = leaf.get("precision_attempts", [])
                expected_attempt_bits = [p for p in PRECISIONS if p < bits]
                if [a.get("quant_bits") for a in attempts] != expected_attempt_bits:
                    raise ValueError(f"precision upgrade witness sequence mismatch at {key!r}")
                if attempts:
                    upgraded_leaf_count += 1
                previous_interval = None
                for attempt in attempts:
                    attempt_bits = attempt["quant_bits"]
                    attempt_result = evaluate(leaf["node"], roots_expected[rid], point_tail, epsilon_text, attempt_bits,
                                              attempt["value_pairs"], mode)
                    if attempt_result["passed"]:
                        raise ValueError(f"unnecessary precision upgrade at {key!r}: {attempt_bits}")
                    lower_precision_counts[str(attempt_bits)] += 1
                    # Each finer dyadic outward interval must nest inside the previous enclosure.
                    current_intervals = decode_pairs_for_check(attempt["value_pairs"], attempt_bits)
                    if previous_interval is not None:
                        if any(not (coarse[0] <= fine[0] <= fine[1] <= coarse[1])
                               for coarse, fine in zip(previous_interval, current_intervals)):
                            raise ValueError(f"outward precision intervals are not nested at {key!r}")
                    previous_interval = current_intervals
                result = evaluate(leaf["node"], roots_expected[rid], point_tail, epsilon_text, bits,
                                  leaf["value_pairs"], mode)
                if not result["passed"]:
                    raise ValueError(f"compact local certificate did not separate at {key!r}: {result.get('reason')}")
                if previous_interval is not None:
                    final_intervals = decode_pairs_for_check(leaf["value_pairs"], bits)
                    if any(not (coarse[0] <= fine[0] <= fine[1] <= coarse[1])
                           for coarse, fine in zip(previous_interval, final_intervals)):
                        raise ValueError(f"final outward intervals are not nested at {key!r}")
                modes[mode] = modes.get(mode, 0) + 1
                chunk_modes[mode] = chunk_modes.get(mode, 0) + 1
                quant_bits_counts[str(bits)] = quant_bits_counts.get(str(bits), 0) + 1
                chunk_bits[str(bits)] = chunk_bits.get(str(bits), 0) + 1
                if minimum_delta is None or result["delta_lower"] < minimum_delta:
                    minimum_delta = result["delta_lower"]
                if mode not in minimum_by_mode or result["minimum"] < minimum_by_mode[mode]:
                    minimum_by_mode[mode] = result["minimum"]
                totals["log_checks"] += 1

            for rid in range(root_start, root_stop):
                nodes = root_nodes[rid]
                if "" not in nodes:
                    raise ValueError(f"missing root node for catalog id {rid}")
                for p in nodes:
                    if p and nodes.get(p[:-1]) != "split":
                        raise ValueError(f"missing split ancestor for {(rid,p)!r}")
                for p in split_map[rid]:
                    if nodes.get(p+"0") is None or nodes.get(p+"1") is None:
                        raise ValueError(f"split child missing for {(rid,p)!r}")
                leaf_paths = sorted(p for p, t in nodes.items() if t == "leaf")
                if not leaf_paths or any(b.startswith(a) for a, b in zip(leaf_paths, leaf_paths[1:])):
                    raise ValueError(f"leaf paths not prefix-free for root {rid}")
                depth = max(map(len, leaf_paths))
                if sum(1 << (depth - len(p)) for p in leaf_paths) != 1 << depth:
                    raise ValueError(f"root {rid} is not completely covered")
            if len(packet["leaves"]) != meta["leaves"] or len(packet["splits"]) != meta["splits"]:
                raise ValueError(f"chunk record count mismatch: {path.name}")
            if chunk_modes != meta["curvature_modes"] or chunk_bits != meta["quant_bits"] or chunk_epsilons != meta["epsilon_counts"]:
                raise ValueError(f"chunk mode/precision/epsilon histogram mismatch: {path.name}")
            next_root = root_stop
            totals["roots"] += root_stop - root_start
            totals["leaves"] += len(packet["leaves"])
            totals["splits"] += len(packet["splits"])
            write_json(prefixpath, {"status": "in_progress", "chunks_checked": chunk_index,
                                    "total_chunks": len(chunks), **totals,
                                    "actual_unique_leaf_nodes": len(global_leaf_nodes),
                                    "actual_unique_epsilon_point_pairs": len(global_leaf_cache_keys)})
            del packet, root_nodes, split_map
        if next_root != ROOT_COUNT:
            raise ValueError(f"root coverage ends at {next_root}, expected {ROOT_COUNT}")
        expected_counts = {"roots": ROOT_COUNT, "leaves": source_counts["leaves"], "splits": source_counts["splits"]}
        if {k: totals[k] for k in expected_counts} != expected_counts:
            raise ValueError("explicit package source count mismatch")
        if len(global_leaf_nodes) != totals["leaves"] or len(global_leaf_cache_keys) != totals["leaves"]:
            raise ValueError("global leaf-node or epsilon-point uniqueness failure")
        if len(global_tree_nodes) != totals["leaves"] + totals["splits"]:
            raise ValueError("global tree node uniqueness/count failure")
        if modes != modes_expected or quant_bits_counts != quant_expected:
            raise ValueError("global curvature/precision histogram differs from manifest")
        if epsilon_counts != source_counts.get("epsilon_counts") or epsilon_counts != manifest.get("epsilon_counts"):
            raise ValueError("global per-leaf epsilon histogram differs from manifest/source counts")
        if upgraded_leaf_count != manifest.get("adaptive_upgrade_leaf_count"):
            raise ValueError("precision-upgrade leaf count differs from manifest")

        summary = {
            "status": "PASS_portable_compact_closed_root_packet",
            "packet_manifest_sha256": sha256_file(bundle / "manifest.json"),
            "source_locks": manifest["source_locks"],
            "roots": totals["roots"], "leaves": totals["leaves"], "splits": totals["splits"],
            "actual_unique_leaf_nodes": len(global_leaf_nodes),
            "actual_unique_epsilon_point_pairs": len(global_leaf_cache_keys),
            "actual_unique_tree_nodes": len(global_tree_nodes),
            "curvature_modes": modes,
            "epsilon_counts": epsilon_counts,
            "quantization_bits": quant_bits_counts,
            "lower_precision_attempt_counts": lower_precision_counts,
            "adaptive_precision_upgrades": upgraded_leaf_count,
            "g": str(G), "supported_epsilons": sorted(map(str, SUPPORTED_EPSILONS)),
            "fixed_log_bits": LOG_BITS, "fixed_log_terms": LOG_TERMS,
            "minimum_Delta_lower": str(minimum_delta),
            "minimum_margin_by_mode": {k: str(v) for k, v in minimum_by_mode.items()},
            "scope": "local closed-root certificate collection only; shared historical H balls; no global g or Cn claim",
            "no_production_recomputation": "no H, gradient, kernel, Green, Arb, or special-function calls",
            "checker_sha256": sha256_file(Path(__file__).resolve()),
            "elapsed_seconds": round(time.monotonic() - start, 6),
        }
        write_json(output, summary)
        write_json(prefixpath, {"status": "complete", **totals,
                                "actual_unique_leaf_nodes": len(global_leaf_nodes),
                                "actual_unique_epsilon_point_pairs": len(global_leaf_cache_keys),
                                "quantization_bits": quant_bits_counts})
        registration["status"] = "completed"
        registration["elapsed_seconds"] = round(time.monotonic() - start, 6)
        registration["checker_sha256"] = summary["checker_sha256"]
        write_json(regpath, registration)
        print(json.dumps({k: summary[k] for k in ("status", "roots", "leaves", "splits", "actual_unique_epsilon_point_pairs", "curvature_modes", "quantization_bits", "elapsed_seconds", "scope")}, indent=2))
        return 0
    except BaseException as exc:
        registration["status"] = "failed"
        registration["failure_stage"] = stage
        registration["failure_reason"] = f"{type(exc).__name__}: {exc}"
        registration["elapsed_seconds"] = round(time.monotonic() - start, 6)
        write_json(regpath, registration)
        write_json(prefixpath, {"status": "failed", "stage": stage,
                                "reason": f"{type(exc).__name__}: {exc}", **totals})
        raise


def decode_pairs_for_check(pairs: list, bits: int) -> list[tuple[F, F]]:
    # Kept separate from evaluate so lower-precision witnesses can be compared.
    if bits not in PRECISIONS or len(pairs) != 5:
        raise ValueError("invalid precision attempt")
    scale = 1 << bits
    out = []
    for pair in pairs:
        if len(pair) != 2 or any(not isinstance(x, str) or str(int(x)) != x for x in pair):
            raise ValueError("invalid quantized integer pair")
        lo_i, width_i = map(int, pair)
        if width_i < 0:
            raise ValueError("negative interval width")
        out.append((F(lo_i, scale), F(lo_i + width_i, scale)))
    return out


if __name__ == "__main__":
    sys.exit(main())
