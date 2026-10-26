#!/usr/bin/env python3
"""Independently check a portable, closed-root upper_cover070 packet."""
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


SOURCE_ATANH_SHA256 = "2ad20352bf6a013d9e861f39c18e72903893dd4ee986cb8c740fab38536f9341"
ROOT_COUNT = 52360
LEAF_COUNT = 97960
SPLIT_COUNT = 45600
MAX_CHUNK_BYTES = 25_000_000
MAX_TOTAL_BYTES = 2_000_000_000
SCALE_BITS = 896
SERIES_TERMS = 224
SCALE = 1 << SCALE_BITS
TARGET = F(1, 10000)
EXPECTED_MODES = {
    "Q=I4/4-J4/36": 77453,
    "signed_basis_81": 16341,
    "full_frame_16": 4166,
}
PROGRAM_STARTED = time.monotonic()


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


def ceildiv(a: int, b: int) -> int:
    return -(-a // b)


def log_atanh(z: F) -> tuple[int, int]:
    if not F(0) <= z <= F(1, 3):
        raise ValueError("fixed-point log argument outside [0,1/3]")
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


LOG2_FIXED = log_atanh(F(1, 3))


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
        raise ArithmeticError("fixed-point log normalization failed")
    lo, hi = log_atanh((y - 1) / (y + 1))
    if e >= 0:
        return lo + e * LOG2_FIXED[0], hi + e * LOG2_FIXED[1]
    return lo + e * LOG2_FIXED[1], hi + e * LOG2_FIXED[0]


def rat(s: str) -> F:
    if not isinstance(s, str):
        raise TypeError("exact rational strings required")
    v = F(s)
    if str(v) != s:
        raise ValueError(f"noncanonical rational: {s!r}")
    return v


def box(root: list[int], path: str) -> tuple[list[F], list[F]]:
    if len(root) != 4 or any(type(x) is not int or x < 0 or x > 31 for x in root):
        raise ValueError("root must have four integers in [0,31]")
    if any(root[i] < root[i + 1] for i in range(3)):
        raise ValueError("root coordinates must be nonincreasing")
    if not isinstance(path, str) or any(c not in "01" for c in path):
        raise ValueError("path must contain interleaved binary split bits")
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


def linear_lb(v: tuple[F, ...], g: list[tuple[F, F]]) -> F:
    return sum((min(a * x, b * x) for x, (a, b) in zip(v, g)), F(0))


def q_single(v: tuple[F, ...]) -> F:
    return sum((x * x for x in v), F(0)) / 4 - sum(v, F(0)) ** 2 / 36


def psi_81(v: tuple[F, ...]) -> F:
    return sum((x * x for x in v), F(0)) / 4 - sum((abs(x) for x in v), F(0)) ** 2 / 36


def q_full(v: tuple[F, ...]) -> F:
    return sum((x * x for x in v), F(0)) / 16


def expected_root_catalog() -> list[tuple[int, int, int, int]]:
    """Canonical nested enumeration of all 31>=i>=j>=k>=l>=0 roots."""
    catalog = [
        (i, j, k, l)
        for i in range(32)
        for j in range(i + 1)
        for k in range(j + 1)
        for l in range(k + 1)
    ]
    return catalog


def guard_root_range(root_range: list[int], next_root: int, total_roots: int) -> int:
    if not isinstance(root_range, list) or len(root_range) != 2:
        raise ValueError("root range must be [start, stop]")
    start, stop = root_range
    if type(start) is not int or type(stop) is not int:
        raise ValueError("root range endpoints must be integers")
    if start != next_root:
        raise ValueError(f"noncontiguous/duplicate root range: expected start {next_root}, got {start}")
    if stop <= start or stop > total_roots:
        raise ValueError(f"invalid root range [{start},{stop})")
    return stop


def guard_catalog_entry(root_id: int, actual: list[int], catalog: list[tuple[int, int, int, int]]) -> None:
    if type(root_id) is not int or not 0 <= root_id < len(catalog):
        raise ValueError(f"root id outside canonical catalog: {root_id!r}")
    if not isinstance(actual, list) or len(actual) != 4 or any(type(x) is not int for x in actual):
        raise ValueError(f"malformed root tuple at id {root_id}")
    if tuple(actual) != catalog[root_id]:
        raise ValueError(f"root tuple does not match canonical catalog at id {root_id}")
    if not (31 >= actual[0] >= actual[1] >= actual[2] >= actual[3] >= 0):
        raise ValueError(f"root tuple is not in the nested nonincreasing catalog at id {root_id}")


def self_test_guards() -> None:
    mini = [(0, 0, 0, 0), (1, 0, 0, 0), (1, 1, 0, 0)]
    full = expected_root_catalog()
    assert len(full) == ROOT_COUNT
    assert full[:5] == [
        (0, 0, 0, 0), (1, 0, 0, 0), (1, 1, 0, 0), (1, 1, 1, 0), (1, 1, 1, 1)
    ]
    assert full[-1] == (31, 31, 31, 31)
    assert guard_root_range([0, 2], 0, 3) == 2
    guard_catalog_entry(1, [1, 0, 0, 0], mini)
    for bad in (
        lambda: guard_root_range([0, 2], 2, 3),  # repeated range
        lambda: guard_root_range([1, 2], 0, 3),  # skipped first ID
        lambda: guard_root_range([2, 4], 2, 3),  # stop beyond catalog
        lambda: guard_catalog_entry(1, [0, 0, 0, 0], mini),  # valid tuple, wrong ID
        lambda: guard_catalog_entry(3, [1, 1, 1, 0], mini),  # root ID out of range
        lambda: guard_catalog_entry(2, [1, 0, 1, 0], mini),  # wrong tuple/order
    ):
        try:
            bad()
        except ValueError:
            continue
        raise AssertionError("a negative root/range guard case was accepted")


def check_leaf(leaf: dict, root: list[int]) -> tuple[str, F, F]:
    node = leaf["node"]
    if not isinstance(node, list) or len(node) != 2 or type(node[0]) is not int or node[0] < 0:
        raise ValueError("malformed leaf node")
    lo, hi = box(root, node[1])
    cert = leaf["certificate"]
    cache = leaf["cache_record"]
    point_s = cache["point"]
    if len(point_s) != 5 or any(not isinstance(x, str) for x in point_s):
        raise ValueError(f"cache point schema mismatch at {node!r}")
    if cert.get("anchor") != ",".join(point_s):
        raise ValueError(f"anchor/cache key mismatch at {node!r}")
    point = [rat(x) for x in point_s]
    if point[0] != 1 or any(not a <= x <= b for a, x, b in zip(lo, point[1:], hi)):
        raise ValueError(f"anchor outside reconstructed box at {node!r}")
    vals_s = cache["values"]
    if len(vals_s) != 5 or any(len(p) != 2 for p in vals_s):
        raise ValueError(f"cache must contain five exact interval pairs at {node!r}")
    vals = [tuple(rat(x) for x in pair) for pair in vals_s]
    if any(a > b for a, b in vals):
        raise ValueError(f"reversed stored interval at {node!r}")
    HL, HU = vals[0]
    if not F(0) < HL <= HU or HL <= TARGET:
        raise ValueError(f"invalid H interval or target gap at {node!r}")
    if type(cache.get("bits")) is not int or cache["bits"] < 1:
        raise ValueError(f"invalid saved precision at {node!r}")

    gl_s = cert["loggradient_directed"]
    if len(gl_s) != 4 or any(len(p) != 2 for p in gl_s):
        raise ValueError(f"four directed physical-gradient intervals required at {node!r}")
    gl = [tuple(rat(x) for x in pair) for pair in gl_s]
    for (glo, ghi), (dlo, dhi) in zip(gl, vals[1:]):
        if glo > ghi:
            raise ValueError(f"reversed directed log-gradient interval at {node!r}")
        ratios = (dlo / HL, dlo / HU, dhi / HL, dhi / HU)
        if glo > min(ratios) or max(ratios) > ghi:
            raise ValueError(f"directed log-gradient interval misses physical derivative/H range at {node!r}")

    delta = rat(cert["Delta_lower"])
    if delta <= 0:
        raise ValueError(f"nonpositive saved Delta at {node!r}")
    loglo, _ = log_fixed(HL / TARGET)
    if delta.numerator * SCALE > loglo * delta.denominator:
        raise ValueError(f"saved Delta exceeds independent fixed-point log lower endpoint at {node!r}")

    lower = tuple(a - x for a, x in zip(lo, point[1:]))
    upper = tuple(b - x for b, x in zip(hi, point[1:]))
    bits_corners = list(product((0, 1), repeat=4))
    mode = cert.get("mechanism", cert.get("quadratic"))
    if mode == "Q=I4/4-J4/36":
        if cert.get("all16_corners") is not True or cert.get("quadratic") != mode:
            raise ValueError(f"old single-basis certificate schema mismatch at {node!r}")
        margins = []
        for bits in bits_corners:
            v = tuple(upper[j] if bits[j] else lower[j] for j in range(4))
            margins.append((delta + linear_lb(v, gl) - 4096 * q_single(v), list(bits)))
        minimum, worst_bits = min(margins)
        if minimum != rat(cert["minimum_log_margin"]) or worst_bits != cert["worst_corner_bits"]:
            raise ValueError(f"old 16-corner minimum/worst bits mismatch at {node!r}")
        if cert.get("passed") is not True or minimum <= 0:
            raise ValueError(f"nonaccepted old single-basis leaf in passed cover: {node!r}")
    elif mode == "signed_basis_81":
        if cert.get("quadratic") != "psi(v)=sum(v^2)/4-sum(abs(v))^2/36" or cert.get("checked_points") != 81:
            raise ValueError(f"81-point certificate schema mismatch at {node!r}")
        axes = [sorted({lower[j], F(0), upper[j]}) for j in range(4)]
        if any(len(a) != 3 for a in axes):
            raise ValueError(f"81-point box has a degenerate coordinate at {node!r}")
        vals81 = []
        offsets = []
        for v in product(*axes):
            vals81.append(delta + linear_lb(v, gl) - 4096 * psi_81(v))
            offsets.append(v)
        minimum = min(vals81)
        worst = offsets[vals81.index(minimum)]
        saved_worst = tuple(rat(x) for x in cert["worst_offset"])
        if minimum != rat(cert["minimum_log_margin"]) or worst != saved_worst:
            raise ValueError(f"81-point minimum/worst offset mismatch at {node!r}")
        if cert.get("passed") is not True or minimum <= 0:
            raise ValueError(f"nonaccepted signed-basis leaf in passed cover: {node!r}")
    elif mode == "full_frame_16":
        if cert.get("quadratic") != "Q_full=I4/16" or cert.get("checked_points") != 16:
            raise ValueError(f"full-frame certificate schema mismatch at {node!r}")
        margins = []
        for bits in bits_corners:
            v = tuple(upper[j] if bits[j] else lower[j] for j in range(4))
            margins.append((delta + linear_lb(v, gl) - 4096 * q_full(v), list(bits)))
        minimum, worst_bits = min(margins)
        if minimum != rat(cert["minimum_log_margin"]) or worst_bits != cert["worst_corner_bits"]:
            raise ValueError(f"full-frame 16-corner minimum/worst bits mismatch at {node!r}")
        if cert.get("passed") is not True or minimum <= 0:
            raise ValueError(f"nonaccepted full-frame leaf in passed cover: {node!r}")
    else:
        raise ValueError(f"unknown certificate mode {mode!r} at {node!r}")
    return mode, delta, minimum


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bundle-dir", type=Path)
    ap.add_argument("--output", type=Path, help="new summary JSON path; sidecars share its stem")
    ap.add_argument("--self-test-guards", action="store_true", help="run only tiny root/range negative cases")
    args = ap.parse_args()
    if args.self_test_guards:
        self_test_guards()
        print("root/range negative guard checks passed")
        return 0
    if args.bundle_dir is None or args.output is None:
        ap.error("--bundle-dir and --output are required unless --self-test-guards is selected")
    began = PROGRAM_STARTED
    bundle_dir = args.bundle_dir.expanduser().resolve(strict=True)
    output = args.output.expanduser().resolve()
    regpath = output.with_name(output.name + ".registration.json")
    prefixpath = output.with_name(output.name + ".prefix.json")
    for p in (output, regpath, prefixpath):
        if p.exists():
            raise SystemExit(f"refusing to overwrite existing output path: {p}")
    output.parent.mkdir(parents=True, exist_ok=True)
    reg = {
        "task": "upper_portable084_independent_check",
        "status": "registered_before_packet_read",
        "registered_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "pid": os.getpid(),
        "bundle_dir": str(bundle_dir),
        "checker_sha256": sha256_file(Path(__file__).resolve()),
        "fixed_log_bits": SCALE_BITS,
        "fixed_log_terms": SERIES_TERMS,
        "fixed_log_source_sha256": SOURCE_ATANH_SHA256,
        "scope": "portable local closed-root certificates; shared historical H balls; no global g or Cn claim",
    }
    write_json(regpath, reg)
    chunk_reports = []
    totals = {"roots": 0, "leaves": 0, "splits": 0, "fixed_log_checks": 0}
    modes = {k: 0 for k in EXPECTED_MODES}
    minimum_by_mode: dict[str, F] = {}
    minimum_delta: F | None = None
    global_leaf_nodes: set[tuple[int, str]] = set()
    global_leaf_anchors: set[str] = set()
    global_tree_nodes: set[tuple[int, str]] = set()
    catalog = expected_root_catalog()
    stage = "read_bundle_manifest"
    try:
        manifest_path = bundle_dir / "manifest.json"
        manifest = json.load(manifest_path.open(encoding="utf-8"))
        if manifest.get("task") != "upper_portable084" or manifest.get("status") != "extracted_closed_root_certificate_package":
            raise ValueError("unexpected bundle manifest task/status")
        expected_source_counts = {
            "roots": ROOT_COUNT,
            "leaves": LEAF_COUNT,
            "splits": SPLIT_COUNT,
            "source_cache_entries": 180040,
            "unique_leaf_anchors": LEAF_COUNT,
        }
        if manifest.get("source_counts") != expected_source_counts:
            raise ValueError("package manifest source-count inventory mismatch")
        if manifest.get("certificate_modes") != EXPECTED_MODES:
            raise ValueError("certificate-mode inventory mismatch")
        ranges = manifest.get("partition", {})
        if ranges.get("max_chunk_bytes") != MAX_CHUNK_BYTES or ranges.get("total_chunk_bytes", MAX_TOTAL_BYTES) > MAX_TOTAL_BYTES:
            raise ValueError("chunk size or total package size exceeds registered limits")
        expected_chunks = manifest.get("chunks")
        if not isinstance(expected_chunks, list) or not expected_chunks:
            raise ValueError("manifest has no chunk list")
        expected_chunk_count = (ROOT_COUNT + 2047) // 2048
        if ranges.get("roots_per_chunk") != 2048 or ranges.get("chunk_count") != expected_chunk_count or len(expected_chunks) != expected_chunk_count:
            raise ValueError("manifest chunk partition counts mismatch")
        if len(catalog) != ROOT_COUNT or len(set(catalog)) != ROOT_COUNT:
            raise ValueError("independently generated root catalog has an incorrect size or duplicates")
        next_root = 0

        for chunk_i, chunk_meta in enumerate(expected_chunks, start=1):
            stage = f"chunk_{chunk_i:02d}_hash_load_tree_math"
            path = bundle_dir / chunk_meta["file"]
            if path.stat().st_size != chunk_meta["bytes"] or path.stat().st_size > 25_000_000:
                raise ValueError(f"chunk size mismatch or over limit: {path.name}")
            if sha256_file(path) != chunk_meta["sha256"]:
                raise ValueError(f"chunk SHA-256 mismatch: {path.name}")
            with path.open("rb") as f:
                header = f.read(8)
            if len(header) < 8 or header[:2] != b"\x1f\x8b" or header[4:8] != b"\0\0\0\0":
                raise ValueError(f"invalid gzip header/mtime: {path.name}")
            with gzip.open(path, "rt", encoding="utf-8") as f:
                packet = json.load(f)
            if packet.get("format") != "upper_portable084_closed_root_chunk_v1" or packet.get("task") != "upper_cover070":
                raise ValueError(f"unexpected chunk schema: {path.name}")
            if packet.get("epsilon") != "1/64" or packet.get("target") != "1/10000":
                raise ValueError(f"epsilon/target mismatch: {path.name}")
            root_start, root_stop = packet["root_range"]
            next_root = guard_root_range([root_start, root_stop], next_root, ROOT_COUNT)
            if [root_start, root_stop] != [min((chunk_i - 1) * 2048, ROOT_COUNT), min(chunk_i * 2048, ROOT_COUNT)]:
                raise ValueError(f"chunk root range does not match canonical partition at {path.name}")
            if [root_start, root_stop] != chunk_meta["root_range"] or len(packet["roots"]) != root_stop - root_start:
                raise ValueError(f"root range mismatch: {path.name}")

            roots = packet["roots"]
            for offset, actual_root in enumerate(roots):
                guard_catalog_entry(root_start + offset, actual_root, catalog)
            node_type: dict[tuple[int, str], str] = {}
            split_by_root: dict[int, dict[str, dict]] = {r: {} for r in range(root_start, root_stop)}
            leaf_by_root: dict[int, dict[str, dict]] = {r: {} for r in range(root_start, root_stop)}
            chunk_modes = {k: 0 for k in EXPECTED_MODES}
            for split in packet["splits"]:
                node = split["node"]
                rid, p = node
                if type(rid) is not int or not root_start <= rid < root_stop:
                    raise ValueError(f"split outside chunk root range: {node!r}")
                if node_type.get((rid, p)) is not None:
                    raise ValueError(f"duplicate split/leaf node: {node!r}")
                global_node = (rid, p)
                if global_node in global_tree_nodes:
                    raise ValueError(f"duplicate tree node across chunks: {node!r}")
                global_tree_nodes.add(global_node)
                if split["axis"] != len(p) % 4 or split["children"] != [[rid, p + "0"], [rid, p + "1"]]:
                    raise ValueError(f"split axis/children mismatch: {node!r}")
                node_type[(rid, p)] = "split"
                split_by_root[rid][p] = split
            for leaf in packet["leaves"]:
                node = leaf["node"]
                rid, p = node
                if type(rid) is not int or not root_start <= rid < root_stop:
                    raise ValueError(f"leaf outside chunk root range: {node!r}")
                if node_type.get((rid, p)) is not None:
                    raise ValueError(f"duplicate split/leaf node: {node!r}")
                global_node = (rid, p)
                if global_node in global_tree_nodes:
                    raise ValueError(f"duplicate tree node across chunks: {node!r}")
                global_tree_nodes.add(global_node)
                leaf_key = (rid, p)
                if leaf_key in global_leaf_nodes:
                    raise ValueError(f"duplicate leaf node across chunks: {node!r}")
                global_leaf_nodes.add(leaf_key)
                anchor = leaf.get("certificate", {}).get("anchor")
                if not isinstance(anchor, str) or anchor in global_leaf_anchors:
                    raise ValueError(f"duplicate or malformed leaf anchor across chunks: {node!r}")
                global_leaf_anchors.add(anchor)
                node_type[(rid, p)] = "leaf"
                leaf_by_root[rid][p] = leaf
            if len(roots) != chunk_meta["roots"] or len(packet["splits"]) != chunk_meta["splits"] or len(packet["leaves"]) != chunk_meta["leaves"]:
                raise ValueError(f"chunk row count mismatch: {path.name}")
            if len(packet["leaves"]) != len(roots) + len(packet["splits"]):
                raise ValueError(f"binary tree leaf/internal count identity fails in {path.name}")

            for offset, rid in enumerate(range(root_start, root_stop)):
                root = roots[offset]
                # Validate root shape/order even when a very shallow tree leaf is the only node.
                box(root, "")
                root_nodes = set(split_by_root[rid]) | set(leaf_by_root[rid])
                if "" not in root_nodes:
                    raise ValueError(f"root node missing for root {rid}")
                for p in root_nodes:
                    if p and p[:-1] not in split_by_root[rid]:
                        raise ValueError(f"missing internal ancestor for root/path {(rid, p)!r}")
                for p, split in split_by_root[rid].items():
                    if p + "0" not in root_nodes or p + "1" not in root_nodes:
                        raise ValueError(f"split child missing for root/path {(rid, p)!r}")

                leaf_paths = sorted(leaf_by_root[rid])
                if not leaf_paths or any(b.startswith(a) for a, b in zip(leaf_paths, leaf_paths[1:])):
                    raise ValueError(f"leaf paths are not a nonempty prefix-free set for root {rid}")
                depth = max(map(len, leaf_paths))
                kraft = sum(1 << (depth - len(p)) for p in leaf_paths)
                if kraft != (1 << depth):
                    raise ValueError(f"leaf paths do not close root {rid}: Kraft sum {kraft}")
                for p, leaf in leaf_by_root[rid].items():
                    mode, delta, minimum = check_leaf(leaf, root)
                    chunk_modes[mode] += 1
                    modes[mode] += 1
                    if minimum_by_mode.get(mode) is None or minimum < minimum_by_mode[mode]:
                        minimum_by_mode[mode] = minimum
                    if minimum_delta is None or delta < minimum_delta:
                        minimum_delta = delta
                    totals["fixed_log_checks"] += 1
            totals["roots"] += len(roots)
            totals["leaves"] += len(packet["leaves"])
            totals["splits"] += len(packet["splits"])
            nonzero_chunk_modes = {k: v for k, v in chunk_modes.items() if v}
            if nonzero_chunk_modes != chunk_meta["certificate_modes"]:
                raise ValueError(f"chunk certificate-mode counts mismatch: {path.name}")
            chunk_reports.append({
                "file": path.name,
                "sha256": chunk_meta["sha256"],
                "roots": len(roots),
                "leaves": len(packet["leaves"]),
                "splits": len(packet["splits"]),
                "certificate_modes": nonzero_chunk_modes,
                "tree_closed": True,
                "local_certificates_checked": len(packet["leaves"]),
            })
            write_json(prefixpath, {
                "status": "in_progress",
                "checked_chunks": chunk_i,
                "total_chunks": len(expected_chunks),
                **totals,
                "certificate_modes": modes,
                "last_chunk": path.name,
            })
            del packet, node_type, split_by_root, leaf_by_root

        stage = "final_totals"
        if next_root != ROOT_COUNT:
            raise ValueError(f"root ranges ended at {next_root}, expected {ROOT_COUNT}")
        if totals["roots"] != ROOT_COUNT or totals["leaves"] != LEAF_COUNT or totals["splits"] != SPLIT_COUNT or totals["fixed_log_checks"] != LEAF_COUNT:
            raise ValueError(f"package totals do not close: {totals!r}")
        if len(global_leaf_nodes) != LEAF_COUNT or len(global_leaf_anchors) != LEAF_COUNT:
            raise ValueError("actual global leaf-node or anchor unique-set cardinality mismatch")
        if len(global_tree_nodes) != LEAF_COUNT + SPLIT_COUNT:
            raise ValueError("actual global tree-node unique-set cardinality mismatch")
        if modes != EXPECTED_MODES:
            raise ValueError(f"full certificate mode totals mismatch: {modes!r}")
        summary = {
            "status": "PASS_portable_closed_root_certificate_package",
            "bundle_manifest_sha256": sha256_file(manifest_path),
            "roots": totals["roots"],
            "leaves": totals["leaves"],
            "internal_splits": totals["splits"],
            "unique_nodes": totals["leaves"] + totals["splits"],
            "unique_leaf_nodes_actual_set": len(global_leaf_nodes),
            "unique_leaf_anchors_actual_set": len(global_leaf_anchors),
            "unique_tree_nodes_actual_set": len(global_tree_nodes),
            "canonical_root_catalog_count": len(catalog),
            "root_ranges_contiguous_and_exact": True,
            "roots_with_exact_closed_tree": totals["roots"],
            "certificate_modes": modes,
            "fixed_log_independent_checks": totals["fixed_log_checks"],
            "fixed_log_bits": SCALE_BITS,
            "fixed_log_terms": SERIES_TERMS,
            "fixed_log_algorithm_source_sha256": SOURCE_ATANH_SHA256,
            "minimum_Delta_lower_endpoint": str(minimum_delta),
            "minimum_saved_margin_by_mode": {k: str(v) for k, v in minimum_by_mode.items()},
            "chunks_checked": chunk_reports,
            "checker_sha256": reg["checker_sha256"],
            "runtime_seconds": round(time.monotonic() - PROGRAM_STARTED, 6),
            "scope": "closed-root local certificate collection only; shared historical H intervals; no global g or Cn claim",
            "no_recomputation": "no H, physical gradient, kernel, Arb, Green function, or special-function evaluations; exact arithmetic on stored intervals only",
        }
        write_json(output, summary)
        write_json(prefixpath, {"status": "complete", **totals, "certificate_modes": modes, "chunks_checked": len(chunk_reports)})
        reg["status"] = "completed"
        reg["elapsed_seconds"] = round(time.monotonic() - PROGRAM_STARTED, 6)
        write_json(regpath, reg)
        print(json.dumps({k: summary[k] for k in ("status", "roots", "leaves", "internal_splits", "certificate_modes", "fixed_log_independent_checks", "runtime_seconds", "scope")}, indent=2))
        return 0
    except BaseException as exc:
        reg["status"] = "failed"
        reg["failure_stage"] = stage
        reg["failure_reason"] = f"{type(exc).__name__}: {exc}"
        reg["elapsed_seconds"] = round(time.monotonic() - PROGRAM_STARTED, 6)
        write_json(regpath, reg)
        write_json(prefixpath, {"status": "failed", "stage": stage,
                               "reason": f"{type(exc).__name__}: {exc}", **totals,
                               "certificate_modes": modes})
        raise


if __name__ == "__main__":
    sys.exit(main())
