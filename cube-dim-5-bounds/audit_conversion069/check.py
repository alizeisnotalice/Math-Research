#!/usr/bin/env python3
"""Independent exact audit of frozen 066 batch0001 logical conversion.

Uses Python standard library only. No author module is imported. Source cache
is streamed once; 066 first-block logical deltas are independently reconstructed.
"""
import argparse
from collections import defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import product, zip_longest
import json
from math import lcm
import os
from pathlib import Path
import resource
import sys
import time

HERE = Path(__file__).resolve().parent
SOURCE_FIELDS = ("Delta_lower", "anchor", "loggradient_directed")


def require(test, message):
    if not test:
        raise AssertionError(message)


def canonical(obj):
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()


def digest(obj):
    return sha256(canonical(obj)).hexdigest()


def file_sha(path):
    h = sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda: f.read(8 << 20), b""):
            h.update(b)
    return h.hexdigest()


def no_duplicates(pairs):
    out = {}
    for k, v in pairs:
        require(k not in out, "duplicate JSON object key")
        out[k] = v
    return out


class Stream:
    """Small-buffer JSON reader; large top-level collections are never loaded."""
    def __init__(self, path):
        self.f = Path(path).open("r", encoding="utf-8")
        self.buf = ""
        self.pos = 0
        self.eof = False
        self.decoder = json.JSONDecoder(object_pairs_hook=no_duplicates)
        self.first_key = True
        self.previous_key = None
        self.char("{")

    def fill(self):
        if self.pos:
            self.buf = self.buf[self.pos:]
            self.pos = 0
        s = self.f.read(1 << 20)
        self.eof = not s
        self.buf += s

    def peek(self):
        while True:
            while self.pos < len(self.buf) and self.buf[self.pos].isspace():
                self.pos += 1
            if self.pos < len(self.buf):
                return self.buf[self.pos]
            if self.eof:
                return ""
            self.fill()

    def char(self, wanted):
        require(self.peek() == wanted, "JSON delimiter mismatch: " + wanted)
        self.pos += 1

    def value(self):
        self.peek()
        while True:
            try:
                value, end = self.decoder.raw_decode(self.buf, self.pos)
            except json.JSONDecodeError:
                require(not self.eof, "truncated or invalid JSON")
                self.fill()
                continue
            # A number could have been cut at a chunk boundary. Wait until its
            # delimiter is present, even when raw_decode already accepts it.
            if end == len(self.buf) and not self.eof:
                self.fill()
                continue
            require(end == len(self.buf) or self.buf[end] in " \n\r\t,]}:", "invalid JSON value delimiter")
            self.pos = end
            return value

    def key(self):
        if self.peek() == "}":
            self.char("}")
            return None
        if not self.first_key:
            self.char(",")
        self.first_key = False
        k = self.value()
        require(isinstance(k, str), "non-string top-level JSON key")
        require(self.previous_key is None or self.previous_key < k, "top-level keys not uniquely sorted")
        self.previous_key = k
        self.char(":")
        return k

    def items(self, object_mode=False):
        begin, end = ("{", "}") if object_mode else ("[", "]")
        self.char(begin)
        first = True
        previous = None
        while self.peek() != end:
            if not first:
                self.char(",")
            first = False
            if object_mode:
                k = self.value()
                require(isinstance(k, str), "non-string cache key")
                require(previous is None or previous < k, "cache keys not uniquely sorted")
                previous = k
                self.char(":")
                yield k, self.value()
            else:
                yield self.value()
        self.char(end)

    def close(self):
        require(self.peek() == "", "trailing JSON content")
        self.f.close()


def node(raw):
    require(isinstance(raw, list) and len(raw) == 2, "bad node shape")
    root, path = raw
    require(type(root) is int and root >= 0 and isinstance(path, str), "bad node types")
    require(all(b in "01" for b in path), "non-binary path")
    return root, path


def physical_box(root_coordinates, n):
    """Interleave bits by coordinate; does not simulate the author's splitter."""
    r, path = n
    require(r < len(root_coordinates), "root index out of range")
    lo, hi = [], []
    for axis, coordinate in enumerate(root_coordinates[r]):
        digits = path[axis::4]
        scale = 1 << len(digits)
        offset = int(digits, 2) if digits else 0
        numerator = coordinate * scale + offset
        lo.append(F(numerator, 32 * scale))
        hi.append(F(numerator + 1, 32 * scale))
    return lo, hi


def covering_parent(n, parents, allow_equal):
    root, path = n
    end = len(path) + (1 if allow_equal else 0)
    hits = [(root, path[:k]) for k in range(end) if (root, path[:k]) in parents]
    require(len(hits) <= 1, "overlapping promoted ancestors")
    return hits[0] if hits else None


def validate_tree(roots, frontier, splits, label):
    require(len(frontier) == len(set(frontier)), label + " duplicate leaf")
    by_root = defaultdict(list)
    internal = set()
    for r, path in frontier:
        require(0 <= r < len(roots), label + " leaf root index")
        by_root[r].append(path)
        internal.update((r, path[:k]) for k in range(len(path)))
    require(set(by_root) == set(range(len(roots))), label + " missing root")
    maximum_depth = 0
    for r, paths in by_root.items():
        depth = max(map(len, paths))
        maximum_depth = max(maximum_depth, depth)
        intervals = sorted(((int(p, 2) if p else 0) << (depth - len(p)),
                            ((int(p, 2) if p else 0) + 1) << (depth - len(p))) for p in paths)
        cursor = 0
        for low, high in intervals:
            require(low == cursor and high > low, label + " prefix overlap or Kraft gap at root " + str(r))
            cursor = high
        require(cursor == 1 << depth, label + " Kraft sum not exactly one")
    require(internal == splits, label + " active splits differ from proper leaf prefixes")
    require(len(frontier) == len(roots) + len(splits), label + " binary tree count identity")
    return {"roots_with_exact_unit_interval_tiling": len(by_root), "leaf_count": len(frontier),
            "active_split_count": len(splits), "maximum_path_depth": maximum_depth,
            "prefix_free": True, "every_root_Kraft_sum_exactly_one": True}


def verify_certificate(row, source, cache, roots):
    n = node(row["node"])
    cert, old_cert = row["certificate"], source["failed_certificate"]
    require(digest(old_cert) == row["source_failed_certificate_sha256"], "old certificate digest mismatch")
    require(old_cert["passed"] is False and F(old_cert["minimum_log_margin"]) <= 0, "source split was not recorded failed")
    require(row["old_minimum_log_margin"] == old_cert["minimum_log_margin"], "old margin changed in ledger")
    for field in SOURCE_FIELDS:
        require(cert[field] == old_cert[field], "inherited certificate field changed: " + field)
    require(type(cert["passed"]) is bool and cert["mechanism"] == "signed_basis_81", "wrong new certificate mechanism")
    require(cert["checked_points"] == 81, "new certificate does not claim 81 points")
    anchor_key = cert["anchor"]
    record = cache[anchor_key]
    require(",".join(record["point"]) == anchor_key, "cache anchor key mismatch")
    point = list(map(F, record["point"]))
    require(len(point) == 5 and point[0] == 1, "anchor is not in x1=1 plane")
    lo, hi = physical_box(roots, n)
    require(list(map(F, row["box_lower"])) == lo and list(map(F, row["box_upper"])) == hi, "ledger physical box mismatch")
    require(row["anchor_inside"] is True, "ledger denies inside anchor")
    require(all(a <= p <= b for a, p, b in zip(lo, point[1:], hi)), "anchor outside reconstructed closed box")
    values = [[F(a), F(b)] for a, b in record["values"]]
    require(len(values) == 5 and all(a <= b for a, b in values), "unordered cache intervals")
    hl, hu = values[0]
    require(hl > 0, "cached H lower endpoint is not positive")
    gradients = [[F(a), F(b)] for a, b in cert["loggradient_directed"]]
    require(len(gradients) == 4 and all(a <= b for a, b in gradients), "unordered directed gradient")
    cross_checks = 0
    for (gl, gu), (dl, du) in zip(gradients, values[1:]):
        for h in (hl, hu):
            for d in (dl, du):
                require(gl * h <= d <= gu * h, "physical dH/H not enclosed")
                cross_checks += 1
    delta = F(cert["Delta_lower"])
    offsets = [(a - p, F(0), b - p) for a, p, b in zip(lo, point[1:], hi)]
    vden = lcm(*(x.denominator for triple in offsets for x in triple))
    gden = lcm(delta.denominator, *(x.denominator for pair in gradients for x in pair))
    integer_offsets = [tuple(int(x * vden) for x in triple) for triple in offsets]
    integer_gradient = [tuple(int(x * gden) for x in pair) for pair in gradients]
    di = int(delta * gden)
    denominator = 36 * gden * vden * vden
    best = None
    minima = set()
    count = 0
    for z in product(*integer_offsets):
        linear = sum(min(a * v, b * v) for (a, b), v in zip(integer_gradient, z))
        curvature = 9 * sum(v * v for v in z) - sum(abs(v) for v in z) ** 2
        require(curvature >= 0, "negative curvature numerator")
        numerator = 36 * di * vden * vden + 36 * vden * linear - 4096 * gden * curvature
        if best is None or numerator < best:
            best, minima = numerator, {z}
        elif numerator == best:
            minima.add(z)
        count += 1
    require(count == 81, "81-point grid count mismatch")
    require(cert["passed"] == (best > 0), "saved verdict differs from exact margin sign")
    if best > 0:
        require(hl > F(1, 10000), "accepted cached H lower endpoint not above target")
    margin = F(best, denominator)
    require(F(cert["minimum_log_margin"]) == margin, "new minimum margin mismatch")
    worst = tuple(F(x) * vden for x in cert["worst_offset"])
    require(len(worst) == 4 and all(x.denominator == 1 for x in worst), "worst offset outside integer grid")
    require(tuple(map(int, worst)) in minima, "declared worst offset not a minimizer")
    return {"node": list(n), "source_split_index": row["source_split_index"], "passed": best > 0, "minimum_log_margin": str(margin), "checked_points": count,
            "gradient_endpoint_cross_checks": cross_checks, "anchor_closed_box": True,
            "anchor_strictly_inside": all(a < p < b for a, p, b in zip(lo, point[1:], hi)),
            "physical_box_lower": list(map(str, lo)), "physical_box_upper": list(map(str, hi))}


def read_source(path, rows, reg):
    selected = {row["source_split_index"] for row in rows}
    require(len(selected) == len(rows), "duplicate source test index")
    anchors = {row["certificate"]["anchor"] for row in rows}
    arrays = {"passed": [], "queue": [], "splits": []}
    selected_splits, cache, metadata = {}, {}, {}
    roots = None
    count = 0
    cache_hash = sha256()
    stream = Stream(path)
    while True:
        field = stream.key()
        if field is None:
            break
        if field == "cache":
            cache_hash.update(b"{")
            first = True
            for key, record in stream.items(True):
                if not first:
                    cache_hash.update(b",")
                first = False
                cache_hash.update(canonical(key)); cache_hash.update(b":"); cache_hash.update(canonical(record))
                count += 1
                if key in anchors:
                    cache[key] = record
            cache_hash.update(b"}")
        elif field in arrays:
            for index, record in enumerate(stream.items()):
                n = node(record if field == "queue" else record["node"])
                arrays[field].append(n)
                if field == "splits":
                    require(record["axis"] == len(n[1]) % 4, "source split axis geometry mismatch")
                    require(record["children"] == [[n[0], n[1] + "0"], [n[0], n[1] + "1"]], "source split child geometry mismatch")
                    if index in selected:
                        selected_splits[index] = record
        else:
            obj = stream.value()
            if field == "roots":
                roots = obj
            elif field in ("epsilon", "target", "model", "active", "inflight", "resume_jobs", "fingerprint061", "status061", "total_kernel_attempts_completed", "total_new_anchor_starts"):
                metadata[field] = obj
    stream.close()
    require(set(cache) == anchors, "tested anchor missing in immutable source cache")
    require(set(selected_splits) == selected, "tested source split index not found")
    require(count == reg["source_counts"]["cache"], "source cache count mismatch")
    require(cache_hash.hexdigest() == reg["source_cache_canonical_sha256"], "source cache canonical hash mismatch")
    require(roots == [[i, j, k, l] for i in range(32) for j in range(i + 1) for k in range(j + 1) for l in range(k + 1)], "source roots differ from independent enumeration")
    require(len(roots) == reg["source_counts"]["roots"], "source root count mismatch")
    for field in arrays:
        require(len(arrays[field]) == reg["source_counts"][field], "source array count mismatch: " + field)
        require(len(arrays[field]) == len(set(arrays[field])), "source node duplicates: " + field)
    require(metadata["epsilon"] == "1/64" and metadata["target"] == "1/10000" and metadata["model"] == "alpha1", "wrong model/epsilon/target")
    require(metadata["active"] is None and metadata["inflight"] == {} and metadata["resume_jobs"] == {}, "source contains live/resume jobs")
    return {"roots": roots, "arrays": arrays, "selected_splits": selected_splits, "cache": cache,
            "metadata": metadata, "cache_entries": count, "cache_canonical_sha256": cache_hash.hexdigest()}


def logical_digest(source_sha, active, accepted, tested):
    return digest({"source_sha256": source_sha, "active_source_indices": active,
                   "accepted_certificates": accepted, "tested_source_split_indices": sorted(tested)})


def audit_selection(source, block_reg, rows):
    arrays = source["arrays"]
    qcount = defaultdict(int)
    for root, path in arrays["queue"]:
        for length in range(len(path)):
            qcount[(root, path[:length])] += 1
    eligible = [(index, n) for index, n in enumerate(arrays["splits"]) if qcount[n] > 0]
    eligible.sort(key=lambda pair: (len(pair[1][1]), -qcount[pair[1]], pair[1][0], pair[1][1]))
    expected = [{"source_split_index": i, "node": list(n), "initial_queue_descendants": qcount[n]}
                for i, n in eligible[:2048]]
    require(block_reg["selected"] == expected, "registered selection is not the deterministic eligible source prefix")
    require(len(rows) == len(expected) == 2048, "first block did not test the registered 2048-prefix")
    for row, selected in zip(rows, expected):
        require(row["node"] == selected["node"] and row["source_split_index"] == selected["source_split_index"], "tested rows do not follow the frozen selection")
        require(node(source["selected_splits"][row["source_split_index"]]["node"]) == node(row["node"]), "tested source index/node mismatch")
    return {"initial_eligible_source_splits": len(eligible), "deterministic_prefix_exact": True,
            "selected": len(expected), "first_selected": expected[0], "last_selected": expected[-1],
            "initial_queue_descendant_counts": sorted(set(x["initial_queue_descendants"] for x in expected))}


def audit_delta(source, rows, retired, reg, block_reg, checkpoint, summary):
    arrays = source["arrays"]
    accepted = [{"node": row["node"], "certificate": row["certificate"]} for row in rows if row["certificate"]["passed"]]
    parents = {node(row["node"]) for row in accepted}
    require(len(parents) == len(accepted), "duplicate accepted parent")
    for p in parents:
        require(covering_parent(p, parents, False) is None, "accepted parent antichain fails")
    old_frontier = arrays["passed"] + arrays["queue"]
    initial_tree = validate_tree(source["roots"], old_frontier, set(arrays["splits"]), "source061")
    initial_active = {k: list(range(len(v))) for k, v in arrays.items()}
    initial_hash = logical_digest(reg["source_sha256"], initial_active, [], [])
    require(initial_hash == block_reg["parent_frontier_sha256"], "first-block logical parent hash mismatch")
    require(block_reg["parent_checkpoint_sha256"] is None and checkpoint["parent_checkpoint_sha256"] is None, "first block has a non-null parent checkpoint")
    active, removed_counts = {}, {}
    for field, nodes in arrays.items():
        expected_removed, keep = [], []
        for index, n in enumerate(nodes):
            parent = covering_parent(n, parents, field == "splits")
            if parent is None:
                keep.append(index)
            else:
                expected_removed.append({"source_index": index, "node": list(n), "unique_new_parent": list(parent)})
        require(retired[field] == expected_removed, "retirement ledger does not exactly equal covered source records: " + field)
        active[field] = keep
        removed_counts[field] = len(expected_removed)
    frontier = [arrays["passed"][i] for i in active["passed"]] + [node(r["node"]) for r in accepted] + [arrays["queue"][i] for i in active["queue"]]
    splits = {arrays["splits"][i] for i in active["splits"]}
    new_tree = validate_tree(source["roots"], frontier, splits, "first_block")
    tested = [r["source_split_index"] for r in rows]
    new_hash = logical_digest(reg["source_sha256"], active, accepted, tested)
    require(new_hash == checkpoint["frontier_sha256"] == summary["frontier_sha256"], "independent logical frontier hash mismatch")
    counts = {"passed": len(active["passed"]) + len(accepted), "queue": len(active["queue"]), "splits": len(active["splits"])}
    require(checkpoint["cumulative_tested"] == summary["cumulative_tested"] == summary["tested"] == summary["selected"] == summary["cursor"] == len(rows), "tested/cursor counters mismatch")
    require(checkpoint["cumulative_accepted"] == summary["cumulative_accepted"] == summary["accepted"] == len(accepted), "accepted counters mismatch")
    require(summary["failed_uncertified"] == len(rows) - len(accepted) and summary["retired"] == removed_counts, "summary retired/failed counters mismatch")
    require(summary["passed"] == counts["passed"] and summary["queue"] == counts["queue"] and summary["split"] == counts["splits"], "summary resulting frontier counters mismatch")
    require(summary["cache_canonical_sha256"] == source["cache_canonical_sha256"] and summary["cache_entries"] == source["cache_entries"], "cache reference differs from immutable source")
    require(summary["kernel_calls"] == 0 and summary["qualification_seconds"] <= 120, "summary exceeds zero-kernel finite limits")
    require(checkpoint["stage"] == "complete" and checkpoint["block"] == summary["block"] == block_reg["block"] == 1, "wrong block stage/identity")
    require(checkpoint["pid"] == summary["pid"], "author PID receipt mismatch")
    require(summary["points"] == 81 * len(rows), "summary point count mismatch")
    return {"source_tree": initial_tree, "first_block_tree": new_tree, "new_counts": counts,
            "retired_counts": removed_counts, "accepted": len(accepted), "failed_uncertified": len(rows) - len(accepted),
            "initial_logical_frontier_sha256": initial_hash, "new_logical_frontier_sha256": new_hash,
            "source_cache_reference_verified": True, "active_source_index_counts": {k: len(v) for k, v in active.items()},
            "logical_retained_source_fields": "inherited by immutable source index; final full-state materialization not read or certified here"}


def write_json(path, obj):
    temp = Path(str(path) + ".tmp")
    temp.write_text(json.dumps(obj, sort_keys=True, indent=2) + "\n")
    temp.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="store_true")
    args = parser.parse_args()
    reg = json.loads((HERE / "registration.json").read_text())
    if not args.run:
        print(json.dumps({"task": reg["task"], "block": reg["block"], "requires_run": True}))
        return
    started = time.monotonic()
    pid = os.getpid()
    print(json.dumps({"pid": pid, "stage": "immutable_input_locks"}), flush=True)
    block = Path(reg["block"])
    inputs = {str(block / k): v for k, v in reg["block_locks"].items()}
    inputs.update({reg["source"]: reg["source_sha256"], str(block.parent / "registration.json"): reg["run_registration_sha256"], reg["author_code"]: reg["author_code_sha256"]})
    protected = {str(HERE.parent / "audit_transition062" / k): v for k, v in reg["protected_062_locks"].items()}
    for path, expected in {**inputs, **protected}.items():
        require(file_sha(path) == expected, "input/protected062 hash mismatch: " + path)
    documents = {name: json.loads((block / name).read_text(), object_pairs_hook=no_duplicates) for name in reg["block_locks"]}
    block_reg, checkpoint, summary = (documents[k] for k in ("registration.json", "checkpoint.json", "summary.json"))
    run_reg = json.loads((block.parent / "registration.json").read_text())
    for name in ("certificates.json", "retired.json", "summary.json", "registration.json"):
        require(checkpoint[name.removesuffix(".json") + "_sha256"] == reg["block_locks"][name], "checkpoint artifact link mismatch: " + name)
    require(block_reg["run_registration_sha256"] == reg["run_registration_sha256"] and block_reg["code_sha256"] == reg["author_code_sha256"], "block code/run registration link mismatch")
    require(block_reg["source_state"] == run_reg["source"] == reg["source"] and block_reg["source_state_sha256"] == run_reg["source_sha256"] == reg["source_sha256"], "source reference mismatch")
    require(block_reg["max_tested"] == 2048 and block_reg["qualification_seconds"] == 120 and run_reg["new_kernel_chain_authorized"] is False, "registered finite zero-kernel scope mismatch")
    rows = documents["certificates.json"]
    require(0 < len(rows) <= 2048, "row count exceeds finite parent cap")
    print(json.dumps({"pid": pid, "stage": "source_stream_and_logical_reconstruction"}), flush=True)
    source_started = time.monotonic()
    source = read_source(reg["source"], rows, reg)
    selection = audit_selection(source, block_reg, rows)
    delta = audit_delta(source, rows, documents["retired.json"], reg, block_reg, checkpoint, summary)
    source_seconds = time.monotonic() - source_started
    result = {"task": "audit_conversion069", "scope": reg["scope"], "pid": pid, "status": "checking_certificates", "inputs": inputs,
              "code_sha256": file_sha(Path(__file__)), "registration_sha256": file_sha(HERE / "registration.json"),
              "source_and_tree_seconds": source_seconds, "source_metadata": source["metadata"],
              "source_cache_canonical_sha256": source["cache_canonical_sha256"], "source_cache_entries": source["cache_entries"],
              "selection": selection, "first_block_delta": delta, "certificate_results": [],
              "zero_new_H_gradient_kernel_Green_log_exp_Gamma_Arb_calls": True, "MATLAB_calls": 0,
              "author_functions_imported": False, "protected062_locks": protected, "limits": reg["limits"]}
    print(json.dumps({"pid": pid, "stage": "all_2048_integer_81_point_certificates"}), flush=True)
    loop_started = time.monotonic()
    for row in rows:
        if time.monotonic() - loop_started >= 120:
            break
        source_split = source["selected_splits"][row["source_split_index"]]
        item = verify_certificate(row, source_split, source["cache"], source["roots"])
        result["certificate_results"].append(item)
        if len(result["certificate_results"]) % 128 == 0:
            write_json(HERE / "result.json", result)
    result["certificate_loop_seconds"] = time.monotonic() - loop_started
    result["tests_verified"] = len(result["certificate_results"])
    result["points_verified"] = sum(r["checked_points"] for r in result["certificate_results"])
    result["accepted_verified"] = sum(r["passed"] for r in result["certificate_results"])
    result["uncertified_verified"] = result["tests_verified"] - result["accepted_verified"]
    positives = [F(r["minimum_log_margin"]) for r in result["certificate_results"] if r["passed"]]
    negatives = [F(r["minimum_log_margin"]) for r in result["certificate_results"] if not r["passed"]]
    result["minimum_accepted_exact_margin"] = str(min(positives)) if positives else None
    result["maximum_uncertified_exact_margin"] = str(max(negatives)) if negatives else None
    result["all_anchors_strictly_inside"] = all(r["anchor_strictly_inside"] for r in result["certificate_results"])
    result["status"] = "passed_first_block_only" if result["tests_verified"] == len(rows) else "certificate_prefix_only"
    if result["tests_verified"] == len(rows):
        require(result["accepted_verified"] == delta["accepted"], "independent accepted count mismatch")
    for path, expected in {**inputs, **protected}.items():
        require(file_sha(path) == expected, "input or protected062 changed during audit: " + path)
    result["input_and_protected062_hashes_unchanged_after_run"] = True
    result["wall_seconds"] = time.monotonic() - started
    result["peak_rss_bytes"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    write_json(HERE / "result.json", result)
    print(json.dumps({k: result[k] for k in ("pid", "status", "tests_verified", "points_verified", "accepted_verified", "uncertified_verified", "source_and_tree_seconds", "certificate_loop_seconds", "wall_seconds", "peak_rss_bytes", "minimum_accepted_exact_margin")}), flush=True)


if __name__ == "__main__":
    main()


