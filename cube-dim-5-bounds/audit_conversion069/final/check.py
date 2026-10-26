#!/usr/bin/env python3
"""Independent exact audit of all frozen 066 blocks and final materialization.

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
            check_budget()
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
        check_budget()
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
            elif field in ("epsilon", "target", "model", "active", "inflight", "resume_jobs", "fingerprint", "fingerprint061", "status061", "total_kernel_attempts_completed", "total_new_anchor_starts"):
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


def write_json(path, obj):
    temp = Path(str(path) + ".tmp")
    temp.write_text(json.dumps(obj, sort_keys=True, indent=2) + "\n")
    temp.replace(path)


RUN_STARTED = None
RECORD = None


class BudgetExpired(Exception):
    pass


def check_budget():
    if RUN_STARTED is not None and time.monotonic() - RUN_STARTED >= 160:
        raise BudgetExpired("checkpoint before180-second total wall limit")


def checkpoint():
    if RECORD is not None:
        RECORD["wall_seconds_at_checkpoint"] = time.monotonic() - RUN_STARTED
        write_json(HERE / "result.json", RECORD)


def final_object_comparison(source_path, final_path, active, accepted, expected_conversion):
    old, new = Stream(source_path), Stream(final_path)
    ak, bk = old.key(), new.key()
    additions = {"task066", "status066", "fingerprint066", "conversion066"}
    extra, direct = {}, {"cache": 0, "passed": 0, "queue": 0, "splits": 0, "old_top_fields": 0, "new_passed": 0}
    missing = object()
    while ak is not None or bk is not None:
        check_budget()
        if bk in additions and (ak is None or bk < ak):
            extra[bk] = new.value(); bk = new.key(); continue
        require(ak == bk and ak is not None, "final state has an unapproved top-level key change")
        field = ak
        if field == "cache":
            for a, b in zip_longest(old.items(True), new.items(True), fillvalue=missing):
                require(a is not missing and b is not missing and a == b, "final cache entry changed, added or missing")
                direct["cache"] += 1
        elif field in active:
            target = new.items()
            for index, record in enumerate(old.items()):
                if index not in active[field]:
                    continue
                other = next(target, missing)
                require(other is not missing and record == other, "retained final record/order changed: " + field)
                direct[field] += 1
            if field == "passed":
                for expected in accepted:
                    require(next(target, missing) == expected, "new final passed record differs from replayed certificate")
                    direct["new_passed"] += 1
            require(next(target, missing) is missing, "extra final array record: " + field)
        else:
            require(old.value() == new.value(), "unchanged source top-level value changed: " + field)
            direct["old_top_fields"] += 1
        ak, bk = old.key(), new.key()
    old.close(); new.close()
    require(set(extra) == additions, "final state additions differ from authorized schema")
    require(extra["task066"] == "upper_cover066" and extra["status066"] == "zero_kernel_conversion_complete_waiting_main_review", "wrong final task/status")
    require(extra["conversion066"] == expected_conversion, "final conversion metadata differs from independently replayed chain")
    expected_fingerprint = digest({"registration": expected_conversion["run_registration_sha256"], "last_checkpoint": expected_conversion["last_checkpoint_sha256"]})
    require(extra["fingerprint066"] == expected_fingerprint, "final fingerprint formula mismatch")
    return {"direct_object_comparisons": direct, "final_new_fields": extra, "all_old_top_fields_preserved": True}


def replay_chain(source, blocks, reg, first_receipt, old_progress=None):
    arrays = source["arrays"]
    active = {k: set(range(len(v))) for k, v in arrays.items()}
    accepted, tested = [], set()
    accepted_by_root = defaultdict(list)
    qcounts = defaultdict(int)
    for r, p in arrays["queue"]:
        for depth in range(len(p)):
            qcounts[(r, p[:depth])] += 1
    order = sorted((i for i, n in enumerate(arrays["splits"]) if qcounts[n] > 0),
                   key=lambda i: (len(arrays["splits"][i][1]), -qcounts[arrays["splits"][i]], *arrays["splits"][i]))
    initial = validate_tree(source["roots"], arrays["passed"] + arrays["queue"], set(arrays["splits"]), "source061")
    prevcp = None
    prevfrontier = logical_digest(reg["source_sha256"], {k: sorted(v) for k, v in active.items()}, [], [])
    summaries, margin_min = [], None
    old_completed = {x["block"]: x for x in (old_progress or {}).get("block_results", [])}
    for bundle in blocks:
        check_budget()
        number, br, cp, rows, retired, summary = (bundle[k] for k in ("number", "registration", "checkpoint", "rows", "retired", "summary"))
        print(json.dumps({"pid": os.getpid(), "stage": "block", "block": number, "tested_rows": len(rows)}), flush=True)
        require(br["block"] == cp["block"] == summary["block"] == number, "block identity mismatch")
        require(br["parent_checkpoint_sha256"] == cp["parent_checkpoint_sha256"] == prevcp, "checkpoint parent chain mismatch")
        require(br["parent_frontier_sha256"] == prevfrontier, "parent logical frontier mismatch")
        current_q = defaultdict(int)
        for index in active["queue"]:
            r, p = arrays["queue"][index]
            for depth in range(len(p)):
                current_q[(r, p[:depth])] += 1
        eligible = [i for i in order if i in active["splits"] and i not in tested]
        require(all(current_q[arrays["splits"][i]] > 0 for i in eligible), "live untested candidate has no remaining queue descendant")
        selected_indices = eligible[:2048]
        selected = [{"source_split_index": i, "node": list(arrays["splits"][i]), "initial_queue_descendants": qcounts[arrays["splits"][i]]} for i in selected_indices]
        require(br["selected"] == selected, "registered block prefix differs from independent selection")
        require(summary["selected"] == len(selected) and summary["cursor"] == len(selected), "block did not finish its registered prefix")
        require(0 < len(rows) <= len(selected) <= 2048, "block test/selection cap exceeded")
        local_parents = set()
        row_cursor = 0
        for chosen in selected_indices:
            n = arrays["splits"][chosen]
            if covering_parent(n, local_parents, True) is not None:
                continue
            require(row_cursor < len(rows), "a non-skipped selected row is missing")
            row = rows[row_cursor]; row_cursor += 1
            require(row["source_split_index"] == chosen and node(row["node"]) == n, "row differs from registered prefix with accepted-ancestor skips")
            require(chosen in active["splits"] and chosen not in tested, "tested index is inactive or repeated")
            tested.add(chosen)
            if row["certificate"]["passed"]:
                require(all(not (n[1].startswith(p) or p.startswith(n[1])) for p in accepted_by_root[n[0]]), "accepted parent overlaps an earlier accepted parent")
                accepted_by_root[n[0]].append(n[1]); local_parents.add(n)
                accepted.append({"node": row["node"], "certificate": row["certificate"]})
        require(row_cursor == len(rows), "extra unselected certificate rows")
        exact_results = []
        loop_start = time.monotonic()
        reused = False
        if number == 1:
            require(first_receipt["status"] == "passed_first_block_only", "first-block receipt is not passed")
            exact_results = first_receipt["certificate_results"]
            require(len(exact_results) == len(rows), "first-block result count mismatch")
            reused = True
        elif number in old_completed:
            previous = old_completed[number]
            require(previous["certificate_input_sha256"] == bundle["hashes"]["certificates.json"], "resume certificate artifact changed")
            require(previous["status"] == "verified", "resume block was not fully verified")
            reused = True
        else:
            partial = (old_progress or {}).get("partial_block", {})
            if partial.get("block") == number:
                require(partial["certificate_input_sha256"] == bundle["hashes"]["certificates.json"], "resume partial input changed")
                exact_results = partial["results"]
            for row in rows[len(exact_results):]:
                check_budget()
                require(time.monotonic() - loop_start < 120, "per-block120-second certificate budget exceeded")
                item = verify_certificate(row, source["selected_splits"][row["source_split_index"]], source["cache"], source["roots"])
                exact_results.append(item)
                if len(exact_results) % 128 == 0:
                    RECORD["partial_block"] = {"block": number, "certificate_input_sha256": bundle["hashes"]["certificates.json"], "results": exact_results}
                    checkpoint()
        loop_seconds = time.monotonic() - loop_start
        if number in old_completed and number != 1:
            checked = dict(old_completed[number])
            checked["reused_verified_receipt"] = True
            checked["this_run_certificate_loop_seconds"] = loop_seconds
            require(checked["tested"] == len(rows), "resume tested count mismatch")
        else:
            require(len(exact_results) == len(rows), "incomplete exact results")
            for row, item in zip(rows, exact_results):
                require(row["node"] == item["node"] and row["source_split_index"] == item["source_split_index"], "verified result row identity mismatch")
                require(row["certificate"]["passed"] == item["passed"] and F(row["certificate"]["minimum_log_margin"]) == F(item["minimum_log_margin"]), "verified receipt verdict/margin changed")
                require(item["checked_points"] == 81, "verified receipt has wrong point count")
            positive = [F(x["minimum_log_margin"]) for x in exact_results if x["passed"]]
            checked = {"block": number, "status": "verified", "certificate_input_sha256": bundle["hashes"]["certificates.json"],
                       "tested": len(rows), "accepted": len(positive), "uncertified": len(rows) - len(positive),
                       "points": 81 * len(rows), "results_sha256": digest(exact_results),
                       "minimum_accepted_exact_margin": str(min(positive)) if positive else None,
                       "certificate_loop_seconds": loop_seconds, "reused_first_block_receipt": number == 1}
        require(checked["accepted"] == len(local_parents), "independently verified accepted count mismatch")
        if checked["minimum_accepted_exact_margin"] is not None:
            value = F(checked["minimum_accepted_exact_margin"])
            require(value > 0, "accepted margin is not strict")
            margin_min = value if margin_min is None else min(margin_min, value)
        retired_counts = {}
        for field, nodes in arrays.items():
            expected_removed, remove = [], set()
            for index in sorted(active[field]):
                parent = covering_parent(nodes[index], local_parents, field == "splits")
                if parent is not None:
                    expected_removed.append({"source_index": index, "node": list(nodes[index]), "unique_new_parent": list(parent)})
                    remove.add(index)
            require(retired[field] == expected_removed, "per-block retirement differs from exact covered active source records: " + field)
            active[field].difference_update(remove); retired_counts[field] = len(remove)
        frontier = [arrays["passed"][i] for i in sorted(active["passed"])] + [node(r["node"]) for r in accepted] + [arrays["queue"][i] for i in sorted(active["queue"])]
        tree = validate_tree(source["roots"], frontier, {arrays["splits"][i] for i in active["splits"]}, "block" + str(number))
        current_hash = logical_digest(reg["source_sha256"], {k: sorted(v) for k, v in active.items()}, accepted, tested)
        require(current_hash == cp["frontier_sha256"] == summary["frontier_sha256"], "logical post-block frontier hash mismatch")
        require(cp["cumulative_tested"] == summary["cumulative_tested"] == len(tested) and cp["cumulative_accepted"] == summary["cumulative_accepted"] == len(accepted), "cumulative counter mismatch")
        require(summary["tested"] == len(rows) and summary["accepted"] == len(local_parents) and summary["failed_uncertified"] == len(rows) - len(local_parents), "block verdict counters mismatch")
        require(summary["retired"] == retired_counts and summary["passed"] == len(active["passed"]) + len(accepted) and summary["queue"] == len(active["queue"]) and summary["split"] == len(active["splits"]), "block frontier counters mismatch")
        require(summary["points"] == len(rows) * 81 and summary["cache_entries"] == source["cache_entries"] and summary["cache_canonical_sha256"] == source["cache_canonical_sha256"], "block points/cache reference mismatch")
        require(summary["kernel_calls"] == 0 and summary["qualification_seconds"] <= 120 and cp["stage"] == "complete", "block finite zero-kernel qualification receipt mismatch")
        checked.update({"selected": len(selected), "skipped_covered_selected": len(selected) - len(rows), "retired_counts": retired_counts,
                        "tree": tree, "logical_frontier_sha256": current_hash, "checkpoint_sha256": bundle["hashes"]["checkpoint.json"],
                        "new_counts": {"passed": len(active["passed"]) + len(accepted), "splits": len(active["splits"]), "queue": len(active["queue"])}})
        summaries.append(checked)
        prevfrontier, prevcp = current_hash, bundle["hashes"]["checkpoint.json"]
        RECORD["block_results"] = summaries
        RECORD.pop("partial_block", None)
        RECORD["stage"] = "verified_through_block_" + str(number)
        checkpoint()
    remaining = [i for i in order if i in active["splits"] and i not in tested]
    require(not remaining, "chain stopped with untested eligible source splits")
    return {"active": active, "accepted": accepted, "tested": tested, "source_tree": initial,
            "last_checkpoint_sha256": prevcp, "final_logical_frontier_sha256": prevfrontier,
            "minimum_accepted_exact_margin": str(margin_min), "remaining_untested_eligible": len(remaining)}


def main():
    global RUN_STARTED, RECORD
    import signal
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--resume-json")
    parser.add_argument("--resume-sha256")
    args = parser.parse_args()
    reg = json.loads((HERE / "registration.json").read_text())
    if not args.run:
        print(json.dumps({"task": reg["task"], "total_wall_cap_seconds": 180, "blocks": 17, "requires_run": True})); return
    RUN_STARTED = time.monotonic()
    signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(BudgetExpired("174-second emergency checkpoint")))
    signal.setitimer(signal.ITIMER_REAL, 174)
    RECORD = {"task": reg["task"], "pid": os.getpid(), "status": "verified_prefix_only", "stage": "input_locks",
              "code_sha256": file_sha(Path(__file__)), "registration_sha256": file_sha(HERE / "registration.json"),
              "manifest_sha256": reg["manifest_sha256"], "block_results": [], "program_total_wall_cap_seconds": 180,
              "first_block_81_point_tests_recomputed": False, "zero_new_H_gradient_kernel_Green_log_exp_Gamma_Arb_calls": True,
              "MATLAB_calls": 0, "author_functions_imported": False, "limits": reg["limits"]}
    print(json.dumps({"pid": os.getpid(), "stage": "input_locks", "total_wall_cap_seconds": 180}), flush=True)
    try:
        require(file_sha(reg["manifest"]) == reg["manifest_sha256"], "final manifest hash mismatch")
        manifest = json.loads(Path(reg["manifest"]).read_text())
        inputs = {reg["manifest"]: reg["manifest_sha256"], reg["source"]: reg["source_sha256"], reg["final_state"]: reg["final_state_sha256"], reg["author_code"]: reg["author_code_sha256"]}
        inputs.update(manifest["external_files"])
        require(inputs[reg["final_state"]] == reg["final_state_sha256"] and manifest["source_sha256"] == reg["source_sha256"] and manifest["code_sha256"] == reg["author_code_sha256"], "manifest source/final/code links mismatch")
        protected = {str(HERE.parent / k): v for k, v in reg["protected_first_block069"].items()}
        protected.update({str(HERE.parent.parent / "audit_transition062" / k): v for k, v in reg["protected_062"].items()})
        bundles, all_rows = [], []
        for index, item in enumerate(manifest["blocks"], 1):
            require(item["block"] == index, "manifest block order mismatch")
            directory = Path(item["directory"])
            hashes = {name + ".json": item[name + "_sha256"] for name in ("registration", "certificates", "retired", "checkpoint", "summary")}
            for name, value in hashes.items(): inputs[str(directory / name)] = value
        for path, expected in {**inputs, **protected}.items():
            check_budget(); require(file_sha(path) == expected, "frozen input/protected file hash mismatch: " + path)
        require(len(manifest["blocks"]) == 17, "wrong frozen block count")
        RECORD["inputs"] = inputs; RECORD["protected_inputs"] = protected
        for index, item in enumerate(manifest["blocks"], 1):
            directory = Path(item["directory"])
            docs = {name: json.loads((directory / (name + ".json")).read_text(), object_pairs_hook=no_duplicates) for name in ("registration", "certificates", "retired", "checkpoint", "summary")}
            hashes = {name + ".json": item[name + "_sha256"] for name in docs}
            cp, br = docs["checkpoint"], docs["registration"]
            for name in ("registration", "certificates", "retired", "summary"):
                require(cp[name + "_sha256"] == hashes[name + ".json"], "checkpoint artifact link mismatch")
            require(cp == {k: v for k, v in item.items() if k not in ("directory", "checkpoint_sha256")}, "manifest checkpoint fields mismatch")
            require(br["source_state"] == reg["source"] and br["source_state_sha256"] == reg["source_sha256"] and br["run_registration_sha256"] == reg["run_registration_sha256"] and br["code_sha256"] == reg["author_code_sha256"], "block frozen source/code/registration link mismatch")
            require(br["max_tested"] == 2048 and br["qualification_seconds"] == 120, "block finite limits mismatch")
            bundle = {"number": index, "registration": br, "checkpoint": cp, "rows": docs["certificates"], "retired": docs["retired"], "summary": docs["summary"], "hashes": hashes}
            bundles.append(bundle); all_rows.extend(bundle["rows"])
        first_receipt = json.loads((HERE.parent / "result.json").read_text())
        require(first_receipt["inputs"][reg["source"]] == reg["source_sha256"], "first receipt source link mismatch")
        prior = None
        if args.resume_json:
            require(args.resume_sha256 and file_sha(args.resume_json) == args.resume_sha256, "resume receipt hash missing or mismatched")
            prior = json.loads(Path(args.resume_json).read_text())
            require(prior["code_sha256"] == RECORD["code_sha256"] and prior["registration_sha256"] == RECORD["registration_sha256"] and prior["inputs"] == inputs, "resume code/registration/input identity changed")
        RECORD["stage"] = "stream_source"; checkpoint()
        print(json.dumps({"pid": os.getpid(), "stage": "stream_source", "total_rows": len(all_rows)}), flush=True)
        source = read_source(reg["source"], all_rows, reg)
        RECORD["source_cache_entries"] = source["cache_entries"]; RECORD["source_cache_canonical_sha256"] = source["cache_canonical_sha256"]
        RECORD["stage"] = "replay_chain"
        chain = replay_chain(source, bundles, reg, first_receipt, prior)
        RECORD["stage"] = "final_direct_object_comparison"; checkpoint()
        print(json.dumps({"pid": os.getpid(), "stage": "final_direct_object_comparison"}), flush=True)
        expected_conversion = {"source_state": reg["source"], "source_state_sha256": reg["source_sha256"], "source_fingerprint": source["metadata"]["fingerprint"],
            "run_registration_sha256": reg["run_registration_sha256"], "last_block": 17, "last_checkpoint_sha256": chain["last_checkpoint_sha256"],
            "frontier_sha256": chain["final_logical_frontier_sha256"], "tested": len(chain["tested"]), "accepted": len(chain["accepted"]),
            "cache_canonical_sha256": source["cache_canonical_sha256"], "new_kernel_chain_started": False}
        materialization = final_object_comparison(reg["source"], reg["final_state"], chain["active"], chain["accepted"], expected_conversion)
        require(materialization["direct_object_comparisons"]["cache"] == source["cache_entries"], "final exact cache comparison count mismatch")
        final_summary = json.loads(Path(reg["final_state"]).with_name("run_receipt.json").read_text())
        require(final_summary["state_sha256"] == reg["final_state_sha256"] and final_summary["source_sha256_after"] == reg["source_sha256"] and final_summary["last_checkpoint_sha256"] == chain["last_checkpoint_sha256"], "final receipt source/state/chain links mismatch")
        require(final_summary["tested"] == len(chain["tested"]) and final_summary["accepted"] == len(chain["accepted"]) and final_summary["failed_uncertified"] == len(chain["tested"]) - len(chain["accepted"]), "final receipt verdict counts mismatch")
        require(final_summary["remaining_untested_eligible"] == 0 and final_summary["kernel_calls"] == 0 and final_summary["global_g_or_C5_claim"] is False, "final receipt scope/exhaustion mismatch")
        final_counts = RECORD["block_results"][-1]["new_counts"]
        require(final_summary["passed"] == final_counts["passed"] and final_summary["split"] == final_counts["splits"] and final_summary["queue"] == final_counts["queue"], "final receipt frontier counts mismatch")
        RECORD["materialization"] = materialization
        RECORD["source_tree"] = chain["source_tree"]
        RECORD["tested"] = len(chain["tested"]); RECORD["accepted"] = len(chain["accepted"])
        RECORD["uncertified"] = RECORD["tested"] - RECORD["accepted"]
        RECORD["total_points_verified_including_first_receipt"] = 81 * RECORD["tested"]
        RECORD["this_run_new_certificate_tests"] = sum(x["tested"] for x in RECORD["block_results"] if not x.get("reused_first_block_receipt") and not x.get("reused_verified_receipt"))
        RECORD["first_block_reused_tests"] = len(first_receipt["certificate_results"])
        RECORD["final_counts"] = RECORD["block_results"][-1]["new_counts"]
        RECORD["minimum_accepted_exact_margin"] = chain["minimum_accepted_exact_margin"]
        RECORD["remaining_untested_eligible"] = 0; RECORD["last_checkpoint_sha256"] = chain["last_checkpoint_sha256"]
        RECORD["final_logical_frontier_sha256"] = chain["final_logical_frontier_sha256"]
        RECORD["stage"] = "post_hash_locks"; checkpoint()
        for path, expected in {**inputs, **protected}.items():
            check_budget(); require(file_sha(path) == expected, "input/protected file changed during audit: " + path)
        RECORD["input_and_protected_locks_unchanged_after_run"] = True
        RECORD["status"] = "passed_all17_blocks_and_final_materialization"
        RECORD["stage"] = "complete"
    except BudgetExpired as exc:
        signal.setitimer(signal.ITIMER_REAL, 0)
        RECORD["status"] = "verified_prefix_only"; RECORD["deadline_reason"] = str(exc)
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        RECORD["wall_seconds"] = time.monotonic() - RUN_STARTED
        RECORD["peak_rss_bytes"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        checkpoint()
    print(json.dumps({k: RECORD.get(k) for k in ("pid", "status", "stage", "tested", "accepted", "uncertified", "final_counts", "this_run_new_certificate_tests", "wall_seconds", "peak_rss_bytes")}), flush=True)
    print(json.dumps({"result_sha256": sha256((HERE / "result.json").read_bytes()).hexdigest(), "completed_blocks": len(RECORD["block_results"])}), flush=True)


if __name__ == "__main__":
    main()


