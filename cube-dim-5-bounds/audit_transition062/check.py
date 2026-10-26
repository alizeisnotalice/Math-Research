#!/usr/bin/env python3
"""Independent exact audit of the frozen 052 -> 061 finite conversion.

Uses Python standard library only. No author module is imported. The simultaneous
JSON streams compare unchanged objects directly, including all cache records.
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
NEW_FIELDS = {"task061", "fingerprint061", "status061", "conversion061"}
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
    require(cert["passed"] is True and cert["mechanism"] == "signed_basis_81", "wrong new certificate mechanism")
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
    require(hl > F(1, 10000), "cached H lower endpoint not above target")
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
    require(count == 81 and best > 0, "81-point exact margin is not strictly positive")
    margin = F(best, denominator)
    require(F(cert["minimum_log_margin"]) == margin, "new minimum margin mismatch")
    worst = tuple(F(x) * vden for x in cert["worst_offset"])
    require(len(worst) == 4 and all(x.denominator == 1 for x in worst), "worst offset outside integer grid")
    require(tuple(map(int, worst)) in minima, "declared worst offset not a minimizer")
    return {"node": list(n), "minimum_log_margin": str(margin), "checked_points": count,
            "gradient_endpoint_cross_checks": cross_checks, "anchor_closed_box": True,
            "anchor_strictly_inside": all(a < p < b for a, p, b in zip(lo, point[1:], hi)),
            "physical_box_lower": list(map(str, lo)), "physical_box_upper": list(map(str, hi))}


def audit_streams(old_path, new_path, certificates, retired):
    parents = {node(row["node"]) for row in certificates}
    require(len(parents) == len(certificates) <= 2048 and len(parents) > 0, "bad new parent count")
    for n in parents:
        require(covering_parent(n, parents, False) is None, "new parents not prefix-free")
    anchors = {r["certificate"]["anchor"] for r in certificates}
    cache = {}
    old_frontier, new_frontier, old_splits, new_splits = [], [], set(), set()
    old_split_sources = {}
    new_fields, old_metadata = {}, {}
    counts = {"old": {}, "new": {}, "retired": {}}
    direct = {"cache": 0, "passed": 0, "splits": 0, "queue": 0, "old_top_fields": 0}
    a, b = Stream(old_path), Stream(new_path)
    ak, bk = a.key(), b.key()
    missing = object()
    roots = None
    cache_digest = sha256()
    while ak is not None or bk is not None:
        if bk in NEW_FIELDS and (ak is None or bk < ak):
            new_fields[bk] = b.value()
            bk = b.key()
            continue
        require(ak == bk and ak is not None, "unapproved top-level field added/removed")
        field = ak
        if field == "cache":
            cache_digest.update(b"{")
            first = True
            for x, y in zip_longest(a.items(True), b.items(True), fillvalue=missing):
                require(x is not missing and y is not missing and x == y, "cache entry changed or missing")
                k, record = x
                if not first:
                    cache_digest.update(b",")
                first = False
                cache_digest.update(canonical(k)); cache_digest.update(b":"); cache_digest.update(canonical(record))
                if k in anchors:
                    cache[k] = record
                direct["cache"] += 1
            cache_digest.update(b"}")
            require(set(cache) == anchors, "selected certificate anchor absent from unchanged cache")
        elif field in ("passed", "splits", "queue"):
            ait, bit = a.items(), b.items()
            expected_retired = []
            old_count, new_count = 0, 0
            for index, record in enumerate(ait):
                n = node(record if field == "queue" else record["node"])
                old_count += 1
                if field == "splits":
                    require(n not in old_splits, "duplicate old split")
                    old_splits.add(n)
                    require(record["axis"] == len(n[1]) % 4, "split axis differs from interleaved path geometry")
                    require(record["children"] == [[n[0], n[1] + "0"], [n[0], n[1] + "1"]], "split children mismatch")
                    if n in parents:
                        old_split_sources[n] = {"index": index, **record}
                else:
                    old_frontier.append(n)
                parent = covering_parent(n, parents, field == "splits")
                if parent is not None:
                    expected_retired.append({"source_index": index, "node": list(n), "unique_new_parent": list(parent)})
                    continue
                other = next(bit, missing)
                require(other is not missing and record == other, "surviving " + field + " record or order changed")
                direct[field] += 1
                new_count += 1
                if field == "splits":
                    new_splits.add(n)
                else:
                    new_frontier.append(n)
            promoted = []
            for record in bit:
                require(field == "passed", "new queue/split record unexpectedly created")
                n = node(record["node"])
                promoted.append(n)
                require(set(record) == {"node", "certificate"}, "unexpected new passed record fields")
                expected = next((r["certificate"] for r in certificates if node(r["node"]) == n), None)
                require(expected is not None and record["certificate"] == expected, "promoted certificate differs between state and ledger")
                new_frontier.append(n)
                new_count += 1
            if field == "passed":
                require(len(promoted) == len(parents) and set(promoted) == parents, "new passed parents differ from certificate ledger")
            require(retired[field] == expected_retired, "retired " + field + " ledger does not exactly enumerate removed records")
            counts["old"][field], counts["new"][field] = old_count, new_count
            counts["retired"][field] = len(expected_retired)
        else:
            x, y = a.value(), b.value()
            require(x == y, "old top-level value changed: " + field)
            direct["old_top_fields"] += 1
            if field == "roots":
                roots = x
            elif field in ("active", "inflight", "resume_jobs", "epsilon", "target", "model", "fingerprint", "pool_evidence", "status", "total_kernel_attempts_completed", "total_new_anchor_starts"):
                old_metadata[field] = x
        ak, bk = a.key(), b.key()
    a.close(); b.close()
    require(set(new_fields) == NEW_FIELDS, "wrong added top-level fields")
    require(roots == [[i, j, k, k2] for i in range(32) for j in range(i + 1) for k in range(j + 1) for k2 in range(k + 1)], "root list differs from independently enumerated coarse chamber boxes")
    require(old_metadata["active"] is None and old_metadata["inflight"] == {} and old_metadata["resume_jobs"] == {}, "source or destination has live/resumable jobs")
    require(old_metadata["epsilon"] == "1/64" and old_metadata["target"] == "1/10000" and old_metadata["model"] == "alpha1", "wrong fixed model parameters")
    require(set(old_split_sources) == parents, "a promoted parent was not an old active split")
    pool = old_metadata["pool_evidence"]
    require(pool["close_called"] is True and pool["join_returned"] is True and len(pool["exits"]) == 8, "old pool closure evidence incomplete")
    require(all(x["exitcode"] == 0 and x["is_alive"] is False for x in pool["exits"]), "old pool workers not all cleanly exited")
    require(sorted(pool["worker_pids"]) == sorted(x["pid"] for x in pool["exits"]), "pool worker receipt identities mismatch")
    for row in certificates:
        require(row["source_split_index"] == old_split_sources[node(row["node"])]["index"], "certificate source split index mismatch")
    return {"counts": counts, "direct_object_comparisons": direct,
            "cache_canonical_sha256": cache_digest.hexdigest(),
            "old_tree": validate_tree(roots, old_frontier, old_splits, "old"),
            "new_tree": validate_tree(roots, new_frontier, new_splits, "new"),
            "new_top_fields": new_fields, "old_metadata": old_metadata,
            "roots": roots, "cache": cache, "old_split_sources": old_split_sources}


def write_json(path, obj):
    temp = Path(str(path) + ".tmp")
    temp.write_text(json.dumps(obj, sort_keys=True, indent=2) + "\n")
    temp.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    reg = json.loads((HERE / "registration.json").read_text())
    if not args.run:
        print(json.dumps({"registered": True, "new_state": reg["new_state"], "requires_run": True}))
        return
    started = time.monotonic()
    pid = os.getpid()
    print(json.dumps({"pid": pid, "stage": "immutable_input_hashes", "registered_parent_cap": 2048}), flush=True)
    old_path, new_path = Path(reg["old_state"]), Path(reg["new_state"])
    batch = new_path.parent
    inputs = {str(old_path): reg["old_state_sha256"], str(new_path): reg["new_state_sha256"],
              str(batch / "certificates.json"): reg["certificates_sha256"],
              str(batch / "retired.json"): reg["retired_sha256"],
              str(batch / "registration.json"): reg["author_registration_sha256"]}
    manifest = batch / "manifest.json"
    if manifest.exists():
        inputs[str(manifest)] = reg["author_manifest_sha256"]
    for path, expected in inputs.items():
        require(file_sha(path) == expected, "frozen input SHA mismatch: " + path)
    certificates = json.loads((batch / "certificates.json").read_text(), object_pairs_hook=no_duplicates)
    retired = json.loads((batch / "retired.json").read_text(), object_pairs_hook=no_duplicates)
    print(json.dumps({"pid": pid, "stage": "paired_stream_exact_object_comparison"}), flush=True)
    tree_started = time.monotonic()
    tree = audit_streams(old_path, new_path, certificates, retired)
    tree_seconds = time.monotonic() - tree_started
    require(tree["counts"]["old"] == {k: reg["old_counts"][k] for k in ("passed", "splits", "queue")}, "registered old counts mismatch")
    require(tree["direct_object_comparisons"]["cache"] == reg["old_counts"]["cache"], "registered cache count mismatch")
    new_fields = tree["new_top_fields"]
    require(new_fields["task061"] == "upper_cover061", "new task identity mismatch")
    require(new_fields["status061"] == "finite_zero_kernel_conversion_complete", "new conversion status mismatch")
    conversion = new_fields["conversion061"]
    require(conversion["source_state_sha256"] == reg["old_state_sha256"] and conversion["source_state"] == reg["old_state"], "new conversion source-state link mismatch")
    require(conversion["source_fingerprint"] == tree["old_metadata"]["fingerprint"], "new conversion source-fingerprint link mismatch")
    for key in ("certificates_sha256", "retired_sha256"):
        require(conversion[key] == reg[key], "new conversion artifact link mismatch: " + key)
    require(conversion["registration_sha256"] == reg["author_registration_sha256"], "new conversion registration link mismatch")
    require(conversion["cache_sha256"] == tree["cache_canonical_sha256"], "new conversion cache digest mismatch")
    require(conversion["accepted_parents"] == conversion["tested_parents"] == conversion["cursor"] == len(certificates) <= 2048, "new conversion parent counters mismatch")
    require(conversion["new_kernel_chain_started"] is False, "new conversion claims a started kernel chain")
    author_registration = json.loads((batch / "registration.json").read_text())
    require(author_registration["max_tested_split_parents"] == 2048 and author_registration["seconds"] == 120, "author finite registration cap mismatch")
    require(author_registration["new_kernel_chain_authorized"] is False, "author registration unexpectedly authorizes a kernel chain")
    require(author_registration["parent_state"] == reg["old_state"] and author_registration["parent_state_sha256"] == reg["old_state_sha256"], "author registered parent link mismatch")
    result = {"task": "audit_transition062", "pid": pid, "status": "certificate_prefix",
              "inputs": inputs, "code_sha256": file_sha(Path(__file__)),
              "registration_sha256": file_sha(HERE / "registration.json"),
              "structural_checks": {k: tree[k] for k in ("counts", "direct_object_comparisons", "cache_canonical_sha256", "old_tree", "new_tree", "new_top_fields", "old_metadata")},
              "tree_and_object_seconds": tree_seconds, "certificate_results": [],
              "zero_new_H_gradient_kernel_Green_log_exp_Gamma_Arb_calls": True,
              "MATLAB_calls": 0, "author_functions_imported": False,
              "inherited_evidence_boundary": reg["limitations"]}
    if args.resume and (HERE / "result.json").exists():
        previous = json.loads((HERE / "result.json").read_text())
        require(previous["inputs"] == inputs and previous["code_sha256"] == result["code_sha256"], "resume receipt changed inputs or code")
        result["certificate_results"] = previous["certificate_results"]
        for row, prior in zip(certificates, result["certificate_results"]):
            require(row["node"] == prior["node"], "resume prefix order mismatch")
    certificate_started = time.monotonic()
    print(json.dumps({"pid": pid, "stage": "integer_81_point_certificates", "resumed_prefix": len(result["certificate_results"])}), flush=True)
    for row in certificates[len(result["certificate_results"]):]:
        if time.monotonic() - certificate_started >= reg["certificate_loop_seconds"]:
            break
        n = node(row["node"])
        item = verify_certificate(row, tree["old_split_sources"][n], tree["cache"], tree["roots"])
        result["certificate_results"].append(item)
        if len(result["certificate_results"]) % 128 == 0:
            write_json(HERE / "result.json", result)
    result["certificate_loop_seconds"] = time.monotonic() - certificate_started
    result["certificates_verified"] = len(result["certificate_results"])
    result["points_verified"] = sum(r["checked_points"] for r in result["certificate_results"])
    result["minimum_exact_margin"] = str(min((F(r["minimum_log_margin"]) for r in result["certificate_results"]), default=F(0)))
    result["all_anchors_strictly_inside"] = all(r["anchor_strictly_inside"] for r in result["certificate_results"])
    result["status"] = "passed" if len(result["certificate_results"]) == len(certificates) else "certificate_prefix"
    for path, expected in inputs.items():
        require(file_sha(path) == expected, "input changed during audit: " + path)
    result["input_hashes_unchanged_after_run"] = True
    result["wall_seconds"] = time.monotonic() - started
    result["peak_rss_bytes"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    write_json(HERE / "result.json", result)
    print(json.dumps({k: result[k] for k in ("pid", "status", "certificates_verified", "points_verified", "minimum_exact_margin", "tree_and_object_seconds", "certificate_loop_seconds", "wall_seconds", "peak_rss_bytes")}), flush=True)


if __name__ == "__main__":
    main()
