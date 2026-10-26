#!/usr/bin/env python3
"""Independent full-frame16corner audit of frozen077 conversion.

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
    cert, old = row["certificate"], source["failed_certificate"]
    require(digest(old) == row["source_failed_certificate_sha256"], "source failed certificate SHA mismatch")
    require(old["passed"] is False and F(old["minimum_log_margin"]) <= 0, "source split is not recorded failed")
    require(row["old_minimum_log_margin"] == old["minimum_log_margin"], "old recorded margin changed")
    for field in SOURCE_FIELDS: require(cert[field] == old[field], "inherited field changed: " + field)
    require(type(cert["passed"]) is bool and cert["mechanism"] == "full_frame_16" and cert["quadratic"] == "Q_full=I4/16" and cert["checked_points"] == 16, "incorrect full-frame certificate schema")
    anchor = cert["anchor"]; record = cache[anchor]
    require(",".join(record["point"]) == anchor, "anchor/cache key mismatch")
    point = list(map(F,record["point"]))
    require(len(point) == 5 and point[0] == 1, "anchor is not physical x1=1")
    lo, hi = physical_box(roots,n)
    require(list(map(F,row["box_lower"])) == lo and list(map(F,row["box_upper"])) == hi, "independent physical box mismatch")
    inside = all(a <= x <= b for a,x,b in zip(lo,point[1:],hi))
    require(row["anchor_inside"] == inside, "incorrect recorded anchor-inside flag")
    values = [[F(a),F(b)] for a,b in record["values"]]
    require(len(values) == 5 and all(a <= b for a,b in values) and values[0][0] > 0, "invalid positive cache intervals")
    grads = [[F(a),F(b)] for a,b in cert["loggradient_directed"]]
    require(len(grads) == 4 and all(a <= b for a,b in grads), "invalid directed gradient")
    checks = 0
    for (gl,gu),(dl,du) in zip(grads,values[1:]):
        for h in values[0]:
            for d in (dl,du):
                require(gl*h <= d <= gu*h, "physical derivative/H enclosure fails"); checks += 1
    delta=F(cert["Delta_lower"])
    offsets=[(a-x,b-x) for a,x,b in zip(lo,point[1:],hi)]
    vd=lcm(*(x.denominator for pair in offsets for x in pair)); gd=lcm(delta.denominator,*(x.denominator for pair in grads for x in pair))
    zz=[tuple(int(x*vd) for x in pair) for pair in offsets]; gg=[tuple(int(x*gd) for x in pair) for pair in grads]; di=int(delta*gd)
    denominator=gd*vd*vd; rows=[]
    for bits in product((0,1),repeat=4):
        z=tuple(zz[j][bits[j]] for j in range(4))
        lin=sum(min(a*v,b*v) for (a,b),v in zip(gg,z))
        numerator=di*vd*vd+vd*lin-256*gd*sum(v*v for v in z)
        rows.append((numerator,bits))
    minimum=min(x[0] for x in rows); margin=F(minimum,denominator)
    require(F(cert["minimum_log_margin"]) == margin and cert["passed"] == (minimum > 0), "16corner exact margin or verdict mismatch")
    worst=tuple(cert["worst_corner_bits"])
    require(len(worst)==4 and all(type(x) is int and x in (0,1) for x in worst) and (minimum,worst) in rows, "saved worst corner is not a minimum")
    require(margin >= F(old["minimum_log_margin"]), "full-frame margin unexpectedly worsens saved old bound")
    if minimum > 0: require(values[0][0] > F(1,10000), "accepted cached H lower endpoint not above target")
    return {"node":list(n),"source_split_index":row["source_split_index"],"passed":minimum>0,"minimum_log_margin":str(margin),"checked_points":16,"gradient_endpoint_cross_checks":checks,"anchor_closed_box":inside,"anchor_strictly_inside":all(a<x<b for a,x,b in zip(lo,point[1:],hi)),"physical_box_lower":list(map(str,lo)),"physical_box_upper":list(map(str,hi))}


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
            elif field in ("epsilon", "target", "model", "active", "inflight", "resume_jobs", "fingerprint", "fingerprint061", "status061", "fingerprint066", "status066", "total_kernel_attempts_completed", "total_new_anchor_starts"):
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
    if RUN_STARTED is not None and time.monotonic() - RUN_STARTED >= 540:
        raise BudgetExpired("checkpoint before600-second total wall limit")


def checkpoint():
    if RECORD is not None:
        RECORD["wall_seconds_at_checkpoint"] = time.monotonic() - RUN_STARTED
        write_json(HERE / "result.json", RECORD)


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
    initial = validate_tree(source["roots"], arrays["passed"] + arrays["queue"], set(arrays["splits"]), "source066")
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
        eligible = [i for i in active["splits"] if i not in tested and current_q[arrays["splits"][i]] > 0]
        eligible.sort(key=lambda i: (len(arrays["splits"][i][1]), -current_q[arrays["splits"][i]], *arrays["splits"][i]))
        require(all(current_q[arrays["splits"][i]] > 0 for i in eligible), "live untested candidate has no remaining queue descendant")
        selected_indices = eligible[:4096]
        selected = [{"source_split_index": i, "node": list(arrays["splits"][i]), "queue_descendants": current_q[arrays["splits"][i]]} for i in selected_indices]
        require(br["selected"] == selected, "registered block prefix differs from independent selection")
        require(summary["selected"] == len(selected) and summary["cursor"] == len(selected), "block did not finish its registered prefix")
        require(0 < len(rows) <= len(selected) <= 4096, "block test/selection cap exceeded")
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
        if number == 1 and first_receipt is not None:
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
                require(time.monotonic() - loop_start < 540, "per-block120-second certificate budget exceeded")
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
                require(item["checked_points"] == 16, "verified receipt has wrong point count")
            positive = [F(x["minimum_log_margin"]) for x in exact_results if x["passed"]]
            checked = {"block": number, "status": "verified", "certificate_input_sha256": bundle["hashes"]["certificates.json"],
                       "tested": len(rows), "accepted": len(positive), "uncertified": len(rows) - len(positive),
                       "points": 16 * len(rows), "results_sha256": digest(exact_results),
                       "minimum_accepted_exact_margin": str(min(positive)) if positive else None,
                       "certificate_loop_seconds": loop_seconds, "reused_first_block_receipt": number == 1 and first_receipt is not None}
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
        require(summary["points"] == len(rows) * 16 and summary["cache_entries"] == source["cache_entries"] and summary["cache_canonical_sha256"] == source["cache_canonical_sha256"], "block points/cache reference mismatch")
        require(summary["kernel_calls"] == 0 and summary["qualification_seconds"] <= 600 and cp["stage"] == "complete", "block finite zero-kernel qualification receipt mismatch")
        checked.update({"selected": len(selected), "skipped_covered_selected": len(selected) - len(rows), "retired_counts": retired_counts,
                        "tree": tree, "logical_frontier_sha256": current_hash, "checkpoint_sha256": bundle["hashes"]["checkpoint.json"],
                        "new_counts": {"passed": len(active["passed"]) + len(accepted), "splits": len(active["splits"]), "queue": len(active["queue"])}})
        summaries.append(checked)
        prevfrontier, prevcp = current_hash, bundle["hashes"]["checkpoint.json"]
        RECORD["block_results"] = summaries
        RECORD.pop("partial_block", None)
        RECORD["stage"] = "verified_through_block_" + str(number)
        checkpoint()
    current_q = defaultdict(int)
    for index in active["queue"]:
        r, p = arrays["queue"][index]
        for depth in range(len(p)): current_q[(r,p[:depth])] += 1
    remaining = [i for i in active["splits"] if i not in tested and current_q[arrays["splits"][i]] > 0]
    if reg["require_exhausted"]: require(not remaining, "chain stopped with untested eligible source splits")
    return {"active": active, "accepted": accepted, "tested": tested, "source_tree": initial,
            "last_checkpoint_sha256": prevcp, "final_logical_frontier_sha256": prevfrontier,
            "minimum_accepted_exact_margin": str(margin_min), "remaining_untested_eligible": len(remaining)}


def main():
    global RUN_STARTED, RECORD
    import signal
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--run',action='store_true');args=parser.parse_args()
    reg=json.loads((HERE/'registration.json').read_text())
    if not args.run: print(json.dumps({'registered':True,'scope':reg['scope'],'wall_cap':600}));return
    RUN_STARTED=time.monotonic()
    signal.signal(signal.SIGALRM,lambda *_: (_ for _ in ()).throw(BudgetExpired('590second emergency checkpoint')));signal.setitimer(signal.ITIMER_REAL,590)
    RECORD={'task':reg['task'],'pid':os.getpid(),'status':'verified_prefix_only','stage':'input_locks','code_sha256':file_sha(Path(__file__)),'registration_sha256':file_sha(HERE/'registration.json'),'block_results':[],'program_total_wall_cap_seconds':600,'author_functions_imported':False,'zero_new_H_gradient_kernel_Green_log_exp_Gamma_Arb_calls':True,'MATLAB_calls':0,'limits':reg['limits']}
    print(json.dumps({'pid':os.getpid(),'stage':'input_locks','total_wall_cap_seconds':600}),flush=True)
    try:
        block=Path(reg['block'])
        inputs={str(block/k):v for k,v in reg['block_locks'].items()}
        inputs.update({reg['source']:reg['source_sha256'],reg['author_code']:reg['author_code_sha256'],str(block.parent/'registration.json'):reg['run_registration_sha256']});inputs.update(reg['theory_locks'])
        for path,h in {**inputs,**reg['protected_old_locks']}.items():require(file_sha(path)==h,'frozen/protected hash mismatch: '+path)
        RECORD['inputs']=inputs;RECORD['protected_old_locks']=reg['protected_old_locks']
        docs={k:json.loads((block/k).read_text(),object_pairs_hook=no_duplicates) for k in reg['block_locks']}
        br,cp,summary=(docs[k] for k in ('registration.json','checkpoint.json','summary.json'))
        for name in ('registration','certificates','retired','summary'):require(cp[name+'_sha256']==reg['block_locks'][name+'.json'],'checkpoint artifact link mismatch')
        require(br['source_state']==reg['source'] and br['source_state_sha256']==reg['source_sha256'] and br['code_sha256']==reg['author_code_sha256'] and br['run_registration_sha256']==reg['run_registration_sha256'],'frozen source/code links mismatch')
        require(br['max_tested_per_process']==4096 and br['total_wall_seconds']==600 and br['normal_checkpoint_seconds']==540 and br['alarm_seconds']==590,'finite process limits mismatch')
        runreg=json.loads((block.parent/'registration.json').read_text());require(runreg['new_kernel_chain_authorized'] is False,'new kernel authorization unexpectedly present')
        rows=docs['certificates.json'];bundle={'number':1,'registration':br,'checkpoint':cp,'rows':rows,'retired':docs['retired.json'],'summary':summary,'hashes':reg['block_locks']}
        RECORD['stage']='source_stream';checkpoint();print(json.dumps({'pid':os.getpid(),'stage':'source_stream','tested_rows':len(rows)}),flush=True)
        source=read_source(reg['source'],rows,reg)
        chain=replay_chain(source,[bundle],reg,None)
        RECORD.update({'source_cache_canonical_sha256':source['cache_canonical_sha256'],'source_cache_entries':source['cache_entries'],'source_tree':chain['source_tree'],'tested':len(chain['tested']),'accepted':len(chain['accepted']),'uncertified':len(chain['tested'])-len(chain['accepted']),'points_verified':16*len(chain['tested']),'minimum_accepted_exact_margin':chain['minimum_accepted_exact_margin'],'remaining_untested_eligible':chain['remaining_untested_eligible'],'new_counts':RECORD['block_results'][-1]['new_counts'],'last_checkpoint_sha256':chain['last_checkpoint_sha256'],'logical_frontier_sha256':chain['final_logical_frontier_sha256']})
        for path,h in {**inputs,**reg['protected_old_locks']}.items():require(file_sha(path)==h,'input/protected changed during audit: '+path)
        RECORD['input_and_protected_locks_unchanged_after_run']=True;RECORD['status']='passed_first_block_only';RECORD['stage']='complete'
    except BudgetExpired as exc:signal.setitimer(signal.ITIMER_REAL,0);RECORD['deadline_reason']=str(exc)
    finally:
        signal.setitimer(signal.ITIMER_REAL,0);RECORD['wall_seconds']=time.monotonic()-RUN_STARTED;RECORD['peak_rss_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss;checkpoint()
    print(json.dumps({k:RECORD.get(k) for k in ('pid','status','tested','accepted','uncertified','points_verified','new_counts','remaining_untested_eligible','wall_seconds','peak_rss_bytes')}),flush=True)
    print(json.dumps({'result_sha256':sha256((HERE/'result.json').read_bytes()).hexdigest()}),flush=True)


if __name__=='__main__':main()


