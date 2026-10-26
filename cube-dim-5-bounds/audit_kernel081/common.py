#!/usr/bin/env python3
"""Independent frozen070 audit primitives; no author module imported.

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


