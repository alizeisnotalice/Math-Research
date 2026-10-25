#!/usr/bin/env python3
"""Independent Arb audit of the 007 Fejer input bit enclosures and error.

This does not import the producer's node generator.  It recomputes the
mathematical time nodes, weights, positivity, interpolation error and tail
from the declared rational parameters, then checks the producer's binary64
bit intervals against the Arb balls.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import struct
from pathlib import Path

from flint import arb, ctx


def from_bits(text: str) -> float:
    return struct.unpack(">d", int(text, 16).to_bytes(8, "big"))[0]


def check_contains(ball: arb, lo: str, hi: str) -> None:
    low, high = arb(from_bits(lo)), arb(from_bits(hi))
    assert low <= ball.lower() and high >= ball.upper() and low <= high


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source_dir", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    src = args.source_dir
    ctx.prec = 256
    declared = json.loads((src / "007_fejer_inputs.json").read_text())
    lines = (src / "007_fejer_nodes.bits").read_text().splitlines()
    assert len(lines) == 175
    rho, pi = arb(49) / 10, arb.pi()
    cpeak = arb(35509) / (arb(764411904) * pi ** 5)
    error, minimum = arb(0), None
    checked = 0
    for line in lines:
        segment, j, tlo, thi, wlo, whi = line.split()
        segment, j = int(segment), int(j)
        assert 0 <= segment < 7 and 0 <= j < 25
        a = arb(1) / (2 ** (7 - segment))
        theta = arb(2 * j + 1) * pi / 50
        t = a * (arb(3) / 2 + theta.cos() / 2)
        w = a / 25
        for k in range(1, 13):
            w -= (a / 25) * 2 * (2 * k * theta).cos() / (4 * k * k - 1)
        assert w.lower() > 0
        check_contains(t, tlo, thi)
        check_contains(w, wlo, whi)
        minimum = w if minimum is None or w.lower() < minimum.lower() else minimum
        if j == 0:
            sigma = a * (arb(3) / 2 - (rho + 1 / rho) / 4)
            error += a * 4 * cpeak * sigma ** -5 * rho ** -24 / (rho - 1)
        checked += 1
    tail = cpeak / 4
    assert declared["nodes_bits_sha256"] == hashlib.sha256((src / "007_fejer_nodes.bits").read_bytes()).hexdigest()
    result = {
        "status": "PASS_independent_arb_input_and_error_audit",
        "precision_bits": ctx.prec,
        "checked_nodes": checked,
        "all_weights_positive": True,
        "minimum_weight_arb": str(minimum),
        "Cpeak_arb": str(cpeak),
        "interpolation_error_arb": str(error),
        "tail_upper_arb": str(tail.upper()),
        "source_nodes_sha256": hashlib.sha256((src / "007_fejer_nodes.bits").read_bytes()).hexdigest(),
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "scope": "Checks only time-input enclosure and analytic error arithmetic; C++ heat arithmetic is separately replayed.",
    }
    raw = json.dumps(result, indent=2, sort_keys=True) + "\n"
    args.output.write_text(raw)
    print(raw, end="")


if __name__ == "__main__":
    main()
